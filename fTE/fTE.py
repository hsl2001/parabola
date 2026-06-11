import os
import sys
import glob
import argparse
import subprocess
import shutil
import collections
import heapq
from Bio import SeqIO, Phylo
from concurrent.futures import ThreadPoolExecutor

def run_cmd(cmd_list, capture_out=False):
    """명령어 실행 및 결과 캡처 헬퍼 함수"""
    result = subprocess.run(cmd_list, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[Error] Command failed: {' '.join(cmd_list)}\n{result.stderr}")
        sys.exit(1)
    return result.stdout if capture_out else None

def stream_windows(fasta_path, window_size):
    step_size = window_size // 2
    with open(fasta_path, 'r') as f:
        current_chrom = None
        buffer = []
        buffer_len = 0
        chrom_pos = 0
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current_chrom = line[1:].split()[0]
                buffer = []
                buffer_len = 0
                chrom_pos = 0
            else:
                buffer.append(line)
                buffer_len += len(line)
                while buffer_len >= window_size:
                    full_str = "".join(buffer)
                    window_seq = full_str[:window_size]
                    yield current_chrom, chrom_pos, window_seq
                    
                    rest = full_str[step_size:]
                    buffer = [rest]
                    buffer_len = len(rest)
                    chrom_pos += step_size

def write_chunk(chunk_file, short_id, seq):
    with open(chunk_file, "w") as f:
        f.write(f">{short_id}\n{seq}\n")

def generate_windows(input_dir, work_dir, window_size=10000, threads=8):
    """윈도우 크기의 절반(step size) 간격으로 유전체를 분할하여 개별 FASTA 파일로 저장 (고성능/저메모리 버전)"""
    print(f"[*] Generating {window_size}bp windows with {window_size//2}bp step...")
    fasta_files = glob.glob(os.path.join(input_dir, "*.fa*"))
    
    chunks_dir = os.path.join(work_dir, "chunks")
    os.makedirs(chunks_dir, exist_ok=True)
    
    id_map = {}
    chunk_files = []
    window_counter = 0
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = []
        for fasta in fasta_files:
            genome_id = os.path.basename(fasta).split('.')[0]
            for chrom, start, seq in stream_windows(fasta, window_size):
                window_counter += 1
                
                short_id = f"W{window_counter:05d}"
                real_id = f"{genome_id}|{chrom}|{start}-{start + window_size}"
                
                chunk_file = os.path.join(chunks_dir, f"{short_id}.fasta")
                id_map[short_id] = real_id
                chunk_files.append(chunk_file)
                
                futures.append(executor.submit(write_chunk, chunk_file, short_id, seq))
                
        # Wait for all file writes to complete
        for future in futures:
            future.result()
                
    print(f"    -> Total {window_counter} valid windows generated.")
    return chunk_files, id_map

def run_parabola(fasta_files, work_dir, id_map, copy_threshold, k=21, scale=1000, threads=8):
    """Parabola sketch 및 triangle을 이용한 거리 행렬 계산"""
    print(f"[*] Running Parabola (k={k}, s={scale}, threads={threads})...")
    
    # 1. Parabola Sketch (배칭 처리하여 인자 제한 우회)
    batch_size = 5000
    for i in range(0, len(fasta_files), batch_size):
        batch = fasta_files[i:i+batch_size]
        cmd_sketch = ["./parabola", "sketch", "-k", str(k), "-s", str(scale), "-p", str(threads)] + batch
        run_cmd(cmd_sketch)
    
    # 2. Parabola Triangle
    sketch_files = [f + ".parabola" for f in fasta_files]
    triangle_out = os.path.join(work_dir, "triangle.tsv")
    
    print("[*] Calculating distance matrix ...")
    with open(triangle_out, "w") as f_out:
        cmd_triangle = ["./parabola", "triangle"] + sketch_files
        res = subprocess.run(cmd_triangle, stdout=f_out, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            print(f"[Error] Command failed: {' '.join(cmd_triangle[:5])} ...\n{res.stderr}")
            sys.exit(1)
            
    # 3. 1차 패스: 동일 유전체 내에서 distance < 1.0인 카운트 계산
    print("[*] Pass 1: Filtering chunks...")
    genome_ids = [id_map[os.path.basename(p).split('.')[0]].split('|')[0] for p in fasta_files]
    counts = [0] * len(fasta_files)
    taxa_names = []
    
    with open(triangle_out, 'r') as f:
        first_line = f.readline()
        if not first_line:
            return ""
        num_taxa = int(first_line.strip())
        
        for i in range(num_taxa):
            line = f.readline()
            parts = line.strip().split('\t')
            taxa_names.append(parts[0])
            g_i = genome_ids[i]
            for j in range(1, len(parts)):
                real_j = j - 1
                if genome_ids[real_j] == g_i:
                    val = float(parts[j].split(',')[0])
                    if 0.0 < val < 1.0:
                        counts[i] += 1
                        counts[real_j] += 1
                        
    keep_indices = [i for i, c in enumerate(counts) if c >= copy_threshold]
    
    if len(keep_indices) < 3:
        print("[*] Warning: Too few chunks kept after filtering. Keeping all chunks.")
        keep_indices = list(range(num_taxa))
    else:
        print(f"[*] Filtered out {num_taxa - len(keep_indices)} unique chunks. Keeping {len(keep_indices)} repetitive chunks.")
        
    # 4. 2차 패스: 살아남은 chunk들만의 거리 행렬 추출
    print("[*] Pass 2: Extracting distance matrix for kept chunks...")
    keep_set = set(keep_indices)
    keep_idx_map = {old_i: new_i for new_i, old_i in enumerate(keep_indices)}
    num_kept = len(keep_indices)
    
    dists_kept = [[0.0] * num_kept for _ in range(num_kept)]
    
    with open(triangle_out, 'r') as f:
        f.readline() # skip num_taxa
        for i in range(num_taxa):
            line = f.readline()
            if i not in keep_set:
                continue
                
            parts = line.strip().split('\t')
            for j in range(1, len(parts)):
                real_j = j - 1
                if real_j in keep_set:
                    val = float(parts[j].split(',')[0])
                    new_i = keep_idx_map[i]
                    new_j = keep_idx_map[real_j]
                    dists_kept[new_i][new_j] = val
                    dists_kept[new_j][new_i] = val
                    
    # 5. FastME 포맷으로 저장 및 dist_dict 생성
    phylip_lines = [str(num_kept)]
    dist_dict = collections.defaultdict(dict)
    for new_i, old_i in enumerate(keep_indices):
        name = taxa_names[old_i]
        row = f"{name} " + " ".join(f"{dists_kept[new_i][new_j]:.6f}" for new_j in range(num_kept))
        phylip_lines.append(row)
        
        for new_j, old_j in enumerate(keep_indices):
            dist_dict[name][taxa_names[old_j]] = dists_kept[new_i][new_j]
            
    matrix_file = os.path.join(work_dir, "distance_matrix.phy")
    with open(matrix_file, "w") as f:
        f.write('\n'.join(phylip_lines) + '\n')
        
    return matrix_file, [fasta_files[i] for i in keep_indices], dist_dict



def calc_silhouette_score(clusters, dist_dict):
    valid_clusters = [c for c in clusters if c]
    if len(valid_clusters) < 2: return -1.0

    total_s = 0.0
    num_leaves = sum(len(c) for c in valid_clusters)
    
    for cidx, comp in enumerate(valid_clusters):
        for leaf in comp:
            if len(comp) == 1:
                total_s += 0.0
            else:
                a_i = sum(dist_dict[leaf][other] for other in comp if other != leaf) / (len(comp) - 1)
                b_i = min((sum(dist_dict[leaf][other] for other in other_comp) / len(other_comp) 
                           for other_cidx, other_comp in enumerate(valid_clusters) if cidx != other_cidx), default=float('inf'))
                
                max_val = max(a_i, b_i)
                total_s += (b_i - a_i) / max_val if max_val > 0 else 0.0
                
    return total_s / num_leaves if num_leaves > 0 else -1.0


def get_connected_components(nodes, edges, return_weighted_adj=False):
    adj = collections.defaultdict(list)
    for u, v, w in edges:
        adj[u].append((v, w) if return_weighted_adj else v)
        adj[v].append((u, w) if return_weighted_adj else u)
        
    visited = set()
    components = []
    for node in nodes:
        if node not in visited:
            comp = set()
            q = collections.deque([node])
            visited.add(node)
            while q:
                curr = q.popleft()
                comp.add(curr)
                for neighbor in adj[curr]:
                    n_node = neighbor[0] if return_weighted_adj else neighbor
                    if n_node not in visited:
                        visited.add(n_node)
                        q.append(n_node)
            components.append(comp)
    return components, adj


def find_founder_candidates(tree_file, id_map, dist_dict, k_clusters="auto", max_k=20):
    """FastME 트리 파싱: 최장 거리 간선들을 제거하여 K개의 클러스터로 분리 및 중심 서열 추출"""
    print(f"[*] Parsing FastME tree for K-clusters cut...")
    tree = Phylo.read(tree_file, "newick")
    
    edges = []
    def traverse(node):
        for child in node.clades:
            edges.append((node, child, max(0.0, child.branch_length or 0.0)))
            traverse(child)
    traverse(tree.root)
    
    all_nodes = set()
    for u, v, w in edges:
        all_nodes.update([u, v])
    if not all_nodes: all_nodes.add(tree.root)
        
    leaves_count = sum(1 for n in all_nodes if n.name and n.is_terminal())
    sorted_edges = sorted(edges, key=lambda x: x[2], reverse=True)
    
    def evaluate_k(k):
        keep_edges = sorted_edges[k-1:] if k > 1 else sorted_edges
        comps, _ = get_connected_components(all_nodes, keep_edges)
        return [[n.name for n in c if n.name and n.is_terminal()] for c in comps]

    if str(k_clusters).lower() == "auto":
        print(f"    -> Finding optimal K using Silhouette Score (testing K=2 to {max_k})...")
        best_k, best_score = 1, -float('inf')
        max_k = min(max_k, leaves_count - 1)
        
        for k in range(2, max_k + 1):
            score = calc_silhouette_score(evaluate_k(k), dist_dict)
            print(f"       [K={k}] Silhouette Score: {score:.4f}")
            if score > best_score:
                best_score, best_k = score, k
                
        if best_k == 1: print("    -> Not enough leaves to cluster. Defaulting to K=1")
        else: print(f"    => Optimal K selected: {best_k} (Score: {best_score:.4f})")
    else:
        try:
            best_k = int(k_clusters)
            print(f"    -> Using specified K={best_k}")
        except ValueError:
            print("[Error] Invalid k_clusters value.")
            sys.exit(1)
            
    # 최적 K 결정 이후 진짜 외군(Basal leaf) 탐색
    keep_edges = sorted_edges[best_k-1:] if best_k > 1 else sorted_edges
    components, _ = get_connected_components(all_nodes, keep_edges)
            
    candidates = []
    all_leaf_names = set(n.name for n in all_nodes if n.name and n.is_terminal())
    
    for idx, comp_nodes in enumerate(components):
        comp_leaves = [n.name for n in comp_nodes if n.name and n.is_terminal()]
        if not comp_leaves: continue
            
        outside_leaves = [l for l in all_leaf_names if l not in comp_leaves]
                
        if not outside_leaves:
            # K=1 이거나 외군이 없는 경우: 클러스터 내부의 중심점(Medoid) 사용
            best_leaf = min(comp_leaves, key=lambda l: 0 if len(comp_leaves) == 1 else sum(dist_dict[l][o] for o in comp_leaves if o != l) / (len(comp_leaves)-1))
        else:
            # K>1: 클러스터 외부 서열들(Outgroup)과의 평균 거리가 가장 짧은 서열을 이 클러스터의 뿌리(Basal)로 지정
            best_leaf = min(comp_leaves, key=lambda l: sum(dist_dict[l][o] for o in outside_leaves) / len(outside_leaves))

        if best_leaf:
            short_id = best_leaf.split('.')[0]
            candidates.append({
                "clade_id": f"Cluster_{idx+1}",
                "founder_short_id": short_id,
                "founder_real_id": id_map.get(short_id, short_id),
                "cluster_size": len(comp_leaves)
            })
            
    return sorted(candidates, key=lambda x: x['cluster_size'], reverse=True)


def save_regions_to_bed(item_list, id_map, results_dir, suffix):
    """주어진 리스트를 기반으로 파싱하여 BED 파일 생성"""
    genome_beds = collections.defaultdict(list)
    for item in item_list:
        short_id = os.path.basename(item).split('.')[0]
        real_id = id_map.get(short_id)
        if real_id and '|' in real_id:
            genome_id, chrom, coords = real_id.split('|')
            start, end = coords.split('-')
            genome_beds[genome_id].append((chrom, int(start), int(end), short_id))
            
    for genome_id, regions in genome_beds.items():
        bed_file = os.path.join(results_dir, f"{genome_id}_{suffix}.bed")
        print(f"[*] Saving BED file: {bed_file}")
        with open(bed_file, "w") as f:
            for chrom, start, end, name in sorted(regions, key=lambda x: (x[0], x[1])):
                f.write(f"{chrom}\t{start}\t{end}\t{name}\n")



def main():
    parser = argparse.ArgumentParser(description="TE Founder Detection with Parabola")
    parser.add_argument("-i", "--input_dir", required=True, help="Target FASTA directory")
    parser.add_argument("-w", "--window", type=int, default=10000, help="Non-overlapping window size (bp)")
    parser.add_argument("-k", "--kmer", type=int, default=21, help="K-mer size for Parabola")
    parser.add_argument("-c", "--scale", type=int, default=1000, help="FracMinHash scale for Parabola")
    parser.add_argument("-p", "--threads", type=int, default=16, help="Number of threads")
    parser.add_argument("-y", "--copy", type=int, default=40, help="Minimum copy number threshold for a TE chunk (default: 3)")
    parser.add_argument("-k_cls", "--k_clusters", type=str, default="auto", help="Number of clusters to cut the tree into, or 'auto' to find optimal K (default: auto)")
    parser.add_argument("-max_k", "--max_k", type=int, default=20, help="Maximum K to test when k_clusters is 'auto' (default: 20)")
    args = parser.parse_args()

    input_name = os.path.basename(args.input_dir.rstrip('/'))
    work_dir = f"te_workspace_{input_name}"
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir, exist_ok=True)
    
    # 1. 윈도우 청크 생성
    chunk_files, id_map = generate_windows(args.input_dir, work_dir, args.window, args.threads)
    
    # 2. Parabola 실행 (Sketch & Triangle)
    matrix_file, kept_files, dist_dict = run_parabola(chunk_files, work_dir, id_map, args.copy, args.kmer, args.scale, args.threads)
    
    # 3. FastME 계통수 구축
    print("[*] Building phylogenetic tree with FastME...")
    tree_file = os.path.join(work_dir, "tree.nwk")
    run_cmd(["fastme", "-i", matrix_file, "-o", tree_file, "-T", "16"])
    
    # 4. 트리 파싱 및 후보 클러스터 추출
    candidates = find_founder_candidates(tree_file, id_map, dist_dict, k_clusters=args.k_clusters, max_k=args.max_k)
    
    print("\n[+] Top TE Clades Detected:")
    for c in candidates:
        print(f" - {c['clade_id']} (Size: {c['cluster_size']}) -> Founder Candidate: {c['founder_real_id']}")
        


    # 6. 최종 결과를 fTE_results/ 디렉토리에 저장
    results_dir = "fTE_results"
    os.makedirs(results_dir, exist_ok=True)

    # 6.1 최종 TE candidate의 위치를 <genome>_fTE.bed 파일로 작성
    save_regions_to_bed([c['founder_short_id'] for c in candidates], id_map, results_dir, "fTE")

    # 6.2 살아남은 window의 위치를 <genome>_window.bed 파일로 작성
    save_regions_to_bed(kept_files, id_map, results_dir, "window")

    final_tree_file = os.path.join(results_dir, f"{input_name}_tree.nwk")
    if os.path.exists(tree_file):
        with open(tree_file, "r") as f:
            tree_content = f.read()
        tree_content = tree_content.replace(f"{work_dir}/chunks/", "")
        tree_content = tree_content.replace(".fasta", "")
        with open(tree_file, "w") as f:
            f.write(tree_content)
        shutil.move(tree_file, final_tree_file)
        print(f"[*] Tree file saved to: {final_tree_file}")
        
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
        print(f"[*] Cleaned up workspace directory: {work_dir}")

if __name__ == "__main__":
    main()