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
                    
                    rest = full_str[window_size:]
                    buffer = [rest]
                    buffer_len = len(rest)
                    chrom_pos += window_size

def write_chunk(chunk_file, short_id, seq):
    with open(chunk_file, "w") as f:
        f.write(f">{short_id}\n{seq}\n")

def generate_non_overlapping_windows(input_dir, work_dir, window_size=10000, threads=8, complexity_threshold=0.5):
    """오버랩 없이 유전체를 분할하여 개별 FASTA 파일로 저장 (고성능/저메모리 버전)"""
    print(f"[*] Generating {window_size}bp non-overlapping windows (complexity >= {complexity_threshold})...")
    fasta_files = glob.glob(os.path.join(input_dir, "*.fa*"))
    
    chunks_dir = os.path.join(work_dir, "chunks")
    os.makedirs(chunks_dir, exist_ok=True)
    
    id_map = {}
    chunk_files = []
    window_counter = 0
    filtered_counter = 0
    
    def get_complexity(seq, k=15):
        if len(seq) < k: return 0.0
        kmers = set(seq[i:i+k] for i in range(len(seq)-k+1))
        return len(kmers) / (len(seq) - k + 1)
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = []
        for fasta in fasta_files:
            genome_id = os.path.basename(fasta).split('.')[0]
            for chrom, start, seq in stream_windows(fasta, window_size):
                window_counter += 1
                
                comp = get_complexity(seq)
                if comp < complexity_threshold:
                    filtered_counter += 1
                    continue
                
                short_id = f"W{window_counter:05d}"
                real_id = f"{genome_id}|{chrom}|{start}-{start + window_size}"
                
                chunk_file = os.path.join(chunks_dir, f"{short_id}.fasta")
                id_map[short_id] = real_id
                chunk_files.append(chunk_file)
                
                futures.append(executor.submit(write_chunk, chunk_file, short_id, seq))
                
        # Wait for all file writes to complete
        for future in futures:
            future.result()
                
    print(f"    -> Total {window_counter} windows scanned, {filtered_counter} TR windows dropped.")
    print(f"    -> {len(chunk_files)} valid windows generated.")
    return chunk_files, id_map

def run_parabola(fasta_files, work_dir, id_map, filter_pct, k=21, scale=1000, threads=8):
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
    cmd_triangle = ["./parabola", "triangle"] + sketch_files
    matrix_out = run_cmd(cmd_triangle, capture_out=True)
    
    # FastME는 full square matrix만 지원하므로 lower-triangular를 square matrix로 변환
    lines = matrix_out.strip().split('\n')
    if not lines:
        return ""
        
    try:
        num_taxa = int(lines[0].strip())
    except ValueError:
        print("[Error] Invalid matrix size from parabola triangle")
        sys.exit(1)
        
    taxa_names = []
    dists = [[0.0] * num_taxa for _ in range(num_taxa)]
    
    for i in range(num_taxa):
        line = lines[i + 1]
        parts = line.split('\t')
        taxa_names.append(parts[0])
        for j in range(1, len(parts)):
            val = parts[j].split(',')[0]
            d = float(val)
            dists[i][j - 1] = d
            dists[j - 1][i] = d
            
    # 3. 한 유전체 내에서 지정된 비율(filter_pct%)의 chunk에서 검출되지 않는 chunk(즉, 동일 유전체 내에서 distance < 1.0인 chunk가 지정 비율 미만인 경우)는 필터링
    keep_indices = []
    genome_to_chunks = {}
    for idx, path in enumerate(fasta_files):
        short_id = os.path.basename(path).split('.')[0]
        genome_id = id_map[short_id].split('|')[0]
        if genome_id not in genome_to_chunks:
            genome_to_chunks[genome_id] = []
        genome_to_chunks[genome_id].append(idx)
        
    for genome_id, indices in genome_to_chunks.items():
        n_g = len(indices)
        threshold = (filter_pct / 100.0) * n_g
        for i in indices:
            detected_count = sum(1 for j in indices if dists[i][j] < 1.0)
            if detected_count >= threshold:
                keep_indices.append(i)
                
    if len(keep_indices) < 3:
        print("[*] Warning: Too few chunks kept after filtering. Keeping all chunks.")
        keep_indices = list(range(num_taxa))
    else:
        print(f"[*] Filtered out {num_taxa - len(keep_indices)} unique chunks. Keeping {len(keep_indices)} repetitive chunks.")
        
    num_kept = len(keep_indices)
    phylip_lines = [str(num_kept)]
    dist_dict = {}
    for i in keep_indices:
        row = f"{taxa_names[i]} " + " ".join(f"{dists[i][j]:.6f}" for j in keep_indices)
        phylip_lines.append(row)
        
        dist_dict[taxa_names[i]] = {}
        for j in keep_indices:
            dist_dict[taxa_names[i]][taxa_names[j]] = dists[i][j]
            
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
    cut_edges = sorted_edges[:best_k-1] if best_k > 1 else []
    components, adj_weighted = get_connected_components(all_nodes, keep_edges, return_weighted_adj=True)
            
    candidates = []
    for idx, comp_nodes in enumerate(components):
        comp_leaves = [n for n in comp_nodes if n.name and n.is_terminal()]
        if not comp_leaves: continue
            
        cut_nodes = []
        for u, v, w in cut_edges:
            if u in comp_nodes and v not in comp_nodes: cut_nodes.append(u)
            elif v in comp_nodes and u not in comp_nodes: cut_nodes.append(v)
                
        best_leaf = None
        if not cut_nodes:
            # K=1 이거나 외부 연결점이 없는 경우 중심점(Medoid) 사용
            leaf_names = [n.name for n in comp_leaves]
            best_leaf = min(leaf_names, key=lambda l: 0 if len(leaf_names) == 1 else sum(dist_dict[l][o] for o in leaf_names if o != l) / (len(leaf_names)-1))
        else:
            # 외부(cut edge)와 연결된 노드들부터 시작하여 BFS/Dijkstra
            pq = [(0, 0.0, id(cn), cn) for cn in cut_nodes]
            shortest_paths = {}
            while pq:
                tdist, bdist, _, curr = heapq.heappop(pq)
                if curr in shortest_paths: continue
                shortest_paths[curr] = (tdist, bdist)
                
                for neighbor, w in adj_weighted[curr]:
                    if neighbor not in shortest_paths:
                        heapq.heappush(pq, (tdist + 1, bdist + w, id(neighbor), neighbor))
                        
            best_leaf_node = min(comp_leaves, key=lambda n: shortest_paths.get(n, (float('inf'), float('inf'))))
            best_leaf = best_leaf_node.name

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
    parser.add_argument("-f", "--filter_pct", type=float, default=0.5, help="Genome coverage filter threshold in percent (default: 0.5)")
    parser.add_argument("-k_cls", "--k_clusters", type=str, default="auto", help="Number of clusters to cut the tree into, or 'auto' to find optimal K (default: auto)")
    parser.add_argument("-max_k", "--max_k", type=int, default=20, help="Maximum K to test when k_clusters is 'auto' (default: 20)")
    parser.add_argument("-x", "--complexity", type=float, default=0.5, help="K-mer complexity threshold to drop tandem repeats (default: 0.5)")
    args = parser.parse_args()

    input_name = os.path.basename(args.input_dir.rstrip('/'))
    work_dir = f"te_workspace_{input_name}"
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir, exist_ok=True)
    
    # 1. 오버랩 없는 윈도우 청크 생성
    chunk_files, id_map = generate_non_overlapping_windows(args.input_dir, work_dir, args.window, args.threads, args.complexity)
    
    # 2. Parabola 실행 (Sketch & Triangle)
    matrix_file, kept_files, dist_dict = run_parabola(chunk_files, work_dir, id_map, args.filter_pct, args.kmer, args.scale, args.threads)
    
    # 3. FastME 계통수 구축
    print("[*] Building phylogenetic tree with FastME...")
    tree_file = os.path.join(work_dir, "tree.nwk")
    run_cmd(["fastme", "-i", matrix_file, "-o", tree_file, "-T", "16", "-s"])
    
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