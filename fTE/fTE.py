import os
import sys
import glob
import argparse
import subprocess
import shutil
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

def generate_non_overlapping_windows(input_dir, work_dir, window_size=10000, threads=8):
    """오버랩 없이 유전체를 분할하여 개별 FASTA 파일로 저장 (고성능/저메모리 버전)"""
    print(f"[*] Generating {window_size}bp non-overlapping windows...")
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
                
    print(f"    -> Total {window_counter} windows generated.")
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
    for i in keep_indices:
        row = f"{taxa_names[i]} " + " ".join(f"{dists[i][j]:.6f}" for j in keep_indices)
        phylip_lines.append(row)
        
    matrix_file = os.path.join(work_dir, "distance_matrix.phy")
    with open(matrix_file, "w") as f:
        f.write('\n'.join(phylip_lines) + '\n')
        
    return matrix_file, [fasta_files[i] for i in keep_indices]



def find_founder_candidates(tree_file, id_map):
    """FastME 트리 파싱: 주요 클레이드(TE 패밀리) 및 중심 서열(Founder 후보) 추출"""
    print("[*] Parsing FastME tree to identify Founder Candidates...")
    tree = Phylo.read(tree_file, "newick")
    candidates = []
    
    major_clades = tree.root.clades
    for idx, clade in enumerate(major_clades):
        if clade.is_terminal(): continue
            
        leaves = clade.get_terminals()
        if len(leaves) < 3: continue # 너무 작은 클러스터는 무시
        
        best_leaf = None
        min_dist = float('inf')
        
        for leaf in leaves:
            dist = tree.distance(clade, leaf)
            if dist < min_dist:
                min_dist = dist
                best_leaf = leaf
                
        if best_leaf:
            short_id = os.path.basename(best_leaf.name).split('.')[0]
            candidates.append({
                "clade_id": f"Clade_{idx+1}",
                "founder_short_id": short_id,
                "founder_real_id": id_map.get(short_id, short_id),
                "cluster_leaves": [l.name for l in leaves],
                "cluster_size": len(leaves)
            })
            
    # 클러스터 크기 순 정렬 후 반환
    return sorted(candidates, key=lambda x: x['cluster_size'], reverse=True)



def main():
    parser = argparse.ArgumentParser(description="TE Founder Detection with Parabola")
    parser.add_argument("-i", "--input_dir", required=True, help="Target FASTA directory")
    parser.add_argument("-w", "--window", type=int, default=10000, help="Non-overlapping window size (bp)")
    parser.add_argument("-k", "--kmer", type=int, default=21, help="K-mer size for Parabola")
    parser.add_argument("-c", "--scale", type=int, default=1000, help="FracMinHash scale for Parabola")
    parser.add_argument("-p", "--threads", type=int, default=16, help="Number of threads")
    parser.add_argument("-f", "--filter_pct", type=float, default=1.0, help="Genome coverage filter threshold in percent (default: 1.0 for 1 percent)")
    args = parser.parse_args()

    input_name = os.path.basename(args.input_dir.rstrip('/'))
    work_dir = f"te_workspace_{input_name}"
    # 이전 불완전한 런의 찌꺼기를 제거하기 위해 시작 시 기존 워크스페이스 디렉토리 삭제 후 생성
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir, exist_ok=True)
    
    # 1. 오버랩 없는 윈도우 청크 생성
    chunk_files, id_map = generate_non_overlapping_windows(args.input_dir, work_dir, args.window, args.threads)
    
    # 2. Parabola 실행 (Sketch & Triangle)
    matrix_file, kept_files = run_parabola(chunk_files, work_dir, id_map, args.filter_pct, args.kmer, args.scale, args.threads)
    
    # 3. FastME 계통수 구축
    print("[*] Building phylogenetic tree with FastME...")
    tree_file = os.path.join(work_dir, "tree.nwk")
    run_cmd(["fastme", "-i", matrix_file, "-o", tree_file, "-T", "16", "-s"])
    
    # 4. 트리 파싱 및 후보 클러스터 추출
    candidates = find_founder_candidates(tree_file, id_map)
    
    print("\n[+] Top TE Clades Detected:")
    for c in candidates:
        print(f" - {c['clade_id']} (Size: {c['cluster_size']}) -> Founder Candidate: {c['founder_real_id']}")
        


    # 6. 최종 결과를 fTE_results/ 디렉토리에 저장
    results_dir = "fTE_results"
    os.makedirs(results_dir, exist_ok=True)

    # 6.1 최종 TE candidate의 위치를 <genome>_fTE.bed 파일로 작성
    genome_beds = {}
    for c in candidates:
        short_id = c['founder_short_id']
        real_id = id_map.get(short_id)
        if real_id and '|' in real_id:
            parts = real_id.split('|')
            genome_id = parts[0]
            chrom = parts[1]
            coords = parts[2]
            start, end = coords.split('-')
            
            if genome_id not in genome_beds:
                genome_beds[genome_id] = []
            genome_beds[genome_id].append((chrom, start, end, short_id))
            
    for genome_id, regions in genome_beds.items():
        bed_file = os.path.join(results_dir, f"{genome_id}_fTE.bed")
        print(f"[*] Saving BED file: {bed_file}")
        with open(bed_file, "w") as f:
            for chrom, start, end, name in regions:
                f.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    # 6.2 살아남은 window의 위치를 <genome>_window.bed 파일로 작성
    genome_windows = {}
    for path in kept_files:
        short_id = os.path.basename(path).split('.')[0]
        real_id = id_map.get(short_id)
        if real_id and '|' in real_id:
            parts = real_id.split('|')
            genome_id = parts[0]
            chrom = parts[1]
            coords = parts[2]
            start, end = coords.split('-')
            
            if genome_id not in genome_windows:
                genome_windows[genome_id] = []
            genome_windows[genome_id].append((chrom, start, end, short_id))
            
    for genome_id, regions in genome_windows.items():
        bed_file = os.path.join(results_dir, f"{genome_id}_window.bed")
        print(f"[*] Saving BED file: {bed_file}")
        with open(bed_file, "w") as f:
            for chrom, start, end, name in sorted(regions, key=lambda x: (x[0], int(x[1]))):
                f.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    final_tree_file = os.path.join(results_dir, "tree.nwk")
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