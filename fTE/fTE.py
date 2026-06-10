import os
import sys
import glob
import argparse
import subprocess
import shutil
from Bio import SeqIO, Phylo

def run_cmd(cmd_list, capture_out=False):
    """명령어 실행 및 결과 캡처 헬퍼 함수"""
    result = subprocess.run(cmd_list, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[Error] Command failed: {' '.join(cmd_list)}\n{result.stderr}")
        sys.exit(1)
    return result.stdout if capture_out else None

def generate_non_overlapping_windows(input_dir, work_dir, window_size=10000):
    """오버랩 없이 유전체를 분할하여 개별 FASTA 파일로 저장"""
    print(f"[*] Generating {window_size}bp non-overlapping windows...")
    fasta_files = glob.glob(os.path.join(input_dir, "*.fa*"))
    
    chunks_dir = os.path.join(work_dir, "chunks")
    os.makedirs(chunks_dir, exist_ok=True)
    
    id_map = {}
    chunk_files = []
    window_counter = 0
    
    for fasta in fasta_files:
        genome_id = os.path.basename(fasta).split('.')[0]
        for record in SeqIO.parse(fasta, "fasta"):
            seq_len = len(record.seq)
            
            # 오버랩 없이(window_size 간격으로) 분할
            for start in range(0, seq_len - window_size + 1, window_size):
                end = start + window_size
                window_seq = record.seq[start:end]
                
                window_counter += 1
                short_id = f"W{window_counter:05d}"
                real_id = f"{genome_id}|{record.id}|{start}-{end}"
                
                chunk_file = os.path.join(chunks_dir, f"{short_id}.fasta")
                with open(chunk_file, "w") as f:
                    f.write(f">{short_id}\n{window_seq}\n")
                
                id_map[short_id] = real_id
                chunk_files.append(chunk_file)
                
    print(f"    -> Total {window_counter} windows generated.")
    return chunk_files, id_map

def run_parabola(fasta_files, work_dir, id_map, filter_pct, k=21, scale=1000, threads=8):
    """Parabola sketch 및 triangle을 이용한 거리 행렬 계산"""
    print(f"[*] Running Parabola (k={k}, s={scale}, threads={threads})...")
    
    # 1. Parabola Sketch
    cmd_sketch = ["./parabola", "sketch", "-k", str(k), "-s", str(scale), "-p", str(threads)] + fasta_files
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

def parse_phylip_matrix_mean(matrix_str):
    """Phylip 행렬 텍스트를 파싱하여 평균 거리 계산 (노이즈 감소 판단용)"""
    lines = matrix_str.strip().split('\n')
    if len(lines) < 2: return 1.0
    
    distances = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) > 1:
            for x in parts[1:]:
                # Comma-separated values (e.g., "0.123,0.456")
                val = x.split(',')[0]
                try:
                    d = float(val)
                    if d > 0:
                        distances.append(d)
                except ValueError:
                    pass
            
    return sum(distances) / len(distances) if distances else 1.0

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

def iterative_boundary_carving(clade_info, work_dir, orig_chunks_dir, k, scale, threads, step_size=100, min_len=2000):
    """특정 클러스터의 시퀀스 양 끝단을 점진적으로 깎아내며 노이즈를 제거"""
    clade_id = clade_info['clade_id']
    leaves = clade_info['cluster_leaves']
    print(f"\n[*] Iterative Boundary Carving for {clade_id} (Size: {len(leaves)})")
    
    carve_dir = os.path.join(work_dir, f"carving_{clade_id}")
    os.makedirs(carve_dir, exist_ok=True)
    
    # 초기 서열 복사
    current_fastas = []
    for leaf in leaves:
        leaf_short_id = os.path.basename(leaf).split('.')[0]
        src = os.path.join(orig_chunks_dir, f"{leaf_short_id}.fasta")
        dst = os.path.join(carve_dir, f"{leaf_short_id}_iter0.fasta")
        shutil.copy(src, dst)
        current_fastas.append(dst)
        
    iteration = 0
    
    # 초기 거리 측정
    cmd_sketch = ["./parabola", "sketch", "-k", str(k), "-s", str(scale), "-p", str(threads)] + current_fastas
    run_cmd(cmd_sketch)
    matrix_out = run_cmd(["./parabola", "triangle"] + [f + ".parabola" for f in current_fastas], capture_out=True)
    best_dist = parse_phylip_matrix_mean(matrix_out)
    
    print(f"  -> Iteration {iteration}: Mean Distance = {best_dist:.6f}")
    
    while True:
        iteration += 1
        trimmed_fastas = []
        valid_trim = False
        
        # 클러스터 내 모든 서열 양 끝단 깎기
        for fasta_file in current_fastas:
            record = SeqIO.read(fasta_file, "fasta")
            if len(record.seq) > min_len + (step_size * 2):
                trimmed_seq = record.seq[step_size:-step_size]
                new_fasta = fasta_file.replace(f"iter{iteration-1}", f"iter{iteration}")
                with open(new_fasta, "w") as f:
                    f.write(f">{record.id}\n{trimmed_seq}\n")
                trimmed_fastas.append(new_fasta)
                valid_trim = True
                
        if not valid_trim:
            print("  -> Reached minimum length limit. Stopping.")
            break
            
        # 깎아낸 서열로 Parabola 재계산
        run_cmd(["./parabola", "sketch", "-k", str(k), "-s", str(scale), "-p", str(threads)] + trimmed_fastas)
        matrix_out = run_cmd(["./parabola", "triangle"] + [f + ".parabola" for f in trimmed_fastas], capture_out=True)
        new_dist = parse_phylip_matrix_mean(matrix_out)
        
        print(f"  -> Iteration {iteration}: Mean Distance = {new_dist:.6f}")
        
        # 평균 거리가 감소(상동성 증가)하면 계속 깎고, 증가하면 중단
        if new_dist < best_dist:
            best_dist = new_dist
            current_fastas = trimmed_fastas
        else:
            print("  -> Distance increased (Hit the TE Core). Stopping carving.")
            print(f"[+] Optimal boundaries found at Iteration {iteration - 1}.")
            break

def main():
    parser = argparse.ArgumentParser(description="TE Founder Detection with Parabola & Iterative Carving")
    parser.add_argument("-i", "--input_dir", required=True, help="Target FASTA directory")
    parser.add_argument("-w", "--window", type=int, default=10000, help="Non-overlapping window size (bp)")
    parser.add_argument("-k", "--kmer", type=int, default=21, help="K-mer size for Parabola")
    parser.add_argument("-c", "--scale", type=int, default=1000, help="FracMinHash scale for Parabola")
    parser.add_argument("-p", "--threads", type=int, default=8, help="Number of threads")
    parser.add_argument("-t", "--trim", type=int, default=100, help="Trimming step size for carving (bp)")
    parser.add_argument("-f", "--filter_pct", type=float, default=1.0, help="Genome coverage filter threshold in percent (default: 1.0 for 1 percent)")
    args = parser.parse_args()

    work_dir = "te_workspace"
    os.makedirs(work_dir, exist_ok=True)
    
    # 1. 오버랩 없는 윈도우 청크 생성
    chunk_files, id_map = generate_non_overlapping_windows(args.input_dir, work_dir, args.window)
    
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
        
    # 5. 모든 발견된 클러스터에 대해 점진적 경계 깎기(Carving) 수행
    for clade in candidates:
        iterative_boundary_carving(
            clade, work_dir, os.path.join(work_dir, "chunks"),
            args.kmer, args.scale, args.threads, step_size=args.trim
        )

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
            genome_beds[genome_id].append((chrom, start, end))
            
    for genome_id, regions in genome_beds.items():
        bed_file = os.path.join(results_dir, f"{genome_id}_fTE.bed")
        print(f"[*] Saving BED file: {bed_file}")
        with open(bed_file, "w") as f:
            for chrom, start, end in regions:
                f.write(f"{chrom}\t{start}\t{end}\n")

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
            genome_windows[genome_id].append((chrom, start, end))
            
    for genome_id, regions in genome_windows.items():
        bed_file = os.path.join(results_dir, f"{genome_id}_window.bed")
        print(f"[*] Saving BED file: {bed_file}")
        with open(bed_file, "w") as f:
            for chrom, start, end in sorted(regions, key=lambda x: (x[0], int(x[1]))):
                f.write(f"{chrom}\t{start}\t{end}\n")

    # 7. te_workspace 제거, tree.nwk 복사 및 nwk 내의 te_workspace/chunks/ 경로 제거
    final_tree_file = os.path.join(results_dir, "tree.nwk")
    if os.path.exists(tree_file):
        with open(tree_file, "r") as f:
            tree_content = f.read()
        tree_content = tree_content.replace("te_workspace/chunks/", "")
        with open(tree_file, "w") as f:
            f.write(tree_content)
        shutil.move(tree_file, final_tree_file)
        print(f"[*] Tree file saved to: {final_tree_file}")
        
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
        print(f"[*] Cleaned up workspace directory: {work_dir}")

if __name__ == "__main__":
    main()