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

def run_parabola(fasta_files, work_dir, k=21, scale=1000, threads=8):
    """Parabola sketch 및 triangle을 이용한 거리 행렬 계산"""
    print(f"[*] Running Parabola (k={k}, s={scale}, threads={threads})...")
    
    # 1. Parabola Sketch
    cmd_sketch = ["./parabola", "sketch", "-k", str(k), "-s", str(scale), "-p", str(threads)] + fasta_files
    run_cmd(cmd_sketch)
    
    # 2. Parabola Triangle
    sketch_files = [f + ".parabola" for f in fasta_files]
    cmd_triangle = ["./parabola", "triangle"] + sketch_files
    matrix_out = run_cmd(cmd_triangle, capture_out=True)
    
    # 결과 행렬 저장
    matrix_file = os.path.join(work_dir, "distance_matrix.phy")
    with open(matrix_file, "w") as f:
        f.write(matrix_out)
        
    return matrix_file

def parse_phylip_matrix_mean(matrix_str):
    """Phylip 행렬 텍스트를 파싱하여 평균 거리 계산 (노이즈 감소 판단용)"""
    lines = matrix_str.strip().split('\n')
    if len(lines) < 2: return 1.0
    
    distances = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) > 1:
            distances.extend([float(x) for x in parts[1:] if float(x) > 0])
            
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
            candidates.append({
                "clade_id": f"Clade_{idx+1}",
                "founder_short_id": best_leaf.name,
                "founder_real_id": id_map.get(best_leaf.name, best_leaf.name),
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
        src = os.path.join(orig_chunks_dir, f"{leaf}.fasta")
        dst = os.path.join(carve_dir, f"{leaf}_iter0.fasta")
        shutil.copy(src, dst)
        current_fastas.append(dst)
        
    iteration = 0
    
    # 초기 거리 측정
    cmd_sketch = ["parabola", "sketch", "-k", str(k), "-s", str(scale), "-p", str(threads)] + current_fastas
    run_cmd(cmd_sketch)
    matrix_out = run_cmd(["parabola", "triangle"] + [f + ".parabola" for f in current_fastas], capture_out=True)
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
        run_cmd(["parabola", "sketch", "-k", str(k), "-s", str(scale), "-p", str(threads)] + trimmed_fastas)
        matrix_out = run_cmd(["parabola", "triangle"] + [f + ".parabola" for f in trimmed_fastas], capture_out=True)
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
    args = parser.parse_args()

    work_dir = "te_workspace"
    os.makedirs(work_dir, exist_ok=True)
    
    # 1. 오버랩 없는 윈도우 청크 생성
    chunk_files, id_map = generate_non_overlapping_windows(args.input_dir, work_dir, args.window)
    
    # 2. Parabola 실행 (Sketch & Triangle)
    matrix_file = run_parabola(chunk_files, work_dir, args.kmer, args.scale, args.threads)
    
    # 3. FastME 계통수 구축
    print("[*] Building phylogenetic tree with FastME...")
    tree_file = os.path.join(work_dir, "tree.nwk")
    run_cmd(["fastme", "-i", matrix_file, "-o", tree_file, "-m", "N", "-s"])
    
    # 4. 트리 파싱 및 후보 클러스터 추출
    candidates = find_founder_candidates(tree_file, id_map)
    
    print("\n[+] Top TE Clades Detected:")
    for c in candidates[:5]: # 상위 5개 클러스터만 출력
        print(f" - {c['clade_id']} (Size: {c['cluster_size']}) -> Founder Candidate: {c['founder_real_id']}")
        
    # 5. 최상위 클러스터에 대해 점진적 경계 깎기(Carving) 수행
    if candidates:
        top_clade = candidates[0]
        iterative_boundary_carving(
            top_clade, work_dir, os.path.join(work_dir, "chunks"),
            args.kmer, args.scale, args.threads, step_size=args.trim
        )

if __name__ == "__main__":
    main()