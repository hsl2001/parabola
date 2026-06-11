import os
import sys
import glob
import argparse
import subprocess
import shutil
import collections
import heapq
import re

from Bio import SeqIO
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
    for record in SeqIO.parse(fasta_path, "fasta-blast"):
        seq = str(record.seq).upper()
        seq_len = len(seq)
        for i in range(0, seq_len - window_size + 1, step_size):
            yield record.id, i, seq[i:i+window_size]

def write_chunk(chunk_file, short_id, seq):
    with open(chunk_file, "w") as f:
        f.write(f">{short_id}\n{seq}\n")

def generate_windows(input_dir, work_dir, window_size=10000, threads=8):
    """윈도우 크기의 절반(step size) 간격으로 유전체를 분할하여 개별 FASTA 파일로 저장 (고성능/저메모리 버전)"""
    print(f"[*] Generating {window_size}bp windows with {window_size//2}bp step...")
    fasta_files = []
    for ext in ("*.fasta", "*.fa", "*.fas", "*.fna"):
        fasta_files.extend(glob.glob(os.path.join(input_dir, ext)))
    
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

def run_parabola(fasta_files, work_dir, id_map, copy_threshold, clip=0.0, k=21, scale=1000, threads=8):
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
    
    genome_to_indices = collections.defaultdict(list)
    for i, gid in enumerate(genome_ids):
        genome_to_indices[gid].append(i)
        
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
            
            for real_j in genome_to_indices[g_i]:
                if real_j >= i:
                    break
                val = float(parts[real_j + 1].split(',')[0])
                if clip < val < 1.0 - clip:
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
    
    kept_list_sorted = sorted(keep_indices)
    
    dists_kept = [[0.0] * num_kept for _ in range(num_kept)]
    
    with open(triangle_out, 'r') as f:
        f.readline() # skip num_taxa
        for i in range(num_taxa):
            line = f.readline()
            if i not in keep_set:
                continue
                
            parts = line.strip().split('\t')
            for real_j in kept_list_sorted:
                if real_j >= i:
                    break
                val = float(parts[real_j + 1].split(',')[0])
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
    parser.add_argument("--clip", type=float, default=0.2, help="Clip threshold for distance filtering (drops dist <= clip and dist >= 1.0 - clip) (default: 0.0)")
    args = parser.parse_args()

    input_name = os.path.basename(args.input_dir.rstrip('/'))
    work_dir = f"te_workspace_{input_name}"
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir, exist_ok=True)
    
    # 1. 윈도우 청크 생성
    chunk_files, id_map = generate_windows(args.input_dir, work_dir, args.window, args.threads)
    
    # 2. Parabola 실행 (Sketch & Triangle)
    matrix_file, kept_files, dist_dict = run_parabola(chunk_files, work_dir, id_map, args.copy, args.clip, args.kmer, args.scale, args.threads)
    
    # 3. FastME 계통수 구축
    print("[*] Building phylogenetic tree with FastME...")
    tree_file = os.path.join(work_dir, "tree.nwk")
    run_cmd(["fastme", "-i", matrix_file, "-o", tree_file, "-T", "16"])

    # 4. 최종 결과를 fTE_results/ 디렉토리에 저장
    results_dir = "fTE_results"
    os.makedirs(results_dir, exist_ok=True)
    
    # 살아남은 window의 위치를 <genome>_window.bed 파일로 작성
    save_regions_to_bed(kept_files, id_map, results_dir, "window")

    # 음수 거리 치환 및 트리 파일 정리 
    print("[*] Formatting tree file...")
    final_tree_file = os.path.join(results_dir, f"{input_name}_tree.nwk")
    
    with open(tree_file, "r") as f:
        tree_content = f.read()
        
    # 1. 음수 거리(Negative branch lengths)를 0.0으로 치환 (예: :-0.00123 -> :0.0)
    tree_content = re.sub(r':-[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?', ':0.0', tree_content)
    
    # 2. 경로 및 확장자 제거 (가독성 향상)
    tree_content = tree_content.replace(f"{work_dir}/chunks/", "")
    tree_content = tree_content.replace(".fasta", "")
    
    with open(final_tree_file, "w") as f:
        f.write(tree_content)
        
    print(f"[*] Tree file saved to: {final_tree_file}")
        
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
        print(f"[*] Cleaned up workspace directory: {work_dir}")

if __name__ == "__main__":
    main()