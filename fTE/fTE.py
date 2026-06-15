import os
import sys
import glob
import argparse
import subprocess
import shutil
import collections
import json

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
    step_size = window_size
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        seq_len = len(seq)
        for i in range(0, seq_len - window_size + 1, step_size):
            yield record.id, i, seq[i:i+window_size]


def write_chunk(chunk_file, short_id, seq):
    with open(chunk_file, "w") as f:
        f.write(f">{short_id}\n{seq}\n")


def generate_windows(input_dir, work_dir, window_size=10000, threads=8):
    """유전체를 윈도우 크기 간격(오버랩 없음)으로 분할하여 개별 FASTA 파일로 저장"""
    print(f"[*] Generating {window_size}bp windows (non-overlapping)...")
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

        for future in futures:
            future.result()

    print(f"    -> Total {window_counter} valid windows generated.")
    return chunk_files, id_map


def run_parabola(fasta_files, work_dir, id_map, copy_threshold, k=21, scale=1000, threads=8):
    """Parabola sketch 및 triangle을 이용한 거리 행렬 계산 및 TE 후보군 선별"""
    print(f"[*] Running Parabola (k={k}, s={scale}, threads={threads})...")

    # 1. Parabola Sketch (배칭 처리)
    batch_size = 5000
    for i in range(0, len(fasta_files), batch_size):
        batch = fasta_files[i:i+batch_size]
        cmd_sketch = ["./parabola", "sketch", "-k", str(k), "-s", str(scale),
                       "-p", str(threads)] + batch
        run_cmd(cmd_sketch)

    # 2. Parabola Triangle
    sketch_files = [f + ".parabola" for f in fasta_files]
    triangle_out = os.path.join(work_dir, "triangle.tsv")

    print("[*] Calculating distance matrix...")
    with open(triangle_out, "w") as f_out:
        cmd_triangle = ["./parabola", "triangle"] + sketch_files
        res = subprocess.run(cmd_triangle, stdout=f_out, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            print(f"[Error] Command failed: {' '.join(cmd_triangle[:5])} ...\n{res.stderr}")
            sys.exit(1)

    # 3. 1차 패스: 동일 유전체 내에서 distance < 1.0인 파트너 수 계산 → TE 후보 선별
    print("[*] Pass 1: Selecting TE candidates...")
    genome_ids = [id_map[os.path.basename(p).split('.')[0]].split('|')[0] for p in fasta_files]
    counts = [0] * len(fasta_files)
    taxa_names = []

    genome_to_indices = collections.defaultdict(list)
    for i, gid in enumerate(genome_ids):
        genome_to_indices[gid].append(i)

    with open(triangle_out, 'r') as f:
        first_line = f.readline()
        if not first_line:
            print("[Error] Empty triangle output.")
            sys.exit(1)
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
                if val < 1.0:
                    counts[i] += 1
                    counts[real_j] += 1

    keep_indices = [i for i, c in enumerate(counts) if c >= copy_threshold]

    if len(keep_indices) < 3:
        print("[*] Warning: Too few chunks kept after filtering. Keeping all chunks.")
        keep_indices = list(range(num_taxa))
    else:
        print(f"[*] Filtered out {num_taxa - len(keep_indices)} unique chunks. "
              f"Keeping {len(keep_indices)} TE candidates.")

    # 4. 2차 패스: 후보 chunk들의 거리 행렬 추출
    print("[*] Pass 2: Extracting distance matrix for TE candidates...")
    keep_set = set(keep_indices)
    keep_idx_map = {old_i: new_i for new_i, old_i in enumerate(keep_indices)}
    num_kept = len(keep_indices)
    kept_list_sorted = sorted(keep_indices)

    dists_kept = [[0.0] * num_kept for _ in range(num_kept)]

    with open(triangle_out, 'r') as f:
        f.readline()  # skip num_taxa
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

    # 5. dist_dict 생성 (taxa_name 기반 거리 사전)
    dist_dict = collections.defaultdict(dict)
    kept_names = []
    for new_i, old_i in enumerate(keep_indices):
        name = taxa_names[old_i]
        kept_names.append(name)
        for new_j, old_j in enumerate(keep_indices):
            dist_dict[name][taxa_names[old_j]] = dists_kept[new_i][new_j]

    return kept_names, [fasta_files[i] for i in keep_indices], dist_dict


def build_ssn(candidates, dist_dict, edge_threshold=1.0):
    """SSN 구축: dist < edge_threshold인 쌍을 엣지로 연결한 인접 리스트 생성"""
    print(f"[*] Building Sequence Similarity Network (edge threshold: dist < {edge_threshold})...")
    adj = collections.defaultdict(set)
    edge_count = 0
    for i, a in enumerate(candidates):
        for j in range(i + 1, len(candidates)):
            b = candidates[j]
            if dist_dict[a][b] < edge_threshold:
                adj[a].add(b)
                adj[b].add(a)
                edge_count += 1
    print(f"    -> SSN: {len(candidates)} nodes, {edge_count} edges.")
    return adj


def find_clusters(candidates, dist_dict, max_intra_dist=0.2):
    """Community Detection: 모든 노드 쌍 거리 < max_intra_dist를 보장하는 클러스터 탐색

    1단계: dist < max_intra_dist 그래프에서 Connected Component 탐색
    2단계: 각 Component 내 모든 쌍 검증, 위반 시 Greedy Clique Cover로 분할
    """
    print(f"[*] Community detection (max intra-cluster dist: {max_intra_dist})...")

    # 1. dist < max_intra_dist 그래프 구축
    adj = collections.defaultdict(set)
    for i, a in enumerate(candidates):
        for j in range(i + 1, len(candidates)):
            b = candidates[j]
            if dist_dict[a][b] < max_intra_dist:
                adj[a].add(b)
                adj[b].add(a)

    # 2. BFS로 Connected Components 탐색
    visited = set()
    components = []
    for node in candidates:
        if node not in visited:
            component = []
            queue = collections.deque([node])
            while queue:
                n = queue.popleft()
                if n in visited:
                    continue
                visited.add(n)
                component.append(n)
                for neighbor in adj[n]:
                    if neighbor not in visited:
                        queue.append(neighbor)
            components.append(component)

    # 3. 각 Component 검증 → 위반 시 Greedy Clique Cover로 분할
    final_clusters = []
    for comp in components:
        if len(comp) <= 1:
            final_clusters.append(comp)
            continue

        # 모든 쌍의 거리 < max_intra_dist 검증
        is_valid = True
        for i_idx, a in enumerate(comp):
            for j_idx in range(i_idx + 1, len(comp)):
                if dist_dict[a][comp[j_idx]] >= max_intra_dist:
                    is_valid = False
                    break
            if not is_valid:
                break

        if is_valid:
            final_clusters.append(comp)
        else:
            # Greedy Clique Cover: 제약 조건을 만족하는 하위 클러스터로 분할
            remaining = set(comp)
            while remaining:
                remaining_list = list(remaining)
                # 평균 거리가 최소인 노드를 seed로 선택
                seed = min(remaining_list, key=lambda n: sum(
                    dist_dict[n].get(m, 1.0) for m in remaining_list if m != n
                ))
                cluster = [seed]
                remaining.remove(seed)
                # 거리순으로 후보 추가: 기존 멤버 전부와 dist < max_intra_dist인 경우에만
                for candidate in sorted(remaining, key=lambda n: dist_dict[seed].get(n, 1.0)):
                    if all(dist_dict[candidate].get(m, 1.0) < max_intra_dist for m in cluster):
                        cluster.append(candidate)
                for n in cluster[1:]:
                    remaining.discard(n)
                final_clusters.append(cluster)

    multi_clusters = [c for c in final_clusters if len(c) >= 2]
    print(f"    -> {len(final_clusters)} total clusters "
          f"({len(multi_clusters)} with >= 2 members).")
    return final_clusters


def identify_founders(clusters, all_candidates, dist_dict, ext_threshold=0.8):
    """각 클러스터에서 Hub(medoid) 식별 및 Founder TE 판정

    Hub: 클러스터 내 평균 거리가 가장 낮은 노드
    Founder: Hub의 클러스터 외부 최근접 거리 > ext_threshold
    """
    print(f"[*] Identifying founder TEs (external dist threshold: > {ext_threshold})...")
    founders = []
    cluster_id = 0

    for cluster in clusters:
        if len(cluster) < 2:
            continue

        cluster_id += 1
        cluster_set = set(cluster)
        external_nodes = [n for n in all_candidates if n not in cluster_set]

        # Hub 식별: 클러스터 내 평균 거리가 가장 낮은 노드 (medoid)
        best_hub = None
        best_avg = float('inf')
        for node in cluster:
            avg = sum(dist_dict[node][m] for m in cluster if m != node) / (len(cluster) - 1)
            if avg < best_avg:
                best_avg = avg
                best_hub = node

        # 클러스터 외부 최근접 거리 계산
        if external_nodes:
            min_ext = min(dist_dict[best_hub].get(n, 1.0) for n in external_nodes)
        else:
            min_ext = float('inf')

        # Founder 판정: 외부 최근접 거리 > ext_threshold
        if min_ext > ext_threshold:
            founders.append({
                'cluster_id': f'C{cluster_id:03d}',
                'hub': best_hub,
                'cluster_size': len(cluster),
                'avg_intra_dist': best_avg,
                'min_ext_dist': min_ext,
                'cluster_members': cluster,
            })

    print(f"    -> {len(founders)} founder TEs identified "
          f"out of {cluster_id} multi-member clusters.")
    return founders


def save_founder_list(founders, id_map, results_dir, input_name):
    """Founder TE 목록을 .txt 파일로 저장"""
    output_file = os.path.join(results_dir, f"{input_name}_founders.txt")

    with open(output_file, "w") as f:
        f.write("# fTE Founder TE Detection Results\n")
        f.write(f"# Input: {input_name}\n")
        f.write("# Cluster_ID\tFounder_ID\tLocation\tCluster_Size\t"
                "Avg_Intra_Dist\tMin_External_Dist\n")

        for founder in founders:
            hub_name = founder['hub']
            short_id = os.path.basename(hub_name).split('.')[0]
            real_id = id_map.get(short_id, short_id)

            # Chr:start..end 형식으로 변환
            location = real_id
            if '|' in real_id:
                try:
                    parts = real_id.split('|')
                    if len(parts) == 3:
                        _, chrom, coords = parts
                        start, end = coords.split('-')
                        location = f"{chrom}:{start}..{end}"
                except Exception:
                    pass

            f.write(f"{founder['cluster_id']}\t{short_id}\t{location}\t"
                    f"{founder['cluster_size']}\t{founder['avg_intra_dist']:.4f}\t"
                    f"{founder['min_ext_dist']:.4f}\n")

    print(f"[*] Founder TE list saved to: {output_file}")
    return output_file


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


def save_network_files(kept_names, dist_dict, edge_threshold, clusters, founders, id_map, results_dir, input_name):
    """Cytoscape용 네트워크 파일들 (Node list TSV, Edge list TSV) 저장"""
    print("[*] Exporting Cytoscape network files (nodes and edges TSV)...")
    
    # 1. 노드 정보 및 역할 분류
    founder_map = {f['hub']: f for f in founders}
    cluster_member_to_id = {}
    cluster_hub_set = set()
    cluster_avg_intra = {}
    
    for c_idx, cluster in enumerate(clusters):
        if len(cluster) < 2:
            continue
        c_id = f"C{c_idx+1:03d}"
        
        # Hub 찾기
        best_hub = None
        best_avg = float('inf')
        for node in cluster:
            cluster_member_to_id[node] = c_id
            avg = sum(dist_dict[node][m] for m in cluster if m != node) / (len(cluster) - 1)
            cluster_avg_intra[node] = avg
            if avg < best_avg:
                best_avg = avg
                best_hub = node
        if best_hub:
            cluster_hub_set.add(best_hub)

    # 2. Node list TSV 저장
    nodes_file = os.path.join(results_dir, f"{input_name}_nodes.tsv")
    with open(nodes_file, "w") as f:
        f.write("Node_ID\tLocation\tRole\tCluster_ID\tAvg_Intra_Dist\tMin_External_Dist\n")
        for name in kept_names:
            short_id = os.path.basename(name).split('.')[0]
            real_id = id_map.get(short_id, short_id)
            
            # Chr:start..end 형식으로 변환
            location = real_id
            if '|' in real_id:
                try:
                    parts = real_id.split('|')
                    if len(parts) == 3:
                        _, chrom, coords = parts
                        start, end = coords.split('-')
                        location = f"{chrom}:{start}..{end}"
                except Exception:
                    pass
            
            # Role 결정
            is_founder = name in founder_map
            is_hub = name in cluster_hub_set
            c_id = cluster_member_to_id.get(name, "Unclustered")
            
            if is_founder:
                role = "Founder"
            elif is_hub:
                role = "Hub"
            elif c_id != "Unclustered":
                role = "Member"
            else:
                role = "Unclustered"
                
            avg_intra = cluster_avg_intra.get(name, 0.0)
            
            if is_founder:
                min_ext = founder_map[name]['min_ext_dist']
            elif is_hub:
                cluster_set = set()
                for cluster in clusters:
                    if name in cluster:
                        cluster_set = set(cluster)
                        break
                external_nodes = [n for n in kept_names if n not in cluster_set]
                min_ext = min(dist_dict[name].get(n, 1.0) for n in external_nodes) if external_nodes else 1.0
            else:
                min_ext = 0.0
                
            f.write(f"{short_id}\t{location}\t{role}\t{c_id}\t{avg_intra:.4f}\t{min_ext:.4f}\n")

    # 3. Edge list TSV 저장
    edges_file = os.path.join(results_dir, f"{input_name}_edges.tsv")
    with open(edges_file, "w") as f:
        f.write("Source\tTarget\tDistance\n")
        for i, a in enumerate(kept_names):
            short_a = os.path.basename(a).split('.')[0]
            for j in range(i + 1, len(kept_names)):
                b = kept_names[j]
                short_b = os.path.basename(b).split('.')[0]
                dist = dist_dict[a][b]
                if dist < edge_threshold:
                    f.write(f"{short_a}\t{short_b}\t{dist:.4f}\n")

    print(f"[*] Node attributes list saved to: {nodes_file}")
    print(f"[*] Edge list saved to: {edges_file}")
    return nodes_file, edges_file


def main():
    parser = argparse.ArgumentParser(
        description="fTE: Founder TE Detection via Sequence Similarity Network"
    )
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Target FASTA directory")
    parser.add_argument("-w", "--window", type=int, default=10000,
                        help="Non-overlapping window size in bp")
    parser.add_argument("-k", "--kmer", type=int, default=21,
                        help="K-mer size for Parabola")
    parser.add_argument("-c", "--scale", type=int, default=100,
                        help="FracMinHash scale for Parabola")
    parser.add_argument("-p", "--threads", type=int, default=16,
                        help="Number of threads ")
    parser.add_argument("-y", "--copy", type=int, default=30,
                        help="Minimum copy number threshold for TE candidate ")
    parser.add_argument("--edge-threshold", type=float, default=1.0,
                        help="SSN edge threshold: connect if dist < value")
    parser.add_argument("--cluster-dist", type=float, default=0.2,
                        help="Max intra-cluster distance")
    parser.add_argument("--founder-dist", type=float, default=0.8,
                        help="Min external distance for founder TE")
    args = parser.parse_args()

    input_name = os.path.basename(args.input_dir.rstrip('/'))
    work_dir = f"te_workspace_{input_name}"
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir, exist_ok=True)

    # 1. 데이터 전처리: 유전체를 윈도우로 분할
    chunk_files, id_map = generate_windows(
        args.input_dir, work_dir, args.window, args.threads
    )

    # 2. 행렬 산출 + 3. 후보군 선별 (Parabola 거리 행렬 → TE 후보 추출)
    kept_names, kept_files, dist_dict = run_parabola(
        chunk_files, work_dir, id_map, args.copy,
        args.kmer, args.scale, args.threads
    )

    # 4. SSN 그래프 구축 (dist < edge_threshold → 엣지)
    ssn_adj = build_ssn(kept_names, dist_dict, args.edge_threshold)

    # SSN에 연결된 노드만 추출 (고립 노드 제외)
    ssn_nodes = [n for n in kept_names if n in ssn_adj]
    print(f"[*] SSN connected nodes: {len(ssn_nodes)} / {len(kept_names)}")

    # 5. Community Detection (intra-cluster dist < cluster_dist 보장)
    clusters = find_clusters(ssn_nodes, dist_dict, args.cluster_dist)

    # 6. Founder TE 특정 (Hub의 외부 최근접 dist > founder_dist)
    founders = identify_founders(clusters, ssn_nodes, dist_dict, args.founder_dist)

    # 7. 결과 저장
    results_dir = "fTE_results"
    os.makedirs(results_dir, exist_ok=True)

    # TE 후보 윈도우 BED 파일
    save_regions_to_bed(kept_files, id_map, results_dir, "window")

    # Founder TE 목록 (.txt)
    save_founder_list(founders, id_map, results_dir, input_name)

    # Cytoscape용 네트워크 파일 저장
    save_network_files(
        kept_names, dist_dict, args.edge_threshold,
        clusters, founders, id_map, results_dir, input_name
    )

    # Cleanup
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
        print(f"[*] Cleaned up workspace directory: {work_dir}")

    print(f"\n[*] Done. {len(founders)} founder TEs identified.")


if __name__ == "__main__":
    main()