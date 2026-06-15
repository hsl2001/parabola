import os
import sys
import glob
import argparse
import subprocess
import shutil
import collections
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO

def run_cmd(cmd_list, capture_out=False):
    res = subprocess.run(cmd_list, capture_output=capture_out, text=True)
    if res.returncode != 0:
        sys.exit(f"[Error] Command failed: {' '.join(cmd_list)}\n{res.stderr if capture_out else ''}")
    return res.stdout if capture_out else None

def format_loc(real_id):
    if real_id and '|' in real_id:
        parts = real_id.split('|')
        if len(parts) == 3:
            return f"{parts[0]}|{parts[1]}:{parts[2].replace('-', '..')}"
    return real_id

def generate_windows(input_dir, work_dir, window_size=10000, threads=8):
    print(f"[*] Generating {window_size}bp windows...")
    fasta_files = []
    for ext in ("*.fasta", "*.fa", "*.fas", "*.fna"):
        fasta_files.extend(glob.glob(os.path.join(input_dir, ext)))

    chunks_dir = os.path.join(work_dir, "chunks")
    os.makedirs(chunks_dir, exist_ok=True)
    
    id_map, chunk_files = {}, []
    
    def write_chunk(chunk_file, short_id, seq):
        with open(chunk_file, "w") as f:
            f.write(f">{short_id}\n{seq}\n")

    with ThreadPoolExecutor(max_workers=threads) as executor:
        for fasta in fasta_files:
            genome_id = os.path.basename(fasta).split('.')[0]
            for record in SeqIO.parse(fasta, "fasta"):
                seq = str(record.seq).upper()
                for start in range(0, len(seq) - window_size + 1, window_size):
                    short_id = f"W{len(chunk_files)+1:05d}"
                    id_map[short_id] = f"{genome_id}|{record.id}|{start}-{start + window_size}"
                    chunk_file = os.path.join(chunks_dir, f"{short_id}.fasta")
                    chunk_files.append(chunk_file)
                    executor.submit(write_chunk, chunk_file, short_id, seq[start:start+window_size])

    print(f"    -> {len(chunk_files)} windows generated.")
    return chunk_files, id_map

def run_parabola(fasta_files, work_dir, id_map, copy_threshold, k, scale, threads):
    print(f"[*] Running Parabola (k={k}, s={scale}, threads={threads})...")
    for i in range(0, len(fasta_files), 5000):
        run_cmd(["./parabola", "sketch", "-k", str(k), "-s", str(scale), "-p", str(threads)] + fasta_files[i:i+5000])

    triangle_out = os.path.join(work_dir, "triangle.tsv")
    print("[*] Calculating distance matrix...")
    with open(triangle_out, "w") as f_out:
        res = subprocess.run(["./parabola", "triangle"] + [f + ".parabola" for f in fasta_files], stdout=f_out, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0: sys.exit(f"[Error] Parabola triangle failed: {res.stderr}")

    print("[*] Filtering TE candidates...")
    genome_ids = [id_map[os.path.basename(p).split('.')[0]].split('|')[0] for p in fasta_files]
    counts, taxa_names = [0] * len(fasta_files), []
    
    g_to_idx = collections.defaultdict(list)
    for i, gid in enumerate(genome_ids): g_to_idx[gid].append(i)

    with open(triangle_out, 'r') as f:
        num_taxa = int(f.readline().strip())
        for i in range(num_taxa):
            parts = f.readline().strip().split('\t')
            taxa_names.append(parts[0])
            for j in g_to_idx[genome_ids[i]]:
                if j >= i: break
                if float(parts[j+1].split(',')[0]) < 1.0:
                    counts[i] += 1
                    counts[j] += 1

    keep_indices = [i for i, c in enumerate(counts) if c >= copy_threshold] or list(range(num_taxa))
    keep_set = set(keep_indices)
    kept_names = [taxa_names[i] for i in keep_indices]
    dist_dict = collections.defaultdict(dict)
    
    print(f"    -> Kept {len(keep_indices)} TE candidates.")
    with open(triangle_out, 'r') as f:
        f.readline()
        for i in range(num_taxa):
            parts = f.readline().strip().split('\t')
            if i not in keep_set: continue
            for j in keep_indices:
                if j >= i: break
                val = float(parts[j + 1].split(',')[0])
                dist_dict[taxa_names[i]][taxa_names[j]] = val
                dist_dict[taxa_names[j]][taxa_names[i]] = val

    return kept_names, [fasta_files[i] for i in keep_indices], dist_dict

def identify_founders(candidates, dist_dict, copy_threshold):
    print("[*] Identifying founder TEs...")
    prelims = []
    for node in candidates:
        dists = [dist_dict[node][other] for other in candidates if node != other]
        over_08, under_02 = sum(1 for d in dists if d > 0.8), sum(1 for d in dists if d <= 0.2)
        if (over_08 + under_02 == len(dists)) and over_08 > 0 and (copy_threshold < under_02 <= copy_threshold * 5):
            prelims.append({'hub': node, 'total': len(dists), 'under_02': under_02, 'over_08': over_08})

    print(f"    -> {len(prelims)} preliminary founders found. Deduplicating...")
    founders_map = {f['hub']: f for f in prelims}
    visited, final_founders = set(), []
    
    for node in founders_map:
        if node in visited: continue
        queue, component = collections.deque([node]), []
        while queue:
            curr = queue.popleft()
            if curr in visited: continue
            visited.add(curr)
            component.append(curr)
            queue.extend(n for n in founders_map if n not in visited and dist_dict[curr][n] <= 0.2)
                    
        best_founder = max(component, key=lambda n: founders_map[n]['under_02'])
        final_founders.append(founders_map[best_founder])

    return final_founders

def save_outputs(kept_names, kept_files, dist_dict, edge_threshold, founders, id_map, results_dir, input_name):
    print(f"[*] Exporting results to {results_dir}...")
    os.makedirs(results_dir, exist_ok=True)
    founder_set = {f['hub'] for f in founders}
    
    # 1. BED File
    genome_beds = collections.defaultdict(list)
    for f in kept_files:
        short_id = os.path.basename(f).split('.')[0]
        real_id = id_map.get(short_id, "")
        if '|' in real_id:
            g_id, chrom, coords = real_id.split('|')
            start, end = coords.split('-')
            genome_beds[g_id].append((chrom, int(start), int(end), short_id))

    for g_id, regions in genome_beds.items():
        with open(os.path.join(results_dir, f"{g_id}_window.bed"), "w") as f:
            for chrom, start, end, name in sorted(regions):
                f.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    # 2. Founders TXT
    with open(os.path.join(results_dir, f"{input_name}_founders.txt"), "w") as f:
        f.write(f"# Input: {input_name}\n# Founder_ID\tLocation\tTotal_Neighbors\tUnder_0.2_Edges\tOver_0.8_Edges\n")
        for fd in founders:
            short_id = os.path.basename(fd['hub']).split('.')[0]
            f.write(f"{short_id}\t{format_loc(id_map.get(short_id))}\t{fd['total']}\t{fd['under_02']}\t{fd['over_08']}\n")

    # 3. Cytoscape TSVs
    with open(os.path.join(results_dir, f"{input_name}_nodes.tsv"), "w") as fn, \
         open(os.path.join(results_dir, f"{input_name}_edges.tsv"), "w") as fe:
         
        fn.write("Node_ID\tLocation\tRole\tTotal_Neighbors\tUnder_0.2_Edges\tOver_0.8_Edges\n")
        fe.write("Source\tTarget\tDistance\n")
        
        for i, a in enumerate(kept_names):
            short_a = os.path.basename(a).split('.')[0]
            dists = [dist_dict[a][b] for b in kept_names if a != b]
            fn.write(f"{short_a}\t{format_loc(id_map.get(short_a))}\t"
                     f"{'Founder' if a in founder_set else 'Node'}\t"
                     f"{len(dists)}\t{sum(1 for d in dists if d <= 0.2)}\t{sum(1 for d in dists if d > 0.8)}\n")
            
            for j in range(i + 1, len(kept_names)):
                b = kept_names[j]
                if dist_dict[a][b] < edge_threshold:
                    fe.write(f"{short_a}\t{os.path.basename(b).split('.')[0]}\t{dist_dict[a][b]:.4f}\n")

def main():
    parser = argparse.ArgumentParser(description="fTE: Founder TE Detection")
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-w", "--window", type=int, default=10000)
    parser.add_argument("-k", "--kmer", type=int, default=21)
    parser.add_argument("-c", "--scale", type=int, default=100)
    parser.add_argument("-p", "--threads", type=int, default=16)
    parser.add_argument("-y", "--copy", type=int, default=30)
    parser.add_argument("--edge-threshold", type=float, default=1.0)
    args = parser.parse_args()

    input_name = os.path.basename(args.input_dir.rstrip('/'))
    work_dir = f"te_workspace_{input_name}"
    if os.path.exists(work_dir): shutil.rmtree(work_dir)

    chunk_files, id_map = generate_windows(args.input_dir, work_dir, args.window, args.threads)
    kept_names, kept_files, dist_dict = run_parabola(chunk_files, work_dir, id_map, args.copy, args.kmer, args.scale, args.threads)
    founders = identify_founders(kept_names, dist_dict, args.copy)
    save_outputs(kept_names, kept_files, dist_dict, args.edge_threshold, founders, id_map, "fTE_results", input_name)

    shutil.rmtree(work_dir)
    print(f"\n[*] Done. {len(founders)} founder TEs identified.")

if __name__ == "__main__":
    main()