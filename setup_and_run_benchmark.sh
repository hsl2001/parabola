#!/usr/bin/env bash
set -e

WORKDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$WORKDIR"

echo "=================================================="
echo "1. Environment Activation & Setup"
echo "=================================================="

echo "Python executable: $(which python3 || echo 'Not found')"

if ! command -v biser >/dev/null 2>&1; then
    echo "Installing BISER via pip..."
    pip install biser || pip install git+https://github.com/0xTCG/biser.git
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "Installing samtools via bioconda..."
    micromamba install -y -c bioconda samtools 2>/dev/null || conda install -y -c bioconda samtools 2>/dev/null || true
fi

echo "BISER path: $(which biser || echo 'Not found')"
echo "Samtools path: $(which samtools || echo 'Not found')"

echo "Building Reverb..."
make clean
make

mkdir -p bench_results

YEAST_DIR="../6-yeast-genomes"
ARAB_DIR="../upload_genome/"

SCER_FA="$YEAST_DIR/Saccharomyces_cerevisiae.fasta"
COL0_FA="$ARAB_DIR/Col-0.fasta"

ARAB10_FILES="$ARAB_DIR/Abd-0.fasta $ARAB_DIR/Altai-5.fasta $ARAB_DIR/Anz-0.fasta $ARAB_DIR/Are-1.fasta $ARAB_DIR/Bay-0.fasta $ARAB_DIR/Bla-1.fasta $ARAB_DIR/Bur-0.fasta $ARAB_DIR/Can-0.fasta $ARAB_DIR/Col-0.fasta $ARAB_DIR/Cvi-0.fasta"
YEAST_ALL_FILES="$YEAST_DIR/*.fasta $YEAST_DIR/*.fna"
ARAB69_FILES="$ARAB_DIR/*.fasta"

if command -v samtools >/dev/null 2>&1; then
    echo "Indexing FASTA files with samtools faidx..."
    for f in "$YEAST_DIR"/*.fasta "$YEAST_DIR"/*.fna "$ARAB_DIR"/*.fasta; do
        if [ -f "$f" ] && [ ! -f "${f}.fai" ]; then
            samtools faidx "$f" 2>/dev/null || true
        fi
    done
fi

TIME_CMD="/usr/bin/time -v"
THREADS=16

# Helper function to execute tool with time monitoring
run_benchmark_case() {
    local label="$1"
    local tool="$2"
    local cmd="$3"
    local log_prefix="bench_results/${tool}_${label}"

    echo "--------------------------------------------------"
    echo "Running [$tool] on [$label]..."
    echo "Command: $cmd"
    eval "$TIME_CMD $cmd > ${log_prefix}.stdout.log 2> ${log_prefix}.time.log" || echo "[$tool] on [$label] returned exit code $?."
}

echo "=================================================="
echo "2. Running Benchmarks (Reverb vs BISER)"
echo "=================================================="

# Case 1: Yeast 1 (S. cerevisiae)
run_benchmark_case "yeast1" "reverb" "./reverb -p $THREADS -o bench_results/reverb_yeast1 $SCER_FA"
run_benchmark_case "yeast1" "biser" "biser -t $THREADS -o bench_results/biser_yeast1 $SCER_FA"

# Case 2: Yeast 6
run_benchmark_case "yeast6" "reverb" "./reverb -p $THREADS -o bench_results/reverb_yeast6 $YEAST_ALL_FILES"
run_benchmark_case "yeast6" "biser" "biser -t $THREADS -o bench_results/biser_yeast6 $YEAST_ALL_FILES"

# Case 3: Arab 1 (Col-0)
run_benchmark_case "arab1" "reverb" "./reverb -p $THREADS -o bench_results/reverb_arab1 $COL0_FA"
run_benchmark_case "arab1" "biser" "biser -t $THREADS -o bench_results/biser_arab1 $COL0_FA"

# Case 4: Arab 10
run_benchmark_case "arab10" "reverb" "./reverb -p $THREADS -o bench_results/reverb_arab10 $ARAB10_FILES"
run_benchmark_case "arab10" "biser" "biser -t $THREADS -o bench_results/biser_arab10 $ARAB10_FILES"

# Case 5: Arab 69
run_benchmark_case "arab69" "reverb" "./reverb -p $THREADS -o bench_results/reverb_arab69 $ARAB69_FILES"
run_benchmark_case "arab69" "biser" "biser -t $THREADS -o bench_results/biser_arab69 $ARAB69_FILES"

echo "=================================================="
echo "3. Evaluating Overlap"
echo "=================================================="

for case_id in yeast1 yeast6 arab1 arab10 arab69; do
    r_bed="bench_results/reverb_${case_id}.dup.bed"
    b_bed="$(ls bench_results/biser_${case_id}*.bed 2>/dev/null | head -n 1)"
    if [ -f "$r_bed" ] && [ -n "$b_bed" ] && [ -f "$b_bed" ]; then
        python3 evaluate_overlap.py "$r_bed" "$b_bed" > "bench_results/${case_id}_overlap.txt" || true
    fi
done

echo "=================================================="
echo "4. Generating Summary Benchmark Report"
echo "=================================================="

python3 - << 'EOF'
import os

def parse_time_log(filepath):
    if not os.path.exists(filepath):
        return {"wall": "N/A", "user": "N/A", "sys": "N/A", "rss_mb": "N/A"}
    wall, user, sys_t, rss = "N/A", "N/A", "N/A", "N/A"
    with open(filepath) as f:
        for line in f:
            if "Elapsed (wall clock) time" in line:
                wall = line.split("):")[-1].strip()
            elif "User time (seconds)" in line:
                user = line.split(":")[-1].strip()
            elif "System time (seconds)" in line:
                sys_t = line.split(":")[-1].strip()
            elif "Maximum resident set size" in line:
                kb = int(line.split(":")[-1].strip())
                rss = f"{kb / 1024:.2f}"
    return {"wall": wall, "user": user, "sys": sys_t, "rss_mb": rss}

cases = [
    ("Yeast-1 (S. cerevisiae - 12.3 Mb)", "yeast1"),
    ("Yeast-6 (6 Yeast Genomes - ~72 Mb)", "yeast6"),
    ("Arab-1 (Arabidopsis Col-0 - ~129 Mb)", "arab1"),
    ("Arab-10 (10 Arabidopsis Genomes - ~1.3 Gb)", "arab10"),
    ("Arab-69 (69 Arabidopsis Genomes - ~9.3 Gb)", "arab69"),
]

report_rows = []

for title, cid in cases:
    r = parse_time_log(f"bench_results/reverb_{cid}.time.log")
    b = parse_time_log(f"bench_results/biser_{cid}.time.log")
    
    overlap_file = f"bench_results/{cid}_overlap.txt"
    overlap_summary = "N/A"
    if os.path.exists(overlap_file):
        lines = [line.strip() for line in open(overlap_file) if line.strip()]
        overlap_summary = "\n".join(lines[-7:])

    report_rows.append(f"### {title}\n")
    report_rows.append("| Tool | Wall Clock Time | User CPU Time | System CPU Time | Peak Memory (Max RSS) |\n|---|---|---|---|---|\n")
    report_rows.append(f"| **Reverb** | {r['wall']} | {r['user']} s | {r['sys']} s | {r['rss_mb']} MB |\n")
    report_rows.append(f"| **BISER** | {b['wall']} | {b['user']} s | {b['sys']} s | {b['rss_mb']} MB |\n\n")
    report_rows.append(f"**Overlap Analysis**:\n```\n{overlap_summary}\n```\n\n---\n\n")

full_report = "# BISER vs Reverb Benchmark Results\n\n" + "".join(report_rows)

with open("benchmark_results.md", "w") as f:
    f.write(full_report)

print("Saved benchmark report to benchmark_results.md")
EOF

echo "All benchmarks completed!"
