#!/usr/bin/env bash
# ==============================================================================
# Segmental Duplication (SD) Pipeline Script using Reverb
# Pipelines for:
#   1) t2t-chm12 (Human T2T CHM13 v2.0)
#   2) t2t-nip   (Rice T2T Nipponbare AGIS1.0)
#   3) col-cen   (Arabidopsis thaliana Col-CEN v1.2 T2T)
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# Defaults
DATA_DIR="${SCRIPT_DIR}/data"
RESULTS_DIR="${SCRIPT_DIR}/results"
THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 8)
TARGET="all"
DOWNLOAD_ONLY=0
SKIP_DOWNLOAD=0
DECOMPRESS=0

# Datasets definition
declare -A URLS_PRIMARY
declare -A URLS_SECONDARY
declare -A FILE_NAMES
declare -A DESCRIPTIONS

# t2t-chm12 (Human T2T CHM13 v2.0 - RefSeq)
URLS_PRIMARY["t2t-chm12"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
URLS_SECONDARY["t2t-chm12"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz"
FILE_NAMES["t2t-chm12"]="t2t_chm13v2.0.fna.gz"
DESCRIPTIONS["t2t-chm12"]="Human T2T (CHM13 v2.0 / RefSeq GCF_009914755.1)"

# t2t-nip (Rice T2T Nipponbare AGIS1.0 - RefSeq)
URLS_PRIMARY["t2t-nip"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_genomic.fna.gz"
URLS_SECONDARY["t2t-nip"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/034/140/825/GCA_034140825.1_AGIS1.0/GCA_034140825.1_AGIS1.0_genomic.fna.gz"
FILE_NAMES["t2t-nip"]="t2t_nip_agis1.0.fna.gz"
DESCRIPTIONS["t2t-nip"]="Rice T2T (Nipponbare AGIS1.0 / RefSeq GCF_034140825.1)"

# col-cen (Arabidopsis thaliana Col-CEN v1.2 T2T)
URLS_PRIMARY["col-cen"]="https://github.com/schatzlab/Col-CEN/raw/main/v1.2/Col-CEN_v1.2.fasta.gz"
URLS_SECONDARY["col-cen"]="https://raw.githubusercontent.com/schatzlab/Col-CEN/main/v1.2/Col-CEN_v1.2.fasta.gz"
FILE_NAMES["col-cen"]="col_cen_v1.2.fasta.gz"
DESCRIPTIONS["col-cen"]="Arabidopsis thaliana T2T (Col-CEN v1.2)"

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Run Segmental Duplication (SD) detection pipeline using Reverb for target genomes.

Options:
  -s, --target TARGET    Genome target to process: t2t-chm12, t2t-nip, col-cen, or all [default: all]
  -t, --threads NUM      Number of CPU threads to use [default: $THREADS]
  -d, --data-dir DIR     Directory to store downloaded genome fasta files [default: $DATA_DIR]
  -r, --results-dir DIR  Directory to store pipeline output results [default: $RESULTS_DIR]
  --download-only        Only download the genome datasets, do not run Reverb
  --skip-download        Skip downloading (assumes genome fasta files exist)
  --decompress           Decompress downloaded .gz files into raw .fa files
  -h, --help             Display this help message and exit

Available targets:
  t2t-chm12   ${DESCRIPTIONS["t2t-chm12"]}
  t2t-nip     ${DESCRIPTIONS["t2t-nip"]}
  col-cen     ${DESCRIPTIONS["col-cen"]}
  all         Run all 3 targets sequentially (t2t-chm12, t2t-nip, col-cen)
EOF
    exit 0
}

# Parse CLI arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--target)
            TARGET="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -d|--data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        -r|--results-dir)
            RESULTS_DIR="$2"
            shift 2
            ;;
        --download-only)
            DOWNLOAD_ONLY=1
            shift
            ;;
        --skip-download)
            SKIP_DOWNLOAD=1
            shift
            ;;
        --decompress)
            DECOMPRESS=1
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "[ERROR] Unknown option: $1"
            usage
            ;;
    esac
done

# Validate target
if [[ "$TARGET" != "all" && "$TARGET" != "t2t-chm12" && "$TARGET" != "t2t-nip" && "$TARGET" != "col-cen" ]]; then
    echo "[ERROR] Invalid target '$TARGET'. Choose from: t2t-chm12, t2t-nip, col-cen, all"
    exit 1
fi

mkdir -p "${DATA_DIR}" "${RESULTS_DIR}"

# Helper function to download file
download_genome() {
    local key="$1"
    local filename="${FILE_NAMES[$key]}"
    local dest_gz="${DATA_DIR}/${filename}"
    local dest_fa="${dest_gz%.gz}"

    if [[ -f "$dest_fa" && -s "$dest_fa" ]]; then
        echo "[INFO] Uncompressed genome file exists: $dest_fa (Skipping download)"
        return 0
    fi

    if [[ -f "$dest_gz" && -s "$dest_gz" ]]; then
        echo "[INFO] Genome compressed archive exists: $dest_gz (Skipping download)"
        if [[ $DECOMPRESS -eq 1 ]]; then
            echo "[INFO] Decompressing $dest_gz -> $dest_fa ..."
            gzip -dc "$dest_gz" > "$dest_fa"
        fi
        return 0
    fi

    echo "[INFO] Downloading ${DESCRIPTIONS[$key]}..."
    local primary_url="${URLS_PRIMARY[$key]}"
    local secondary_url="${URLS_SECONDARY[$key]}"

    local download_success=0

    if command -v wget >/dev/null 2>&1; then
        if wget --quiet --show-progress -O "$dest_gz" "$primary_url"; then
            download_success=1
        elif [[ -n "$secondary_url" ]]; then
            echo "[WARN] Primary URL failed. Trying secondary URL: $secondary_url"
            if wget --quiet --show-progress -O "$dest_gz" "$secondary_url"; then
                download_success=1
            fi
        fi
    elif command -v curl >/dev/null 2>&1; then
        if curl -L --progress-bar -o "$dest_gz" "$primary_url"; then
            download_success=1
        elif [[ -n "$secondary_url" ]]; then
            echo "[WARN] Primary URL failed. Trying secondary URL: $secondary_url"
            if curl -L --progress-bar -o "$dest_gz" "$secondary_url"; then
                download_success=1
            fi
        fi
    else
        echo "[ERROR] Neither wget nor curl found. Please install wget or curl."
        exit 1
    fi

    if [[ $download_success -ne 1 || ! -s "$dest_gz" ]]; then
        echo "[ERROR] Failed to download dataset for $key"
        rm -f "$dest_gz"
        exit 1
    fi

    echo "[INFO] Download completed: $dest_gz"

    if [[ $DECOMPRESS -eq 1 ]]; then
        echo "[INFO] Decompressing $dest_gz -> $dest_fa ..."
        gzip -dc "$dest_gz" > "$dest_fa"
    fi
}

# Helper function to get FASTA path for key
get_fasta_path() {
    local key="$1"
    local filename="${FILE_NAMES[$key]}"
    local dest_gz="${DATA_DIR}/${filename}"
    local dest_fa="${dest_gz%.gz}"

    if [[ -f "$dest_fa" && -s "$dest_fa" ]]; then
        echo "$dest_fa"
    elif [[ -f "$dest_gz" && -s "$dest_gz" ]]; then
        echo "$dest_gz"
    else
        echo ""
    fi
}

# Determine target list
if [[ "$TARGET" == "all" ]]; then
    TARGET_LIST=("t2t-chm12" "t2t-nip" "col-cen")
else
    TARGET_LIST=("$TARGET")
fi

echo "=========================================================="
echo " Segmental Duplication Detection Pipeline (Reverb)"
echo " Target(s): ${TARGET_LIST[*]}"
echo " Threads:   ${THREADS}"
echo " Data Dir:  ${DATA_DIR}"
echo " Out Dir:   ${RESULTS_DIR}"
echo "=========================================================="

# 1. Download Phase
if [[ $SKIP_DOWNLOAD -eq 0 ]]; then
    echo ""
    echo "----------------------------------------------------------"
    echo "[STEP 1/3] Downloading Genomes"
    echo "----------------------------------------------------------"
    for t in "${TARGET_LIST[@]}"; do
        download_genome "$t"
    done
else
    echo "[INFO] Skipping download step as requested (--skip-download)."
fi

if [[ $DOWNLOAD_ONLY -eq 1 ]]; then
    echo "[INFO] Download complete (--download-only mode). Exiting pipeline."
    exit 0
fi

# 2. Build Reverb Phase
echo ""
echo "----------------------------------------------------------"
echo "[STEP 2/3] Building Reverb Binary"
echo "----------------------------------------------------------"
if [[ ! -f "./reverb" ]]; then
    echo "[INFO] Compiling Reverb..."
    make clean && make
else
    echo "[INFO] Reverb binary already exists (./reverb)."
fi

# 3. Execution Phase
echo ""
echo "----------------------------------------------------------"
echo "[STEP 3/3] Running Segmental Duplication Detection"
echo "----------------------------------------------------------"

for t in "${TARGET_LIST[@]}"; do
    fasta_input="$(get_fasta_path "$t")"
    if [[ -z "$fasta_input" || ! -f "$fasta_input" ]]; then
        echo "[ERROR] Input FASTA file for $t not found in ${DATA_DIR}. Run download first."
        exit 1
    fi

    out_prefix="${RESULTS_DIR}/${t}_sd"
    bed_out="${out_prefix}.dup.bed"

    echo ""
    echo ">>> Running Reverb for [$t] (${DESCRIPTIONS[$t]})"
    echo "    Input:   $fasta_input"
    echo "    Output:  $bed_out"
    echo "    Threads: $THREADS"
    
    start_time=$SECONDS
    ./reverb -p "$THREADS" -o "$out_prefix" "$fasta_input"
    elapsed=$(( SECONDS - start_time ))

    if [[ -f "$bed_out" && -s "$bed_out" ]]; then
        lines=$(wc -l < "$bed_out")
        size=$(du -h "$bed_out" | cut -f1)
        echo "[SUCCESS] Finished [$t] in ${elapsed}s. Results saved to $bed_out ($lines SD regions, $size)."
    else
        echo "[WARNING] Reverb finished, but output file $bed_out was not created or is empty."
    fi
done

echo ""
echo "=========================================================="
echo " Pipeline Complete!"
echo " Results directory: ${RESULTS_DIR}"
echo "=========================================================="
