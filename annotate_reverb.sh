#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <reverb.dup.bed> <reference.gff> <reference_prefix>"
    echo "Example: $0 output.dup.bed TAIR10_TE.gff TAIR10-"
    echo "Note: <reference_prefix> is the genome name prefix added by reverb (e.g., 'GenomeA-')"
    exit 1
fi

BED_FILE=$1
GFF_FILE=$2
REF_PREFIX=$3

if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools is not installed or not in PATH."
    exit 1
fi

TMP_DIR=$(mktemp -d)
trap "rm -rf $TMP_DIR" EXIT

# 1. Extract reference lines and strip the prefix to match GFF chromosomes
awk -v pref="$REF_PREFIX" 'BEGIN{OFS="\t"} $1 ~ "^"pref {
    chrom = $1;
    sub("^"pref, "", chrom);
    sub("^-", "", chrom);
    print chrom, $2, $3, $4, $5, $6;
}' "$BED_FILE" > "$TMP_DIR/ref.bed"

# 2. Intersect with GFF
# -wo: write original A and B entries plus the number of base pairs of overlap.
# A has 6 columns. B (GFF) has 9 columns. Overlap length is in the 16th column.
bedtools intersect -a "$TMP_DIR/ref.bed" -b "$GFF_FILE" -wo > "$TMP_DIR/intersected.txt"

# 3. Create a mapping from cluster_id to TE annotation (Max overlap match)
awk -F'\t' 'BEGIN { OFS="\t" }  {
    cluster_id = $4
    attr = $15
    overlap_len = $16
    
    name = "NA"; family = "NA"; class = "NA"
    
    n = split(attr, a, ";")
    for (i=1; i<=n; i++) {
        sub(/^ +/, "", a[i])
        lower_a = tolower(a[i])
        if (lower_a ~ /^name=/) { name = substr(a[i], index(a[i], "=")+1) }
        else if (lower_a ~ /^family_name=/ || lower_a ~ /^locus_biotype=/) { family = substr(a[i], index(a[i], "=")+1) }
        else if (lower_a ~ /^classification=/ || lower_a ~ /^class=/) { class = substr(a[i], index(a[i], "=")+1) }
    }
    
    # Overlap 길이가 가장 큰(best match) 주석으로 갱신
    if (!(cluster_id in max_overlap) || overlap_len > max_overlap[cluster_id]) {
        max_overlap[cluster_id] = overlap_len
        map[cluster_id] = name "\t" family "\t" class
    }
}
END {
    for (k in map) {
        print k, map[k]
    }
}' "$TMP_DIR/intersected.txt" > "$TMP_DIR/cluster_map.txt"

# 4. Annotate the original BED file
echo -e "#chrom\tstart\tend\tcluster_id\tsubcluster_id\tcopy_count\tName\tfamily_name\tclassification"

awk -F'\t' 'BEGIN { OFS="\t" }
NR==FNR {
    map[$1] = $2 "\t" $3 "\t" $4
    next
}
/^#/ { next }
{
    cluster_id = $4
    annot = map[cluster_id]
    if (annot == "") {
        annot = "NA\tNA\tNA"
    }
    print $0, annot
}' "$TMP_DIR/cluster_map.txt" "$BED_FILE"
