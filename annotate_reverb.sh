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
# -wa: write original A entry, -wb: write original B entry
# A has 6 columns. B (GFF) has 9 columns.
bedtools intersect -a "$TMP_DIR/ref.bed" -b "$GFF_FILE" -wa -wb > "$TMP_DIR/intersected.txt"

# 3. Create a mapping from cluster_id to TE annotation
# We use awk to extract Name, family_name, and classification from the 15th column (GFF attributes)
awk -F'\t' 'BEGIN { OFS="\t" } {
    cluster_id = $4
    attr = $15
    
    name = "NA"; family = "NA"; class = "NA"
    
    n = split(attr, a, ";")
    for (i=1; i<=n; i++) {
        sub(/^ +/, "", a[i])
        lower_a = tolower(a[i])
        if (lower_a ~ /^name=/) { name = substr(a[i], 6) }
        else if (lower_a ~ /^family_name=/) { family = substr(a[i], 13) }
        else if (lower_a ~ /^classification=/) { class = substr(a[i], 16) }
    }
    
    # If a cluster overlaps multiple TEs, we keep the first valid one we see
    if (!(cluster_id in map) || map[cluster_id] ~ /^NA\tNA\tNA/) {
        map[cluster_id] = name "\t" family "\t" class
    }
}
END {
    for (k in map) {
        print k, map[k]
    }
}' "$TMP_DIR/intersected.txt" > "$TMP_DIR/cluster_map.txt"

# 4. Annotate the original BED file
# We output the original columns + Name + family_name + classification
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
