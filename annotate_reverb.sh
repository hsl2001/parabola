#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <reverb.dup.bed> <reference.gff> <reference_prefix>"
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

awk -v pref="$REF_PREFIX" 'BEGIN{OFS="\t"} {
    if ($1 ~ "^"pref) {
        chrom = $1
        sub("^"pref, "", chrom)
        sub("^-", "", chrom)
        print chrom, $2, $3, $4
    }
}' "$BED_FILE" > "$TMP_DIR/ref.bed"

bedtools intersect -a "$TMP_DIR/ref.bed" -b "$GFF_FILE" -wo > "$TMP_DIR/intersected.txt"

awk -F'\t' 'BEGIN { OFS="\t" }  $7=="mobile_element" {
    cluster_id = $4
    attr = $13
    overlap_len = $14
    
    if (!(cluster_id in max_overlap) || overlap_len > max_overlap[cluster_id]) {
        max_overlap[cluster_id] = overlap_len
        
        name = "NA"; family = "NA"; class = "NA"
        n = split(attr, a, ";")
        for (i=1; i<=n; i++) {
            sub(/^ +/, "", a[i])
            lower_a = tolower(a[i])
            if (lower_a ~ /^name=/) { name = substr(a[i], index(a[i], "=")+1) }
            else if (lower_a ~ /^family_name=/ || lower_a ~ /^locus_biotype=/) { family = substr(a[i], index(a[i], "=")+1) }
            else if (lower_a ~ /^classification=/ || lower_a ~ /^class=/) { class = substr(a[i], index(a[i], "=")+1) }
        }
        
        map[cluster_id] = name "\t" family "\t" class
    }
}
END {
    for (c in map) {
        print c, map[c]
    }
}' "$TMP_DIR/intersected.txt" > "$TMP_DIR/cluster_map.txt"

awk -F'\t' 'BEGIN { OFS="\t"; print "#chrom\tstart\tend\tcluster_id\tsubcluster_id\tName\tfamily_name\tclassification" }
NR==FNR {
    map[$1] = $2 "\t" $3 "\t" $4
    next
}
/^#/ { next }
{
    cluster_id = $4
    annot = (cluster_id in map) ? map[cluster_id] : "NA\tNA\tNA"
    print $0, annot
}' "$TMP_DIR/cluster_map.txt" "$BED_FILE"
