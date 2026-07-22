#!/usr/bin/env python3
import sys
import os

def merge_intervals(intervals):
    """
    Given a list of (start, end) tuples, merges overlapping/adjacent intervals.
    Returns merged intervals and total unique base pairs.
    """
    if not intervals:
        return [], 0
    intervals.sort(key=lambda x: x[0])
    merged = []
    curr_start, curr_end = intervals[0]
    for s, e in intervals[1:]:
        if s <= curr_end:
            curr_end = max(curr_end, e)
        else:
            merged.append((curr_start, curr_end))
            curr_start, curr_end = s, e
    merged.append((curr_start, curr_end))
    
    total_bp = sum(e - s for s, e in merged)
    return merged, total_bp

def load_bed(bed_path):
    """
    Parses a BED file (Reverb or BISER) into a dictionary of chromosome -> list of (start, end) intervals.
    Handles both standard BED (chrom, start, end) and BEDPE (chrom1, start1, end1, chrom2, start2, end2).
    """
    chrom_intervals = {}
    if not os.path.exists(bed_path):
        print(f"[Warning] File not found: {bed_path}")
        return chrom_intervals

    with open(bed_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                parts = line.split()
            if len(parts) < 3:
                continue

            # First interval
            try:
                chrom1 = parts[0]
                s1 = int(parts[1])
                e1 = int(parts[2])
                if s1 > e1: s1, e1 = e1, s1
                chrom_intervals.setdefault(chrom1, []).append((s1, e1))
            except ValueError:
                continue

            # If BEDPE format (6+ columns with second set of coordinates)
            if len(parts) >= 6:
                try:
                    chrom2 = parts[3]
                    s2 = int(parts[4])
                    e2 = int(parts[5])
                    if s2 > e2: s2, e2 = e2, s2
                    chrom_intervals.setdefault(chrom2, []).append((s2, e2))
                except ValueError:
                    pass

    return chrom_intervals

def compute_overlap(reverb_bed, biser_bed):
    reverb_raw = load_bed(reverb_bed)
    biser_raw = load_bed(biser_bed)

    reverb_merged = {}
    reverb_total_bp = 0
    for chrom, ivs in reverb_raw.items():
        m, bp = merge_intervals(ivs)
        reverb_merged[chrom] = m
        reverb_total_bp += bp

    biser_merged = {}
    biser_total_bp = 0
    for chrom, ivs in biser_raw.items():
        m, bp = merge_intervals(ivs)
        biser_merged[chrom] = m
        biser_total_bp += bp

    # Compute intersection
    intersection_bp = 0
    all_chroms = set(reverb_merged.keys()).union(set(biser_merged.keys()))

    for chrom in all_chroms:
        r_ivs = reverb_merged.get(chrom, [])
        b_ivs = biser_merged.get(chrom, [])
        
        # Calculate overlap between two sorted interval lists
        i, j = 0, 0
        while i < len(r_ivs) and j < len(b_ivs):
            r_start, r_end = r_ivs[i]
            b_start, b_end = b_ivs[j]

            overlap_start = max(r_start, b_start)
            overlap_end = min(r_end, b_end)

            if overlap_start < overlap_end:
                intersection_bp += (overlap_end - overlap_start)

            if r_end < b_end:
                i += 1
            else:
                j += 1

    union_bp = reverb_total_bp + biser_total_bp - intersection_bp
    jaccard = (intersection_bp / union_bp * 100.0) if union_bp > 0 else 0.0
    r_cov = (intersection_bp / reverb_total_bp * 100.0) if reverb_total_bp > 0 else 0.0
    b_cov = (intersection_bp / biser_total_bp * 100.0) if biser_total_bp > 0 else 0.0

    print("=== Genomic Overlap Analysis ===")
    print(f"Reverb BED: {reverb_bed}")
    print(f"  Total unique SD region coverage: {reverb_total_bp:,} bp")
    print(f"BISER BED:  {biser_bed}")
    print(f"  Total unique SD region coverage: {biser_total_bp:,} bp")
    print(f"Intersection (Overlap):           {intersection_bp:,} bp")
    print(f"Union:                          {union_bp:,} bp")
    print(f"Overlap % of Reverb:             {r_cov:.2f}%")
    print(f"Overlap % of BISER:              {b_cov:.2f}%")
    print(f"Jaccard Index:                   {jaccard:.2f}%")

    return {
        "reverb_total_bp": reverb_total_bp,
        "biser_total_bp": biser_total_bp,
        "intersection_bp": intersection_bp,
        "reverb_cov_pct": r_cov,
        "biser_cov_pct": b_cov,
        "jaccard_pct": jaccard
    }

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: evaluate_overlap.py <reverb_dup.bed> <biser.bed>")
        sys.exit(1)
    compute_overlap(sys.argv[1], sys.argv[2])
