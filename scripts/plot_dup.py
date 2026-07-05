#!/usr/bin/env python3
"""
Reverb visualization: plot segmental duplications from BEDPE output.

Usage:
    python plot_dup.py <prefix>

Reads:
    <prefix>.dup.bedpe   — SD pair relationships
    <prefix>.dup.bed     — merged SD regions

Outputs:
    <prefix>.dotplot.png  — self-similarity dot plot
    <prefix>.copynum.png  — copy number track
"""

import argparse
import sys
from collections import defaultdict

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.collections import LineCollection
except ImportError:
    print("Error: matplotlib is required. Install with: pip install matplotlib",
          file=sys.stderr)
    sys.exit(1)


def parse_bedpe(path):
    """Parse BEDPE file, return list of (chrom1, s1, e1, chrom2, s2, e2, family, dist, cc)."""
    records = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            records.append({
                "chrom1": cols[0], "start1": int(cols[1]), "end1": int(cols[2]),
                "chrom2": cols[3], "start2": int(cols[4]), "end2": int(cols[5]),
                "family": cols[6], "distance": float(cols[7]),
                "copy_count": int(cols[8]),
            })
    return records


def parse_bed(path):
    """Parse BED file, return list of (chrom, start, end, family, cc)."""
    records = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 5:
                continue
            records.append({
                "chrom": cols[0], "start": int(cols[1]), "end": int(cols[2]),
                "family": cols[3], "copy_count": int(cols[4]),
            })
    return records


def get_chrom_offsets(records_bedpe, records_bed):
    """Build cumulative chromosome offsets for linear genome coordinate."""
    chrom_max = defaultdict(int)
    for r in records_bedpe:
        chrom_max[r["chrom1"]] = max(chrom_max[r["chrom1"]], r["end1"])
        chrom_max[r["chrom2"]] = max(chrom_max[r["chrom2"]], r["end2"])
    for r in records_bed:
        chrom_max[r["chrom"]] = max(chrom_max[r["chrom"]], r["end"])

    chroms = sorted(chrom_max.keys())
    offsets = {}
    cum = 0
    for ch in chroms:
        offsets[ch] = cum
        cum += chrom_max[ch]
    return offsets, chroms, cum


def plot_dotplot(records, offsets, chroms, total_len, out_path):
    """Self-similarity dot plot."""
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_facecolor("#0a0a1a")
    fig.patch.set_facecolor("#0a0a1a")

    # Draw chromosome boundaries
    for ch in chroms:
        pos = offsets[ch]
        ax.axhline(y=pos, color="#222244", linewidth=0.3)
        ax.axvline(x=pos, color="#222244", linewidth=0.3)

    # Plot SD pairs as points
    xs, ys, cs = [], [], []
    for r in records:
        x = offsets[r["chrom1"]] + (r["start1"] + r["end1"]) / 2
        y = offsets[r["chrom2"]] + (r["start2"] + r["end2"]) / 2
        xs.append(x)
        ys.append(y)
        xs.append(y)
        ys.append(x)
        # Color by distance: closer = brighter
        intensity = max(0.2, 1.0 - r["distance"] * 5)
        cs.append(intensity)
        cs.append(intensity)

    scatter = ax.scatter(xs, ys, c=cs, cmap="plasma", s=1.5, alpha=0.8,
                         vmin=0, vmax=1, edgecolors="none")

    ax.set_xlim(0, total_len)
    ax.set_ylim(0, total_len)
    ax.set_aspect("equal")
    ax.set_xlabel("Genome Position (bp)", color="white", fontsize=11)
    ax.set_ylabel("Genome Position (bp)", color="white", fontsize=11)
    ax.set_title("Reverb: Segmental Duplication Dot Plot", color="white",
                 fontsize=14, fontweight="bold")
    ax.tick_params(colors="white", labelsize=8)

    # Chromosome labels
    for ch in chroms:
        mid = offsets[ch] + (total_len / len(chroms)) / 2 if len(chroms) > 0 else 0
        ax.text(mid, -total_len * 0.02, ch, color="white", fontsize=7,
                ha="center", va="top", rotation=45)

    plt.colorbar(scatter, ax=ax, label="Similarity", fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Dot plot saved: {out_path}")


def plot_copynum(records_bed, offsets, chroms, total_len, out_path):
    """Copy number track along the genome."""
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.set_facecolor("#0a0a1a")
    fig.patch.set_facecolor("#0a0a1a")

    for r in records_bed:
        x_start = offsets[r["chrom"]] + r["start"]
        x_end = offsets[r["chrom"]] + r["end"]
        cc = r["copy_count"]
        color = plt.cm.viridis(min(cc / 20.0, 1.0))
        ax.fill_between([x_start, x_end], 0, cc, color=color, alpha=0.8,
                        linewidth=0)

    # Chromosome boundaries
    for ch in chroms:
        ax.axvline(x=offsets[ch], color="#444466", linewidth=0.5, linestyle="--")
        ax.text(offsets[ch] + 1000, ax.get_ylim()[1] * 0.95, ch,
                color="white", fontsize=7, va="top")

    ax.set_xlim(0, total_len)
    ax.set_xlabel("Genome Position (bp)", color="white", fontsize=11)
    ax.set_ylabel("Copy Count", color="white", fontsize=11)
    ax.set_title("Reverb: Copy Number Track", color="white",
                 fontsize=14, fontweight="bold")
    ax.tick_params(colors="white", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("white")
    ax.spines["bottom"].set_color("white")

    plt.tight_layout()
    plt.savefig(out_path, dpi=200, facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Copy number track saved: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Reverb: Visualize segmental duplication results")
    parser.add_argument("prefix", help="Output prefix from reverb (e.g., 'reverb')")
    args = parser.parse_args()

    bedpe_path = f"{args.prefix}.dup.bedpe"
    bed_path = f"{args.prefix}.dup.bed"

    try:
        records_bedpe = parse_bedpe(bedpe_path)
    except FileNotFoundError:
        print(f"Error: {bedpe_path} not found", file=sys.stderr)
        sys.exit(1)

    try:
        records_bed = parse_bed(bed_path)
    except FileNotFoundError:
        print(f"Error: {bed_path} not found", file=sys.stderr)
        sys.exit(1)

    if not records_bedpe:
        print("No SD pairs found in BEDPE file.", file=sys.stderr)
        sys.exit(0)

    offsets, chroms, total_len = get_chrom_offsets(records_bedpe, records_bed)

    print(f"[plot_dup] {len(records_bedpe)} SD pairs, {len(records_bed)} regions, "
          f"{len(chroms)} chromosomes")

    plot_dotplot(records_bedpe, offsets, chroms, total_len,
                 f"{args.prefix}.dotplot.png")
    plot_copynum(records_bed, offsets, chroms, total_len,
                 f"{args.prefix}.copynum.png")


if __name__ == "__main__":
    main()
