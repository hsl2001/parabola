import sys
from collections import defaultdict

degree = defaultdict(int)
edges = []

with open("reverb.dup.bedpe") as f:
    for line in f:
        if line.startswith("#"): continue
        cols = line.strip().split("\t")
        if len(cols) < 9: continue
        w1 = f"{cols[0]}:{cols[1]}-{cols[2]}"
        w2 = f"{cols[3]}:{cols[4]}-{cols[5]}"
        degree[w1] += 1
        degree[w2] += 1
        edges.append((w1, w2, float(cols[7])))

print(f"Total windows with edges: {len(degree)}")
degrees = list(degree.values())
degrees.sort(reverse=True)
print(f"Top 20 degrees: {degrees[:20]}")
print(f"Number of windows with degree > 50: {sum(1 for d in degrees if d > 50)}")
print(f"Number of windows with degree > 10: {sum(1 for d in degrees if d > 10)}")

# Check distance distribution for high-degree nodes vs low-degree nodes
