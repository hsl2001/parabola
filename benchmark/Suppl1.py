import argparse
import csv
import random
import re
import subprocess
import time
import os
import shutil
from pathlib import Path
from dataclasses import dataclass
from typing import Callable

WORK = Path("1benchmark_run")
GENOME_FASTA = Path("../genome.fasta")
PARABOLA = Path("../parabola")
SEED = 42
THREADS = "8"

C_VALS = [30, 50, 100, 500, 1000]
DISTS = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
N_VALS = [1, 2, 4, 8]
REPLICATES = 1
TIME_C = C_VALS[0]

K_VALS = [21, 23, 25, 27, 29, 31]

MUTATION_MAP = {
    65: [67, 71, 84],
    67: [65, 71, 84],
    71: [65, 67, 84],
    84: [65, 67, 71]
}
BASES = [65, 67, 71, 84]

def get_stem(p_str):
    name = Path(p_str).name
    for ext in [".parabola", ".msh", ".fa", ".fasta"]:
        if name.endswith(ext):
            name = name[: -len(ext)]
    return name

def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def capture(cmd):
    return subprocess.check_output(cmd, stderr=subprocess.DEVNULL).decode()

def get_seq_length(path):
    size = 0
    with open(path, "rb") as f:
        f.readline()
        for line in f:
            size += len(line.strip())
    return size

def build_ref(input_path, rpath):
    with open(input_path, "rt") as fin, open(rpath, "wt") as fref:
        fref.write(">ref\n")
        for line in fin:
            if not line.startswith(">"):
                fref.write(line.strip().upper())
        fref.write("\n")

def mutate(input_fasta, output_fasta, name, rate, snp_only=False):
    with open(input_fasta, 'rb') as fin:
        fin.readline()
        seq = bytearray(fin.read().replace(b'\n', b''))

    total_len = len(seq)
    if total_len == 0:
        return 0.0

    target_muts = int(total_len * rate)

    if target_muts > 0:
        mut_positions = set(random.sample(range(total_len), target_muts))
        new_seq = bytearray()

        for i in range(total_len):
            if i in mut_positions:
                if snp_only or random.random() >= 0.10:
                    old_b = seq[i]
                    choices = MUTATION_MAP.get(old_b, BASES)
                    new_seq.append(random.choice(choices))
                else:
                    if random.random() >= 0.5:
                        new_seq.append(random.choice(BASES))
                        new_seq.append(seq[i])
            else:
                new_seq.append(seq[i])

        seq = new_seq

    with open(output_fasta, 'wb') as fout:
        fout.write(f">{name}\n".encode('ascii'))
        fout.write(seq)
        fout.write(b'\n')

    return target_muts / total_len

def parse_triangle(out):
    matrix = {}
    lines = out.strip().splitlines()
    if not lines:
        return matrix
    names = []
    for line in lines[1:]:
        cols = line.split("\t")
        name = get_stem(cols[0])
        for i, dist_str in enumerate(cols[1:]):
            if dist_str.strip():
                try:
                    matrix[(name, names[i])] = matrix[(names[i], name)] = float(dist_str) if float(dist_str) < 0.3 else 'nan'
                except ValueError:
                    pass
        names.append(name)
    return matrix

def parse_fastani_directional(out):
    direction_ani = {}
    for line in out.strip().splitlines():
        cols = re.split(r"\s+", line.strip())
        if len(cols) < 3:
            continue
        q, r = get_stem(cols[0]), get_stem(cols[1])
        try:
            ani = float(cols[2]) / 100.0
        except ValueError:
            continue
        direction_ani[(q, r)] = max(0.0, min(1.0, ani))
    return direction_ani

def finalize_directional_ani(direction_ani, names=None):
    matrix = {}
    name_set = set(names or [])
    for q, r in direction_ani.keys():
        name_set.add(q)
        name_set.add(r)

    for n in name_set:
        matrix[(n, n)] = 0.0

    for (q, r), ani_qr in direction_ani.items():
        if q == r:
            matrix[(q, q)] = 0.0
            continue
        ani_rq = direction_ani.get((r, q))
        ani = (ani_qr + ani_rq) / 2.0 if ani_rq is not None else ani_qr
        dist = max(0.0, 1.0 - ani)
        matrix[(q, r)] = matrix[(r, q)] = dist

    return matrix

def par_batch(fastas, c, k):
    t0 = time.perf_counter()
    run([
            str(PARABOLA), "sketch", "-k", str(k),
            "-s", str(c), "-p", THREADS
        ] + [str(p) for p in fastas])
    out = capture([str(PARABOLA), "triangle"] + [str(p) + ".parabola" for p in fastas])
    return time.perf_counter() - t0, parse_triangle(out)

def mas_batch(fastas, c, k):
    s = max(1, get_seq_length(fastas[0]) // c)
    t0 = time.perf_counter()
    run([
            "mash", "sketch", "-p", THREADS, "-k", str(k), "-s", str(s),
            "-o", str(WORK / "mall")
        ] + [str(p) for p in fastas])
    out = capture(["mash", "triangle", "-p", THREADS, str(WORK / "mall.msh")])
    return time.perf_counter() - t0, parse_triangle(out)

def ska_batch(fastas, c, k):
    lp = WORK / "skani_list.txt"
    lp.write_text("\n".join(str(p) for p in fastas))
    t0 = time.perf_counter()
    out = capture(["skani", "dist", "--ql", str(lp), "--rl", str(lp), "-c", str(c), "-t", THREADS])
    t1 = time.perf_counter()
    matrix = {}
    for line in out.strip().splitlines()[1:]:
        cols = line.split("\t")
        if len(cols) >= 3:
            r, q = get_stem(cols[0]), get_stem(cols[1])
            try:
                matrix[(r, q)] = matrix[(q, r)] = max(0.0, (100.0 - float(cols[2])) / 100.0)
            except ValueError:
                pass
    return (t1 - t0) / 2, matrix

def fastani_batch(fastas, c, k):
    lp = WORK / "fastani_list.txt"
    out_path = WORK / "fastani.out"
    lp.write_text("\n".join(str(p) for p in fastas))
    t0 = time.perf_counter()
    run(["fastANI", "--ql", str(lp), "--rl", str(lp), "-t", THREADS, "-o", str(out_path)])
    out = out_path.read_text() if out_path.exists() else ""
    t1 = time.perf_counter()
    direction_ani = parse_fastani_directional(out)
    if not direction_ani:
        pass
    names = [get_stem(p) for p in fastas]
    return (t1 - t0) / 2, finalize_directional_ani(direction_ani, names)

@dataclass
class Tool:
    name: str
    run_batch: Callable
    clean: Callable
    run_once_per_replicate: bool = False

def gen_clean_fn(exts):
    def clean():
        for ext in exts:
            for p in WORK.glob(ext):
                try: p.unlink()
                except OSError: pass
            for p in Path(".").glob(ext):
                try: p.unlink()
                except OSError: pass
    return clean

TOOLS = [
    Tool("Skani", ska_batch, gen_clean_fn(["*.txt", "*.sketch"])),
    Tool("FastANI", fastani_batch, gen_clean_fn(["fastani_list.txt", "fastani.out"]), run_once_per_replicate=True),
    Tool("Mash", mas_batch, gen_clean_fn(["*.msh"])),
    Tool("Parabola", par_batch, gen_clean_fn(["*.parabola"])),
]

def build_aug_fastas(rpath):
    fastas = [rpath]
    for i in range(1, max(N_VALS)):
        apath = WORK / f"aug_{i}.fa"
        if apath.exists():
            apath.unlink()
        try:
            os.link(rpath, apath)
        except OSError:
            shutil.copy(rpath, apath)
        fastas.append(apath)
    return fastas

def build_accuracy_fastas(rpath, snp_only):
    fastas = [rpath]
    expected_dists = {}

    for d in DISTS:
        mname = f"mut_d{d:.3f}"
        mpath = WORK / f"{mname}.fa"

        td = mutate(rpath, mpath, mname, d, snp_only)

        fastas.append(mpath)
        expected_dists[mname] = td

    return fastas, expected_dists

def bench_time_once(tool, aug_fastas, k):
    row_times = []
    for n in N_VALS:
        tool.clean()
        elapsed, _ = tool.run_batch(aug_fastas[:n], TIME_C, k)
        row_times.append(f"{elapsed:.4f}")
    tool.clean()
    return row_times

def score_accuracy(matrix, expected_dists):
    row_errs = []

    for d in DISTS:
        mname = f"mut_d{d:.3f}"
        td = expected_dists[mname]

        pred = matrix.get(("ref", mname)) or matrix.get((mname, "ref"))
        if pred is None:
            row_errs.append("nan")
            continue

        try:
            err = abs(td - pred)
            row_errs.append(f"{err:.6f}")
        except TypeError:
            row_errs.append("nan")

    return row_errs

def accuracy_c_values(tool):
    return [TIME_C] if tool.run_once_per_replicate else C_VALS

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--snp", action="store_true")
    args = parser.parse_args()

    WORK.mkdir(exist_ok=True)

    with open("1benchmark.csv", "w", newline="") as cf:
        w = csv.writer(cf)
        time_headers = [f"Time_N={n}" for n in N_VALS]
        err_headers = [f"ERR_d={d:.3f}" for d in DISTS]

        w.writerow(
            ["Genome", "k", "Tool", "C_rate", "Mode"]
            + time_headers
            + err_headers
        )

        mode_str = "SNP_ONLY" if args.snp else "SNP_INDEL"
        print(f"[*] Starting benchmark in {mode_str} mode...")

        for rep in range(1, REPLICATES + 1):
            current_seed = SEED + rep
            random.seed(current_seed)

            rpath = WORK / "ref.fa"
            if rep == 1 or not rpath.exists():
                build_ref(GENOME_FASTA, rpath)

            aug_fastas = build_aug_fastas(rpath)

            acc_fastas, expected_dists = build_accuracy_fastas(rpath, args.snp)

            for tool in TOOLS:
                if tool.name in ["Skani", "FastANI"]:
                    k = "-"
                    row_times = bench_time_once(tool, aug_fastas, k)
                    for c in accuracy_c_values(tool):
                        tool.clean()
                        _, matrix = tool.run_batch(acc_fastas, c, k)
                        print(f"{tool.name} k={k} C={c} Replicate={rep}")
                        row_errs = score_accuracy(matrix, expected_dists)
                        w.writerow(["genome", k, tool.name, c, mode_str] + row_times + row_errs)
                        cf.flush()
                    tool.clean()

                elif tool.name == "Mash" or tool.name == "Parabola":
                    for k in K_VALS:
                        row_times = bench_time_once(tool, aug_fastas, k)
                        for c in accuracy_c_values(tool):
                            tool.clean()
                            _, matrix = tool.run_batch(acc_fastas, c, k)
                            print(f"{tool.name} k={k} C={c} Replicate={rep}")
                            row_errs = score_accuracy(matrix, expected_dists)
                            w.writerow(["genome", k, tool.name, c, mode_str] + row_times + row_errs)
                            cf.flush()
                        tool.clean()

if __name__ == "__main__":
    main()