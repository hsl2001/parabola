import os, sys, subprocess, csv
from pathlib import Path
from glob import glob
import math

WORK = Path("embryophyte_run")
PARABOLA = Path("../parabola/parabola")

S = 1000
THREADS = "48"
K_VALUES = range(21, 32)

FASTAS = list(glob("./*.fasta"))

CLADES = {
    "All": {
        "name": "All",
        "taxa": ["Saccharomyces_cerevisiae", "Chlamydomonas_reinhardtii", "Marchantia_polymorpha", "Amborella_trichopoda", "Zostera_marina", "Oryza_sativa", "Zea_mays", "Aquilegia_coerulea", "Solanum_lycopersicum", "Glycine_max", "Salix_purpurea", "Aethionema_arabicum", "Capsella_rubella", "Arabidopsis_lyrata", "Arabidopsis_thaliana_Col-0", "Arabidopsis_thaliana_Cvi-0"],
        "ref": "(Saccharomyces_cerevisiae,(Chlamydomonas_reinhardtii,(Marchantia_polymorpha,(Amborella_trichopoda,((Zostera_marina,(Oryza_sativa,Zea_mays)),(Aquilegia_coerulea,(Solanum_lycopersicum,((Glycine_max,Salix_purpurea),(Aethionema_arabicum,(Capsella_rubella,(Arabidopsis_lyrata,(Arabidopsis_thaliana_Col-0,Arabidopsis_thaliana_Cvi-0))))))))))));"
    },
    "Viridiplantae": {
        "name": "Viridiplantae",
        "taxa": ["Chlamydomonas_reinhardtii", "Marchantia_polymorpha", "Amborella_trichopoda", "Zostera_marina", "Oryza_sativa", "Zea_mays", "Aquilegia_coerulea", "Solanum_lycopersicum", "Glycine_max", "Salix_purpurea", "Aethionema_arabicum", "Capsella_rubella", "Arabidopsis_lyrata", "Arabidopsis_thaliana_Col-0", "Arabidopsis_thaliana_Cvi-0"],
        "ref": "(Chlamydomonas_reinhardtii,(Marchantia_polymorpha,(Amborella_trichopoda,((Zostera_marina,(Oryza_sativa,Zea_mays)),(Aquilegia_coerulea,(Solanum_lycopersicum,((Glycine_max,Salix_purpurea),(Aethionema_arabicum,(Capsella_rubella,(Arabidopsis_lyrata,(Arabidopsis_thaliana_Col-0,Arabidopsis_thaliana_Cvi-0)))))))))));"
    },
    "Embryophyta": {
        "name": "Embryophyte",
        "taxa": ["Marchantia_polymorpha", "Amborella_trichopoda", "Zostera_marina", "Oryza_sativa", "Zea_mays", "Aquilegia_coerulea", "Solanum_lycopersicum", "Glycine_max", "Salix_purpurea", "Aethionema_arabicum", "Capsella_rubella", "Arabidopsis_lyrata", "Arabidopsis_thaliana_Col-0", "Arabidopsis_thaliana_Cvi-0"],
        "ref": "(Marchantia_polymorpha,(Amborella_trichopoda,((Zostera_marina,(Oryza_sativa,Zea_mays)),(Aquilegia_coerulea,(Solanum_lycopersicum,((Glycine_max,Salix_purpurea),(Aethionema_arabicum,(Capsella_rubella,(Arabidopsis_lyrata,(Arabidopsis_thaliana_Col-0,Arabidopsis_thaliana_Cvi-0))))))))));"
    },
    "Angiosperms": {
        "name": "Angiosperm",
        "taxa": ["Amborella_trichopoda", "Zostera_marina", "Oryza_sativa", "Zea_mays", "Aquilegia_coerulea", "Solanum_lycopersicum", "Glycine_max", "Salix_purpurea", "Aethionema_arabicum", "Capsella_rubella", "Arabidopsis_lyrata", "Arabidopsis_thaliana_Col-0", "Arabidopsis_thaliana_Cvi-0"],
        "ref": "(Amborella_trichopoda,((Zostera_marina,(Oryza_sativa,Zea_mays)),(Aquilegia_coerulea,(Solanum_lycopersicum,((Glycine_max,Salix_purpurea),(Aethionema_arabicum,(Capsella_rubella,(Arabidopsis_lyrata,(Arabidopsis_thaliana_Col-0,Arabidopsis_thaliana_Cvi-0)))))))));"
    },
    "Eudicots": {
        "name": "Eudicot",
        "taxa": ["Aquilegia_coerulea", "Solanum_lycopersicum", "Glycine_max", "Salix_purpurea", "Aethionema_arabicum", "Capsella_rubella", "Arabidopsis_lyrata", "Arabidopsis_thaliana_Col-0", "Arabidopsis_thaliana_Cvi-0"],
        "ref": "(Aquilegia_coerulea,(Solanum_lycopersicum,((Glycine_max,Salix_purpurea),(Aethionema_arabicum,(Capsella_rubella,(Arabidopsis_lyrata,(Arabidopsis_thaliana_Col-0,Arabidopsis_thaliana_Cvi-0)))))));"
    },
    "Rosids": {
        "name": "Rosid",
        "taxa": ["Glycine_max", "Salix_purpurea", "Aethionema_arabicum", "Capsella_rubella", "Arabidopsis_lyrata", "Arabidopsis_thaliana_Col-0", "Arabidopsis_thaliana_Cvi-0"],
        "ref": "((Glycine_max,Salix_purpurea),(Aethionema_arabicum,(Capsella_rubella,(Arabidopsis_lyrata,(Arabidopsis_thaliana_Col-0,Arabidopsis_thaliana_Cvi-0)))));"
    },
    "Brassicaceae": {
        "name": "Brassicaceae",
        "taxa": ["Aethionema_arabicum", "Capsella_rubella", "Arabidopsis_lyrata", "Arabidopsis_thaliana_Col-0", "Arabidopsis_thaliana_Cvi-0"],
        "ref": "(Aethionema_arabicum,(Capsella_rubella,(Arabidopsis_lyrata,(Arabidopsis_thaliana_Col-0,Arabidopsis_thaliana_Cvi-0))));"
    }
}

def get_stem(p_str):
    return Path(p_str).stem.split(".")[0]

def capture(cmd):
    return subprocess.run(cmd, capture_output=True, text=True).stdout

def parse_triangle(out):
    matrix = {}
    lines = out.strip().splitlines()
    if not lines: return matrix
    names = [get_stem(cols[0]) for cols in (line.split("\t") for line in lines[1:])]
    for line in lines[1:]:
        cols = line.split("\t")
        r_name = get_stem(cols[0])
        for c_idx, dist_str in enumerate(cols[1:]):
            if dist_str.strip():
                try:
                    dist_val = float(dist_str)
                except ValueError:
                    dist_val = float('nan')

                c_name = names[c_idx]
                matrix[(r_name, c_name)] = matrix[(c_name, r_name)] = dist_val
    return matrix

def build_fastme_tree(names, dist_matrix, work_dir, prefix="tree"):
    phylip_path = work_dir / f"{prefix}_dist.txt"
    tree_path = work_dir / f"{prefix}_tree.nwk"

    if tree_path.exists():
        tree_path.unlink()

    with open(phylip_path, "w") as f:
        f.write(f"{len(names)}\n")
        for name1 in names:
            f.write(f"{name1} ")
            row_dists = [f"{dist_matrix.get((name1, n2), float('nan')):.6f}" if name1 != n2 else "0.000000" for n2 in names]
            f.write(" ".join(row_dists) + "\n")

    subprocess.run(["fastme", "-i", str(phylip_path), "-o", str(tree_path), "-s"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    if tree_path.exists():
        with open(tree_path, "r") as f:
            return f.read().strip().replace('\n', '').replace('\r', '')
    return None

class ToolWrapper:
    def __init__(self, name):
        self.name = name
    def clean(self): pass
    def run_batch(self, fastas, seed, k): raise NotImplementedError

class MashRunner(ToolWrapper):
    def __init__(self):
        super().__init__("Mash")
    def clean(self):
        for p in WORK.glob("*.msh"):
            try: p.unlink()
            except OSError: pass
    def run_batch(self, fastas, seed, k):
        out_msh = str(WORK / f"mash_{seed}_k{k}")
        sketch_cmd = ["mash", "sketch", "-p", THREADS, "-k", str(k), "-s", str(S), "-S", str(seed), "-o", out_msh] + [str(p) for p in fastas]
        subprocess.run(sketch_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        dist_cmd = ["mash", "triangle", "-p", THREADS, f"{out_msh}.msh"]
        out = capture(dist_cmd)
        return parse_triangle(out)

class ParabolaRunner(ToolWrapper):
    def __init__(self):
        super().__init__("Parabola")
    def clean(self):
        for p in FASTAS:
            try: Path(str(p) + ".parabola").unlink()
            except OSError: pass
    def run_batch(self, fastas, seed, k):
        fasta_strs = [str(p) for p in fastas]
        sketch_cmd = [str(PARABOLA), "sketch", "-k", str(k), "-s", str(S), "-p", THREADS, "-e", str(seed)] + fasta_strs
        subprocess.run(sketch_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        triangle_cmd = [str(PARABOLA), "triangle"] + [f"{p}.parabola" for p in fasta_strs]
        out = capture(triangle_cmd)
        return parse_triangle(out)

def main():
    WORK.mkdir(exist_ok=True)
    available_names = [get_stem(p) for p in FASTAS]

    tools = [MashRunner(), ParabolaRunner()]

    csv_file_path = "nRF_combined.csv"
    with open(csv_file_path, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['Tool', 'Taxon', 'K', 'RF', 'nRF'])

        for target_tool in tools:
            for k in K_VALUES:
                print(f"\n{'='*50}")
                print(f"[*] K = {k} | Tool = {target_tool.name} | Calculating...")
                print(f"{'='*50}")

                dist_matrix = target_tool.run_batch(FASTAS, 42, k)

                for clade_id, data in CLADES.items():
                    name = data["name"]
                    clade_taxa = [t for t in data["taxa"] if t in available_names]

                    if len(clade_taxa) < 3:
                        csv_writer.writerow([target_tool.name.lower(), name, k, "NaN", "NaN"])
                        continue

                    ref_path = WORK / f"ref_{clade_id}.nwk"
                    with open(ref_path, "w") as f:
                        f.write(data["ref"])

                    prefix = f"{target_tool.name.lower()}_{clade_id}_k{k}"
                    build_fastme_tree(clade_taxa, dist_matrix, WORK, prefix=prefix)

                    pred_tree_path = WORK / f"{prefix}_tree.nwk"

                    if not pred_tree_path.exists():
                        csv_writer.writerow([target_tool.name.lower(), name, k, "NaN", "NaN"])
                        continue

                    out = capture(["phykit", "rf", str(ref_path), str(pred_tree_path)])
                    rf_result_raw = out.strip()

                    if rf_result_raw:
                        parts = rf_result_raw.split()
                        rf_val = parts[0] if len(parts) > 0 else "NaN"
                        nrf_val = parts[1] if len(parts) > 1 else "NaN"
                    else:
                        rf_val, nrf_val = "NaN", "NaN"

                    print(f"{name:<10} | {rf_result_raw}")
                    csv_writer.writerow([target_tool.name.lower(), name, k, rf_val, nrf_val])

                target_tool.clean()
                csv_file.flush()

    print(f"\n[*] Results saved to {csv_file_path}")

if __name__ == "__main__":
    main()
