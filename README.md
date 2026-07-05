# Parabola
## Fast and accurate genome sketching in divergent sequences

### Getting started
```bash
# Quick installation by pre-compiled binaries
wget https://github.com/hsl2001/parabola/releases/download/1.0/parabola
chmod a+x parabola
./parabola -h
export PATH=$PATH:`pwd`

# Quick installation from github
git clone https://github.com/hsl2001/parabola.git
cd parabola
make
# or, `make static` for static compile
./parabola -h
export PATH=$PATH:`pwd`

# Simply measure distance between test1 and test2
parabola dist benchmark/test1.fasta benchmark/test2.fasta
# Reference File      : benchmark/test1.fasta
# Query File          : benchmark/test2.fasta
# Shared Hashes       : 4
# ------------------------------------
# Average Containment : 0.020254
# Parabola Distance   : 0.169466

# Compute lower triangular distance metrix
parabola triangle benchmark/test1.fasta benchmark/test2.fasta benchmark/test3.fasta
# 3
# benchmark/test1.fasta
# benchmark/test2.fasta   0.169466
# benchmark/test3.fasta   0.173765        0.182878
```

### Background
Alignment-free methods were developed to avoid limitations of alignment-based methods, such as large computational demands which often reach as infeasible on large genomic datasets.
As one of the alignment-free methods, genome sketching method circumvents complex alignments from large genomic sequences. 
However, existing genome sketching methods often underestimates evolutionary distance (often assumed as genomic dissimilarity) due to sequence homoplasy, which indicates homology arised by random chances.
To resolve these accuracy limitations, Parabola introduces a distance estimation method driven by a 3-bit encoding system called Reverb and Average Containment metric.

### Synopsis
#### Installation using conda
```bash
conda install -c conda-forge -c bioconda parabola
```
Once installed via Conda, the `parabola` command will be available directly in your terminal path, allowing you to seamlessly sketch and compare genomes without manual compilation.

### Implementation details
#### Reverb encoding system
Parabola utilizes a symmetric 3-bit nucleotide encoding system, Reverb, where the bases A, C, G, and T are encoded as $001_{(2)}$, $110_{(2)}$, $011_{(2)}$, and $100_{(2)}$, respectively.
By design, complementary base pairs (A/T and C/G) are bit-reflections of each other.
This structural symmetry allows Parabola to derive the reverse complement hash instantaneously using a 128-bit reversal operation, bypassing the need for an independent, character-wise traversal of the reverse strand.

#### Average Containment distance estimation
Parabola calculates the Average Containment $C_{avg}$ of two sketches $A$ and $B$ as the arithmetic mean of the proportions of shared hashes:
$$C_{avg} = \frac{1}{2} \left( \frac{|A \cap B|}{|A|} + \frac{|A \cap B|}{|B|} \right)$$
The distance is then estimated as $d = 1 - C_{avg}^{1/k}$.

#### FracMinHash
Parabola employs a scaled hashing approach (FracMinHash) to compress genome sequences.
The algorithm establishes a strict hash threshold based on the defined scale factor (e.g., `UINT64_MAX / scale`).
Only generated *k*-mer hashes that fall below this numeric threshold are retained, creating a uniform fractional representation of the genome that remains stable even under high compression rates.

#### Overflowing Max-Heap
To manage the filtered hashes, Parabola buffers valid sequence hashes into a dynamically resizing memory pool.
Once the input stream is fully processed, the pool is finalized by sorting the hashes and deduplicating them.

#### The `three` command for three-way intersection
Parabola features a specific `three` command designed to simultaneously calculate distances between a single reference sketch and two separate query sketches.
The syntax for this operation is `parabola three <ref> <query1> <query2>`. 
Internally, the algorithm concurrently traverses the sorted hash sets of all three sequences to evaluate their intersections in a single pass.
As a result, it outputs three pairwise distance reports at once: reference versus query 1, reference versus query 2, and query 1 versus query 2.

### NGS Support
Parabola (v1.0) is currently only tested in assembled genomes.
We support NGS reads (both short or long) but do not garuntee the quality of results.

#### Count-Min-Sketch
To track *k*-mer frequencies on the fly without consuming large memory, Parabola implements a Count-Min Sketch data structure.
When the `-m` threshold is greater than 1, Parabola generates three distinct hash values for each *k*-mer and increments their corresponding bins in the sketch table.
A *k*-mer is only inserted into the final Parabola sketch pool if its minimum estimated count across the three bins meets the specified threshold requirement.

#### Setting parameters for `read` mode
When processing raw NGS reads, users must apply the `-r` flag to pool all reads from the input files into a single, cohesive sketch.
Because reads contain sequencing errors, it is also highly recommended to pair this with the `-m` (minimum count) parameter (e.g., `-m 2`) to filter out singleton *k*-mers that are likely artifacts of sequencing errors.

### Citing Parabola
If Parabola is useful in your research, please cite:

Lim, H., & Hyun, Y. (2026). Fast and accurate sketching in divergent genomes with Parabola.

### License
Parabola uses wonderful work from [klib](https://attractivechaos.github.io/klib/#About), specifically `kseq.h` for sequence IO, `ketopt.h` for argument parsing and `kthreads.h` for multiprocessing.
The full source codes of Parabola are available under **MIT license**.
