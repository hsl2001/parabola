# Parabola
## Fast and accurate genome sketching in divergent sequences

### Getting started
```bash
# Quick installation by pre-compiled binaries
wget https://github.com/hsl2001/release/~~
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
# Shared Hashes       : 4 / 392
# ------------------------------------
# Jaccard Index       : 0.010204
# Parabola Distance   : 0.197915

# Compute lower triangular distance metrix
parabola triangle benchmark/test1.fasta benchmark/test2.fasta benchmark/test3.fasta
# 3
# benchmark/test1.fasta
# benchmark/test2.fasta   0.197915
# benchmark/test3.fasta   0.180695        0.183343
```

### Background
Alignment-free methods were developed to avoid limitations of alignment-based methods, such as large computational demands which often reach as infeasible on large genomic datasets.
As one of the alignment-free methods, genome sketching method circumvents complex alignments from large genomic sequences. 
However, existing genome sketching methods often underestimates evolutionary distance (often assumed as genomic dissimilarity) due to sequence homoplasy, which indicates homology arised by random chances.
To resolve these accuracy limitations, Parabola introduces a simplicity-aware distance calibration method driven by a 3-bit encoding system called Reverb.
Because low-complexity (high-simplicity) sequences gain homoplasy more easily, Parabola assesses sequence simplicity simultaneously during the hashing phase. 
This allows the algorithm to isolate true phylogenetic signatures while discarding homoplasious noise, resulting in highly accurate distance estimations in sequences exceeding 20% divergence.

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

#### Simplicity calculation with Simpson's index
To measure sequence simplicity without adding computational overhead, Parabola calculates the Simpson index using bitwise operations.
Leveraging the trifold population count of the 3-bit encoded hashes, Parabola dynamically calculates the simplicity score using the paraboloid formula $D = A^2 + C^2 + G^2 + T^2$.
This strand-agnostic scoring system maximizes the score for high-entropy sequences (where base frequencies approach 25%) and heavily penalizes uninformative extremums, such as homopolymers.

#### FracMinHash
Parabola employs a scaled hashing approach (FracMinHash) to compress genome sequences.
The algorithm establishes a strict hash threshold based on the defined scale factor (e.g., `UINT64_MAX / scale`).
Only generated *k*-mer hashes that fall below this numeric threshold are retained, creating a uniform fractional representation of the genome that remains stable even under high compression rates.

#### Overflowing Max-Heap
To manage the filtered hashes, Parabola buffers valid sequence hashes and their calculated simplicity scores into a dynamically resizing memory pool.
Once the input stream is fully processed, the pool is finalized by sorting the hashes and deduplicating them.
During deduplication, if identical hashes are found, Parabola retains only the lowest simplicity score associated with that hash, optimizing the final sketch representation.

#### Isolating orthologous hash intersection with simplicity
Parabola calibrates the final evolutionary distance by treating the intersection of two sketch sets as a mixture distribution of both orthologous hashes and hashes with homoplasy.
By assuming the variance of simplicity for the homoplasy subpopulation converges to zero, Parabola calculates the homoplasious ratio ($x$).
It then applies this ratio to a Poisson mutation model to generate a calibration term ($\alpha$), isolating true orthologous signals and correcting the naive Jaccard distance.

#### Juke-Cantor distance correction
Parabola provides an optional Jukes-Cantor distance correction, which can be enabled by passing the `-j` flag in the command-line interface.
When this option is active, the calculated base distance ($d$) is adjusted using the formula $d_{JC} = -0.75 \times \ln(1.0 - \frac{4}{3}d)$.
To prevent errors from extreme divergence, if the uncorrected distance is greater than or equal to 0.75, the corrected distance is automatically capped at a maximum value of 1.0.

#### The `three` command for three-way intersection
Parabola features a specific `three` command designed to simultaneously calculate distances between a single reference sketch and two separate query sketches.
The syntax for this operation is `parabola three [-j] <ref> <query1> <query2>`. 
Internally, the algorithm concurrently traverses the sorted hash sets of all three sequences to evaluate their intersections and simplicities in a single pass.
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
