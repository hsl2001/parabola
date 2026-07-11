#include <math.h>
#include <stdio.h>
#include <time.h>
#include <zlib.h>

#include "klib/ketopt.h"
#include "klib/kseq.h"
#include "reverb.h"

#define MIX_CONST1 0xff51afd7ed558ccdULL
#define MIX_CONST2 0xc4ceb9fe1a85ec53ULL

/* Generic dynamic-array push: grows `arr` by doubling `cap` as needed. */
#define DA_PUSH(arr, n, cap, val)                                              \
  do {                                                                         \
    if ((n) >= (cap)) {                                                        \
      (cap) = (cap) ? (cap) * 2 : 1024;                                        \
      (arr) = realloc((arr), (cap) * sizeof(*(arr)));                          \
    }                                                                          \
    (arr)[(n)++] = (val);                                                      \
  } while (0)

KSEQ_INIT(gzFile, gzread)

// ==============================================================
// UTILITIES
// ==============================================================

/* Reverb encoding: A = 001, C = 110, G = 011, T = 100 */
static const int8_t BASE_LOOKUP[256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, 1,  -1, 6,  -1, -1, -1, 3,  -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, 4,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 1,  -1, 6,  -1, -1, -1, 3,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 4,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1};

static uint64_t mix_hash(__uint128_t hash_value, uint64_t seed) {
  __uint128_t p = (hash_value ^ MIX_CONST1) * ((__uint128_t)seed ^ MIX_CONST2);
  return (uint64_t)(p ^ (p >> 64));
}

static uint64_t reverse_bits64(uint64_t n) {
#if defined(__aarch64__)
  uint64_t r;
  __asm__("rbit %0, %1" : "=r"(r) : "r"(n));
  return r;
#else
  uint64_t r = __builtin_bswap64(n);
  r = ((r & 0x5555555555555555ULL) << 1) | ((r & 0xAAAAAAAAAAAAAAAAULL) >> 1);
  r = ((r & 0x3333333333333333ULL) << 2) | ((r & 0xCCCCCCCCCCCCCCCCULL) >> 2);
  r = ((r & 0x0F0F0F0F0F0F0F0FULL) << 4) | ((r & 0xF0F0F0F0F0F0F0F0ULL) >> 4);
  return r;
#endif
}

static __uint128_t reverse_bits128(__uint128_t n) {
  return ((__uint128_t)reverse_bits64((uint64_t)n) << 64) |
         reverse_bits64((uint64_t)(n >> 64));
}

// ==============================================================
// INIT
// ==============================================================

void reverb_init(Reverb *r, size_t hash_window) {
  size_t k = hash_window < 42 ? hash_window : 42; /* 3*42=126 bits ≤ 128 */
  uint32_t kmer_bits = 3 * (uint32_t)k;

  __uint128_t remover_mask =
      (kmer_bits > 3) ? (((__uint128_t)1 << (kmer_bits - 3)) - 1) : 0;

  r->hash_window = k;
  r->remover_mask = remover_mask;
  r->kmer_bits = kmer_bits;
  r->rc_shift = (kmer_bits > 0) ? (128 - kmer_bits) : 128;
}

// ==============================================================
// HASH POOL
// ==============================================================

typedef struct {
  size_t size;
  size_t cap;
  uint64_t hash_threshold;
  uint64_t *hashes;
} HashPool;

static int cmp_uint64(const void *a, const void *b) {
  uint64_t x = *(const uint64_t *)a;
  uint64_t y = *(const uint64_t *)b;
  return (x > y) - (x < y);
}

static void pool_init(HashPool *pool, uint64_t threshold) {
  *pool = (HashPool){.hash_threshold = threshold};
}

static void pool_try_insert(HashPool *pool, uint64_t h) {
  if (h >= pool->hash_threshold)
    return;
  DA_PUSH(pool->hashes, pool->size, pool->cap, h);
}

static void pool_finalize(HashPool *pool, uint64_t **out_hashes,
                          size_t *out_size) {
  size_t n = pool->size;
  if (n) {
    qsort(pool->hashes, n, sizeof(uint64_t), cmp_uint64);
    /* Deduplicate in-place */
    size_t u = 0;
    for (size_t i = 0; i < n; i++) {
      if (u == 0 || pool->hashes[i] != pool->hashes[u - 1])
        pool->hashes[u++] = pool->hashes[i];
    }
    n = u;
    uint64_t *nh = realloc(pool->hashes, n * sizeof(uint64_t));
    if (nh)
      pool->hashes = nh;
  }
  *out_size = n;
  *out_hashes = n ? pool->hashes : NULL;
  pool->hashes = NULL;
  pool->size = pool->cap = 0;
}

// ==============================================================
// SKETCH EXTRACTION
// ==============================================================

__attribute__((hot)) static void extract_and_insert(const Reverb *r,
                                                    HashPool *pool,
                                                    const uint8_t *seq,
                                                    size_t len) {
  size_t K = r->hash_window;
  __uint128_t fwd = 0;
  size_t valid = 0;

  for (size_t idx = 0; idx < len; idx++) {
    int8_t lv = BASE_LOOKUP[seq[idx]];
    if (lv < 0) {
      fwd = 0;
      valid = 0;
      continue;
    }

    fwd = ((fwd & r->remover_mask) << 3) | (uint8_t)lv;
    if (valid < K)
      valid++;
    if (valid < K)
      continue;

    __uint128_t rev = reverse_bits128(fwd) >> r->rc_shift;
    __uint128_t canon = fwd < rev ? fwd : rev;
    uint64_t h = mix_hash(canon, r->hash_seed);

    if (h < pool->hash_threshold)
      pool_try_insert(pool, h);
  }
}

void reverb_sketch_free(ReverbSketch *sk) {
  if (sk) {
    free(sk->name);
    free(sk->hashes);
  }
}

// ==============================================================
// DISTANCE CALCULATION
// ==============================================================

ReverbDistResult reverb_dist(const ReverbSketch *ref,
                             const ReverbSketch *query) {
  ReverbDistResult res = {0.0, 1.0, 0};
  if (ref->kmer_size == 0 || query->kmer_size != ref->kmer_size)
    return res;

  size_t shared = 0, i = 0, j = 0;
  while (i < ref->sketch_size && j < query->sketch_size) {
    if (ref->hashes[i] == query->hashes[j]) {
      shared++;
      i++;
      j++;
    } else if (ref->hashes[i] < query->hashes[j]) {
      i++;
    } else {
      j++;
    }
  }

  res.shared_hashes = shared;
  if (ref->sketch_size > 0 && query->sketch_size > 0) {
    double ca = (double)shared / (double)ref->sketch_size;
    double cb = (double)shared / (double)query->sketch_size;
    res.containment = 0.5 * (ca + cb);
    res.distance = 1.0 - pow(res.containment, 1.0 / (double)ref->kmer_size);
  }
  return res;
}

// ==============================================================
// UNION-FIND
// ==============================================================

void uf_init(UnionFind *uf, size_t n) {
  uf->n = n;
  uf->parent = (uint32_t *)malloc(n * sizeof(uint32_t));
  uf->rank = (uint32_t *)calloc(n, sizeof(uint32_t));
  for (size_t i = 0; i < n; i++)
    uf->parent[i] = (uint32_t)i;
}

uint32_t uf_find(UnionFind *uf, uint32_t x) {
  while (uf->parent[x] != x) {
    uf->parent[x] = uf->parent[uf->parent[x]]; /* path splitting */
    x = uf->parent[x];
  }
  return x;
}

void uf_union(UnionFind *uf, uint32_t a, uint32_t b) {
  a = uf_find(uf, a);
  b = uf_find(uf, b);
  if (a == b)
    return;
  if (uf->rank[a] < uf->rank[b]) {
    uint32_t t = a;
    a = b;
    b = t;
  }
  uf->parent[b] = a;
  if (uf->rank[a] == uf->rank[b])
    uf->rank[a]++;
}

void uf_free(UnionFind *uf) {
  free(uf->parent);
  free(uf->rank);
}

// ==============================================================
// PARAMETERS
// ==============================================================

// Params removed, handled directly in cmd_dup

static void print_usage(void) {
  printf("Reverb: Ultra-fast Alignment-free Segmental Duplication Detection\n\n"
         "Usage: reverb [options] fasta1 [fasta2 ...]\n\n"
         "Options:\n"
         "  -k: kmer size (default: 21, max: 42)\n"
         "  -s: scale factor (default: 10)\n"
         "  -w: window size in bp (default: 10000)\n"
         "  -t: step size in bp (default: window/2)\n"
         "  -b: minimum valid bases per window (default: 1000)\n"
         "  -d: maximum distance to consider as copy (default: 0.05)\n"
         "  -m: minimum copy count (default: 2)\n"
         "  -M: maximum copy count to filter ubiquitous repeats (default: 50)\n"
         "  -o: output file prefix (default: reverb)\n"
         "  -p: number of threads (default: 8)\n"
         "  -h, --help: show this help message\n"
         "\n");
}

// ==============================================================
// DUP: WINDOW STREAMING
// ==============================================================

typedef struct {
  char *chrom;
  size_t start;
  size_t end;
} WindowCoord;

/* Inverted hash index entry: maps a hash value to its source window */
typedef struct {
  uint64_t hash;
  uint32_t window_id;
} HashWindowEntry;

static int dup_stream(const char *filename, const Reverb *r, uint64_t scale,
                      size_t window_size, size_t step_size, size_t min_bases,
                      ReverbSketch **sketches, WindowCoord **coords,
                      size_t *num_sketches, size_t *cap_sketches,
                      FILE *bed_fp) {
  gzFile fp = gzopen(filename, "r");
  if (!fp)
    return -1;
  kseq_t *ks = kseq_init(fp);
  if (!ks) {
    gzclose(fp);
    return -1;
  }

  while (kseq_read(ks) >= 0) {
    size_t len = ks->seq.l;
    for (size_t i = 0; i + window_size <= len; i += step_size) {
      /* Count valid bases in this window */
      size_t valid_bases = 0;
      for (size_t j = 0; j < window_size; j++) {
        if (BASE_LOOKUP[(uint8_t)ks->seq.s[i + j]] >= 0)
          valid_bases++;
      }
      if (valid_bases < min_bases)
        continue;

      /* Ensure capacity */
      if (*num_sketches >= *cap_sketches) {
        *cap_sketches = *cap_sketches == 0 ? 256 : *cap_sketches * 2;
        *sketches = realloc(*sketches, *cap_sketches * sizeof(ReverbSketch));
        *coords = realloc(*coords, *cap_sketches * sizeof(WindowCoord));
      }

      ReverbSketch *sk = &(*sketches)[*num_sketches];
      WindowCoord *wc = &(*coords)[*num_sketches];
      memset(sk, 0, sizeof(ReverbSketch));

      char namebuf[512];
      snprintf(namebuf, sizeof(namebuf), "%s_%zu_%zu", ks->name.s, i,
               i + window_size);
      sk->name = strdup(namebuf);
      sk->kmer_size = r->hash_window;
      sk->hash_threshold = UINT64_MAX / scale;

      wc->chrom = strdup(ks->name.s);
      wc->start = i;
      wc->end = i + window_size;

      if (bed_fp)
        fprintf(bed_fp, "%s\t%zu\t%zu\t%s\n", ks->name.s, i, i + window_size,
                sk->name);

      HashPool pool;
      pool_init(&pool, sk->hash_threshold);
      extract_and_insert(r, &pool, (const uint8_t *)ks->seq.s + i, window_size);
      pool_finalize(&pool, &sk->hashes, &sk->sketch_size);

      if (sk->sketch_size > 0) {
        (*num_sketches)++;
      } else {
        free(sk->name);
        free(wc->chrom);
      }
    }
  }
  kseq_destroy(ks);
  gzclose(fp);
  return 0;
}

// ==============================================================
// DUP: INVERTED HASH INDEX
// ==============================================================

static int cmp_hash_window_entry(const void *a, const void *b) {
  const HashWindowEntry *ea = (const HashWindowEntry *)a;
  const HashWindowEntry *eb = (const HashWindowEntry *)b;
  if (ea->hash != eb->hash)
    return (ea->hash > eb->hash) - (ea->hash < eb->hash);
  return (ea->window_id > eb->window_id) - (ea->window_id < eb->window_id);
}

typedef struct {
  uint32_t win_a;
  uint32_t win_b;
  uint32_t shared_count;
} CandidatePair;

static int cmp_candidate_pair(const void *a, const void *b) {
  const CandidatePair *pa = (const CandidatePair *)a;
  const CandidatePair *pb = (const CandidatePair *)b;
  if (pa->win_a != pb->win_a)
    return (pa->win_a < pb->win_a) ? -1 : 1;
  if (pa->win_b != pb->win_b)
    return (pa->win_b < pb->win_b) ? -1 : 1;
  return 0;
}

/* Build inverted hash index and find candidate pairs.
 * Returns edges where distance < max_dist, excluding adjacent windows
 * on the same chromosome. */
static size_t build_candidate_edges(ReverbSketch *sketches, WindowCoord *coords,
                                    size_t n_windows, double max_dist,
                                    size_t window_size,
                                    ReverbDupEdge **out_edges) {
  /* 1. Flatten all (hash, window_id) entries */
  size_t total_entries = 0;
  for (size_t i = 0; i < n_windows; i++)
    total_entries += sketches[i].sketch_size;

  if (total_entries == 0) {
    *out_edges = NULL;
    return 0;
  }

  HashWindowEntry *entries = malloc(total_entries * sizeof(HashWindowEntry));
  size_t idx = 0;
  for (size_t i = 0; i < n_windows; i++) {
    for (size_t j = 0; j < sketches[i].sketch_size; j++) {
      entries[idx++] = (HashWindowEntry){.hash = sketches[i].hashes[j],
                                         .window_id = (uint32_t)i};
    }
  }

  /* 2. Sort by hash value */
  qsort(entries, total_entries, sizeof(HashWindowEntry), cmp_hash_window_entry);

  /* 3. For each group sharing the same hash, emit (win_a, win_b) pairs */
  size_t n_candidates = 0, cap_candidates = 0;
  CandidatePair *candidates = NULL;

  size_t run_start = 0;
  while (run_start < total_entries) {
    size_t run_end = run_start + 1;
    while (run_end < total_entries &&
           entries[run_end].hash == entries[run_start].hash)
      run_end++;

    size_t run_len = run_end - run_start;
    if (run_len >= 2 && run_len <= 100) {
      for (size_t i = run_start; i < run_end; i++) {
        for (size_t j = i + 1; j < run_end; j++) {
          uint32_t a = entries[i].window_id;
          uint32_t b = entries[j].window_id;
          if (a == b)
            continue;
          if (a > b) {
            uint32_t t = a;
            a = b;
            b = t;
          }
          DA_PUSH(candidates, n_candidates, cap_candidates,
                  ((CandidatePair){a, b, 1}));
        }
      }
    }
    run_start = run_end;
  }
  free(entries);

  if (n_candidates == 0) {
    *out_edges = NULL;
    return 0;
  }

  /* 4. Sort and merge duplicates to get shared_count */
  qsort(candidates, n_candidates, sizeof(CandidatePair), cmp_candidate_pair);

  size_t n_unique = 0;
  for (size_t i = 0; i < n_candidates; i++) {
    if (n_unique > 0 && candidates[n_unique - 1].win_a == candidates[i].win_a &&
        candidates[n_unique - 1].win_b == candidates[i].win_b) {
      candidates[n_unique - 1].shared_count++;
    } else {
      candidates[n_unique++] = candidates[i];
    }
  }

  /* 5. For pairs with enough shared hashes, compute full distance */
  size_t n_edges = 0, cap_edges = 0;
  ReverbDupEdge *edges = NULL;

  for (size_t i = 0; i < n_unique; i++) {
    uint32_t a = candidates[i].win_a;
    uint32_t b = candidates[i].win_b;

    /* Filter out adjacent/overlapping windows on the same chromosome */
    if (strcmp(coords[a].chrom, coords[b].chrom) == 0) {
      size_t dist_bp = (coords[a].start < coords[b].start)
                           ? coords[b].start - coords[a].start
                           : coords[a].start - coords[b].start;
      if (dist_bp < window_size)
        continue;
    }

    if (candidates[i].shared_count < 2)
      continue;

    ReverbDistResult d = reverb_dist(&sketches[a], &sketches[b]);
    if (d.distance < max_dist)
      DA_PUSH(edges, n_edges, cap_edges, ((ReverbDupEdge){a, b, d.distance}));
  }

  free(candidates);
  *out_edges = edges;
  return n_edges;
}

// ==============================================================
// DUP: SEGMENT MERGE
// ==============================================================

static int cmp_dup_region(const void *a, const void *b) {
  const ReverbDupRegion *ra = (const ReverbDupRegion *)a;
  const ReverbDupRegion *rb = (const ReverbDupRegion *)b;
  int c = strcmp(ra->chrom, rb->chrom);
  if (c != 0)
    return c;
  return (ra->start > rb->start) - (ra->start < rb->start);
}

/* Merge adjacent/overlapping regions in the same SD family.
 * Returns the new count of merged regions. */
static size_t merge_dup_regions(ReverbDupRegion *regions, size_t n) {
  if (n <= 1)
    return n;

  qsort(regions, n, sizeof(ReverbDupRegion), cmp_dup_region);

  size_t out = 0;
  for (size_t i = 1; i < n; i++) {
    if (regions[i].family_id == regions[out].family_id &&
        strcmp(regions[i].chrom, regions[out].chrom) == 0 &&
        regions[i].start <= regions[out].end) {
      if (regions[i].end > regions[out].end)
        regions[out].end = regions[i].end;
      regions[out].avg_distance =
          (regions[out].avg_distance + regions[i].avg_distance) / 2.0;
    } else {
      out++;
      if (out != i)
        regions[out] = regions[i];
    }
  }
  return out + 1;
}

// ==============================================================
// DUP: OUTPUT HELPERS
// ==============================================================

static int cmp_edge_for_bedpe(const void *p1, const void *p2) {
  const ReverbDupEdge *e1 = (const ReverbDupEdge *)p1;
  const ReverbDupEdge *e2 = (const ReverbDupEdge *)p2;
  if (e1->win_a != e2->win_a)
    return (e1->win_a < e2->win_a) ? -1 : 1;
  if (e1->distance < e2->distance)
    return -1;
  if (e1->distance > e2->distance)
    return 1;
  return 0;
}

static size_t write_bedpe_output(const char *path, ReverbDupEdge *edges,
                                 size_t n_edges, WindowCoord *coords,
                                 UnionFind *uf, uint32_t *copy_counts,
                                 char **hub_label, int min_copy, int max_copy) {
  FILE *fp = fopen(path, "w");
  if (!fp) {
    fprintf(stderr, "Error: cannot open %s\n", path);
    return 0;
  }

  if (n_edges > 0) {
    qsort(edges, n_edges, sizeof(ReverbDupEdge), cmp_edge_for_bedpe);
  }

  fprintf(fp, "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tfamily\t"
              "distance\tcopy_count\n");
  size_t reported = 0;
  uint32_t last_win_a = (uint32_t)-1;
  for (size_t i = 0; i < n_edges; i++) {
    uint32_t a = edges[i].win_a, b = edges[i].win_b;
    if (a == last_win_a)
      continue;
    last_win_a = a;

    uint32_t fam = uf_find(uf, a);
    uint32_t cc = copy_counts[fam];
    if ((int)cc < min_copy || (max_copy > 0 && (int)cc > max_copy))
      continue;
    const char *label = hub_label[fam] ? hub_label[fam] : "unknown";
    fprintf(fp, "%s\t%zu\t%zu\t%s\t%zu\t%zu\t%s\t%.6f\t%u\n", coords[a].chrom,
            coords[a].start, coords[a].end, coords[b].chrom, coords[b].start,
            coords[b].end, label, edges[i].distance, cc);
    reported++;
  }
  fclose(fp);
  return reported;
}

typedef struct {
  size_t n_output;
  size_t n_families;
} BedOutputStats;

static BedOutputStats write_bed_output(const char *path,
                                       ReverbDupRegion *regions,
                                       size_t n_merged, char **hub_label,
                                       int min_copy, int max_copy) {
  BedOutputStats stats = {0, 0};
  FILE *fp = fopen(path, "w");
  if (!fp) {
    fprintf(stderr, "Error: cannot open %s\n", path);
    return stats;
  }

  fprintf(fp, "#chrom\tstart\tend\tfamily\tcopy_count\n");
  uint32_t last_fam = UINT32_MAX;
  for (size_t i = 0; i < n_merged; i++) {
    if ((int)regions[i].copy_count < min_copy ||
        (max_copy > 0 && (int)regions[i].copy_count > max_copy))
      continue;
    const char *label = hub_label[regions[i].family_id]
                            ? hub_label[regions[i].family_id]
                            : "unknown";
    fprintf(fp, "%s\t%zu\t%zu\t%s\t%u\n", regions[i].chrom, regions[i].start,
            regions[i].end, label, regions[i].copy_count);
    if (regions[i].family_id != last_fam) {
      stats.n_families++;
      last_fam = regions[i].family_id;
    }
    stats.n_output++;
  }
  fclose(fp);
  return stats;
}

static void print_summary(int num_files, size_t num_sketches, size_t n_edges,
                          size_t reported_edges, size_t n_families,
                          size_t n_merged, double elapsed,
                          const char *out_prefix) {
  fprintf(stderr,
          "\n=== Reverb Summary ===\n"
          "Input files        : %d\n"
          "Total windows      : %zu\n"
          "Candidate edges    : %zu\n"
          "Reported SD pairs  : %zu\n"
          "SD families        : %zu\n"
          "Merged segments    : %zu\n"
          "Elapsed time       : %.1fs\n"
          "Output prefix      : %s\n"
          "=======================\n",
          num_files, num_sketches, n_edges, reported_edges, n_families,
          n_merged, elapsed, out_prefix);
}

// ==============================================================
// CMD: DUP
// ==============================================================

int cmd_dup(int argc, char **argv) {
  uint32_t def_kmer_size = 21;
  uint64_t def_scale = 10;
  uint64_t def_hash_seed = 42;
  size_t window_size = 5000;
  size_t step_size = 0; /* 0 = auto (window/2) */
  size_t min_bases = 1000;
  double max_dist = 0.05;
  int min_copy = 3;
  int max_copy = 50;
  const char *out_prefix = "reverb";
  int n_threads = 8;

  ketopt_t opt = KETOPT_INIT;
  int c;
  while ((c = ketopt(&opt, argc, argv, 1, "k:s:e:w:t:b:d:m:M:o:p:h", 0)) >= 0) {
    if (c == 'h') {
      print_usage();
      return 0;
    } else if (c == 'k')
      def_kmer_size = (uint32_t)atoi(opt.arg);
    else if (c == 's')
      def_scale = (uint64_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'e')
      def_hash_seed = strtoull(opt.arg, NULL, 0);
    else if (c == 'w')
      window_size = (size_t)strtoull(opt.arg, NULL, 10);
    else if (c == 't')
      step_size = (size_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'b')
      min_bases = (size_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'd')
      max_dist = atof(opt.arg);
    else if (c == 'm')
      min_copy = atoi(opt.arg);
    else if (c == 'M')
      max_copy = atoi(opt.arg);
    else if (c == 'o')
      out_prefix = opt.arg;
    else if (c == 'p')
      n_threads = atoi(opt.arg) < 1 ? 1 : atoi(opt.arg);
    else
      return 1;
  }

  if (step_size == 0)
    step_size = window_size / 2;
  (void)n_threads;

  int num_files = argc - opt.ind;
  if (num_files < 1) {
    fprintf(stderr, "Error: missing input FASTA files\n");
    return 1;
  }
  char **in_files = argv + opt.ind;

  struct timespec t_start, t_end;
  clock_gettime(CLOCK_MONOTONIC, &t_start);

  Reverb r;
  reverb_init(&r, def_kmer_size);
  r.hash_seed = def_hash_seed;

  /* Shared state that needs cleanup */
  int ret = 0;
  ReverbSketch *sketches = NULL;
  WindowCoord *coords = NULL;
  size_t num_sketches = 0, cap_sketches = 0;
  ReverbDupEdge *edges = NULL;
  size_t n_edges = 0;
  ReverbDupRegion *dup_regions = NULL;
  uint32_t *comp_size = NULL;
  char **hub_label = NULL;
  UnionFind uf = {0};

  /* Phase 1: Extract sliding window sketches */
  char path_buf[PATH_MAX];
  snprintf(path_buf, sizeof(path_buf), "%s.window.bed", out_prefix);
  FILE *bed_fp = fopen(path_buf, "w");
  if (!bed_fp) {
    fprintf(stderr, "Error: cannot open %s for writing\n", path_buf);
    ret = 1;
    goto cleanup;
  }

  fprintf(stderr, "[reverb] Extracting windows (w=%zu, step=%zu) ...\n",
          window_size, step_size);
  for (int i = 0; i < num_files; i++)
    dup_stream(in_files[i], &r, def_scale, window_size, step_size, min_bases,
               &sketches, &coords, &num_sketches, &cap_sketches, bed_fp);
  fclose(bed_fp);

  if (num_sketches == 0) {
    fprintf(stderr, "Error: no valid windows extracted.\n");
    ret = 1;
    goto cleanup;
  }
  fprintf(stderr, "[reverb] Total windows: %zu\n", num_sketches);

  /* Phase 2: Build inverted hash index and find candidate edges */
  fprintf(stderr, "[reverb] Building hash index and finding candidates ...\n");
  n_edges = build_candidate_edges(sketches, coords, num_sketches, max_dist,
                                  window_size, &edges);
  fprintf(stderr, "[reverb] Duplicated edges: %zu\n", n_edges);

  if (n_edges == 0) {
    fprintf(stderr, "[reverb] No segmental duplications found.\n");
    goto cleanup;
  }

  /* Phase 3: Union-Find clustering */
  fprintf(stderr, "[reverb] Clustering SD families ...\n");
  uf_init(&uf, num_sketches);
  for (size_t i = 0; i < n_edges; i++)
    uf_union(&uf, edges[i].win_a, edges[i].win_b);

  comp_size = calloc(num_sketches, sizeof(uint32_t));
  for (size_t i = 0; i < num_sketches; i++)
    comp_size[uf_find(&uf, (uint32_t)i)]++;

  /* Find hub window (highest degree) per family */
  uint32_t *win_degree = calloc(num_sketches, sizeof(uint32_t));
  for (size_t i = 0; i < n_edges; i++) {
    win_degree[edges[i].win_a]++;
    win_degree[edges[i].win_b]++;
  }

  uint32_t *hub_window = malloc(num_sketches * sizeof(uint32_t));
  uint32_t *hub_degree = calloc(num_sketches, sizeof(uint32_t));
  for (size_t i = 0; i < num_sketches; i++)
    hub_window[i] = (uint32_t)i;
  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t fam = uf_find(&uf, (uint32_t)i);
    if (win_degree[i] > hub_degree[fam]) {
      hub_degree[fam] = win_degree[i];
      hub_window[fam] = (uint32_t)i;
    }
  }

  /* Build hub labels: "chrom@start-end" */
  hub_label = calloc(num_sketches, sizeof(char *));
  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t fam = uf_find(&uf, (uint32_t)i);
    if (fam == i && comp_size[fam] >= 2) {
      uint32_t hw = hub_window[fam];
      char buf[512];
      snprintf(buf, sizeof(buf), "%s@%zu-%zu", coords[hw].chrom,
               coords[hw].start, coords[hw].end);
      hub_label[fam] = strdup(buf);
    }
  }
  free(win_degree);
  free(hub_degree);
  free(hub_window);

  /* Collect windows that belong to multi-member families */
  size_t n_dup_regions = 0, cap_dup_regions = 0;
  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t fam = uf_find(&uf, (uint32_t)i);
    if (comp_size[fam] < 2)
      continue;
    DA_PUSH(dup_regions, n_dup_regions, cap_dup_regions,
            ((ReverbDupRegion){.chrom = coords[i].chrom,
                               .start = coords[i].start,
                               .end = coords[i].end,
                               .family_id = fam,
                               .copy_count = comp_size[fam],
                               .avg_distance = 0.0}));
  }

  size_t n_merged = merge_dup_regions(dup_regions, n_dup_regions);

  /* Recount per-family copy counts after merge, reusing comp_size */
  memset(comp_size, 0, num_sketches * sizeof(uint32_t));
  for (size_t i = 0; i < n_merged; i++)
    comp_size[dup_regions[i].family_id]++;
  for (size_t i = 0; i < n_merged; i++)
    dup_regions[i].copy_count = comp_size[dup_regions[i].family_id];

  /* Phase 4: Write BEDPE output (pairs) */
  snprintf(path_buf, sizeof(path_buf), "%s.dup.bedpe", out_prefix);
  size_t reported_edges =
      write_bedpe_output(path_buf, edges, n_edges, coords, &uf, comp_size,
                         hub_label, min_copy, max_copy);

  /* Phase 5: Write BED output */
  snprintf(path_buf, sizeof(path_buf), "%s.dup.bed", out_prefix);
  BedOutputStats bed_stats = write_bed_output(path_buf, dup_regions, n_merged,
                                              hub_label, min_copy, max_copy);

  clock_gettime(CLOCK_MONOTONIC, &t_end);
  double elapsed =
      (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

  print_summary(num_files, num_sketches, n_edges, reported_edges,
                bed_stats.n_families, bed_stats.n_output, elapsed, out_prefix);

cleanup:
  free(edges);
  free(dup_regions);
  free(comp_size);
  if (hub_label) {
    for (size_t i = 0; i < num_sketches; i++)
      free(hub_label[i]);
    free(hub_label);
  }
  if (uf.parent)
    uf_free(&uf);
  for (size_t i = 0; i < num_sketches; i++) {
    reverb_sketch_free(&sketches[i]);
    free(coords[i].chrom);
  }
  free(sketches);
  free(coords);

  return ret;
}

// ==============================================================
// MAIN
// ==============================================================

int main(int argc, char **argv) {
  if (argc < 2) {
    print_usage();
    return 1;
  }

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--help") == 0) {
      print_usage();
      return 0;
    }
  }

  return cmd_dup(argc, argv);
}
