#include <math.h>
#include <stdio.h>
#include <zlib.h>

#include "klib/ketopt.h"
#include "klib/khash.h"
#include "klib/kseq.h"
#include "reverb.h"

void kt_for(int n_threads, void (*func)(void *, long, int), void *data, long n);

#define MIX_CONST1 0xff51afd7ed558ccdULL
#define MIX_CONST2 0xc4ceb9fe1a85ec53ULL

/* Generic dynamic-array push: grows `arr` by doubling `cap` as needed. */
#define DA_PUSH(arr, n, cap, val)                                              \
  do {                                                                         \
    /* If exceed capacity */                                                   \
    if ((n) >= (cap)) {                                                        \
      (cap) = (cap) ? (cap) * 2 : 1024;                                        \
      (arr) = realloc((arr), (cap) * sizeof(*(arr)));                          \
    }                                                                          \
    (arr)[(n)++] = (val);                                                      \
  } while (0)

#define CMP(a, b) (((a) > (b)) - ((a) < (b)))
#define SWAP(type, a, b)                                                       \
  do {                                                                         \
    type _t = (a);                                                             \
    (a) = (b);                                                                 \
    (b) = _t;                                                                  \
  } while (0)

#define ABS_DIFF(a, b) ((a) > (b) ? (a) - (b) : (b) - (a))

/* Reader initiation */
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
  /* reverse bits (rbits) asm command for ARM chips */
  uint64_t r;
  __asm__("rbit %0, %1" : "=r"(r) : "r"(n));
  return r;
#else
  /* manually reverse bits for x86_64 chips */
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

void init_reverb(Reverb *r, size_t hash_window) {
  size_t k = hash_window < 42 ? hash_window : 42; /* 3*42=126 bits ≤ 128 */
  uint32_t kmer_bits = 3 * (uint32_t)k;

  __uint128_t remover_mask =
      (kmer_bits > 3) ? (((__uint128_t)1 << (kmer_bits - 3)) - 1) : 0;
  /* remover mask to forget previous base */

  r->hash_window = k;
  r->remover_mask = remover_mask;
  r->kmer_bits = kmer_bits;
  r->rc_shift =
      (kmer_bits > 0) ? (128 - kmer_bits) : 128; /* reverse_complement shift */
}

// ==============================================================
// HASH POOL
// ==============================================================

typedef struct {
  size_t size;             /* ?????????? */
  size_t cap;              /* ?????????? */
  uint64_t hash_threshold; /* FracMinHash threshold */
  uint64_t *hashes;
} HashPool;

/* Wrapper for macro to use in `qsort`*/
static int compare_uint64(const void *a, const void *b) {
  return CMP(*(const uint64_t *)a, *(const uint64_t *)b);
}

static void init_hash_pool(HashPool *pool, uint64_t threshold) {
  *pool = (HashPool){.hash_threshold = threshold};
}

static void insert_hash_pool(HashPool *pool, uint64_t h) {
  if (h >= pool->hash_threshold)
    return;
  DA_PUSH(pool->hashes, pool->size, pool->cap, h);
}

static void finalize_hash_pool(HashPool *pool, uint64_t **out_hashes,
                               size_t *out_size) {
  size_t n = pool->size;
  if (n) {
    qsort(pool->hashes, n, sizeof(uint64_t), compare_uint64);
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

// ==============================================================
// SKETCH EXTRACTION
// ==============================================================

/* Extract reverb hash and insert in HashPool */
/* Hot spot code */
__attribute__((hot)) static void extract_hash(const Reverb *r, HashPool *pool,
                                              const uint8_t *seq, size_t len) {
  size_t K = r->hash_window;
  __uint128_t fwd = 0;
  size_t valid = 0;

  for (size_t idx = 0; idx < len; idx++) {
    /* Convert to numerical values */
    int8_t lv = BASE_LOOKUP[seq[idx]];
    if (lv < 0) {
      fwd = 0;
      valid = 0;
      continue;
    }

    /* Concat for fwd hash */
    fwd = ((fwd & r->remover_mask) << 3) | (uint8_t)lv;
    if (valid < K)
      valid++;
    if (valid < K)
      continue;

    /* Reverse bits for reverse complement */
    __uint128_t rev = reverse_bits128(fwd) >> r->rc_shift;
    /* Min operation to canonicalize */
    __uint128_t canon = fwd < rev ? fwd : rev;
    /* Mix hash to avoid collision */
    uint64_t h = mix_hash(canon, r->hash_seed);

    /* FracMinHash */
    if (h < pool->hash_threshold)
      insert_hash_pool(pool, h);
  }
}

// ==============================================================
// DISTANCE CALCULATION
// ==============================================================

/* Calculate distance between two sketch sets */
ReverbDistResult calculate_reverb_dist(const ReverbSketch *ref,
                                       const ReverbSketch *query,
                                       uint32_t kmer_size) {
  ReverbDistResult res = {0.0, 1.0, 0};

  size_t shared = 0, i = 0, j = 0;

  /* Set opperations */
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

  /* Calculate distance with containment method */
  res.shared_hashes = shared;
  if (ref->sketch_size > 0 && query->sketch_size > 0) {
    res.containment =
        0.5 * shared * (1.0 / ref->sketch_size + 1.0 / query->sketch_size);
    res.distance = 1.0 - pow(res.containment, 1.0 / (double)kmer_size);
  }
  return res;
}

// ==============================================================
// UNION-FIND
// Union-find algorithm to determine two nodes are in same set or not
// ==============================================================

void init_unionfind(UnionFind *uf, size_t n) {
  uf->n = n;
  uf->parent = (uint32_t *)malloc(n * sizeof(uint32_t));
  uf->rank = (uint32_t *)calloc(n, sizeof(uint32_t));
  for (size_t i = 0; i < n; i++)
    uf->parent[i] = (uint32_t)i;
}

uint32_t find_unionfind(UnionFind *uf, uint32_t x) {
  while (uf->parent[x] != x) {
    uf->parent[x] = uf->parent[uf->parent[x]]; /* path splitting */
    x = uf->parent[x];
  }
  return x;
}

void union_unionfind(UnionFind *uf, uint32_t a, uint32_t b) {
  a = find_unionfind(uf, a);
  b = find_unionfind(uf, b);
  if (a == b)
    return;
  if (uf->rank[a] < uf->rank[b])
    SWAP(uint32_t, a, b);
  uf->parent[b] = a;
  if (uf->rank[a] == uf->rank[b])
    uf->rank[a]++;
}

void free_unionfind(UnionFind *uf) {
  free(uf->parent);
  free(uf->rank);
}

// ==============================================================
// PARAMETERS
// ==============================================================

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
// WINDOW STREAMING
// ==============================================================

typedef struct {
  uint32_t seq_id;
  size_t start;
  size_t end;
  size_t sketch_offset; /* byte offset in SketchStore */
  uint32_t sketch_size; /* number of hashes */
} WindowCoord;

/* Inverted hash index entry: maps a hash value to its source window */
typedef struct {
  uint64_t hash;
  uint32_t window_id;
} HashWindowEntry;

// ==============================================================
// INVERTED HASH INDEX
// ==============================================================

static int compare_hash_entry(const void *a, const void *b) {
  const HashWindowEntry *ea = (const HashWindowEntry *)a,
                        *eb = (const HashWindowEntry *)b;
  return ea->hash != eb->hash ? CMP(ea->hash, eb->hash)
                              : CMP(ea->window_id, eb->window_id);
}

typedef struct {
  const uint64_t *all_hashes;
  WindowCoord *coords;
  HashWindowEntry *entries;
  size_t total_entries;
  size_t n_windows;
  double max_dist;
  size_t window_size;
  uint32_t kmer_size;
  ReverbDupEdge **t_edges;
  size_t *t_n_edges;
  size_t *t_cap_edges;
} EdgeWorkerData;

KHASH_MAP_INIT_INT(u32, uint16_t)

/* Helper: compute distance between two windows using in-memory hashes */
static ReverbDistResult calculate_window_dist(const uint64_t *all_hashes,
                                              const WindowCoord *wa,
                                              const WindowCoord *wb,
                                              uint32_t kmer_size) {
  ReverbSketch sa = {.sketch_size = wa->sketch_size,
                     .hashes = (uint64_t *)(all_hashes + wa->sketch_offset)};
  ReverbSketch sb = {.sketch_size = wb->sketch_size,
                     .hashes = (uint64_t *)(all_hashes + wb->sketch_offset)};
  return calculate_reverb_dist(&sa, &sb, kmer_size);
}

static void process_edge(void *data, long i, int tid) {
  EdgeWorkerData *w = (EdgeWorkerData *)data;
  uint32_t a = (uint32_t)i;

  const uint64_t *a_hashes = w->all_hashes + w->coords[a].sketch_offset;

  khash_t(u32) *counts = kh_init(u32);
  int ret;

  for (uint32_t k = 0; k < w->coords[a].sketch_size; k++) {
    uint64_t hash = a_hashes[k];

    // Binary search in entries
    size_t left = 0, right = w->total_entries;
    while (left < right) {
      size_t mid = left + (right - left) / 2;
      if (w->entries[mid].hash < hash)
        left = mid + 1;
      else
        right = mid;
    }

    if (left < w->total_entries && w->entries[left].hash == hash) {
      size_t run_end = left + 1;
      while (run_end < w->total_entries && w->entries[run_end].hash == hash)
        run_end++;
      size_t run_len = run_end - left;

      if (run_len >= 2 && run_len <= 100) {
        for (size_t idx = left; idx < run_end; idx++) {
          uint32_t b = w->entries[idx].window_id;
          if (b > a) {
            khint_t ki = kh_put(u32, counts, b, &ret);
            if (ret) /* new key */
              kh_val(counts, ki) = 1;
            else
              kh_val(counts, ki)++;
          }
        }
      }
    }
  }

  khint_t ki;
  for (ki = kh_begin(counts); ki != kh_end(counts); ++ki) {
    if (!kh_exist(counts, ki))
      continue;
    uint32_t b = kh_key(counts, ki);
    uint16_t cnt = kh_val(counts, ki);
    if (cnt >= 2) {
      if (w->coords[a].seq_id != w->coords[b].seq_id ||
          ABS_DIFF(w->coords[a].start, w->coords[b].start) >= w->window_size) {
        ReverbDistResult d = calculate_window_dist(w->all_hashes, &w->coords[a],
                                                   &w->coords[b], w->kmer_size);
        if (d.distance < w->max_dist) {
          DA_PUSH(w->t_edges[tid], w->t_n_edges[tid], w->t_cap_edges[tid],
                  ((ReverbDupEdge){a, b, d.distance}));
        }
      }
    }
  }

  kh_destroy(u32, counts);
}

/* Build inverted hash index and find candidate pairs.
 * Uses fully in-memory arrays for straightforward computation.
 * Returns edges where distance < max_dist, excluding adjacent windows
 * on the same chromosome. */
static size_t build_candidate_edges(const uint64_t *all_hashes,
                                    WindowCoord *coords, size_t n_windows,
                                    double max_dist, size_t window_size,
                                    int n_threads, uint32_t kmer_size,
                                    ReverbDupEdge **out_edges) {
  /* 1. Flatten all (hash, window_id) entries */
  size_t total_entries = 0;
  for (size_t i = 0; i < n_windows; i++)
    total_entries += coords[i].sketch_size;

  if (total_entries == 0) {
    *out_edges = NULL;
    return 0;
  }

  /* Allocate entries in memory directly */
  size_t entries_bytes = total_entries * sizeof(HashWindowEntry);
  HashWindowEntry *entries = malloc(entries_bytes);

  size_t idx = 0;
  for (size_t i = 0; i < n_windows; i++) {
    const uint64_t *hashes = all_hashes + coords[i].sketch_offset;
    for (uint32_t j = 0; j < coords[i].sketch_size; j++) {
      entries[idx++] =
          (HashWindowEntry){.hash = hashes[j], .window_id = (uint32_t)i};
    }
  }

  /* 2. Sort by hash value */
  qsort(entries, total_entries, sizeof(HashWindowEntry), compare_hash_entry);

  /* 3. Parallel distance computation using query-driven search */
  EdgeWorkerData w;
  w.all_hashes = all_hashes;
  w.coords = coords;
  w.entries = entries;
  w.total_entries = total_entries;
  w.n_windows = n_windows;
  w.max_dist = max_dist;
  w.window_size = window_size;
  w.kmer_size = kmer_size;
  w.t_edges = calloc(n_threads, sizeof(ReverbDupEdge *));
  w.t_n_edges = calloc(n_threads, sizeof(size_t));
  w.t_cap_edges = calloc(n_threads, sizeof(size_t));

  kt_for(n_threads, process_edge, &w, n_windows);
  free(entries);

  size_t n_edges = 0;
  for (int t = 0; t < n_threads; t++)
    n_edges += w.t_n_edges[t];

  ReverbDupEdge *edges = malloc(n_edges * sizeof(ReverbDupEdge));
  size_t offset = 0;
  for (int t = 0; t < n_threads; t++) {
    if (w.t_n_edges[t] > 0) {
      memcpy(edges + offset, w.t_edges[t],
             w.t_n_edges[t] * sizeof(ReverbDupEdge));
      offset += w.t_n_edges[t];
      free(w.t_edges[t]);
    }
  }
  free(w.t_edges);
  free(w.t_n_edges);
  free(w.t_cap_edges);

  if (n_edges == 0) {
    free(edges);
    edges = NULL;
  }

  *out_edges = edges;
  return n_edges;
}

// ==============================================================
// DUP: SEGMENT MERGE
// ==============================================================

static int compare_dup_region(const void *a, const void *b) {
  const ReverbDupRegion *ra = (const ReverbDupRegion *)a,
                        *rb = (const ReverbDupRegion *)b;
  int c = strcmp(ra->chrom, rb->chrom);
  return c ? c : CMP(ra->start, rb->start);
}

/* Merge adjacent/overlapping regions in the same SD family.
 * Returns the new count of merged regions. */
static size_t merge_dup_regions(ReverbDupRegion *regions, size_t n) {
  if (n <= 1)
    return n;

  qsort(regions, n, sizeof(ReverbDupRegion), compare_dup_region);

  size_t out = 0;
  for (size_t i = 1; i < n; i++) {
    if (strcmp(regions[i].cluster_id, regions[out].cluster_id) == 0 &&
        strcmp(regions[i].chrom, regions[out].chrom) == 0 &&
        regions[i].start <= regions[out].end) {
      if (regions[i].end > regions[out].end)
        regions[out].end = regions[i].end;
      free(regions[i].cluster_id);
    } else {
      out++;
      if (out != i)
        regions[out] = regions[i];
    }
  }
  return out + 1;
}

// ==============================================================
// CMD: DUP
// ==============================================================

static void get_basename(const char *filename, char *basename, size_t size) {
  const char *slash = strrchr(filename, '/');
  const char *base = slash ? slash + 1 : filename;
  strncpy(basename, base, size - 1);
  basename[size - 1] = '\0';
  char *dot = strchr(basename, '.');
  if (dot)
    *dot = '\0';
}

typedef struct {
  char *genome;
  char *seq;
} GenomeSeqLen;

typedef struct {
  const char *filename;
  char bname[256];
  const Reverb *r;
  uint64_t scale;
  size_t window_size;
  size_t step_size;
  size_t min_bases;

  uint64_t *all_hashes;
  size_t num_all_hashes;
  size_t cap_all_hashes;

  WindowCoord *coords;
  size_t num_sketches;
  size_t cap_sketches;

  GenomeSeqLen *seq_lens;
  size_t num_seqs;
  size_t cap_seqs;
} StreamWorkerData;

static void stream_pangenome_worker(void *data, long i, int tid) {
  (void)tid;
  StreamWorkerData *w = &((StreamWorkerData *)data)[i];

  gzFile fp = gzopen(w->filename, "r");
  if (!fp)
    return;
  kseq_t *ks = kseq_init(fp);
  if (!ks) {
    gzclose(fp);
    return;
  }

  while (kseq_read(ks) >= 0) {
    size_t len = ks->seq.l;

    if (w->num_seqs >= w->cap_seqs) {
      w->cap_seqs = w->cap_seqs == 0 ? 16 : w->cap_seqs * 2;
      w->seq_lens = realloc(w->seq_lens, w->cap_seqs * sizeof(GenomeSeqLen));
    }
    w->seq_lens[w->num_seqs].genome = strdup(w->bname);
    w->seq_lens[w->num_seqs].seq = strdup(ks->name.s);
    w->num_seqs++;

    for (size_t idx = 0; idx + w->window_size <= len; idx += w->step_size) {
      size_t valid_bases = 0;
      for (size_t j = 0; j < w->window_size; j++) {
        if (BASE_LOOKUP[(uint8_t)ks->seq.s[idx + j]] >= 0)
          valid_bases++;
      }
      if (valid_bases < w->min_bases)
        continue;

      if (w->num_sketches >= w->cap_sketches) {
        w->cap_sketches = w->cap_sketches == 0 ? 256 : w->cap_sketches * 2;
        w->coords = realloc(w->coords, w->cap_sketches * sizeof(WindowCoord));
      }

      WindowCoord *wc = &w->coords[w->num_sketches];

      wc->seq_id = (uint32_t)(w->num_seqs - 1);
      wc->start = idx;
      wc->end = idx + w->window_size;

      /* Extract sketch, write to disk, free immediately */
      HashPool pool;
      init_hash_pool(&pool, UINT64_MAX / w->scale);
      extract_hash(w->r, &pool, (const uint8_t *)ks->seq.s + idx,
                   w->window_size);

      uint64_t *hashes = NULL;
      size_t sketch_size = 0;
      finalize_hash_pool(&pool, &hashes, &sketch_size);

      if (sketch_size > 0) {
        size_t hash_idx = w->num_all_hashes;
        while (w->num_all_hashes + sketch_size > w->cap_all_hashes) {
          w->cap_all_hashes =
              (w->cap_all_hashes == 0) ? 1048576 : (w->cap_all_hashes * 2);
          w->all_hashes =
              realloc(w->all_hashes, (w->cap_all_hashes) * sizeof(uint64_t));
        }
        memcpy(w->all_hashes + hash_idx, hashes,
               sketch_size * sizeof(uint64_t));
        w->num_all_hashes += sketch_size;

        wc->sketch_offset = hash_idx;
        wc->sketch_size = (uint32_t)sketch_size;
        free(hashes);
        w->num_sketches++;
      } else {
        if (hashes)
          free(hashes);
      }
    }
  }
  kseq_destroy(ks);
  gzclose(fp);
}
static void extract_flankings(char **files, int num_files, const Reverb *r,
                              uint64_t scale, ReverbDupRegion *regions,
                              size_t n_regions, size_t flank_size) {
  for (int f = 0; f < num_files; f++) {
    char bname[256];
    get_basename(files[f], bname, sizeof(bname));

    gzFile fp = gzopen(files[f], "r");
    if (!fp)
      continue;
    kseq_t *ks = kseq_init(fp);
    if (!ks) {
      gzclose(fp);
      continue;
    }

    int chr_idx = 1;
    while (kseq_read(ks) >= 0) {
      char chr_name[512];
      snprintf(chr_name, sizeof(chr_name), "%s-Chr%d", bname, chr_idx++);

      for (size_t i = 0; i < n_regions; i++) {
        if (strcmp(regions[i].chrom, chr_name) == 0) {
          size_t start = regions[i].start;
          size_t end = regions[i].end;
          size_t left_start = start > flank_size ? start - flank_size : 0;
          size_t right_end =
              end + flank_size > ks->seq.l ? ks->seq.l : end + flank_size;

          size_t left_len = start - left_start;
          size_t right_len = right_end - end;

          uint8_t *flank_seq = malloc(left_len + right_len);
          if (left_len > 0)
            memcpy(flank_seq, ks->seq.s + left_start, left_len);
          if (right_len > 0)
            memcpy(flank_seq + left_len, ks->seq.s + end, right_len);

          // flank_sketch fields removed for memory optimization

          HashPool pool;
          init_hash_pool(&pool, UINT64_MAX / scale);
          extract_hash(r, &pool, flank_seq, left_len + right_len);
          finalize_hash_pool(&pool, &regions[i].flank_sketch.hashes,
                             &regions[i].flank_sketch.sketch_size);

          free(flank_seq);
        }
      }
    }
    kseq_destroy(ks);
    gzclose(fp);
  }
}

typedef struct {
  uint32_t i, j;
} SubclusterPair;

typedef struct {
  ReverbDupRegion *regions;
  double max_dist;
  size_t n_merged;
  SubclusterPair **t_pairs;
  size_t *t_n_pairs;
  size_t *t_cap_pairs;
  uint32_t kmer_size;
} SubclusterData;

static void process_subcluster(void *data, long i, int tid) {
  SubclusterData *w = (SubclusterData *)data;
  if (w->regions[i].flank_sketch.sketch_size == 0)
    return;
  for (size_t j = i + 1; j < w->n_merged; j++) {
    if (strcmp(w->regions[i].cluster_id, w->regions[j].cluster_id) != 0)
      continue; // must be same cluster
    if (w->regions[j].flank_sketch.sketch_size == 0)
      continue;

    ReverbDistResult d = calculate_reverb_dist(
        &w->regions[i].flank_sketch, &w->regions[j].flank_sketch, w->kmer_size);
    if (d.distance < w->max_dist) {
      DA_PUSH(w->t_pairs[tid], w->t_n_pairs[tid], w->t_cap_pairs[tid],
              ((SubclusterPair){(uint32_t)i, (uint32_t)j}));
    }
  }
}

static void perform_subclustering(ReverbDupRegion *regions, size_t n_merged,
                                  double max_dist, int n_threads,
                                  uint32_t kmer_size) {
  UnionFind sub_uf;
  init_unionfind(&sub_uf, n_merged);

  SubclusterData w;
  w.regions = regions;
  w.max_dist = max_dist;
  w.n_merged = n_merged;
  w.kmer_size = kmer_size;
  w.t_pairs = calloc(n_threads, sizeof(SubclusterPair *));
  w.t_n_pairs = calloc(n_threads, sizeof(size_t));
  w.t_cap_pairs = calloc(n_threads, sizeof(size_t));

  kt_for(n_threads, process_subcluster, &w, n_merged);

  for (int t = 0; t < n_threads; t++) {
    for (size_t k = 0; k < w.t_n_pairs[t]; k++) {
      union_unionfind(&sub_uf, w.t_pairs[t][k].i, w.t_pairs[t][k].j);
    }
    if (w.t_pairs[t])
      free(w.t_pairs[t]);
  }
  free(w.t_pairs);
  free(w.t_n_pairs);
  free(w.t_cap_pairs);

  // Map union-find parents to sub_cluster_id
  uint32_t *mapping = calloc(n_merged, sizeof(uint32_t));
  uint32_t current_id = 1;

  for (size_t i = 0; i < n_merged; i++) {
    uint32_t p = find_unionfind(&sub_uf, i);
    if (mapping[p] == 0) {
      mapping[p] = current_id++;
    }
    regions[i].subcluster_id = mapping[p];
  }
  free(mapping);
  free_unionfind(&sub_uf);
}

int run_pangenome(int num_files, char **files, size_t flank_size,
                  const Reverb *r, uint64_t scale, size_t window_size,
                  size_t step_size, size_t min_bases, double max_dist,
                  int min_copy, int max_copy, const char *out_prefix,
                  int n_threads) {

  uint64_t *all_hashes = NULL;

  WindowCoord *coords = NULL;
  size_t num_sketches = 0;
  ReverbDupEdge *edges = NULL;

  char path_buf[PATH_MAX];
  snprintf(path_buf, sizeof(path_buf), "%s.window.bed", out_prefix);
  FILE *bed_fp = fopen(path_buf, "w");

  fprintf(stderr, "[reverb] Extracting windows across pangenome...\\n");
  GenomeSeqLen *seq_lens = NULL;

  StreamWorkerData *workers = calloc(num_files, sizeof(StreamWorkerData));
  for (int i = 0; i < num_files; i++) {
    workers[i].filename = files[i];
    get_basename(files[i], workers[i].bname, sizeof(workers[i].bname));
    workers[i].r = r;
    workers[i].scale = scale;
    workers[i].window_size = window_size;
    workers[i].step_size = step_size;
    workers[i].min_bases = min_bases;
  }

  kt_for(n_threads, stream_pangenome_worker, workers, num_files);

  size_t total_hashes = 0, total_sketches = 0, total_seqs = 0;
  for (int i = 0; i < num_files; i++) {
    total_hashes += workers[i].num_all_hashes;
    total_sketches += workers[i].num_sketches;
    total_seqs += workers[i].num_seqs;
  }

  all_hashes = malloc(total_hashes * sizeof(uint64_t));
  coords = malloc(total_sketches * sizeof(WindowCoord));
  seq_lens = malloc(total_seqs * sizeof(GenomeSeqLen));

  size_t g_hash_offset = 0;
  size_t g_sketch_offset = 0;
  size_t g_seq_offset = 0;

  for (int i = 0; i < num_files; i++) {
    StreamWorkerData *w = &workers[i];

    // Copy all_hashes
    if (w->num_all_hashes > 0) {
      memcpy(all_hashes + g_hash_offset, w->all_hashes,
             w->num_all_hashes * sizeof(uint64_t));
    }

    // Copy seq_lens
    if (w->num_seqs > 0) {
      memcpy(seq_lens + g_seq_offset, w->seq_lens,
             w->num_seqs * sizeof(GenomeSeqLen));
    }

    // Copy coords and write to bed
    for (size_t j = 0; j < w->num_sketches; j++) {
      WindowCoord c = w->coords[j];
      c.sketch_offset += g_hash_offset;
      c.seq_id += g_seq_offset;
      coords[g_sketch_offset + j] = c;

      if (bed_fp) {
        char chr_name[512];
        snprintf(chr_name, sizeof(chr_name), "%s-%s", seq_lens[c.seq_id].genome,
                 seq_lens[c.seq_id].seq);
        fprintf(bed_fp, "%s\t%zu\t%zu\t%s_%zu_%zu\n", chr_name, c.start, c.end,
                chr_name, c.start, c.end);
      }
    }

    g_hash_offset += w->num_all_hashes;
    g_sketch_offset += w->num_sketches;
    g_seq_offset += w->num_seqs;

    free(w->all_hashes);
    free(w->coords);
    free(w->seq_lens);
  }
  free(workers);

  num_sketches = total_sketches;
  fclose(bed_fp);

  size_t n_edges =
      build_candidate_edges(all_hashes, coords, num_sketches, max_dist,
                            window_size, n_threads, r->hash_window, &edges);

  uint32_t *genome_id = calloc(num_sketches, sizeof(uint32_t));
  for (size_t i = 0; i < num_sketches; i++) {
    for (int f = 0; f < num_files; f++) {
      char bname[256];
      get_basename(files[f], bname, sizeof(bname));
      if (strcmp(seq_lens[coords[i].seq_id].genome, bname) == 0) {
        genome_id[i] = f;
        break;
      }
    }
  }

  // Single Global UnionFind
  UnionFind uf;
  init_unionfind(&uf, num_sketches);
  for (size_t i = 0; i < n_edges; i++) {
    union_unionfind(&uf, edges[i].win_a, edges[i].win_b);
  }

  // Count instances per genome per family
  uint32_t *max_intra_copy = calloc(num_sketches, sizeof(uint32_t));
  uint32_t *counts = calloc(num_sketches * num_files, sizeof(uint32_t));
  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t fam = find_unionfind(&uf, (uint32_t)i);
    uint32_t g_id = genome_id[i];
    counts[fam * num_files + g_id]++;
    if (counts[fam * num_files + g_id] > max_intra_copy[fam]) {
      max_intra_copy[fam] = counts[fam * num_files + g_id];
    }
  }
  free(counts);

  uint8_t *final_is_sd = calloc(num_sketches, sizeof(uint8_t));
  char **hub_label = calloc(num_sketches, sizeof(char *));
  uint32_t next_cluster_id = 1;

  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t fam = find_unionfind(&uf, (uint32_t)i);
    if (max_intra_copy[fam] >= (uint32_t)min_copy &&
        (max_copy <= 0 || max_intra_copy[fam] <= (uint32_t)max_copy)) {
      final_is_sd[fam] = 1;
      if (!hub_label[fam]) {
        char buf[64];
        snprintf(buf, sizeof(buf), "%u", next_cluster_id++);
        hub_label[fam] = strdup(buf);
      }
    }
  }
  free(genome_id);
  free(max_intra_copy);

  uint32_t *comp_size = calloc(num_sketches, sizeof(uint32_t));
  for (size_t i = 0; i < num_sketches; i++)
    comp_size[find_unionfind(&uf, (uint32_t)i)]++;

  size_t n_dup_regions = 0, cap_dup_regions = 0;
  ReverbDupRegion *dup_regions = NULL;
  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t fam = find_unionfind(&uf, (uint32_t)i);
    if (!final_is_sd[fam])
      continue;

    const char *label = hub_label[fam] ? hub_label[fam] : "unknown";

    char chrom_name[512];
    snprintf(chrom_name, sizeof(chrom_name), "%s-%s",
             seq_lens[coords[i].seq_id].genome, seq_lens[coords[i].seq_id].seq);

    DA_PUSH(dup_regions, n_dup_regions, cap_dup_regions,
            ((ReverbDupRegion){.chrom = strdup(chrom_name),
                               .start = coords[i].start,
                               .end = coords[i].end,
                               .cluster_id = strdup(label),
                               .copy_count = comp_size[fam],
                               .subcluster_id = 0,
                               .flank_sketch = {0}}));
  }

  size_t n_merged = merge_dup_regions(dup_regions, n_dup_regions);
  free(final_is_sd);

  // Output dup.bedpe
  snprintf(path_buf, sizeof(path_buf), "%s.dup.bedpe", out_prefix);
  FILE *out_bedpe = fopen(path_buf, "w");
  fprintf(
      out_bedpe,
      "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tcluster_id\tdistance\n");
  for (size_t i = 0; i < n_edges; i++) {
    uint32_t a = edges[i].win_a;
    uint32_t b = edges[i].win_b;
    uint32_t fam = find_unionfind(&uf, a);
    if (!hub_label[fam])
      continue;

    char chrom_a[512], chrom_b[512];
    snprintf(chrom_a, sizeof(chrom_a), "%s-%s",
             seq_lens[coords[a].seq_id].genome, seq_lens[coords[a].seq_id].seq);
    snprintf(chrom_b, sizeof(chrom_b), "%s-%s",
             seq_lens[coords[b].seq_id].genome, seq_lens[coords[b].seq_id].seq);

    fprintf(out_bedpe, "%s\t%zu\t%zu\t%s\t%zu\t%zu\t%s\t%.6f\n", chrom_a,
            coords[a].start, coords[a].end, chrom_b, coords[b].start,
            coords[b].end, hub_label[fam], edges[i].distance);
  }
  fclose(out_bedpe);

  fprintf(stderr,
          "[INFO] Extracting flanking sequences for sub-clustering...\n");
  extract_flankings(files, num_files, r, scale, dup_regions, n_merged,
                    flank_size == 0 ? window_size / 5 : flank_size);

  fprintf(stderr, "[INFO] Sub-clustering based on flanking similarities...\n");
  perform_subclustering(dup_regions, n_merged, max_dist, n_threads,
                        r->hash_window);

  // Output dup.bed
  snprintf(path_buf, sizeof(path_buf), "%s.dup.bed", out_prefix);
  FILE *out_bed = fopen(path_buf, "w");

  fprintf(out_bed,
          "#chrom\tstart\tend\tcluster_id\tsubcluster_id\tcopy_count\n");

  uint32_t max_subcluster = 0;
  for (size_t i = 0; i < n_merged; i++) {
    fprintf(out_bed, "%s\t%zu\t%zu\t%s\t%u\t%u\n", dup_regions[i].chrom,
            dup_regions[i].start, dup_regions[i].end, dup_regions[i].cluster_id,
            dup_regions[i].subcluster_id, dup_regions[i].copy_count);

    if (dup_regions[i].subcluster_id > max_subcluster)
      max_subcluster = dup_regions[i].subcluster_id;
  }
  fclose(out_bed);

  for (size_t i = 0; i < num_sketches; i++) {
    if (hub_label[i])
      free(hub_label[i]);
  }
  free(hub_label);

  free(all_hashes);
  return 0;
}

int run_dup(int argc, char **argv) {

  // ==============================================================
  // CLI defaults
  // ==============================================================

  uint32_t def_kmer_size = 21;
  uint64_t def_scale = 10;
  uint64_t def_hash_seed = 42;
  size_t window_size = 5000;
  size_t step_size = 0; /* 0 = auto (window/2) */
  size_t min_bases = 1000;
  double max_dist = 0.03;
  int min_copy = 3;
  int max_copy = 30;
  const char *out_prefix = "reverb";
  size_t flank_size = 0; /* 0 = auto (window/5) */
  int n_threads = 8;

  ketopt_t opt = KETOPT_INIT;
  int c;
  while ((c = ketopt(&opt, argc, argv, 1, "k:s:e:w:t:b:d:m:M:o:p:f:h", 0)) >=
         0) {
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
    else if (c == 'f')
      flank_size = (size_t)strtoull(opt.arg, NULL, 10);
    else
      return 1;
  }

  if (step_size == 0)
    step_size = window_size / 2;

  if (opt.ind == argc) {
    fprintf(stderr, "[ERROR] Input FASTA files are required.\n");
    return 1;
  }

  int num_files = argc - opt.ind;
  char **files = &argv[opt.ind];

  Reverb r;
  init_reverb(&r, def_kmer_size);
  r.hash_seed = def_hash_seed;
  return run_pangenome(num_files, files, flank_size, &r, def_scale, window_size,
                       step_size, min_bases, max_dist, min_copy, max_copy,
                       out_prefix, n_threads);
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

  return run_dup(argc, argv);
}
