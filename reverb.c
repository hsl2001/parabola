#include <math.h>
#include <stdio.h>
#include <time.h>
#include <zlib.h>

#include "klib/ketopt.h"
#include "klib/kseq.h"
#include "reverb.h"

void kt_for(int n_threads, void (*func)(void *, long, int), void *data, long n);

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

#define CMP(a, b) (((a) > (b)) - ((a) < (b)))
#define SWAP(type, a, b)                                                       \
  do {                                                                         \
    type _t = (a);                                                             \
    (a) = (b);                                                                 \
    (b) = _t;                                                                  \
  } while (0)

#define ABS_DIFF(a, b) ((a) > (b) ? (a) - (b) : (b) - (a))

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
  return CMP(*(const uint64_t *)a, *(const uint64_t *)b);
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
    res.containment =
        0.5 * shared * (1.0 / ref->sketch_size + 1.0 / query->sketch_size);
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
  if (uf->rank[a] < uf->rank[b])
    SWAP(uint32_t, a, b);
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

// ==============================================================
// DUP: INVERTED HASH INDEX
// ==============================================================

static int cmp_hash_window_entry(const void *a, const void *b) {
  const HashWindowEntry *ea = (const HashWindowEntry *)a,
                        *eb = (const HashWindowEntry *)b;
  return ea->hash != eb->hash ? CMP(ea->hash, eb->hash)
                              : CMP(ea->window_id, eb->window_id);
}

typedef struct {
  ReverbSketch *sketches;
  WindowCoord *coords;
  HashWindowEntry *entries;
  size_t total_entries;
  size_t n_windows;
  double max_dist;
  size_t window_size;
  ReverbDupEdge **t_edges;
  size_t *t_n_edges;
  size_t *t_cap_edges;
} EdgeWorkerData;

static void edge_worker(void *data, long i, int tid) {
  EdgeWorkerData *w = (EdgeWorkerData *)data;
  uint32_t a = (uint32_t)i;

  uint16_t *counts = calloc(w->n_windows, sizeof(uint16_t));
  uint32_t *touched = malloc(w->n_windows * sizeof(uint32_t));
  size_t n_touched = 0;

  for (size_t k = 0; k < w->sketches[a].sketch_size; k++) {
    uint64_t hash = w->sketches[a].hashes[k];

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
            if (counts[b] == 0)
              touched[n_touched++] = b;
            counts[b]++;
          }
        }
      }
    }
  }

  for (size_t t = 0; t < n_touched; t++) {
    uint32_t b = touched[t];
    if (counts[b] >= 2) {
      if (strcmp(w->coords[a].chrom, w->coords[b].chrom) != 0 ||
          ABS_DIFF(w->coords[a].start, w->coords[b].start) >= w->window_size) {
        ReverbDistResult d = reverb_dist(&w->sketches[a], &w->sketches[b]);
        if (d.distance < w->max_dist) {
          DA_PUSH(w->t_edges[tid], w->t_n_edges[tid], w->t_cap_edges[tid],
                  ((ReverbDupEdge){a, b, d.distance}));
        }
      }
    }
    counts[b] = 0;
  }

  free(counts);
  free(touched);
}

/* Build inverted hash index and find candidate pairs.
 * Returns edges where distance < max_dist, excluding adjacent windows
 * on the same chromosome. */
static size_t build_candidate_edges(ReverbSketch *sketches, WindowCoord *coords,
                                    size_t n_windows, double max_dist,
                                    size_t window_size, int n_threads,
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

  /* 3. Parallel distance computation using query-driven search */
  EdgeWorkerData w;
  w.sketches = sketches;
  w.coords = coords;
  w.entries = entries;
  w.total_entries = total_entries;
  w.n_windows = n_windows;
  w.max_dist = max_dist;
  w.window_size = window_size;
  w.t_edges = calloc(n_threads, sizeof(ReverbDupEdge *));
  w.t_n_edges = calloc(n_threads, sizeof(size_t));
  w.t_cap_edges = calloc(n_threads, sizeof(size_t));

  kt_for(n_threads, edge_worker, &w, n_windows);
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

static int cmp_dup_region(const void *a, const void *b) {
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

  qsort(regions, n, sizeof(ReverbDupRegion), cmp_dup_region);

  size_t out = 0;
  for (size_t i = 1; i < n; i++) {
    if (strcmp(regions[i].cluster_id, regions[out].cluster_id) == 0 &&
        strcmp(regions[i].chrom, regions[out].chrom) == 0 &&
        regions[i].start <= regions[out].end) {
      if (regions[i].end > regions[out].end)
        regions[out].end = regions[i].end;
      regions[out].avg_distance =
          (regions[out].avg_distance + regions[i].avg_distance) / 2.0;
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

#include <ctype.h>
#include <dirent.h>
#include <stdint.h>
#include <sys/stat.h>

static int cmp_natural(const void *p1, const void *p2) {
  const char *a = *(const char **)p1;
  const char *b = *(const char **)p2;
  while (*a && *b) {
    if (isdigit((unsigned char)*a) && isdigit((unsigned char)*b)) {
      char *ea, *eb;
      long va = strtol(a, &ea, 10);
      long vb = strtol(b, &eb, 10);
      if (va != vb)
        return va - vb;
      a = ea;
      b = eb;
    } else {
      int ca = tolower((unsigned char)*a);
      int cb = tolower((unsigned char)*b);
      if (ca != cb)
        return ca - cb;
      a++;
      b++;
    }
  }
  return *a - *b;
}

static void hsv2rgb(float h, float s, float v, float *r, float *g, float *b) {
  int i = (int)(h * 6);
  float f = h * 6 - i;
  float p = v * (1 - s);
  float q = v * (1 - f * s);
  float t = v * (1 - (1 - f) * s);
  switch (i % 6) {
  case 0:
    *r = v, *g = t, *b = p;
    break;
  case 1:
    *r = q, *g = v, *b = p;
    break;
  case 2:
    *r = p, *g = v, *b = t;
    break;
  case 3:
    *r = p, *g = q, *b = v;
    break;
  case 4:
    *r = t, *g = p, *b = v;
    break;
  case 5:
    *r = v, *g = p, *b = q;
    break;
  default:
    *r = 0, *g = 0, *b = 0;
    break;
  }
}

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
  size_t length;
} GenomeSeqLen;

static int dup_stream_pangenome(const char *filename, const char *bname,
                                const Reverb *r, uint64_t scale,
                                size_t window_size, size_t step_size,
                                size_t min_bases, ReverbSketch **sketches,
                                WindowCoord **coords, size_t *num_sketches,
                                size_t *cap_sketches, FILE *bed_fp,
                                GenomeSeqLen **seq_lens, size_t *num_seqs,
                                size_t *cap_seqs) {
  gzFile fp = gzopen(filename, "r");
  if (!fp)
    return -1;
  kseq_t *ks = kseq_init(fp);
  if (!ks) {
    gzclose(fp);
    return -1;
  }

  while (kseq_read(ks) >= 0) {
    char chr_name[512];
    snprintf(chr_name, sizeof(chr_name), "%s-%s", bname, ks->name.s);
    size_t len = ks->seq.l;

    if (*num_seqs >= *cap_seqs) {
      *cap_seqs = *cap_seqs == 0 ? 16 : *cap_seqs * 2;
      *seq_lens = realloc(*seq_lens, *cap_seqs * sizeof(GenomeSeqLen));
    }
    (*seq_lens)[*num_seqs].genome = strdup(bname);
    (*seq_lens)[*num_seqs].seq = strdup(ks->name.s);
    (*seq_lens)[*num_seqs].length = len;
    (*num_seqs)++;

    for (size_t i = 0; i + window_size <= len; i += step_size) {
      size_t valid_bases = 0;
      for (size_t j = 0; j < window_size; j++) {
        if (BASE_LOOKUP[(uint8_t)ks->seq.s[i + j]] >= 0)
          valid_bases++;
      }
      if (valid_bases < min_bases)
        continue;

      if (*num_sketches >= *cap_sketches) {
        *cap_sketches = *cap_sketches == 0 ? 256 : *cap_sketches * 2;
        *sketches = realloc(*sketches, *cap_sketches * sizeof(ReverbSketch));
        *coords = realloc(*coords, *cap_sketches * sizeof(WindowCoord));
      }

      ReverbSketch *sk = &(*sketches)[*num_sketches];
      WindowCoord *wc = &(*coords)[*num_sketches];
      memset(sk, 0, sizeof(ReverbSketch));

      size_t name_len = strlen(chr_name) + 64;
      char *namebuf = malloc(name_len);
      snprintf(namebuf, name_len, "%s_%zu_%zu", chr_name, i, i + window_size);
      sk->name = namebuf;
      sk->kmer_size = r->hash_window;
      sk->hash_threshold = UINT64_MAX / scale;

      wc->chrom = strdup(chr_name);
      wc->start = i;
      wc->end = i + window_size;

      if (bed_fp)
        fprintf(bed_fp, "%s\t%zu\t%zu\t%s\n", chr_name, i, i + window_size,
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

static void do_pass2(char **files, int num_files, const Reverb *r,
                     uint64_t scale, ReverbDupRegion *regions, size_t n_regions,
                     size_t flank_size) {
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

          regions[i].flank_sketch.kmer_size = r->hash_window;
          regions[i].flank_sketch.hash_threshold = UINT64_MAX / scale;
          regions[i].flank_sketch.name = strdup(chr_name); // placeholder

          HashPool pool;
          pool_init(&pool, regions[i].flank_sketch.hash_threshold);
          extract_and_insert(r, &pool, flank_seq, left_len + right_len);
          pool_finalize(&pool, &regions[i].flank_sketch.hashes,
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
} SubclusterData;

static void subcluster_worker(void *data, long i, int tid) {
  SubclusterData *w = (SubclusterData *)data;
  if (w->regions[i].flank_sketch.sketch_size == 0)
    return;
  for (size_t j = i + 1; j < w->n_merged; j++) {
    if (strcmp(w->regions[i].cluster_id, w->regions[j].cluster_id) != 0)
      continue; // must be same cluster
    if (w->regions[j].flank_sketch.sketch_size == 0)
      continue;

    ReverbDistResult d =
        reverb_dist(&w->regions[i].flank_sketch, &w->regions[j].flank_sketch);
    if (d.distance < w->max_dist) {
      DA_PUSH(w->t_pairs[tid], w->t_n_pairs[tid], w->t_cap_pairs[tid],
              ((SubclusterPair){(uint32_t)i, (uint32_t)j}));
    }
  }
}

static void do_subclustering(ReverbDupRegion *regions, size_t n_merged,
                             double max_dist, int n_threads) {
  UnionFind sub_uf;
  uf_init(&sub_uf, n_merged);

  SubclusterData w;
  w.regions = regions;
  w.max_dist = max_dist;
  w.n_merged = n_merged;
  w.t_pairs = calloc(n_threads, sizeof(SubclusterPair *));
  w.t_n_pairs = calloc(n_threads, sizeof(size_t));
  w.t_cap_pairs = calloc(n_threads, sizeof(size_t));

  kt_for(n_threads, subcluster_worker, &w, n_merged);

  for (int t = 0; t < n_threads; t++) {
    for (size_t k = 0; k < w.t_n_pairs[t]; k++) {
      uf_union(&sub_uf, w.t_pairs[t][k].i, w.t_pairs[t][k].j);
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
    uint32_t p = uf_find(&sub_uf, i);
    if (mapping[p] == 0) {
      mapping[p] = current_id++;
    }
    regions[i].subcluster_id = mapping[p];
  }
  free(mapping);
  uf_free(&sub_uf);
}

typedef struct {
  char name[256];
  uint8_t *vector;
} GenomeVector;

static void natural_sort_genomes(GenomeVector *gv, int num_genomes,
                                 int num_subclusters) {
  // Nearest-neighbor TSP starting from 0
  if (num_genomes <= 1)
    return;

  int *visited = calloc(num_genomes, sizeof(int));
  visited[0] = 1; // start with rep

  GenomeVector *ordered = malloc(num_genomes * sizeof(GenomeVector));
  ordered[0] = gv[0];

  int last_idx = 0;
  for (int i = 1; i < num_genomes; i++) {
    int best_j = -1;
    int min_dist = 1e9;

    for (int j = 1; j < num_genomes; j++) {
      if (visited[j])
        continue;
      int dist = 0;
      for (int k = 1; k <= num_subclusters; k++) {
        if (gv[last_idx].vector[k] != gv[j].vector[k])
          dist++;
      }
      if (dist < min_dist) {
        min_dist = dist;
        best_j = j;
      }
    }

    visited[best_j] = 1;
    ordered[i] = gv[best_j];
    last_idx = best_j;
  }

  for (int i = 0; i < num_genomes; i++) {
    gv[i] = ordered[i];
  }
  free(ordered);
  free(visited);
}

static void generate_svg(const char *filename, ReverbDupRegion *regions,
                         size_t n_merged, const char *target_cluster,
                         GenomeVector *gv, int num_genomes,
                         const char **chr_suffixes, int num_chrs,
                         int *chr_lengths) {
  FILE *fp = fopen(filename, "w");
  if (!fp)
    return;

  int col_width = 150; // Fixed column width to avoid squishing
  int width = 200 + num_chrs * col_width;
  int height = 50 * num_genomes + 100;
  fprintf(
      fp,
      "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\">\n",
      width, height);
  fprintf(fp, "<style>text { font-family: Arial; font-size: 14px; }</style>\n");
  fprintf(fp, "<rect width=\"100%%\" height=\"100%%\" fill=\"white\" />\n");

  char title[256];
  get_basename(filename, title, sizeof(title));
  fprintf(fp,
          "<text x=\"10\" y=\"20\" font-weight=\"bold\" "
          "font-size=\"16px\">%s</text>\n",
          title);

  // Draw Y axis labels
  for (int i = 0; i < num_genomes; i++) {
    int y = 50 + i * 50;
    fprintf(fp, "<text x=\"10\" y=\"%d\">%s</text>\n", y + 5, gv[i].name);
  }

  // Draw chromosomes
  for (int c = 0; c < num_chrs; c++) {
    int x_offset = 200 + c * col_width;
    fprintf(fp, "<text x=\"%d\" y=\"30\" fill=\"#888888\">%s</text>\n",
            x_offset, chr_suffixes[c]);

    for (int i = 0; i < num_genomes; i++) {
      int y = 50 + i * 50;
      int max_len =
          chr_lengths[i * num_chrs + c]; // length of this chr in this genome
      if (max_len > 0) {
        // scale max_len to col_width
        double scale = (double)(col_width - 20) / max_len;
        fprintf(fp,
                "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" "
                "stroke=\"#e0e0e0\" stroke-width=\"20\" />\n",
                x_offset, y, x_offset + (int)(max_len * scale), y);
      }
    }
  }

  // Draw family members
  for (size_t r = 0; r < n_merged; r++) {
    if (target_cluster && strcmp(regions[r].cluster_id, target_cluster) != 0)
      continue;

    int y_idx = -1;
    int bname_len = 0;
    for (int i = 0; i < num_genomes; i++) {
      int len = strlen(gv[i].name);
      if (strncmp(regions[r].chrom, gv[i].name, len) == 0 &&
          regions[r].chrom[len] == '-') {
        y_idx = i;
        bname_len = len;
        break;
      }
    }
    if (y_idx == -1)
      continue;

    // Find chr column
    const char *suffix = regions[r].chrom + bname_len + 1;
    int c_idx = -1;
    for (int c = 0; c < num_chrs; c++) {
      if (strcmp(chr_suffixes[c], suffix) == 0) {
        c_idx = c;
        break;
      }
    }
    if (c_idx == -1)
      continue;

    int y = 50 + y_idx * 50;
    int x_offset = 200 + c_idx * col_width;
    int max_len = chr_lengths[y_idx * num_chrs + c_idx];
    if (max_len == 0)
      continue;
    double scale = (double)(col_width - 20) / max_len;

    double p_start = regions[r].start;
    double p_end = regions[r].end;

    int x1 = x_offset + (int)(p_start * scale);
    int x2 = x_offset + (int)(p_end * scale);
    if (x2 - x1 < 5)
      x2 = x1 + 5;

    // Clamp coordinates to column boundaries
    if (x1 < x_offset)
      x1 = x_offset;
    if (x2 > x_offset + col_width - 20)
      x2 = x_offset + col_width - 20;
    if (x1 > x2)
      continue; // In case it gets completely outside

    float r_c, g_c, b_c;
    // color based on subcluster_id
    uint32_t sc = regions[r].subcluster_id;
    hsv2rgb((sc * 0.618033988749895) - (int)(sc * 0.618033988749895), 0.8, 0.9,
            &r_c, &g_c, &b_c);

    fprintf(fp,
            "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" "
            "stroke=\"rgb(%d,%d,%d)\" stroke-width=\"20\" opacity=\"0.9\" />\n",
            x1, y, x2, y, (int)(r_c * 255), (int)(g_c * 255), (int)(b_c * 255));
  }

  fprintf(fp, "</svg>\n");
  fclose(fp);
}

static void write_bin(const char *prefix, ReverbDupRegion *regions,
                      size_t n_merged, GenomeVector *gv, int num_genomes,
                      const char **chr_suffixes, int num_chrs,
                      int *chr_lengths) {
  char path[PATH_MAX];
  snprintf(path, sizeof(path), "%s.bin", prefix);
  FILE *fp = fopen(path, "wb");
  if (!fp)
    return;

  // Write magic
  fwrite("REV2", 1, 4, fp);

  // Write genomes
  fwrite(&num_genomes, sizeof(int), 1, fp);
  for (int i = 0; i < num_genomes; i++) {
    int len = strlen(gv[i].name);
    fwrite(&len, sizeof(int), 1, fp);
    fwrite(gv[i].name, 1, len, fp);
  }

  // Write chrs
  fwrite(&num_chrs, sizeof(int), 1, fp);
  for (int i = 0; i < num_chrs; i++) {
    int len = strlen(chr_suffixes[i]);
    fwrite(&len, sizeof(int), 1, fp);
    fwrite(chr_suffixes[i], 1, len, fp);
  }
  fwrite(chr_lengths, sizeof(int), num_genomes * num_chrs, fp);

  // Write regions
  fwrite(&n_merged, sizeof(size_t), 1, fp);
  for (size_t i = 0; i < n_merged; i++) {
    int chrom_len = strlen(regions[i].chrom);
    fwrite(&chrom_len, sizeof(int), 1, fp);
    fwrite(regions[i].chrom, 1, chrom_len, fp);
    fwrite(&regions[i].start, sizeof(size_t), 1, fp);
    fwrite(&regions[i].end, sizeof(size_t), 1, fp);

    int cluster_id_len = strlen(regions[i].cluster_id);
    fwrite(&cluster_id_len, sizeof(int), 1, fp);
    fwrite(regions[i].cluster_id, 1, cluster_id_len, fp);

    fwrite(&regions[i].subcluster_id, sizeof(uint32_t), 1, fp);
  }
  fclose(fp);
}

int cmd_vis(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: reverb vis <prefix.bin> <cluster_id|all>\n");
    return 1;
  }

  FILE *fp = fopen(argv[1], "rb");
  if (!fp) {
    fprintf(stderr, "Cannot open %s\n", argv[1]);
    return 1;
  }

  char magic[5] = {0};
  fread(magic, 1, 4, fp);
  if (strcmp(magic, "REV2") != 0) {
    fprintf(stderr, "Invalid bin file\n");
    return 1;
  }

  int num_genomes;
  fread(&num_genomes, sizeof(int), 1, fp);
  GenomeVector *gv = calloc(num_genomes, sizeof(GenomeVector));
  for (int i = 0; i < num_genomes; i++) {
    int len;
    fread(&len, sizeof(int), 1, fp);
    fread(gv[i].name, 1, len, fp);
    gv[i].name[len] = '\0';
  }

  int num_chrs;
  fread(&num_chrs, sizeof(int), 1, fp);
  char **chr_suffixes = malloc(num_chrs * sizeof(char *));
  for (int i = 0; i < num_chrs; i++) {
    int len;
    fread(&len, sizeof(int), 1, fp);
    char buf[256] = {0};
    fread(buf, 1, len, fp);
    chr_suffixes[i] = strdup(buf);
  }

  int *chr_lengths = malloc(num_genomes * num_chrs * sizeof(int));
  fread(chr_lengths, sizeof(int), num_genomes * num_chrs, fp);

  size_t n_merged;
  fread(&n_merged, sizeof(size_t), 1, fp);
  ReverbDupRegion *regions = calloc(n_merged, sizeof(ReverbDupRegion));
  for (size_t i = 0; i < n_merged; i++) {
    int chrom_len;
    fread(&chrom_len, sizeof(int), 1, fp);
    char buf[256] = {0};
    fread(buf, 1, chrom_len, fp);
    regions[i].chrom = strdup(buf);
    fread(&regions[i].start, sizeof(size_t), 1, fp);
    fread(&regions[i].end, sizeof(size_t), 1, fp);

    int cluster_id_len;
    fread(&cluster_id_len, sizeof(int), 1, fp);
    char cbuf[512] = {0};
    fread(cbuf, 1, cluster_id_len, fp);
    regions[i].cluster_id = strdup(cbuf);

    fread(&regions[i].subcluster_id, sizeof(uint32_t), 1, fp);
  }
  fclose(fp);
  char *target_cluster = NULL;
  char svg_path[PATH_MAX] = {0};

  if (argc >= 3 && strcmp(argv[2], "all") != 0) {
    target_cluster = argv[2];
  }

  if (target_cluster) {
    snprintf(svg_path, sizeof(svg_path), "%s.svg", target_cluster);
    generate_svg(svg_path, regions, n_merged, target_cluster, gv, num_genomes,
                 (const char **)chr_suffixes, num_chrs, chr_lengths);
    fprintf(stderr, "Wrote %s\n", svg_path);
  } else {
    // all
    char **seen = calloc(n_merged, sizeof(char *));
    size_t n_seen = 0;

    for (size_t i = 0; i < n_merged; i++) {
      char *cluster = regions[i].cluster_id;
      int found = 0;
      for (size_t j = 0; j < n_seen; j++) {
        if (strcmp(seen[j], cluster) == 0) {
          found = 1;
          break;
        }
      }
      if (!found) {
        seen[n_seen++] = cluster;
        snprintf(svg_path, sizeof(svg_path), "%s.svg", cluster);
        generate_svg(svg_path, regions, n_merged, cluster, gv, num_genomes,
                     (const char **)chr_suffixes, num_chrs, chr_lengths);
      }
    }
    free(seen);
    fprintf(stderr, "Wrote SVGs for all families\n");
  }

  for (int i = 0; i < num_genomes; i++) {
    free(gv[i].vector);
  }
  free(gv);
  for (int i = 0; i < num_chrs; i++) {
    free(chr_suffixes[i]);
  }
  free(chr_suffixes);
  free(chr_lengths);
  for (size_t i = 0; i < n_merged; i++) {
    free(regions[i].chrom);
    free(regions[i].cluster_id);
  }
  free(regions);
  return 0;
}

// ==============================================================
// TE ANNOTATION
// ==============================================================

typedef struct {
  char chrom[64];
  size_t start;
  size_t end;
  char name[128];
  char family[128];
  char classification[128];
} TE;

static TE *load_tes(const char *gff_file, size_t *num_tes) {
  FILE *f = fopen(gff_file, "r");
  if (!f) {
    fprintf(stderr, "[reverb] Warning: Could not open GFF file: %s\n",
            gff_file);
    *num_tes = 0;
    return NULL;
  }

  size_t cap = 10000;
  TE *tes = malloc(cap * sizeof(TE));
  size_t count = 0;

  char line[4096];
  while (fgets(line, sizeof(line), f)) {
    if (line[0] == '#')
      continue;

    char *parts[9];
    int n_parts = 0;
    char *p = line;
    char *tab;
    while ((tab = strchr(p, '\t')) != NULL && n_parts < 8) {
      *tab = '\0';
      parts[n_parts++] = p;
      p = tab + 1;
    }
    parts[n_parts++] = p; // info column

    if (n_parts < 9)
      continue;

    const char *ftype = parts[2];
    if (strcmp(ftype, "mobile_element") == 0 ||
        strcmp(ftype, "transposable_element") == 0 ||
        strcmp(ftype, "repeat_region") == 0) {

      TE te;
      memset(&te, 0, sizeof(TE));
      size_t chrom_len = strlen(parts[0]);
      if (chrom_len >= sizeof(te.chrom)) chrom_len = sizeof(te.chrom) - 1;
      memcpy(te.chrom, parts[0], chrom_len);
      te.chrom[chrom_len] = '\0';
      te.start = strtoull(parts[3], NULL, 10);
      te.end = strtoull(parts[4], NULL, 10);

      strcpy(te.name, "NA");
      strcpy(te.family, "NA");
      strcpy(te.classification, "NA");

      char *info = parts[8];
      char *item = strtok(info, ";\n");
      while (item) {
        if (strncmp(item, "Name=", 5) == 0) {
          strncpy(te.name, item + 5, sizeof(te.name) - 1);
        } else if (strncmp(item, "family_name=", 12) == 0) {
          strncpy(te.family, item + 12, sizeof(te.family) - 1);
        } else if (strncmp(item, "classification=", 15) == 0) {
          strncpy(te.classification, item + 15, sizeof(te.classification) - 1);
        }
        item = strtok(NULL, ";\n");
      }

      if (count >= cap) {
        cap *= 2;
        tes = realloc(tes, cap * sizeof(TE));
      }
      tes[count++] = te;
    }
  }
  fclose(f);

  *num_tes = count;
  fprintf(stderr, "[reverb] Loaded %zu TEs from %s\n", count, gff_file);
  return tes;
}

static TE *annotate_cluster(const char *cluster_id, TE *tes, size_t num_tes) {
  if (!tes || num_tes == 0)
    return NULL;

  char chrom[64];
  size_t c_start = 0, c_end = 0;

  // Parse cluster_id: e.g., TAIR12_1Feb26-Chr1@1465000-1470000
  const char *at = strchr(cluster_id, '@');
  if (!at)
    return NULL;

  const char *dash = at;
  while (dash > cluster_id && *dash != '-')
    dash--;
  if (*dash != '-')
    return NULL;

  size_t chrom_len = at - dash - 1;
  if (chrom_len >= sizeof(chrom))
    chrom_len = sizeof(chrom) - 1;
  strncpy(chrom, dash + 1, chrom_len);
  chrom[chrom_len] = '\0';

  sscanf(at + 1, "%zu-%zu", &c_start, &c_end);

  size_t max_overlap = 0;
  TE *best_te = NULL;

  for (size_t i = 0; i < num_tes; i++) {
    if (strcmp(tes[i].chrom, chrom) == 0) {
      if (tes[i].start <= c_end && tes[i].end >= c_start) {
        size_t o_start = (c_start > tes[i].start) ? c_start : tes[i].start;
        size_t o_end = (c_end < tes[i].end) ? c_end : tes[i].end;
        size_t overlap = o_end - o_start + 1;
        if (overlap > max_overlap) {
          max_overlap = overlap;
          best_te = &tes[i];
        }
      }
    }
  }

  return best_te;
}

int cmd_pangenome(int argc, char **argv, const char *pangenome_dir,
                  const char *rep_fasta, size_t flank_size, const Reverb *r,
                  uint64_t scale, size_t window_size, size_t step_size,
                  size_t min_bases, double max_dist, int min_copy, int max_copy,
                  const char *out_prefix, const char *gff_file, int n_threads) {
  (void)argc;
  (void)argv;

  // 1. Directory parsing
  char **files = NULL;
  int num_files = 0;

  // Add representative genome at index 0
  files = realloc(files, (num_files + 1) * sizeof(char *));
  files[num_files++] = strdup(rep_fasta);

  if (pangenome_dir) {
    DIR *d = opendir(pangenome_dir);
    if (!d) {
      fprintf(stderr, "Cannot open %s\n", pangenome_dir);
      return 1;
    }
    struct dirent *ent;

    char rep_bname[256];
    get_basename(rep_fasta, rep_bname, sizeof(rep_bname));

    while ((ent = readdir(d)) != NULL) {
      if (strstr(ent->d_name, ".fasta") || strstr(ent->d_name, ".fa")) {
        if (strstr(ent->d_name, rep_bname))
          continue;
        char path[1024];
        snprintf(path, sizeof(path), "%s/%s", pangenome_dir, ent->d_name);
        files = realloc(files, (num_files + 1) * sizeof(char *));
        files[num_files++] = strdup(path);
      }
    }
    closedir(d);

    if (num_files > 1) {
      qsort(files + 1, num_files - 1, sizeof(char *), cmp_natural);
    }
  }

  ReverbSketch *sketches = NULL;
  WindowCoord *coords = NULL;
  size_t num_sketches = 0, cap_sketches = 0;
  ReverbDupEdge *edges = NULL;

  char path_buf[PATH_MAX];
  snprintf(path_buf, sizeof(path_buf), "%s.window.bed", out_prefix);
  FILE *bed_fp = fopen(path_buf, "w");

  fprintf(stderr, "[reverb] Extracting windows across pangenome...\n");
  GenomeSeqLen *seq_lens = NULL;
  size_t num_seqs = 0, cap_seqs = 0;

  for (int i = 0; i < num_files; i++) {
    char bname[256];
    get_basename(files[i], bname, sizeof(bname));
    dup_stream_pangenome(files[i], bname, r, scale, window_size, step_size,
                         min_bases, &sketches, &coords, &num_sketches,
                         &cap_sketches, bed_fp, &seq_lens, &num_seqs,
                         &cap_seqs);
  }
  fclose(bed_fp);

  size_t n_edges = build_candidate_edges(
      sketches, coords, num_sketches, max_dist, window_size, n_threads, &edges);

  uint32_t *genome_id = calloc(num_sketches, sizeof(uint32_t));
  for (size_t i = 0; i < num_sketches; i++) {
    for (int f = 0; f < num_files; f++) {
      char bname[256];
      get_basename(files[f], bname, sizeof(bname));
      int blen = strlen(bname);
      if (strncmp(coords[i].chrom, bname, blen) == 0 &&
          coords[i].chrom[blen] == '-') {
        genome_id[i] = f;
        break;
      }
    }
  }

  // Step 1: Intra-genome UnionFind for Representative Genome (genome_id == 0)
  // ONLY
  UnionFind uf_intra;
  uf_init(&uf_intra, num_sketches);
  for (size_t i = 0; i < n_edges; i++) {
    if (genome_id[edges[i].win_a] == 0 && genome_id[edges[i].win_b] == 0) {
      uf_union(&uf_intra, edges[i].win_a, edges[i].win_b);
    }
  }

  uint32_t *comp_size_intra = calloc(num_sketches, sizeof(uint32_t));
  for (size_t i = 0; i < num_sketches; i++) {
    comp_size_intra[uf_find(&uf_intra, (uint32_t)i)]++;
  }

  uint8_t *is_sd = calloc(num_sketches, sizeof(uint8_t));
  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t fam = uf_find(&uf_intra, (uint32_t)i);
    if (comp_size_intra[fam] >= (uint32_t)min_copy &&
        (max_copy <= 0 || comp_size_intra[fam] <= (uint32_t)max_copy)) {
      is_sd[i] = 1;
    }
  }
  free(comp_size_intra);
  uf_free(&uf_intra);

  // Step 2: Final UnionFind
  UnionFind uf;
  uf_init(&uf, num_sketches);
  for (size_t i = 0; i < n_edges; i++) {
    uint32_t a = edges[i].win_a;
    uint32_t b = edges[i].win_b;
    if (genome_id[a] == genome_id[b]) {
      uf_union(&uf, a, b);
    } else if (is_sd[a] || is_sd[b]) {
      uf_union(&uf, a, b);
    }
  }

  uint8_t *final_is_sd = calloc(num_sketches, sizeof(uint8_t));
  for (size_t i = 0; i < num_sketches; i++) {
    if (is_sd[i])
      final_is_sd[uf_find(&uf, (uint32_t)i)] = 1;
  }
  free(is_sd);

  uint32_t *comp_size = calloc(num_sketches, sizeof(uint32_t));
  for (size_t i = 0; i < num_sketches; i++)
    comp_size[uf_find(&uf, (uint32_t)i)]++;

  char **hub_label = calloc(num_sketches, sizeof(char *));
  for (size_t i = 0; i < num_sketches; i++) {
    if (genome_id[i] == 0) {
      uint32_t fam = uf_find(&uf, (uint32_t)i);
      if (!hub_label[fam]) {
        char buf[512];
        snprintf(buf, sizeof(buf), "%s@%zu-%zu", coords[i].chrom,
                 coords[i].start, coords[i].end);
        hub_label[fam] = strdup(buf);
      }
    }
  }
  free(genome_id);

  size_t n_dup_regions = 0, cap_dup_regions = 0;
  ReverbDupRegion *dup_regions = NULL;
  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t fam = uf_find(&uf, (uint32_t)i);
    if (!final_is_sd[fam])
      continue;

    const char *label = hub_label[fam] ? hub_label[fam] : "unknown";

    DA_PUSH(dup_regions, n_dup_regions, cap_dup_regions,
            ((ReverbDupRegion){.chrom = coords[i].chrom,
                               .start = coords[i].start,
                               .end = coords[i].end,
                               .cluster_id = strdup(label),
                               .copy_count = comp_size[fam],
                               .avg_distance = 0.0,
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
    uint32_t fam = uf_find(&uf, a);
    if (!hub_label[fam])
      continue;
    fprintf(out_bedpe, "%s\t%zu\t%zu\t%s\t%zu\t%zu\t%s\t%.6f\n",
            coords[a].chrom, coords[a].start, coords[a].end, coords[b].chrom,
            coords[b].start, coords[b].end, hub_label[fam], edges[i].distance);
  }
  fclose(out_bedpe);

  if (!pangenome_dir) {
    snprintf(path_buf, sizeof(path_buf), "%s.dup.bed", out_prefix);
    FILE *out_bed = fopen(path_buf, "w");
    fprintf(out_bed, "#chrom\tstart\tend\tcluster_id\tcopy_count\n");
    for (size_t i = 0; i < n_merged; i++) {
      fprintf(out_bed, "%s\t%zu\t%zu\t%s\t%u\n", dup_regions[i].chrom,
              dup_regions[i].start, dup_regions[i].end,
              dup_regions[i].cluster_id, dup_regions[i].copy_count);
    }
    fclose(out_bed);

    for (size_t i = 0; i < num_sketches; i++) {
      if (hub_label[i])
        free(hub_label[i]);
    }
    free(hub_label);

    return 0;
  }

  fprintf(stderr, "[reverb] Pass 2: Extracting flanking sequences...\n");
  do_pass2(files, num_files, r, scale, dup_regions, n_merged,
           flank_size == 0 ? window_size : flank_size);

  fprintf(stderr, "[reverb] Pass 2: Sub-clustering flanking sequences...\n");
  do_subclustering(dup_regions, n_merged, max_dist, n_threads);

  // Output dup.bed
  snprintf(path_buf, sizeof(path_buf), "%s.dup.bed", out_prefix);
  FILE *out_bed = fopen(path_buf, "w");

  size_t num_tes = 0;
  TE *tes = NULL;
  if (gff_file) {
    tes = load_tes(gff_file, &num_tes);
    fprintf(out_bed, "#chrom\tstart\tend\tcluster_id\tsubcluster_"
                     "id\tName\tfamily_name\tclassification\n");
  } else {
    fprintf(out_bed, "#chrom\tstart\tend\tcluster_id\tsubcluster_id\n");
  }

  char rep_bname[256];
  get_basename(rep_fasta, rep_bname, sizeof(rep_bname));
  size_t rep_len = strlen(rep_bname);

  uint32_t max_subcluster = 0;
  for (size_t i = 0; i < n_merged; i++) {
    if (gff_file) {
      TE *best_te = annotate_cluster(dup_regions[i].cluster_id, tes, num_tes);
      if (best_te) {
        int is_rep = (strncmp(dup_regions[i].chrom, rep_bname, rep_len) == 0 &&
                      (dup_regions[i].chrom[rep_len] == '-' ||
                       dup_regions[i].chrom[rep_len] == '\0'));
        const char *out_name = is_rep ? best_te->name : "NA";

        fprintf(out_bed, "%s\t%zu\t%zu\t%s\t%u\t%s\t%s\t%s\n",
                dup_regions[i].chrom, dup_regions[i].start, dup_regions[i].end,
                dup_regions[i].cluster_id, dup_regions[i].subcluster_id,
                out_name, best_te->family, best_te->classification);
      } else {
        fprintf(out_bed, "%s\t%zu\t%zu\t%s\t%u\tNA\tNA\tNA\n",
                dup_regions[i].chrom, dup_regions[i].start, dup_regions[i].end,
                dup_regions[i].cluster_id, dup_regions[i].subcluster_id);
      }
    } else {
      fprintf(out_bed, "%s\t%zu\t%zu\t%s\t%u\n", dup_regions[i].chrom,
              dup_regions[i].start, dup_regions[i].end,
              dup_regions[i].cluster_id, dup_regions[i].subcluster_id);
    }

    if (dup_regions[i].subcluster_id > max_subcluster)
      max_subcluster = dup_regions[i].subcluster_id;
  }
  fclose(out_bed);
  if (tes)
    free(tes);

  // Sort genomes for SVG
  GenomeVector *gv = calloc(num_files, sizeof(GenomeVector));
  for (int i = 0; i < num_files; i++) {
    char bname[256];
    get_basename(files[i], bname, sizeof(bname));
    snprintf(gv[i].name, sizeof(gv[i].name), "%s", bname);
    gv[i].vector = calloc(max_subcluster + 1, sizeof(uint8_t));
  }
  for (size_t i = 0; i < n_merged; i++) {
    int g_idx = -1;
    for (int g = 0; g < num_files; g++) {
      int len = strlen(gv[g].name);
      if (strncmp(dup_regions[i].chrom, gv[g].name, len) == 0 &&
          dup_regions[i].chrom[len] == '-') {
        g_idx = g;
        break;
      }
    }
    if (g_idx != -1) {
      gv[g_idx].vector[dup_regions[i].subcluster_id] = 1;
    }
  }

  natural_sort_genomes(gv, num_files, max_subcluster);

  char **chr_suffixes = NULL;
  int num_chrs = 0;
  int cap_chrs = 0;
  for (size_t i = 0; i < num_seqs; i++) {
    int found = 0;
    for (int j = 0; j < num_chrs; j++) {
      if (strcmp(chr_suffixes[j], seq_lens[i].seq) == 0) {
        found = 1;
        break;
      }
    }
    if (!found) {
      if (num_chrs >= cap_chrs) {
        cap_chrs = cap_chrs == 0 ? 8 : cap_chrs * 2;
        chr_suffixes = realloc(chr_suffixes, cap_chrs * sizeof(char *));
      }
      chr_suffixes[num_chrs++] = strdup(seq_lens[i].seq);
    }
  }

  qsort(chr_suffixes, num_chrs, sizeof(char *), cmp_natural);

  int *chr_lengths = calloc(num_files * num_chrs, sizeof(int));
  for (size_t i = 0; i < num_seqs; i++) {
    int y_idx = -1;
    for (int g = 0; g < num_files; g++) {
      if (strcmp(gv[g].name, seq_lens[i].genome) == 0) {
        y_idx = g;
        break;
      }
    }
    int c_idx = -1;
    for (int c = 0; c < num_chrs; c++) {
      if (strcmp(chr_suffixes[c], seq_lens[i].seq) == 0) {
        c_idx = c;
        break;
      }
    }
    if (y_idx != -1 && c_idx != -1) {
      chr_lengths[y_idx * num_chrs + c_idx] = seq_lens[i].length;
    }
  }

  fprintf(stderr, "[reverb] Generating SVG visualizations...\n");
  char dir_buf[PATH_MAX];
  snprintf(dir_buf, sizeof(dir_buf), "%s_vis", out_prefix);
  mkdir(dir_buf, 0755);

  // Get unique families
  char **seen = calloc(n_merged, sizeof(char *));
  size_t n_seen = 0;

  for (size_t i = 0; i < n_merged; i++) {
    char *cluster = dup_regions[i].cluster_id;
    int found = 0;
    for (size_t j = 0; j < n_seen; j++) {
      if (strcmp(seen[j], cluster) == 0) {
        found = 1;
        break;
      }
    }
    if (!found) {
      seen[n_seen++] = cluster;
      size_t svg_len = strlen(dir_buf) + strlen(cluster) + 6;
      char *svg_path = malloc(svg_len);
      snprintf(svg_path, svg_len, "%s/%s.svg", dir_buf, cluster);
      generate_svg(svg_path, dup_regions, n_merged, cluster, gv, num_files,
                   (const char **)chr_suffixes, num_chrs, chr_lengths);
      free(svg_path);
    }
  }
  free(seen);

  write_bin(out_prefix, dup_regions, n_merged, gv, num_files,
            (const char **)chr_suffixes, num_chrs, chr_lengths);

  for (size_t i = 0; i < num_sketches; i++) {
    if (hub_label[i])
      free(hub_label[i]);
  }
  free(hub_label);

  return 0;
}

int cmd_dup(int argc, char **argv) {

  // ==============================================================
  // CLI defaults
  // ==============================================================

  uint32_t def_kmer_size = 21;
  uint64_t def_scale = 10;
  uint64_t def_hash_seed = 42;
  size_t window_size = 5000;
  size_t step_size = 0; /* 0 = auto (window/2) */
  size_t min_bases = 1000;
  double max_dist = 0.05;
  int min_copy = 3;
  int max_copy = 30;
  const char *out_prefix = "reverb";
  const char *pangenome_dir = NULL;
  const char *rep_fasta = NULL;
  const char *gff_file = NULL;
  size_t flank_size = 0;
  int n_threads = 8;

  ketopt_t opt = KETOPT_INIT;
  int c;
  while ((c = ketopt(&opt, argc, argv, 1,
                     "k:s:e:w:t:b:d:m:M:o:p:hI:i:f:g:", 0)) >= 0) {
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
    else if (c == 'I')
      pangenome_dir = opt.arg;
    else if (c == 'i')
      rep_fasta = opt.arg;
    else if (c == 'f')
      flank_size = (size_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'g')
      gff_file = opt.arg;
    else
      return 1;
  }

  if (step_size == 0)
    step_size = window_size / 2;

  if (!rep_fasta) {
    fprintf(stderr, "Error: -i <fasta> is required.\n");
    return 1;
  }

  Reverb r;
  reverb_init(&r, def_kmer_size);
  r.hash_seed = def_hash_seed;
  return cmd_pangenome(argc, argv, pangenome_dir, rep_fasta, flank_size, &r,
                       def_scale, window_size, step_size, min_bases, max_dist,
                       min_copy, max_copy, out_prefix, gff_file, n_threads);
}

// ==============================================================
// MAIN
// ==============================================================

int main(int argc, char **argv) {
  if (argc > 1 && strcmp(argv[1], "vis") == 0) {
    return cmd_vis(argc - 1, argv + 1);
  }

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
