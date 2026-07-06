#include <math.h>
#include <stdio.h>
#include <time.h>
#include <zlib.h>

#include "klib/ketopt.h"
#include "klib/kseq.h"
#include "reverb.h"

#define MIX_CONST1 0xff51afd7ed558ccdULL
#define MIX_CONST2 0xc4ceb9fe1a85ec53ULL

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

static size_t min_size_t(size_t a, size_t b) { return a < b ? a : b; }

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
  size_t k = min_size_t(hash_window, 42); /* 3*42 = 126 bits fits __uint128_t */
  uint32_t kmer_bits = 3 * (uint32_t)k;

  __uint128_t remover_mask =
      (kmer_bits > 3) ? (((__uint128_t)1 << (kmer_bits - 3)) - 1) : 0;
  __uint128_t layer_mask = 0;
  for (uint32_t j = 0; j < (uint32_t)k; j++)
    layer_mask |= (__uint128_t)1 << (3 * j);

  r->hash_window = k;
  r->remover_mask = remover_mask;
  r->layer_mask = layer_mask;
  r->kmer_bits = kmer_bits;
  r->rc_shift = (kmer_bits > 0) ? (128 - kmer_bits) : 128;
}

// ==============================================================
// HASH POOL (OVERFLOWING BUFFER)
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
  pool->size = 0;
  pool->cap = 0;
  pool->hash_threshold = threshold;
  pool->hashes = NULL;
}

static void pool_try_insert(HashPool *pool, uint64_t h) {
  if (h >= pool->hash_threshold)
    return;
  if (pool->size == pool->cap) {
    size_t nc = pool->cap ? pool->cap * 2 : 64;
    uint64_t *nh = realloc(pool->hashes, nc * sizeof(uint64_t));
    if (!nh)
      return;
    pool->hashes = nh;
    pool->cap = nc;
  }
  pool->hashes[pool->size++] = h;
}

static void pool_finalize(HashPool *pool, uint64_t **out_hashes,
                          size_t *out_size) {
  size_t n = pool->size;
  if (n) {
    qsort(pool->hashes, n, sizeof(uint64_t), cmp_uint64);
    size_t u = 0;
    for (size_t i = 0; i < n; i++) {
      if (u == 0 || pool->hashes[i] != pool->hashes[u - 1]) {
        pool->hashes[u++] = pool->hashes[i];
      }
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
__attribute__((hot)) static void
extract_and_insert(const Reverb *r, HashPool *pool, const uint8_t *seq,
                   size_t len, int min_count, uint8_t *cms_table,
                   size_t cms_mask) {
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

    if (valid >= K) {
      __uint128_t rev = reverse_bits128(fwd) >> r->rc_shift;
      __uint128_t canon = fwd < rev ? fwd : rev;
      uint64_t h = mix_hash(canon, r->hash_seed);

      if (h >= pool->hash_threshold)
        continue;

      if (min_count > 1 && cms_table) {
        uint64_t k1 = h, k2 = h * MIX_CONST1, k3 = h * MIX_CONST2;
        size_t h1 = (size_t)((k1 ^ (k1 >> 32)) & cms_mask);
        size_t h2 = (size_t)((k2 ^ (k2 >> 32)) & cms_mask);
        size_t h3 = (size_t)((k3 ^ (k3 >> 32)) & cms_mask);

        uint8_t c1 = cms_table[h1], c2 = cms_table[h2], c3 = cms_table[h3];
        uint8_t min_c =
            (c1 < c2) ? ((c1 < c3) ? c1 : c3) : ((c2 < c3) ? c2 : c3);
        uint8_t thr =
            (uint8_t)(min_count > (int)UINT8_MAX ? UINT8_MAX : min_count);

        if (min_c < thr) {
          if (c1 == min_c && c1 < UINT8_MAX)
            cms_table[h1]++;
          if (c2 == min_c && c2 < UINT8_MAX)
            cms_table[h2]++;
          if (c3 == min_c && c3 < UINT8_MAX)
            cms_table[h3]++;
          if ((uint8_t)(min_c + 1) != thr)
            continue;
        } else {
          continue;
        }
      }

      pool_try_insert(pool, h);
    }
  }
}

// ==============================================================
// FILE STREAMING
// ==============================================================
void reverb_stream(const char *filename, const Reverb *r, HashPool *pool,
                   int min_count, uint8_t *cms_table, size_t cms_mask) {
  gzFile fp = gzopen(filename, "r");
  if (!fp)
    return;
  kseq_t *ks = kseq_init(fp);
  if (!ks) {
    gzclose(fp);
    return;
  }

  while (kseq_read(ks) >= 0) {
    if (ks->seq.l < r->hash_window)
      continue;
    extract_and_insert(r, pool, (const uint8_t *)ks->seq.s, ks->seq.l,
                       min_count, cms_table, cms_mask);
  }
  kseq_destroy(ks);
  gzclose(fp);
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

ReverbTripleDistResult reverb_dist_three(const ReverbSketch *ref,
                                         const ReverbSketch *q1,
                                         const ReverbSketch *q2) {
  ReverbTripleDistResult t;
  t.d01 = t.d02 = t.d12 = (ReverbDistResult){0.0, 1.0, 0};

  if (ref->kmer_size == 0 || q1->kmer_size != ref->kmer_size ||
      q2->kmer_size != ref->kmer_size) {
    return t;
  }

  size_t nR = ref->sketch_size, n1 = q1->sketch_size, n2 = q2->sketch_size;
  const uint64_t *R = ref->hashes, *Q1 = q1->hashes, *Q2 = q2->hashes;

  size_t sh01 = 0, sh02 = 0, sh12 = 0;

  size_t i = 0, j = 0, k = 0;
  while (i < nR || j < n1 || k < n2) {
    uint64_t hr = (i < nR) ? R[i] : UINT64_MAX;
    uint64_t h1 = (j < n1) ? Q1[j] : UINT64_MAX;
    uint64_t h2 = (k < n2) ? Q2[k] : UINT64_MAX;

    uint64_t hm = hr < h1 ? (hr < h2 ? hr : h2) : (h1 < h2 ? h1 : h2);
    if (hm == UINT64_MAX)
      break;

    int inR = (hr == hm);
    int inQ1 = (h1 == hm);
    int inQ2 = (h2 == hm);

    if (inR && inQ1)
      sh01++;
    if (inR && inQ2)
      sh02++;
    if (inQ1 && inQ2)
      sh12++;

    if (inR)
      i++;
    if (inQ1)
      j++;
    if (inQ2)
      k++;
  }

  double inv_k = 1.0 / (double)ref->kmer_size;

#define FILL(r, sh, a, b)                                                      \
  do {                                                                         \
    (r).shared_hashes = (sh);                                                  \
    if ((a)->sketch_size > 0 && (b)->sketch_size > 0) {                        \
      double ca = (double)(sh) / (double)(a)->sketch_size;                     \
      double cb = (double)(sh) / (double)(b)->sketch_size;                     \
      (r).containment = 0.5 * (ca + cb);                                       \
      (r).distance = 1.0 - pow((r).containment, inv_k);                        \
    }                                                                          \
  } while (0)

  FILL(t.d01, sh01, ref, q1);
  FILL(t.d02, sh02, ref, q2);
  FILL(t.d12, sh12, q1, q2);
#undef FILL

  return t;
}

void reverb_info(const ReverbSketch *sk) {
  if (!sk)
    return;
  printf("--- Reverb Sketch Info ---\n");
  printf("Name                 : %s\n", sk->name ? sk->name : "N/A");
  printf("K-mer Size (K)       : %u\n", sk->kmer_size);
  printf("Total Sketch Size    : %zu\n", sk->sketch_size);
  printf("--------------------------\n");
}

// ==============================================================
// SERIALIZATION
// ==============================================================
int reverb_sketch_save(const ReverbSketch *sk, const char *filepath) {
  FILE *f = fopen(filepath, "wb");
  if (!f)
    return -1;

  uint32_t name_len = sk->name ? (uint32_t)strlen(sk->name) : 0;
  uint64_t sz = sk->sketch_size;

  fwrite("RVRB", 1, 4, f);
  fwrite(&sk->kmer_size, sizeof(uint32_t), 1, f);
  fwrite(&sk->hash_threshold, sizeof(uint64_t), 1, f);
  fwrite(&sz, sizeof(uint64_t), 1, f);
  fwrite(&name_len, sizeof(uint32_t), 1, f);
  if (name_len > 0)
    fwrite(sk->name, 1, name_len, f);
  if (sk->sketch_size > 0) {
    fwrite(sk->hashes, sizeof(uint64_t), sk->sketch_size, f);
  }
  fclose(f);
  return 0;
}

int reverb_sketch_load(ReverbSketch *sk, const char *filepath) {
  FILE *f = fopen(filepath, "rb");
  if (!f)
    return -1;

  char magic[4];
  if (fread(magic, 1, 4, f) != 4 || memcmp(magic, "RVRB", 4) != 0) {
    fclose(f);
    return -1;
  }

  uint32_t name_len = 0;
  uint64_t sz = 0;
  sk->name = NULL;
  sk->hashes = NULL;

  if (fread(&sk->kmer_size, sizeof(uint32_t), 1, f) != 1)
    goto err;
  if (fread(&sk->hash_threshold, sizeof(uint64_t), 1, f) != 1)
    goto err;
  if (fread(&sz, sizeof(uint64_t), 1, f) != 1)
    goto err;
  if (fread(&name_len, sizeof(uint32_t), 1, f) != 1)
    goto err;

  sk->sketch_size = (size_t)sz;

  if (name_len > 0) {
    sk->name = (char *)malloc(name_len + 1);
    if (!sk->name || fread(sk->name, 1, name_len, f) != name_len)
      goto err;
    sk->name[name_len] = '\0';
  }

  if (sk->sketch_size > 0) {
    sk->hashes = (uint64_t *)malloc(sk->sketch_size * sizeof(uint64_t));
    if (!sk->hashes)
      goto err;
    if (fread(sk->hashes, sizeof(uint64_t), sk->sketch_size, f) !=
        sk->sketch_size)
      goto err;
  }

  fclose(f);
  return 0;

err:
  free(sk->name);
  free(sk->hashes);
  sk->name = NULL;
  sk->hashes = NULL;
  fclose(f);
  return -1;
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
    uf->parent[x] = uf->parent[uf->parent[x]]; /* path compression */
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
typedef struct {
  uint32_t kmer_size;
  uint64_t scale;
  uint64_t hash_seed;
  uint8_t num_threads;
  uint8_t min_count;
  uint8_t read_pooling;
} SketchBuildParams;

static SketchBuildParams sketch_build_defaults(void) {
  return (SketchBuildParams){.kmer_size = 21,
                             .scale = 1000,
                             .hash_seed = MIX_CONST1,
                             .num_threads = 1,
                             .min_count = 0,
                             .read_pooling = 0};
}

static int build_sketch_from_fasta(ReverbSketch *sk, const char *path,
                                   const SketchBuildParams *params) {
  memset(sk, 0, sizeof(ReverbSketch));
  sk->name = strdup(path);
  sk->kmer_size = params->kmer_size;
  sk->hash_threshold = UINT64_MAX / params->scale;

  if (!params->scale)
    return -1;

  Reverb r;
  reverb_init(&r, params->kmer_size);
  r.hash_seed = params->hash_seed;

  HashPool pool;
  pool_init(&pool, sk->hash_threshold);
  reverb_stream(path, &r, &pool, params->min_count, NULL, 0);
  pool_finalize(&pool, &sk->hashes, &sk->sketch_size);

  return sk->sketch_size > 0 ? 0 : -1;
}

// ==============================================================
// IO
// ==============================================================
static int has_extension(const char *filename, const char *ext) {
  size_t len = strlen(filename);
  size_t elen = strlen(ext);
  if (len < elen)
    return 0;
  return strcasecmp(filename + len - elen, ext) == 0;
}

static int is_sequence_file(const char *filename) {
  return has_extension(filename, ".fa") || has_extension(filename, ".fasta") ||
         has_extension(filename, ".fq") || has_extension(filename, ".fastq") ||
         has_extension(filename, ".fa.gz") ||
         has_extension(filename, ".fasta.gz") ||
         has_extension(filename, ".fq.gz") ||
         has_extension(filename, ".fastq.gz");
}

static int prepare_input_sketches(ReverbSketch *sketches, int *ready,
                                  char *const *paths, int num_files,
                                  SketchBuildParams *params, int *failed) {
  if (failed)
    *failed = -1;
  int params_seeded = 0;

  for (int i = 0; i < num_files; i++) {
    ready[i] = 0;

    if (!is_sequence_file(paths[i])) {
      if (reverb_sketch_load(&sketches[i], paths[i]) == 0) {
        ready[i] = 1;
        if (!params_seeded) {
          params->kmer_size = sketches[i].kmer_size;
          params_seeded = 1;
        }
      } else if (has_extension(paths[i], ".reverb")) {
        if (failed)
          *failed = i;
        return -1;
      }
    }
  }

  for (int i = 0; i < num_files; i++) {
    if (ready[i])
      continue;

    if (build_sketch_from_fasta(&sketches[i], paths[i], params) != 0) {
      if (failed)
        *failed = i;
      return -1;
    }
    ready[i] = 1;
  }

  return 0;
}
static void free_ready_sketches(ReverbSketch *sketches, const int *ready,
                                int num_files) {
  for (int i = 0; i < num_files; i++)
    if (ready[i])
      reverb_sketch_free(&sketches[i]);
}

// ==============================================================
// PRINT REPORTS
// ==============================================================

static void print_distance_report(const ReverbSketch *ref,
                                  const ReverbSketch *query,
                                  const ReverbDistResult *dist) {
  printf("Reference File      : %s\n", ref->name);
  printf("Query File          : %s\n", query->name);
  printf("Shared Hashes       : %zu\n", dist->shared_hashes);
  printf("------------------------------------\n");
  printf("Average Containment : %f\n", dist->containment);
  printf("Reverb Distance     : %f\n", dist->distance);
}

static void print_usage(void) {
  printf(
      "Reverb: Ultra-fast Alignment-free Segmental Duplication Detection\n\n"
      "Usage: reverb [command] [options] [arguments]\n\n"
      "Default (no command): runs 'dup' on the given FASTA files.\n\n"
      "Commands:\n"
      "  help | -h | --help\n"
      "  dup    [-k K] [-s S] [-w win] [-t step] [-b min_bases] [-d max_dist]\n"
      "         [-m min_copy] [-M max_copy] [-o prefix] [-p threads] fasta1 "
      "[fasta2 ...]\n"
      "         -k: kmer size (default: 21, max: 42)\n"
      "         -s: scale factor (default: 10)\n"
      "         -w: window size in bp (default: 10000)\n"
      "         -t: step size in bp (default: window/2)\n"
      "         -b: minimum valid bases per window (default: 1000)\n"
      "         -d: maximum distance to consider as copy (default: 0.2)\n"
      "         -m: minimum copy count (default: 2)\n"
      "         -M: maximum copy count to filter ubiquitous repeats (default: "
      "unlimited)\n"
      "         -o: output file prefix (default: reverb)\n"
      "         -p: number of threads (default: 8)\n"
      "  sketch [-k K] [-s S] [-e E] [-p threads] [-r] [-m min_count] "
      "[-o out_file] fasta1 [fasta2 ...]\n"
      "         -k: kmer size (default: 21, max: 42)\n"
      "         -s: scale factor (keep 1/S of k-mers, default: 1000)\n"
      "         -e: hash seed\n"
      "         -p: number of threads (default: 1)\n"
      "         -r: pool all reads from input files into one sketch\n"
      "         -m: minimum k-mer count filter (used with -r, default: 2)\n"
      "         -o: output file path\n"
      "  dist   <sketch|fasta1> <sketch|fasta2>\n"
      "  three  <ref> <query1> <query2>\n"
      "  triangle <sketch|fasta1> ... <sketch|fastaN>\n"
      "  info   <sketch>\n"
      "\n");
}

// ==============================================================
// SKETCH JOBS
// ==============================================================
typedef struct {
  const char *in_file;
  const char *out_file;
  uint32_t kmer_size;
  uint64_t hash_threshold;
  uint64_t hash_seed;
  uint8_t min_count;
  int result;
} SketchJob;

static void sketch_worker_kt(void *data, long i, int _unused) {
  (void)_unused;
  SketchJob *job = &((SketchJob *)data)[i];

  ReverbSketch sk = {0};
  sk.name = job->in_file ? strdup(job->in_file) : NULL;
  sk.kmer_size = job->kmer_size;
  sk.hash_threshold = job->hash_threshold;

  if (job->hash_threshold) {
    Reverb r;
    reverb_init(&r, job->kmer_size);
    r.hash_seed = job->hash_seed;

    HashPool pool;
    pool_init(&pool, job->hash_threshold);
    reverb_stream(job->in_file, &r, &pool, job->min_count, NULL, 0);
    pool_finalize(&pool, &sk.hashes, &sk.sketch_size);
  }

  job->result = reverb_sketch_save(&sk, job->out_file);
  reverb_sketch_free(&sk);
}

// ==============================================================
// CMD: SKETCH
// ==============================================================
static int setup_dist_cmds(int argc, char **argv, int min_files,
                           ReverbSketch **sk, int **ready) {
  ketopt_t opt = KETOPT_INIT;
  int c;

  while ((c = ketopt(&opt, argc - 1, argv + 1, 1, "", 0)) >= 0) {
    return -1;
  }

  int num_files = argc - (opt.ind + 1);
  if (num_files < min_files) {
    fprintf(stderr, "Error: requires at least %d input file(s)\n", min_files);
    return -1;
  }

  *sk = (ReverbSketch *)calloc(num_files, sizeof(ReverbSketch));
  *ready = (int *)calloc(num_files, sizeof(int));
  SketchBuildParams params = sketch_build_defaults();
  int failed = -1;

  if (prepare_input_sketches(*sk, *ready, argv + opt.ind + 1, num_files,
                             &params, &failed) != 0) {
    fprintf(stderr, "Error: failed to load/build sketch at index %d\n", failed);
    free_ready_sketches(*sk, *ready, num_files);
    free(*sk);
    free(*ready);
    return -1;
  }
  return num_files;
}

int cmd_sketch(int argc, char **argv) {
  SketchBuildParams def = sketch_build_defaults();
  const char *out_opt = NULL;

  ketopt_t opt = KETOPT_INIT;
  int c;
  while ((c = ketopt(&opt, argc - 1, argv + 1, 1, "k:s:e:p:r:m:o:", 0)) >= 0) {
    if (c == 'k')
      def.kmer_size = (uint32_t)atoi(opt.arg);
    else if (c == 's')
      def.scale = (uint64_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'e') {
      uint64_t hs = strtoull(opt.arg, NULL, 0);
      def.hash_seed = hs % 2 ? hs + 1 : hs;
    } else if (c == 'p')
      def.num_threads = atoi(opt.arg) < 1 ? 1 : atoi(opt.arg);
    else if (c == 'r')
      def.read_pooling = 1;
    else if (c == 'm')
      def.min_count = atoi(opt.arg);
    else if (c == 'o')
      out_opt = opt.arg;
    else
      return 1;
  }

  int num_files = argc - (opt.ind + 1);
  if (num_files < 1)
    return 1;
  char **in_files = argv + opt.ind + 1;

  SketchJob *jobs = (SketchJob *)calloc(num_files, sizeof(SketchJob));
  for (int i = 0; i < num_files; i++) {
    size_t len = out_opt ? strlen(out_opt) + 1 : strlen(in_files[i]) + 16;
    jobs[i] = (SketchJob){.in_file = in_files[i],
                          .kmer_size = def.kmer_size,
                          .hash_threshold = UINT64_MAX / def.scale,
                          .hash_seed = def.hash_seed,
                          .min_count = def.min_count,
                          .result = 0,
                          .out_file = (char *)malloc(len)};
    if (!jobs[i].out_file)
      return 1;
    out_opt ? (long)strcpy((char *)jobs[i].out_file, out_opt)
            : snprintf((char *)jobs[i].out_file, len, "%s.reverb", in_files[i]);
  }

  kt_for(def.num_threads, sketch_worker_kt, jobs, num_files);

  int ret = 0;
  for (int i = 0; i < num_files; i++) {
    if (jobs[i].result != 0)
      ret = 1;
    free((char *)jobs[i].out_file);
  }
  free(jobs);
  return ret;
}

// ==============================================================
// CMD: DIST
// ==============================================================
int cmd_dist(int argc, char **argv) {
  ReverbSketch *sk;
  int *ready;

  if (setup_dist_cmds(argc, argv, 2, &sk, &ready) < 2)
    return 1;

  ReverbDistResult dist = reverb_dist(&sk[0], &sk[1]);
  print_distance_report(&sk[0], &sk[1], &dist);

  free_ready_sketches(sk, ready, 2);
  free(sk);
  free(ready);
  return 0;
}

// ==============================================================
// CMD: THREE
// ==============================================================
int cmd_three(int argc, char **argv) {
  ReverbSketch *sk;
  int *ready;

  if (setup_dist_cmds(argc, argv, 3, &sk, &ready) < 3)
    return 1;

  ReverbTripleDistResult tri = reverb_dist_three(&sk[0], &sk[1], &sk[2]);
  print_distance_report(&sk[0], &sk[1], &tri.d01);
  printf("\n");
  print_distance_report(&sk[0], &sk[2], &tri.d02);
  printf("\n");
  print_distance_report(&sk[1], &sk[2], &tri.d12);

  free_ready_sketches(sk, ready, 3);
  free(sk);
  free(ready);
  return 0;
}

// ==============================================================
// CMD: TRIANGLE
// ==============================================================
int cmd_triangle(int argc, char **argv) {
  ReverbSketch *sk;
  int *ready;

  int n = setup_dist_cmds(argc, argv, 1, &sk, &ready);
  if (n < 1)
    return 1;

  printf("%d\n", n);
  for (int i = 0; i < n; i++) {
    printf("%s", sk[i].name ? sk[i].name : "N/A");
    for (int j = 0; j < i; j++) {
      ReverbDistResult d = reverb_dist(&sk[i], &sk[j]);
      printf("\t%f", d.distance);
    }
    printf("\n");
  }

  free_ready_sketches(sk, ready, n);
  free(sk);
  free(ready);
  return 0;
}

// ==============================================================
// CMD: INFO
// ==============================================================
int cmd_info(int argc, char **argv) {
  if (argc < 3)
    return 1;

  ReverbSketch sk;
  if (reverb_sketch_load(&sk, argv[2]) != 0)
    return 1;

  reverb_info(&sk);
  reverb_sketch_free(&sk);
  return 0;
}

// ==============================================================
// DUP: WINDOW STREAMING
// ==============================================================

/* Window metadata for tracking genomic coordinates */
typedef struct {
  char *chrom;
  size_t start;
  size_t end;
} WindowCoord;

static int dup_stream(const char *filename, const Reverb *r,
                      SketchBuildParams *params, size_t window_size,
                      size_t step_size, size_t min_bases,
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
      size_t chunk_len = window_size;

      size_t valid_bases = 0;
      for (size_t j = 0; j < chunk_len; j++) {
        int8_t lv = BASE_LOOKUP[(uint8_t)ks->seq.s[i + j]];
        if (lv >= 0)
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

      char namebuf[512];
      snprintf(namebuf, sizeof(namebuf), "%s_%zu_%zu", ks->name.s, i,
               i + chunk_len);
      sk->name = strdup(namebuf);
      sk->kmer_size = r->hash_window;
      sk->hash_threshold = UINT64_MAX / params->scale;

      wc->chrom = strdup(ks->name.s);
      wc->start = i;
      wc->end = i + chunk_len;

      if (bed_fp) {
        fprintf(bed_fp, "%s\t%zu\t%zu\t%s\n", ks->name.s, i, i + chunk_len,
                sk->name);
      }

      HashPool pool;
      pool_init(&pool, sk->hash_threshold);
      extract_and_insert(r, &pool, (const uint8_t *)ks->seq.s + i, chunk_len, 0,
                         NULL, 0);
      pool_finalize(&pool, &sk->hashes, &sk->sketch_size);

      if (sk->sketch_size > 0) {
        (*num_sketches)++;
      } else {
        free(sk->name);
        free(wc->chrom);
      }
    }

    /* Handle trailing region if it's large enough */
    size_t last_start =
        (len >= window_size)
            ? ((len - window_size) / step_size) * step_size + step_size
            : 0;
    if (last_start + window_size <= len)
      ; /* already covered */
    else if (len > window_size && last_start < len - window_size + step_size) {
      /* trailing chunk: last window_size bases of the contig */
      size_t trail_start = len - window_size;
      if (trail_start > 0 && trail_start != last_start - step_size) {
        /* skip if it overlaps too much with the previous window */
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
  if (ea->hash < eb->hash)
    return -1;
  if (ea->hash > eb->hash)
    return 1;
  if (ea->window_id < eb->window_id)
    return -1;
  if (ea->window_id > eb->window_id)
    return 1;
  return 0;
}

/* A candidate pair with accumulated shared hash count */
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
      entries[idx].hash = sketches[i].hashes[j];
      entries[idx].window_id = (uint32_t)i;
      idx++;
    }
  }

  /* 2. Sort by hash value */
  qsort(entries, total_entries, sizeof(HashWindowEntry), cmp_hash_window_entry);

  /* 3. For each group of entries with the same hash, collect window pairs.
   *    Use a hash map approximation: accumulate (win_a, win_b) -> shared_count
   *    via a sorted array of CandidatePairs. */
  size_t n_candidates = 0;
  size_t cap_candidates = 0;
  CandidatePair *candidates = NULL;

  size_t run_start = 0;
  while (run_start < total_entries) {
    size_t run_end = run_start + 1;
    while (run_end < total_entries &&
           entries[run_end].hash == entries[run_start].hash)
      run_end++;

    size_t run_len = run_end - run_start;
    /* Skip hashes that appear in too many windows (likely repetitive noise) */
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

          if (n_candidates >= cap_candidates) {
            cap_candidates = cap_candidates ? cap_candidates * 2 : 4096;
            candidates =
                realloc(candidates, cap_candidates * sizeof(CandidatePair));
          }
          candidates[n_candidates++] = (CandidatePair){a, b, 1};
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

  /* 4. Sort candidates and merge duplicates to get shared_count */
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
  size_t n_edges = 0;
  size_t cap_edges = 1024;
  ReverbDupEdge *edges = malloc(cap_edges * sizeof(ReverbDupEdge));

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

    /* Require at least 2 shared hashes to proceed */
    if (candidates[i].shared_count < 2)
      continue;

    ReverbDistResult d = reverb_dist(&sketches[a], &sketches[b]);
    if (d.distance < max_dist) {
      if (n_edges >= cap_edges) {
        cap_edges *= 2;
        edges = realloc(edges, cap_edges * sizeof(ReverbDupEdge));
      }
      edges[n_edges++] = (ReverbDupEdge){a, b, d.distance};
    }
  }

  free(candidates);
  *out_edges = edges;
  return n_edges;
}

// ==============================================================
// DUP: SEGMENT MERGE
// ==============================================================

/* Compare ReverbDupRegion by (chrom, start) for sorting */
static int cmp_dup_region(const void *a, const void *b) {
  const ReverbDupRegion *ra = (const ReverbDupRegion *)a;
  const ReverbDupRegion *rb = (const ReverbDupRegion *)b;
  int c = strcmp(ra->chrom, rb->chrom);
  if (c != 0)
    return c;
  if (ra->start < rb->start)
    return -1;
  if (ra->start > rb->start)
    return 1;
  return 0;
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
      /* Merge: extend end, update avg distance */
      if (regions[i].end > regions[out].end)
        regions[out].end = regions[i].end;
      regions[out].avg_distance =
          (regions[out].avg_distance + regions[i].avg_distance) / 2.0;
    } else {
      out++;
      if (out != i) {
        regions[out] = regions[i];
      }
    }
  }
  return out + 1;
}

// ==============================================================
// CMD: DUP (core feature)
// ==============================================================
int cmd_dup(int argc, char **argv) {
  SketchBuildParams def = sketch_build_defaults();
  def.scale = 10;
  def.hash_seed = 42;
  size_t window_size = 10000;
  size_t step_size = 0; /* 0 = auto (window/2) */
  size_t min_bases = 1000;
  double max_dist = 0.05;
  int min_copy = 2;
  int max_copy = 50; /* 0 = unlimited */
  const char *out_prefix = "reverb";
  int n_threads = 8;

  ketopt_t opt = KETOPT_INIT;
  int c;
  while ((c = ketopt(&opt, argc - 1, argv + 1, 1,
                     "k:s:e:w:t:b:d:m:M:o:p:", 0)) >= 0) {
    if (c == 'k')
      def.kmer_size = (uint32_t)atoi(opt.arg);
    else if (c == 's')
      def.scale = (uint64_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'e')
      def.hash_seed = strtoull(opt.arg, NULL, 0);
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
  (void)n_threads; /* reserved for future parallelization */

  int num_files = argc - (opt.ind + 1);
  if (num_files < 1) {
    fprintf(stderr, "Error: missing input FASTA files\n");
    return 1;
  }
  char **in_files = argv + opt.ind + 1;

  struct timespec t_start, t_end;
  clock_gettime(CLOCK_MONOTONIC, &t_start);

  Reverb r;
  reverb_init(&r, def.kmer_size);
  r.hash_seed = def.hash_seed;

  ReverbSketch *sketches = NULL;
  WindowCoord *coords = NULL;
  size_t num_sketches = 0;
  size_t cap_sketches = 0;

  /* Open window BED output */
  char path_buf[PATH_MAX];
  snprintf(path_buf, sizeof(path_buf), "%s.window.bed", out_prefix);
  FILE *bed_fp = fopen(path_buf, "w");
  if (!bed_fp) {
    fprintf(stderr, "Error: cannot open %s for writing\n", path_buf);
    return 1;
  }

  /* Phase 1: Extract sliding window sketches */
  fprintf(stderr, "[reverb] Extracting windows (w=%zu, step=%zu) ...\n",
          window_size, step_size);
  for (int i = 0; i < num_files; i++) {
    dup_stream(in_files[i], &r, &def, window_size, step_size, min_bases,
               &sketches, &coords, &num_sketches, &cap_sketches, bed_fp);
  }
  fclose(bed_fp);

  if (num_sketches == 0) {
    fprintf(stderr, "Error: no valid windows extracted.\n");
    free(sketches);
    free(coords);
    return 1;
  }
  fprintf(stderr, "[reverb] Total windows: %zu\n", num_sketches);

  /* Phase 2: Build inverted hash index and find candidate edges */
  fprintf(stderr, "[reverb] Building hash index and finding candidates ...\n");
  ReverbDupEdge *edges = NULL;
  size_t n_edges = build_candidate_edges(sketches, coords, num_sketches,
                                         max_dist, window_size, &edges);
  fprintf(stderr, "[reverb] Duplicated edges: %zu\n", n_edges);

  if (n_edges == 0) {
    fprintf(stderr, "[reverb] No segmental duplications found.\n");
    goto cleanup;
  }

  /* Phase 3: Union-Find clustering */
  fprintf(stderr, "[reverb] Clustering SD families ...\n");
  UnionFind uf;
  uf_init(&uf, num_sketches);

  for (size_t i = 0; i < n_edges; i++) {
    uf_union(&uf, edges[i].win_a, edges[i].win_b);
  }

  /* Count component sizes */
  uint32_t *comp_size = calloc(num_sketches, sizeof(uint32_t));
  for (size_t i = 0; i < num_sketches; i++)
    comp_size[uf_find(&uf, (uint32_t)i)]++;

  /* Phase 4: Write BEDPE output (pairs) */
  snprintf(path_buf, sizeof(path_buf), "%s.dup.bedpe", out_prefix);
  FILE *bedpe_fp = fopen(path_buf, "w");
  if (!bedpe_fp) {
    fprintf(stderr, "Error: cannot open %s\n", path_buf);
    uf_free(&uf);
    free(comp_size);
    goto cleanup;
  }

  fprintf(bedpe_fp, "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tfamily\tdista"
                    "nce\tcopy_count\n");
  size_t reported_edges = 0;
  for (size_t i = 0; i < n_edges; i++) {
    uint32_t a = edges[i].win_a;
    uint32_t b = edges[i].win_b;
    uint32_t fam = uf_find(&uf, a);
    uint32_t cc = comp_size[fam];
    if ((int)cc < min_copy)
      continue;
    if (max_copy > 0 && (int)cc > max_copy)
      continue;

    fprintf(bedpe_fp, "%s\t%zu\t%zu\t%s\t%zu\t%zu\tSD_%u\t%.6f\t%u\n",
            coords[a].chrom, coords[a].start, coords[a].end, coords[b].chrom,
            coords[b].start, coords[b].end, fam, edges[i].distance, cc);
    reported_edges++;
  }
  fclose(bedpe_fp);

  /* Phase 5: Build merged BED regions */
  snprintf(path_buf, sizeof(path_buf), "%s.dup.bed", out_prefix);
  FILE *dup_bed_fp = fopen(path_buf, "w");
  if (!dup_bed_fp) {
    fprintf(stderr, "Error: cannot open %s\n", path_buf);
    uf_free(&uf);
    free(comp_size);
    goto cleanup;
  }

  /* Collect all windows that belong to families >= min_copy */
  size_t n_dup_regions = 0;
  size_t cap_dup_regions = 256;
  ReverbDupRegion *dup_regions =
      malloc(cap_dup_regions * sizeof(ReverbDupRegion));

  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t fam = uf_find(&uf, (uint32_t)i);
    if ((int)comp_size[fam] < min_copy)
      continue;
    if (max_copy > 0 && (int)comp_size[fam] > max_copy)
      continue;

    if (n_dup_regions >= cap_dup_regions) {
      cap_dup_regions *= 2;
      dup_regions =
          realloc(dup_regions, cap_dup_regions * sizeof(ReverbDupRegion));
    }
    dup_regions[n_dup_regions++] =
        (ReverbDupRegion){.chrom = coords[i].chrom, /* borrowed pointer */
                          .start = coords[i].start,
                          .end = coords[i].end,
                          .family_id = fam,
                          .copy_count = comp_size[fam],
                          .avg_distance = 0.0};
  }

  size_t n_merged = merge_dup_regions(dup_regions, n_dup_regions);

  fprintf(dup_bed_fp, "#chrom\tstart\tend\tfamily\tcopy_count\n");
  size_t n_families = 0;
  uint32_t last_fam = UINT32_MAX;
  for (size_t i = 0; i < n_merged; i++) {
    fprintf(dup_bed_fp, "%s\t%zu\t%zu\tSD_%u\t%u\n", dup_regions[i].chrom,
            dup_regions[i].start, dup_regions[i].end, dup_regions[i].family_id,
            dup_regions[i].copy_count);
    if (dup_regions[i].family_id != last_fam) {
      n_families++;
      last_fam = dup_regions[i].family_id;
    }
  }
  fclose(dup_bed_fp);

  clock_gettime(CLOCK_MONOTONIC, &t_end);
  double elapsed =
      (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;

  /* Summary */
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

  free(dup_regions);
  uf_free(&uf);
  free(comp_size);

cleanup:
  free(edges);
  for (size_t i = 0; i < num_sketches; i++) {
    reverb_sketch_free(&sketches[i]);
    free(coords[i].chrom);
  }
  free(sketches);
  free(coords);

  return 0;
}

// ==============================================================
// MAIN
// ==============================================================

static int is_known_command(const char *cmd) {
  return strcmp(cmd, "dup") == 0 || strcmp(cmd, "sketch") == 0 ||
         strcmp(cmd, "dist") == 0 || strcmp(cmd, "three") == 0 ||
         strcmp(cmd, "triangle") == 0 || strcmp(cmd, "info") == 0 ||
         strcmp(cmd, "help") == 0;
}

int main(int argc, char **argv) {
  if (argc < 2 || strcmp(argv[1], "help") == 0 || strcmp(argv[1], "-h") == 0 ||
      strcmp(argv[1], "--help") == 0) {
    print_usage();
    return argc < 2 ? 1 : 0;
  }

  const char *cmd = argv[1];

  /* If the first argument is not a known command, treat everything
   * as arguments to the default 'dup' command. This handles both:
   *   reverb genome.fasta          (file as first arg)
   *   reverb -o prefix genome.fa   (options as first arg) */
  if (!is_known_command(cmd)) {
    return cmd_dup(argc, argv);
  }

  if (strcmp(cmd, "dup") == 0)
    return cmd_dup(argc, argv);
  if (strcmp(cmd, "sketch") == 0)
    return cmd_sketch(argc, argv);
  if (strcmp(cmd, "dist") == 0)
    return cmd_dist(argc, argv);
  if (strcmp(cmd, "three") == 0)
    return cmd_three(argc, argv);
  if (strcmp(cmd, "triangle") == 0)
    return cmd_triangle(argc, argv);
  if (strcmp(cmd, "info") == 0)
    return cmd_info(argc, argv);

  fprintf(stderr, "Error: unknown command '%s'\n", cmd);
  print_usage();
  return 1;
}
