#include <math.h>
#include <stdio.h>
#include <zlib.h>

#include "klib/ketopt.h"
#include "klib/kseq.h"
#include "parabola.h"

#define MIX_CONST1 0xff51afd7ed558ccdULL
#define MIX_CONST2 0xc4ceb9fe1a85ec53ULL

KSEQ_INIT(gzFile, gzread)

// ==============================================================
// UTILITIES
// ==============================================================
/* A = 001, C = 110, G = 011, T = 100 */
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
// UTILITIES
// ==============================================================
void parabola_init(Parabola *p, size_t hash_window) {
  size_t k = min_size_t(hash_window, 42); /* 3*42 = 126 bits fits __uint128_t */
  uint32_t kmer_bits = 3 * (uint32_t)k;

  __uint128_t remover_mask =
      (kmer_bits > 3) ? (((__uint128_t)1 << (kmer_bits - 3)) - 1) : 0;
  __uint128_t layer_mask = 0;
  for (uint32_t j = 0; j < (uint32_t)k; j++)
    layer_mask |= (__uint128_t)1 << (3 * j);

  p->hash_window = k;
  p->remover_mask = remover_mask;
  p->layer_mask = layer_mask;
  p->kmer_bits = kmer_bits;
  p->rc_shift = (kmer_bits > 0) ? (128 - kmer_bits) : 128;
}

// ==============================================================
// OVERFLOWING MAX HEAP
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

static void pool_finalize(HashPool *pool, uint64_t **out_hashes, size_t *out_size) {
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
__attribute__((hot))
static void extract_and_insert(const Parabola *p, HashPool *pool,
                               const uint8_t *seq, size_t len, int min_count,
                               uint8_t *cms_table, size_t cms_mask) {
  size_t K = p->hash_window;
  __uint128_t fwd = 0;
  size_t valid = 0;

  for (size_t idx = 0; idx < len; idx++) {
    int8_t lv = BASE_LOOKUP[seq[idx]];
    if (lv < 0) {
      fwd = 0;
      valid = 0;
      continue;
    }

    fwd = ((fwd & p->remover_mask) << 3) | (uint8_t)lv;
    if (valid < K)
      valid++;

    if (valid >= K) {
      __uint128_t rev = reverse_bits128(fwd) >> p->rc_shift;
      __uint128_t canon = fwd < rev ? fwd : rev;
      uint64_t h = mix_hash(canon, p->hash_seed);

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
void parabola_stream(const char *filename, const Parabola *p, HashPool *pool,
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
    if (ks->seq.l < p->hash_window)
      continue;
    extract_and_insert(p, pool, (const uint8_t *)ks->seq.s, ks->seq.l,
                       min_count, cms_table, cms_mask);
  }
  kseq_destroy(ks);
  gzclose(fp);
}

void parabola_sketch_free(ParabolaSketch *sk) {
  if (sk) {
    free(sk->name);
    free(sk->hashes);
  }
}

// ==============================================================
// DISTANCE CALCULATION
// ==============================================================
ParabolaDistResult parabola_dist(const ParabolaSketch *ref,
                                 const ParabolaSketch *query) {
  ParabolaDistResult res = {0.0, 1.0, 0};
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

ParabolaTripleDistResult parabola_dist_three(const ParabolaSketch *ref,
                                             const ParabolaSketch *q1,
                                             const ParabolaSketch *q2) {
  ParabolaTripleDistResult t;
  t.d01 = t.d02 = t.d12 = (ParabolaDistResult){0.0, 1.0, 0};

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

    if (inR && inQ1) sh01++;
    if (inR && inQ2) sh02++;
    if (inQ1 && inQ2) sh12++;

    if (inR) i++;
    if (inQ1) j++;
    if (inQ2) k++;
  }

  double inv_k = 1.0 / (double)ref->kmer_size;

  #define FILL(r, sh, a, b) do { \
    (r).shared_hashes = (sh); \
    if ((a)->sketch_size > 0 && (b)->sketch_size > 0) { \
      double ca = (double)(sh) / (double)(a)->sketch_size; \
      double cb = (double)(sh) / (double)(b)->sketch_size; \
      (r).containment = 0.5 * (ca + cb); \
      (r).distance = 1.0 - pow((r).containment, inv_k); \
    } \
  } while(0)

  FILL(t.d01, sh01, ref, q1);
  FILL(t.d02, sh02, ref, q2);
  FILL(t.d12, sh12, q1, q2);
  #undef FILL

  return t;
}

void parabola_info(const ParabolaSketch *sk) {
  if (!sk)
    return;
  printf("--- Parabola Sketch Info ---\n");
  printf("Name                 : %s\n", sk->name ? sk->name : "N/A");
  printf("K-mer Size (K)       : %u\n", sk->kmer_size);
  printf("Total Sketch Size    : %zu\n", sk->sketch_size);
  printf("----------------------------\n");
}

// ==============================================================
// SERIALIZATION
// ==============================================================
int parabola_sketch_save(const ParabolaSketch *sk, const char *filepath) {
  FILE *f = fopen(filepath, "wb");
  if (!f)
    return -1;

  uint32_t name_len = sk->name ? (uint32_t)strlen(sk->name) : 0;
  uint64_t sz = sk->sketch_size;

  fwrite("PARA", 1, 4, f);
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

int parabola_sketch_load(ParabolaSketch *sk, const char *filepath) {
  FILE *f = fopen(filepath, "rb");
  if (!f)
    return -1;

  char magic[4];
  if (fread(magic, 1, 4, f) != 4 || memcmp(magic, "PARA", 4) != 0) {
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

static int build_sketch_from_fasta(ParabolaSketch *sk, const char *path,
                                   const SketchBuildParams *params) {
  memset(sk, 0, sizeof(ParabolaSketch));
  sk->name = strdup(path);
  sk->kmer_size = params->kmer_size;
  sk->hash_threshold = UINT64_MAX / params->scale;

  if (!params->scale)
    return -1;

  Parabola p;
  parabola_init(&p, params->kmer_size);
  p.hash_seed = params->hash_seed;

  HashPool pool;
  pool_init(&pool, sk->hash_threshold);
  parabola_stream(path, &p, &pool, params->min_count, NULL, 0);
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

static int prepare_input_sketches(ParabolaSketch *sketches, int *ready,
                                  char *const *paths, int num_files,
                                  SketchBuildParams *params, int *failed) {
  if (failed)
    *failed = -1;
  int params_seeded = 0;

  for (int i = 0; i < num_files; i++) {
    ready[i] = 0;

    if (!is_sequence_file(paths[i])) {
      if (parabola_sketch_load(&sketches[i], paths[i]) == 0) {
        ready[i] = 1;
        if (!params_seeded) {
          params->kmer_size = sketches[i].kmer_size;
          params_seeded = 1;
        }
      } else if (has_extension(paths[i], ".parabola")) {
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
static void free_ready_sketches(ParabolaSketch *sketches, const int *ready,
                                int num_files) {
  for (int i = 0; i < num_files; i++)
    if (ready[i])
      parabola_sketch_free(&sketches[i]);
}

// ==============================================================
// PRINT REPORTS
// ==============================================================

static void print_distance_report(const ParabolaSketch *ref,
                                  const ParabolaSketch *query,
                                  const ParabolaDistResult *dist) {
  printf("Reference File      : %s\n", ref->name);
  printf("Query File          : %s\n", query->name);
  printf("Shared Hashes       : %zu\n", dist->shared_hashes);
  printf("------------------------------------\n");
  printf("Average Containment : %f\n", dist->containment);
  printf("Parabola Distance   : %f\n", dist->distance);
}

static void print_usage(void) {
  printf("Parabola: Fast Genomic Sketching Tool\n\n"
         "Usage: parabola <command> [options] [arguments]\n"
         "Commands:\n"
         "  help | -h | --help\n"
         "  sketch [-k K] [-s S|N/M] [-e E] [-p threads] [-r] [-m min_count] "
         "[-o out_file] fasta1 [fasta2 ...]\n"
         "         -k: kmer size (default: 21, max: 42)\n"
         "         -s: scale factor (keep 1/S of k-mers, default: 1000)\n"
         "         -e: hash seed\n"
         "         -p: number of threads (default: 1)\n"
         "         -r: pool all reads from input files into one sketch\n"
         "         -m: minimum k-mer count filter (used with -r, default: 2)\n"
         "         -o: output file path (required with -r; optional for single "
         "input)\n"
         "  dist   <sketch|fasta1> <sketch|fasta2>\n"
         "  three  <ref> <query1> <query2>\n"
         "  triangle <sketch|fasta1> ... <sketch|fastaN>\n"
         "  dup    [-k K] [-s S] [-w win] [-b min_bases] [-d identity] [-y "
         "min_copy] "
         "fasta1 [fasta2 ...]\n"
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

  ParabolaSketch sk = {0};
  sk.name = job->in_file ? strdup(job->in_file) : NULL;
  sk.kmer_size = job->kmer_size;
  sk.hash_threshold = job->hash_threshold;

  if (job->hash_threshold) {
    Parabola p;
    parabola_init(&p, job->kmer_size);
    p.hash_seed = job->hash_seed;

    HashPool pool;
    pool_init(&pool, job->hash_threshold);
    parabola_stream(job->in_file, &p, &pool, job->min_count, NULL, 0);
    pool_finalize(&pool, &sk.hashes, &sk.sketch_size);
  }

  job->result = parabola_sketch_save(&sk, job->out_file);
  parabola_sketch_free(&sk);
}

// ==============================================================
// CMD
// ==============================================================
static int setup_dist_cmds(int argc, char **argv, int min_files,
                           ParabolaSketch **sk, int **ready) {
  ketopt_t opt = KETOPT_INIT;
  int c;

  while ((c = ketopt(&opt, argc - 1, argv + 1, 1, "", 0)) >= 0) {
    return -1; // unknown or error
  }

  int num_files = argc - (opt.ind + 1);
  if (num_files < min_files) {
    fprintf(stderr, "Error: requires at least %d input file(s)\n", min_files);
    return -1;
  }

  *sk = (ParabolaSketch *)calloc(num_files, sizeof(ParabolaSketch));
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

// ==============================================================
// CMD
// ==============================================================
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
    out_opt
        ? (long)strcpy((char *)jobs[i].out_file, out_opt)
        : snprintf((char *)jobs[i].out_file, len, "%s.parabola", in_files[i]);
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

int cmd_dist(int argc, char **argv) {
  ParabolaSketch *sk;
  int *ready;

  if (setup_dist_cmds(argc, argv, 2, &sk, &ready) < 2)
    return 1;

  ParabolaDistResult dist = parabola_dist(&sk[0], &sk[1]);
  print_distance_report(&sk[0], &sk[1], &dist);

  free_ready_sketches(sk, ready, 2);
  free(sk);
  free(ready);
  return 0;
}

int cmd_three(int argc, char **argv) {
  ParabolaSketch *sk;
  int *ready;

  if (setup_dist_cmds(argc, argv, 3, &sk, &ready) < 3)
    return 1;

  ParabolaTripleDistResult tri = parabola_dist_three(&sk[0], &sk[1], &sk[2]);
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

static int dup_stream(const char *filename, const Parabola *p,
                      SketchBuildParams *params, size_t window_size,
                      size_t min_bases, ParabolaSketch **sketches,
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
    for (size_t i = 0; i < len; i += window_size) {
      size_t chunk_len = window_size;
      if (i + chunk_len > len)
        chunk_len = len - i;

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
        *sketches = realloc(*sketches, *cap_sketches * sizeof(ParabolaSketch));
      }

      ParabolaSketch *sk = &(*sketches)[*num_sketches];
      memset(sk, 0, sizeof(ParabolaSketch));

      char namebuf[512];
      snprintf(namebuf, sizeof(namebuf), "%s_%zu_%zu", ks->name.s, i,
               i + chunk_len);
      sk->name = strdup(namebuf);
      sk->kmer_size = p->hash_window;
      sk->hash_threshold = UINT64_MAX / params->scale;

      if (bed_fp) {
        fprintf(bed_fp, "%s\t%zu\t%zu\t%s\n", ks->name.s, i, i + chunk_len,
                sk->name);
      }

      HashPool pool;
      pool_init(&pool, sk->hash_threshold);
      extract_and_insert(p, &pool, (const uint8_t *)ks->seq.s + i, chunk_len, 0,
                         NULL, 0);
      pool_finalize(&pool, &sk->hashes, &sk->sketch_size);

      if (sk->sketch_size > 0) {
        (*num_sketches)++;
      } else {
        parabola_sketch_free(sk);
      }
    }
  }
  kseq_destroy(ks);
  gzclose(fp);
  return 0;
}

typedef struct {
  ParabolaSketch *sketches;
  int *is_dup;
  int *copy_count;
  int **copy_indices;
  int *notcopy_count;
  double L;
  double H;
} DupWorkerData;

static void dup_dist_worker(void *data, long i, int _unused) {
  (void)_unused;
  DupWorkerData *wd = (DupWorkerData *)data;
  for (size_t j = 0; j < (size_t)i; j++) {
    if (!__atomic_load_n(&wd->is_dup[i], __ATOMIC_RELAXED) &&
        !__atomic_load_n(&wd->is_dup[j], __ATOMIC_RELAXED))
      continue;
    ParabolaDistResult d = parabola_dist(&wd->sketches[i], &wd->sketches[j]);
    double dist = d.distance;
    if (dist >= wd->L && dist <= wd->H) {
      __atomic_store_n(&wd->is_dup[i], 0, __ATOMIC_RELAXED);
      __atomic_store_n(&wd->is_dup[j], 0, __ATOMIC_RELAXED);
    } else if (dist < wd->L) {
      int idx_i = __sync_fetch_and_add(&wd->copy_count[i], 1);
      wd->copy_indices[i][idx_i] = j;
      int idx_j = __sync_fetch_and_add(&wd->copy_count[j], 1);
      wd->copy_indices[j][idx_j] = i;
    } else if (dist > wd->H) {
      __sync_fetch_and_add(&wd->notcopy_count[i], 1);
      __sync_fetch_and_add(&wd->notcopy_count[j], 1);
    }
  }
}

int cmd_dup(int argc, char **argv) {
  SketchBuildParams def = sketch_build_defaults();
  def.scale = 10;
  def.hash_seed = 42;
  size_t window_size = 10000;
  size_t min_bases = 1000;
  double identity_th = 0.8;
  int min_copy = 3;

  ketopt_t opt = KETOPT_INIT;
  int c;
  while ((c = ketopt(&opt, argc - 1, argv + 1, 1, "k:s:e:w:b:d:y:", 0)) >= 0) {
    if (c == 'k')
      def.kmer_size = (uint32_t)atoi(opt.arg);
    else if (c == 's')
      def.scale = (uint64_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'e')
      def.hash_seed = strtoull(opt.arg, NULL, 0);
    else if (c == 'w')
      window_size = (size_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'b')
      min_bases = (size_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'd')
      identity_th = atof(opt.arg);
    else if (c == 'y')
      min_copy = atoi(opt.arg);
    else
      return 1;
  }

  int num_files = argc - (opt.ind + 1);
  if (num_files < 1) {
    fprintf(stderr, "Error: missing input fasta files\n");
    return 1;
  }
  char **in_files = argv + opt.ind + 1;

  double H = identity_th;
  double L = 1.0 - identity_th;

  Parabola p;
  parabola_init(&p, def.kmer_size);
  p.hash_seed = def.hash_seed;

  ParabolaSketch *sketches = NULL;
  size_t num_sketches = 0;
  size_t cap_sketches = 0;

  FILE *bed_fp = fopen("window.bed", "w");
  if (!bed_fp) {
    fprintf(stderr, "Error: cannot open window.bed for writing\n");
    return 1;
  }

  for (int i = 0; i < num_files; i++) {
    dup_stream(in_files[i], &p, &def, window_size, min_bases, &sketches,
               &num_sketches, &cap_sketches, bed_fp);
  }
  fclose(bed_fp);

  if (num_sketches == 0) {
    fprintf(stderr, "Error: no valid windows extracted.\n");
    free(sketches);
    return 1;
  }

  int *is_dup = calloc(num_sketches, sizeof(int));
  int *copy_count = calloc(num_sketches, sizeof(int));
  int **copy_indices = malloc(num_sketches * sizeof(int *));
  for (size_t i = 0; i < num_sketches; i++) {
    copy_indices[i] = malloc(num_sketches * sizeof(int));
  }
  int *notcopy_count = calloc(num_sketches, sizeof(int));
  for (size_t i = 0; i < num_sketches; i++)
    is_dup[i] = 1;

  int n_threads =
      def.num_threads > 1 ? def.num_threads : 8; // Default to 8 threads for dup
  DupWorkerData wd = {sketches,      is_dup, copy_count, copy_indices,
                      notcopy_count, L,      H};
  kt_for(n_threads, dup_dist_worker, &wd, num_sketches);

  FILE *dup_fp = fopen("dup.txt", "w");
  if (!dup_fp) {
    fprintf(stderr, "Error: cannot open dup.txt for writing\n");
    for (size_t i = 0; i < num_sketches; i++)
      parabola_sketch_free(&sketches[i]);
    free(sketches);
    free(is_dup);
    free(copy_count);
    free(notcopy_count);
    return 1;
  }

  for (size_t i = 0; i < num_sketches; i++) {
    if (is_dup[i] && copy_count[i] >= min_copy && notcopy_count[i] > 0) {
      fprintf(dup_fp, "%s\t%d\t", sketches[i].name, copy_count[i]);
      for (int c = 0; c < copy_count[i]; c++) {
        fprintf(dup_fp, "%s%s", sketches[copy_indices[i][c]].name,
                c == copy_count[i] - 1 ? "" : ",");
      }
      fprintf(dup_fp, "\n");
    }
  }
  fclose(dup_fp);

  for (size_t i = 0; i < num_sketches; i++) {
    parabola_sketch_free(&sketches[i]);
    free(copy_indices[i]);
  }
  free(copy_indices);
  free(sketches);
  free(is_dup);
  free(copy_count);
  free(notcopy_count);

  return 0;
}

int cmd_triangle(int argc, char **argv) {
  ParabolaSketch *sk;
  int *ready;

  int n = setup_dist_cmds(argc, argv, 1, &sk, &ready);
  if (n < 1)
    return 1;

  printf("%d\n", n);
  for (int i = 0; i < n; i++) {
    printf("%s", sk[i].name ? sk[i].name : "N/A");
    for (int j = 0; j < i; j++) {
      ParabolaDistResult d = parabola_dist(&sk[i], &sk[j]);
      printf("\t%f", d.distance);
    }
    printf("\n");
  }

  free_ready_sketches(sk, ready, n);
  free(sk);
  free(ready);
  return 0;
}

int cmd_info(int argc, char **argv) {
  if (argc < 3)
    return 1;

  ParabolaSketch sk;
  if (parabola_sketch_load(&sk, argv[2]) != 0)
    return 1;

  parabola_info(&sk);
  parabola_sketch_free(&sk);
  return 0;
}

// ==============================================================
// MAIN
// ==============================================================
int main(int argc, char **argv) {
  if (argc < 2 || strcmp(argv[1], "help") == 0 || strcmp(argv[1], "-h") == 0 ||
      strcmp(argv[1], "--help") == 0) {
    print_usage();
    return argc < 2 ? 1 : 0;
  }
  if (argv[1][0] == '-') {
    fprintf(stderr, "Error: unknown option '%s'\n", argv[1]);
    print_usage();
    return 1;
  }

  const char *cmd = argv[1];

  if (strcmp(cmd, "sketch") == 0)
    return cmd_sketch(argc, argv);
  if (strcmp(cmd, "dist") == 0)
    return cmd_dist(argc, argv);
  if (strcmp(cmd, "three") == 0)
    return cmd_three(argc, argv);
  if (strcmp(cmd, "triangle") == 0)
    return cmd_triangle(argc, argv);
  if (strcmp(cmd, "dup") == 0)
    return cmd_dup(argc, argv);
  if (strcmp(cmd, "info") == 0)
    return cmd_info(argc, argv);

  fprintf(stderr, "Error: unknown command '%s'\n", cmd);
  print_usage();
  return 1;
}