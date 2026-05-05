#ifndef PARABOLA_H
#define PARABOLA_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  uint32_t hash_window;
  __uint128_t remover_mask;
  __uint128_t layer_mask;
  uint32_t kmer_bits;
  uint32_t rc_shift;
  uint64_t hash_seed;
} Parabola;

typedef struct {
  uint64_t hash;
  uint32_t simplicity;
} HSPair;

typedef struct {
  size_t size;
  size_t cap;
  uint64_t hash_threshold;
  uint64_t *hashes;
  uint32_t *simps;
} HashPool;

typedef struct {
  const uint8_t *read_sequence;
  size_t sequence_length;
  const Parabola *parameters;
  size_t current_index;
  __uint128_t forward_reverb;
  size_t valid_length;
} ParabolaIterator;

typedef struct {
  char *name;
  uint32_t kmer_size;
  size_t sketch_size;
  uint64_t hash_threshold;
  uint64_t *hashes;
  uint32_t *simplicities;
} ParabolaSketch;

typedef struct {
  double distance;
  double jaccard;
  size_t shared_hashes;
  size_t total_hashes;
} ParabolaDistResult;

typedef struct {
  ParabolaDistResult d01;
  ParabolaDistResult d02;
  ParabolaDistResult d12;
} ParabolaTripleDistResult;

void parabola_init(Parabola *p, size_t hash_window);
void kt_for(int n_threads, void (*func)(void *, long, int), void *data, long n);
void parabola_sketch_create(ParabolaSketch *sk, const char *name,
                            uint32_t kmer_size, size_t target_size,
                            uint64_t hash_seed, const char **files,
                            int num_files, int min_count, int num_threads);
void parabola_stream(const char *filename, const Parabola *p, HashPool *pool,
                     int min_count, uint8_t *cms_table, size_t cms_mask);
void parabola_sketch_free(ParabolaSketch *sk);
ParabolaDistResult parabola_dist(const ParabolaSketch *ref,
                                 const ParabolaSketch *query, int use_jc);
ParabolaTripleDistResult parabola_dist_three(const ParabolaSketch *ref,
                                             const ParabolaSketch *q1,
                                             const ParabolaSketch *q2,
                                             int use_jc);
void parabola_info(const ParabolaSketch *sk);
int parabola_sketch_save(const ParabolaSketch *sk, const char *filepath);
int parabola_sketch_load(ParabolaSketch *sk, const char *filepath);

#ifdef __cplusplus
}
#endif

#endif