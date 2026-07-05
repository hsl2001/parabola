#ifndef PARABOLA_H
#define PARABOLA_H

#include <stddef.h>
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
  char *name;
  uint32_t kmer_size;
  size_t sketch_size;
  uint64_t hash_threshold;
  uint64_t *hashes;
} ParabolaSketch;

typedef struct {
  double containment;
  double distance;
  size_t shared_hashes;
} ParabolaDistResult;

typedef struct {
  ParabolaDistResult d01;
  ParabolaDistResult d02;
  ParabolaDistResult d12;
} ParabolaTripleDistResult;

void parabola_init(Parabola *p, size_t hash_window);
void kt_for(int n_threads, void (*func)(void *, long, int), void *data, long n);
void parabola_sketch_free(ParabolaSketch *sk);
ParabolaDistResult parabola_dist(const ParabolaSketch *ref,
                                 const ParabolaSketch *query);
ParabolaTripleDistResult parabola_dist_three(const ParabolaSketch *ref,
                                             const ParabolaSketch *q1,
                                             const ParabolaSketch *q2);
void parabola_info(const ParabolaSketch *sk);
int parabola_sketch_save(const ParabolaSketch *sk, const char *filepath);
int parabola_sketch_load(ParabolaSketch *sk, const char *filepath);

#ifdef __cplusplus
}
#endif

#endif