#ifndef REVERB_H
#define REVERB_H

#include <stddef.h>
#include <stdint.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  uint32_t hash_window;
  __uint128_t remover_mask;
  uint32_t kmer_bits;
  uint32_t rc_shift;
  uint64_t hash_seed;
} Reverb;

typedef struct {
  size_t sketch_size;
  uint64_t *hashes;
} ReverbSketch;

typedef struct {
  double containment;
  double distance;
  size_t shared_hashes;
} ReverbDistResult;

/* (HashWindowEntry is defined internally in reverb.c) */

/* Edge in the duplication graph */
typedef struct {
  uint32_t win_a;
  uint32_t win_b;
  double distance;
} ReverbDupEdge;

typedef struct {
  uint32_t *parent;
  uint32_t *rank;
  size_t n;
} UnionFind;

/* Merged SD region */
typedef struct {
  char *chrom;
  size_t start;
  size_t end;
  char *cluster_id;
  uint32_t copy_count;
  double avg_distance;
  uint32_t subcluster_id;
  ReverbSketch flank_sketch;
} ReverbDupRegion;

void reverb_init(Reverb *r, size_t hash_window);
void reverb_sketch_free(ReverbSketch *sk);
ReverbDistResult reverb_dist(const ReverbSketch *ref,
                             const ReverbSketch *query, uint32_t kmer_size);

/* Union-Find operations */
void uf_init(UnionFind *uf, size_t n);
uint32_t uf_find(UnionFind *uf, uint32_t x);
void uf_union(UnionFind *uf, uint32_t a, uint32_t b);
void uf_free(UnionFind *uf);

#ifdef __cplusplus
}
#endif

#endif
