#ifndef DP_H

#include "polygon.h"

typedef struct {
	int begin;
	int end;
	char is_problem;
} segment_t;

typedef struct {
	segment_t *segs;
	int num_segs;
} reduced_ring_t;

mpoly_t compute_reduced_pointset(mpoly_t *in_mpoly, double tolerance);
reduced_ring_t compute_reduced_ring(ring_t *orig_string, double res);
void fix_topology(mpoly_t *mpoly, reduced_ring_t *reduced_rings);
mpoly_t reduction_to_mpoly(mpoly_t *in_mpoly, reduced_ring_t *reduced_rings);

#endif // ifndef DP_H
