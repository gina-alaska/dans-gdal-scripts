/*
Copyright (c) 2007, Regents of the University of Alaska

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the name of the Geographic Information Network of Alaska nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/



#include "common.h"
#include "polygon.h"

#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a)<(b)?(a):(b))

#ifdef MAX
#undef MAX
#endif
#define MAX(a,b) ((a)>(b)?(a):(b))

#ifdef SGN
#undef SGN
#endif
#define SGN(a) ((a)<0?-1:(a)>0?1:0)

typedef struct {
	double x, y;
	int ring_idx, vert_idx;
} vertsort_entry;

static int compare_verts(const void *ap, const void *bp) {
	vertsort_entry *a = (vertsort_entry *)ap;
	vertsort_entry *b = (vertsort_entry *)bp;
	return 
		a->x < b->x ? -1 :
		a->x > b->x ?  1 :
		a->y < b->y ? -1 :
		a->y > b->y ?  1 :
		0;
}

// This function is only meant to be called on polygons
// that have orthogonal sides on an integer lattice.
void bevel_self_intersections(mpoly_t *mp, double amount) {
	int total_pts = 0;
	for(int i=0; i<mp->num_rings; i++) {
		total_pts += mp->rings[i].npts;
	}

	printf("allocating %d megs for beveler\n",
		(int)(total_pts*sizeof(vertsort_entry)) >> 20);
	vertsort_entry *entries = (vertsort_entry *)malloc_or_die(
		sizeof(vertsort_entry) * total_pts);
	int entry_idx = 0;
	for(int r_idx=0; r_idx<mp->num_rings; r_idx++) {
		ring_t *ring = mp->rings + r_idx;
		for(int v_idx=0; v_idx<ring->npts; v_idx++) {
			vertsort_entry *entry = entries + (entry_idx++);
			entry->x = ring->pts[v_idx].x;
			entry->y = ring->pts[v_idx].y;
			entry->ring_idx = r_idx;
			entry->vert_idx = v_idx;
		}
	}

	printf("qsort\n");
	qsort(entries, total_pts, sizeof(vertsort_entry), compare_verts);

	printf("collecting\n");
	int total_num_touch = 0;
	int prev_was_same = 0;
	for(int i=0; i<total_pts-1; i++) {
		//printf("%lf %lf %d %d\n",
		//	entries[i].x, entries[i].y,
		//	entries[i].ring_idx, entries[i].vert_idx);
		if(
			entries[i].x == entries[i+1].x &&
			entries[i].y == entries[i+1].y
		) {
			if(prev_was_same) {
				fatal_error("should not have triple intersections in beveler");
			}
			entries[total_num_touch++] = entries[i];
			prev_was_same = 1;
		} else {
			prev_was_same = 0;
		}
	}

	printf("found %d self-intersections\n", total_num_touch);
	if(!total_num_touch) return;

	entries = (vertsort_entry *)realloc_or_die(entries, sizeof(vertsort_entry) * total_num_touch);
	int *ring_num_touch = (int *)malloc_or_die(sizeof(int) * mp->num_rings);
	for(int i=0; i<mp->num_rings; i++) ring_num_touch[i] = 0;
	for(int i=0; i<total_num_touch; i++) {
		ring_num_touch[entries[i].ring_idx]++;
	}

	for(int r_idx=0; r_idx<mp->num_rings; r_idx++) {
		if(!ring_num_touch[r_idx]) continue;
		ring_t *ring = mp->rings + r_idx;

		if(VERBOSE >= 2) printf("ring %d: num_touch=%d\n", r_idx, ring_num_touch[r_idx]);

		char *touch_mask = (char *)malloc_or_die(ring->npts);
		memset(touch_mask, 0, ring->npts);
		for(int i=0; i<total_num_touch; i++) {
			if(entries[i].ring_idx == r_idx) {
				touch_mask[entries[i].vert_idx] = 1;
			}
		}

		int new_numpts = ring->npts + ring_num_touch[r_idx];
		vertex_t *new_pts = (vertex_t *)malloc_or_die(sizeof(vertex_t) * new_numpts);
		int vout_idx = 0;
		for(int v_idx=0; v_idx<ring->npts; v_idx++) {
			if(touch_mask[v_idx]) {
				vertex_t this_v = ring->pts[v_idx];
				vertex_t prev_v = ring->pts[(v_idx+ring->npts-1) % ring->npts];
				vertex_t next_v = ring->pts[(v_idx+1) % ring->npts];
				double sx_prev = SGN(prev_v.x - this_v.x);
				double sy_prev = SGN(prev_v.y - this_v.y);
				double sx_next = SGN(next_v.x - this_v.x);
				double sy_next = SGN(next_v.y - this_v.y);
				new_pts[vout_idx++] = (vertex_t){ this_v.x + sx_prev*amount, this_v.y + sy_prev*amount };
				new_pts[vout_idx++] = (vertex_t){ this_v.x + sx_next*amount, this_v.y + sy_next*amount };
			} else {
				new_pts[vout_idx++] = ring->pts[v_idx];
			}
		}
		if(vout_idx != new_numpts) fatal_error("wrong number of points in beveled ring");
		free(ring->pts);
		ring->npts = new_numpts;
		ring->pts = new_pts;
		free(touch_mask);
	}

	free(ring_num_touch);
	free(entries);

	printf("beveler finish\n");
}
