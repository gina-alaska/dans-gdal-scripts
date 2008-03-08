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

#if 0

/*
void pinch_self_intersections(mpoly_t *mp) {
	int r_idx;
	for(r_idx=0; r_idx<mp->num_rings; r_idx++) {
		ring_t *ring = &mp->rings[r_idx];
		int v1in_idx, v1out_idx, v2_idx;
		for(v1in_idx=0, v1out_idx=0; v1in_idx<ring->npts; v1in_idx++) {
			vertex_t *v1 = &ring->pts[v1in_idx];
			for(v2_idx=0; v2_idx<v1out_idx; v2_idx++) {
				vertex_t *v2 = &ring->pts[v2_idx];
				int touches = (v1->x == v2->x) && (v1->y == v2->y);
				if(touches) {
					fprintf("touch at ring %d vert %d vs %d : xy=(%f,%f)\n", r_idx, v1out_idx, v2_idx, v1->x, v1->y);
					mp->rings = (ring_t *)realloc_or_die(mp->rings, sizeof(ring_t)*(mp->num_rings+1));
					ring = &mp->rings[r_idx]; // must recompute this since location of mp->rings changed
					ring_t *newring = &mp->rings[mp->num_rings++];
					newring->npts = v1out_idx - v2_idx;
					newring->pts = (vertex_t *)malloc_or_die(sizeof(vertex_t) * newring->npts);
					memcpy(newring->pts, ring->pts+v2_idx, sizeof(vertex_t) * newring->npts);
					v1out_idx = v2_idx;
				}
			}
			if(v1out_idx != v1in_idx) {
				ring->pts[v1out_idx] = *v1;
			}
			v1out_idx++;
		}
		ring->npts = v1out_idx;
	}
}
*/

#define BORDER_TOUCH_HASH_SQRTSIZE 5000
#define BORDER_TOUCH_HASH_SIZE (BORDER_TOUCH_HASH_SQRTSIZE*BORDER_TOUCH_HASH_SQRTSIZE)

static inline int get_touch_hash_key(vertex_t *v) {
	int key = (int)(v->x + v->y * (double)BORDER_TOUCH_HASH_SQRTSIZE);
	key = ((key % BORDER_TOUCH_HASH_SIZE) + BORDER_TOUCH_HASH_SIZE) % BORDER_TOUCH_HASH_SIZE;
	if(key<0 || key>=BORDER_TOUCH_HASH_SIZE) fatal_error("hash key out of range");
	return key;
}

static char *mpoly_border_touch_create_hashtable(mpoly_t *mp) {
	char *table = (char *)malloc_or_die(BORDER_TOUCH_HASH_SIZE);
	memset(table, 0, BORDER_TOUCH_HASH_SIZE);
	int r_idx;
	for(r_idx=0; r_idx<mp->num_rings; r_idx++) {
		ring_t *ring = &mp->rings[r_idx];
		int v_idx;
		for(v_idx=0; v_idx<ring->npts; v_idx++) {
			vertex_t *v = &ring->pts[v_idx];
			int key = get_touch_hash_key(v);
			if(table[key] < 2) table[key]++;
		}
	}
	return table;
}

static inline int mpoly_border_touches_point(char *table, mpoly_t *mp, int r1_idx, int v1_idx, bbox_t *bboxes) {
	vertex_t *v1 = &mp->rings[r1_idx].pts[v1_idx];
	double v1x = v1->x;
	double v1y = v1->y;

	int key = get_touch_hash_key(v1);
	if(table[key] < 2) return 0;

	if(VERBOSE >= 2) printf("hash hit for %d,%d\n", r1_idx, v1_idx);

	int r2_idx;
	for(r2_idx=0; r2_idx<mp->num_rings; r2_idx++) {
		ring_t *ring = &mp->rings[r2_idx];
		if(ring->npts == 0) continue;

		bbox_t *bbox = bboxes + r2_idx;
		if(v1x < bbox->min_x || v1x > bbox->max_x) continue;
		if(v1y < bbox->min_y || v1y > bbox->max_y) continue;

		if(VERBOSE >= 2) printf("full search for %d,%d\n", r1_idx, v1_idx);

		int r2npts = ring->npts;
		vertex_t *v2 = ring->pts;
		int v2_idx;
		for(v2_idx=0; v2_idx<r2npts; v2_idx++) {
			int same = (r1_idx == r2_idx) && (v1_idx == v2_idx);
			int touches = (v1x == v2->x) && (v1y == v2->y);
			if(touches && !same) {
				//if(VERBOSE) printf("touches for %d,%d vs. %d,%d\n", r1_idx, v1_idx, r2_idx, v2_idx);
				return 1;
			}
			v2++;
		}
	}

	//if(VERBOSE) printf("no touch for %d,%d\n", r1_idx, v1_idx);

	return 0;
}

// This function is only meant to be called on polygons
// that have orthogonal sides on an integer lattice.
static void bevel_self_intersections_old(mpoly_t *mp, double amount) {
	char *table = mpoly_border_touch_create_hashtable(mp);

	int total_pts = 0;
	int r_idx;
	for(r_idx=0; r_idx<mp->num_rings; r_idx++) {
		ring_t *ring = &mp->rings[r_idx];
		total_pts += ring->npts;
	}

	if(VERBOSE) printf("bevel with amount=%lf, total_pts=%d\n", amount, total_pts);

	bbox_t *bboxes = (bbox_t *)malloc_or_die(sizeof(bbox_t) * mp->num_rings);
	for(r_idx=0; r_idx<mp->num_rings; r_idx++) {
		bboxes[r_idx] = make_bbox(mp->rings + r_idx);
	}

	printf("Beveling: ");
	int done_pts = 0;
	for(r_idx=0; r_idx<mp->num_rings; r_idx++) {
		ring_t *ring = &mp->rings[r_idx];
		int num_touch = 0;
		char *touch_mask = NULL;
		int v_idx;
		for(v_idx=0; v_idx<ring->npts; v_idx++) {
			GDALTermProgress((double)(done_pts*2+v_idx)/(double)(total_pts*2), NULL, NULL);
			int touches = mpoly_border_touches_point(table, mp, r_idx, v_idx, bboxes);
			if(touches) {
				if(!touch_mask) {
					touch_mask = (char *)malloc_or_die(ring->npts);
					memset(touch_mask, 0, ring->npts);
				}
				touch_mask[v_idx] = 1;
				num_touch++;
			}
		}
		int old_npts = ring->npts;
		if(num_touch) {
			if(VERBOSE >= 2) printf("ring %d: num_touch=%d\n", r_idx, num_touch);
			int new_numpts = ring->npts + num_touch;
			vertex_t *new_pts = (vertex_t *)malloc_or_die(sizeof(vertex_t) * new_numpts);
			int vout_idx = 0;
			for(v_idx=0; v_idx<ring->npts; v_idx++) {
				GDALTermProgress((double)(done_pts*2+ring->npts+v_idx)/(double)(total_pts*2), NULL, NULL);
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
			ring->npts = new_numpts;
			ring->pts = new_pts;
		}
		done_pts += old_npts;
		if(touch_mask) free(touch_mask);
	}

	GDALTermProgress(1, NULL, NULL);

	free(table);
	free(bboxes);
}

#endif

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
