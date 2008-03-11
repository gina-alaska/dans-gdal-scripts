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

#define EPSILON 10E-10

typedef struct {
	int begin;
	int end;
	char is_problem;
} segment_t;

typedef struct {
	segment_t *segs;
	int num_segs;
} reduced_ring_t;

reduced_ring_t reduce_linestring_detail(ring_t *orig_string, double res);
double get_dist_to_seg(double seg_vec_x, double seg_vec_y, 
	vertex_t *seg_vert1, vertex_t *seg_vert2, vertex_t *test_vert);
ring_t make_ring_from_segs(ring_t *c_in, reduced_ring_t *r_in);
mpoly_t reduction_to_mpoly(mpoly_t *in_mpoly, reduced_ring_t *reduced_rings);
void fix_topology(mpoly_t *mpoly, reduced_ring_t *reduced_rings);
char segs_cross(ring_t *c1, segment_t *s1, ring_t *c2, segment_t *s2);

mpoly_t compute_reduced_pointset(mpoly_t *in_mpoly, double tolerance) {
	if(VERBOSE) printf("reducing...\n");

	if(!in_mpoly->num_rings) {
		return empty_polygon();
	}

	reduced_ring_t *reduced_rings = (reduced_ring_t *)
		malloc_or_die(sizeof(reduced_ring_t) * in_mpoly->num_rings);

	int c_idx;
	for(c_idx=0; c_idx<in_mpoly->num_rings; c_idx++) {
		reduced_rings[c_idx] = reduce_linestring_detail(&in_mpoly->rings[c_idx], tolerance);
	}

	fix_topology(in_mpoly, reduced_rings);

	return reduction_to_mpoly(in_mpoly, reduced_rings);
}

// Implementation of Douglas-Peucker polyline reduction algorithm
// rewrite of code from http://www.3dsoftware.com/Cartography/Programming/PolyLineReduction
// (and adapted from src/linework/dp.c in SwathViewer)

#define VECLEN(x,y) sqrt((x)*(x)+(y)*(y))

reduced_ring_t reduce_linestring_detail(ring_t *orig_string, double res) {
//printf("enter dp\n");

	int num_in = orig_string->npts;
	vertex_t *pts_in = orig_string->pts;
	double tolerance = res;
	int i;

	segment_t *stack = (segment_t *)malloc_or_die(sizeof(segment_t) * num_in);
	int stack_ptr = 0;

	segment_t *keep_segs = (segment_t *)malloc_or_die(sizeof(segment_t) * num_in);
	int num_keep_segs = 0;

	// must keep closure segment
	keep_segs[num_keep_segs].begin = num_in-1;
	keep_segs[num_keep_segs].end = 0;
	num_keep_segs++;

	stack[stack_ptr].begin = 0;
	stack[stack_ptr].end = num_in-1;
	stack_ptr++;
	while(stack_ptr) {
		stack_ptr--;
		int seg_begin = stack[stack_ptr].begin;
		int seg_end = stack[stack_ptr].end;
//printf("stack_ptr=%d / range=[%d,%d]\n", stack_ptr, seg_begin, seg_end);

		double max_dist = -1.0;
		int idx_of_max = -1; // to prevent compiler warning

		double seg_vec_x = pts_in[seg_end].x - pts_in[seg_begin].x;
		double seg_vec_y = pts_in[seg_end].y - pts_in[seg_begin].y;
		double seg_vec_len = VECLEN(seg_vec_x, seg_vec_y);
		if(seg_vec_len > 0.0) {
			// normalize vector
			seg_vec_x /= seg_vec_len;
			seg_vec_y /= seg_vec_len;
			for(i=seg_begin+1; i<seg_end; i++) {
				double dist_to_seg = get_dist_to_seg(seg_vec_x, seg_vec_y,
					pts_in+seg_begin, pts_in+seg_end, pts_in+i);
				if(dist_to_seg < 0.0) fatal_error("dist_to_seg < 0.0");
				if(isnan(dist_to_seg)) fatal_error("dist_to_seg == NaN");

				if(dist_to_seg > max_dist) {
					max_dist = dist_to_seg;
					idx_of_max = i;
				}
			}
		} else {
			// Segment is length zero, so we can't use get_dist_to_seg.
			// Instead, just use cartesian distance
			for(i=seg_begin+1; i<seg_end; i++) {
				double dx = pts_in[i].x - pts_in[seg_begin].x;
				double dy = pts_in[i].y - pts_in[seg_begin].y;
				double dist_to_seg = VECLEN(dx, dy);

				if(dist_to_seg > max_dist) {
					max_dist = dist_to_seg;
					idx_of_max = i;
				}
			}
		}

//printf("max=%.15f, toler=%.15f, idx=%i\n", max_dist, tolerance, idx_of_max);
		if(max_dist >= tolerance) {
			if(idx_of_max <= 0) fatal_error(
				"idx_of_max out of range (perhaps it wasn't set?)");

			// add point and recursively divide subsegments
			stack[stack_ptr].begin = seg_begin;
			stack[stack_ptr].end = idx_of_max;
			stack_ptr++;
			if(stack_ptr >= num_in) fatal_error("stack overflow in dp.c");

			stack[stack_ptr].begin = idx_of_max;
			stack[stack_ptr].end = seg_end;
			stack_ptr++;
			if(stack_ptr >= num_in) fatal_error("stack overflow in dp.c");
		} else {
			// segment doesn't need subdivision - tag
			// endpoint for inclusion
			if(num_keep_segs == num_in) fatal_error("output stack overflow in dp.c");
			keep_segs[num_keep_segs].begin = seg_begin;
			keep_segs[num_keep_segs].end = seg_end;
			num_keep_segs++;
		}
	}

//printf("keeping %d of %d\n", num_keep_segs, num_in);
//printf("exit dp\n");

	free(stack);

	reduced_ring_t ret;
	ret.segs = keep_segs;
	ret.num_segs = num_keep_segs;

	return ret;
}

inline double get_dist_to_seg(double seg_vec_x, double seg_vec_y, 
vertex_t *seg_vert1, vertex_t *seg_vert2, vertex_t *test_vert) {
	double vert_vec_x = test_vert->x - seg_vert1->x;
	double vert_vec_y = test_vert->y - seg_vert1->y;
	double scalar_prod = vert_vec_x*seg_vec_x + vert_vec_y*seg_vec_y;
	if(scalar_prod < 0.0) {
		// past beginning of segment - distance from segment
		// is equal to distance from first endpoint of segment
		return VECLEN(vert_vec_x, vert_vec_y);
	} else {
		vert_vec_x = seg_vert2->x - test_vert->x;
		vert_vec_y = seg_vert2->y - test_vert->y;
		scalar_prod = vert_vec_x*seg_vec_x + vert_vec_y*seg_vec_y;
		if(scalar_prod < 0.0) {
			// past end of segment - distance from segment
			// is equal to distance from second endpoint of segment
			return VECLEN(vert_vec_x, vert_vec_y);
		} else {
			double c_squared = vert_vec_x*vert_vec_x + vert_vec_y*vert_vec_y;
			double a_squared = scalar_prod * scalar_prod;
			double b_squared = c_squared - a_squared;
			if(b_squared < 0) {
				if(fabs(b_squared) < a_squared * EPSILON) {
					b_squared = 0; // correct for round-off error
				} else {
					fatal_error("a_squared > c_squared");
				}
			}
			// use pythagorean theorem to find distance from
			// vertex to segment
			return sqrt(b_squared);
		}
	}
}

ring_t make_ring_from_segs(ring_t *c_in, reduced_ring_t *r_in) {
	int i;
	char *keep_pts = (char *)malloc_or_die(c_in->npts);
	for(i=0; i<c_in->npts; i++) keep_pts[i] = 0;

	for(i=0; i<r_in->num_segs; i++) {
		keep_pts[r_in->segs[i].begin]++;
		keep_pts[r_in->segs[i].end]++;
	}

	int num_to_keep = 0;
	for(i=0; i<c_in->npts; i++) {
		if(keep_pts[i] && keep_pts[i] != 2) {
			fatal_error("point must be in 0 or 2 segs");
		}
		if(keep_pts[i]) num_to_keep++;
	}
	// number of points == number of line segments
	if(num_to_keep != r_in->num_segs) fatal_error("num_to_keep != r_in->num_segs");

//printf("keeping %d of %d\n", num_to_keep, num_in);

	ring_t new_string = *c_in; // copy parent_id, is_hole, etc.
	new_string.npts = num_to_keep;
	new_string.pts = (vertex_t *)malloc_or_die(sizeof(vertex_t) * num_to_keep);
	vertex_t *pts_out = new_string.pts;
	int idx_out = 0;
	for(i=0; i<c_in->npts; i++) {
		if(keep_pts[i]) {
			pts_out[idx_out++] = c_in->pts[i];
		}
	}
	if(idx_out != new_string.npts) {
		fatal_error("count mismatch after point copy in dp.c");
	}

	free(keep_pts);

	return new_string;
}

mpoly_t reduction_to_mpoly(mpoly_t *in_mpoly, reduced_ring_t *reduced_rings) {
	// This function takes only rings with more than two points.
	// (a polygon with two or less points has no area)
	// Care is taken to make sure that the containment info
	// (specifically the parent_id pointer) is still proper.
	// It is difficult to come up with a test case so I just
	// have to hope that is works.

	// First, find which rings have at least three points.
	char *keep_rings = (char *)malloc_or_die(in_mpoly->num_rings);
	int total_npts_in=0, total_npts_out=0;
	int c_idx;
	for(c_idx=0; c_idx<in_mpoly->num_rings; c_idx++) {
		ring_t *c_in = &in_mpoly->rings[c_idx];
		reduced_ring_t *r_in = &reduced_rings[c_idx];

		if(VERBOSE >= 2) printf("ring %d: %d => %d pts\n", c_idx, c_in->npts, r_in->num_segs);
		total_npts_in += c_in->npts;

		if(r_in->num_segs > 2) {
			total_npts_out += r_in->num_segs;
			keep_rings[c_idx] = 1;
		} else {
			keep_rings[c_idx] = 0;
		}
	}

	// If an outer ring has been removed, remove its holes
	// also.  Extremely unlikely but just in case...
	for(c_idx=0; c_idx<in_mpoly->num_rings; c_idx++) {
		ring_t *c_in = &in_mpoly->rings[c_idx];
		if(c_in->parent_id >= 0) {
			if(!keep_rings[c_in->parent_id]) {
				keep_rings[c_idx] = 0;
			}
		}
	}

	// Map index of old rings to index of new rings.
	int *new_idx_map = (int *)malloc_or_die(sizeof(int) * in_mpoly->num_rings);
	int out_idx = 0;
	for(c_idx=0; c_idx<in_mpoly->num_rings; c_idx++) {
		if(keep_rings[c_idx]) {
			//printf("map ring %d => %d\n", c_idx, out_idx);
			new_idx_map[c_idx] = out_idx++;
		} else {
			new_idx_map[c_idx] = -1;
		}
	}
	int num_out_rings = out_idx;

	if(!num_out_rings) return empty_polygon();

	// Now create the new rings and update parent_id using the
	// remapped index values.
	ring_t *out_rings = (ring_t *)malloc_or_die(sizeof(ring_t)*num_out_rings);
	for(c_idx=0; c_idx<in_mpoly->num_rings; c_idx++) {
		ring_t *c_in = &in_mpoly->rings[c_idx];
		reduced_ring_t *r_in = &reduced_rings[c_idx];

		out_idx = new_idx_map[c_idx];
		if(out_idx < 0) continue;

		ring_t new_string = make_ring_from_segs(c_in, r_in);

		int old_parent_id = in_mpoly->rings[c_idx].parent_id;
		if(old_parent_id >= 0) {
			int new_parent_id = new_idx_map[old_parent_id];
			//printf("map parent_id %d => %d\n", old_parent_id, new_parent_id);
			new_string.parent_id = new_parent_id;
		}

		out_rings[out_idx] = new_string;
	}

	if(VERBOSE) printf("reduced %d => %d rings, %d => %d pts\n",
		in_mpoly->num_rings, num_out_rings, total_npts_in, total_npts_out);

	return (mpoly_t){ num_out_rings, out_rings };
}

void fix_topology(mpoly_t *mpoly, reduced_ring_t *reduced_rings) {
	int r1_idx, r2_idx;
	int seg1_idx, seg2_idx;

	printf("Analyzing topology: ");
	fflush(stdout);

	// clear problem flags
	for(r1_idx=0; r1_idx < mpoly->num_rings; r1_idx++) {
		reduced_ring_t *r1 = &reduced_rings[r1_idx];
		for(seg1_idx=0; seg1_idx < r1->num_segs; seg1_idx++) {
			r1->segs[seg1_idx].is_problem = 0;
		}
	}

	bbox_t *bboxes = make_bboxes(mpoly);

	// flag segments that cross
	int have_problems = 0;
	for(r1_idx=0; r1_idx < mpoly->num_rings; r1_idx++) {
		GDALTermProgress(pow((double)r1_idx / (double)mpoly->num_rings, 2), NULL, NULL);
		ring_t *c1 = &mpoly->rings[r1_idx];
		reduced_ring_t *r1 = &reduced_rings[r1_idx];
		bbox_t bbox1 = bboxes[r1_idx];
		if(bbox1.empty) continue;
		for(r2_idx=0; r2_idx < mpoly->num_rings; r2_idx++) {
			if(r2_idx > r1_idx) continue; // symmetry optimization

			bbox_t bbox2 = bboxes[r2_idx];
			if(bbox2.empty) continue;

			if(
				bbox1.min_x > bbox2.max_x ||
				bbox1.min_y > bbox2.max_y ||
				bbox2.min_x > bbox1.max_x ||
				bbox2.min_y > bbox1.max_y
			) continue;

			ring_t *c2 = &mpoly->rings[r2_idx];
			reduced_ring_t *r2 = &reduced_rings[r2_idx];

			for(seg1_idx=0; seg1_idx < r1->num_segs; seg1_idx++) {
				for(seg2_idx=0; seg2_idx < r2->num_segs; seg2_idx++) {
					if(r2_idx == r1_idx && seg2_idx > seg1_idx) continue; // symmetry optimization

					char crosses = segs_cross(c1, &r1->segs[seg1_idx], c2, &r2->segs[seg2_idx]);
					if(crosses) {
						//printf("found a crossing: %d,%d,%d,%d\n",
						//	r1_idx, seg1_idx, r2_idx, seg2_idx);
						r1->segs[seg1_idx].is_problem = 1;
						r2->segs[seg2_idx].is_problem = 1;
						have_problems += 2;
					}
				} // seg loop
			} // seg loop
		} // ring loop
	} // ring loop
	GDALTermProgress(1, NULL, NULL);

	free(bboxes);

	if(have_problems) {
		if(VERBOSE) printf("fixing %d crossed segments from reduction\n", have_problems/2);
	}

	int did_something = 1;
	while(have_problems && did_something) {
		//printf("%d crossings to fix\n", have_problems/2);
		did_something = 0;
		// subdivide problem segments
		for(r1_idx=0; r1_idx < mpoly->num_rings; r1_idx++) {
			reduced_ring_t *r1 = &reduced_rings[r1_idx];
			int orig_num_segs = r1->num_segs; // this number will change as we go, so copy it
			for(seg1_idx=0; seg1_idx < orig_num_segs; seg1_idx++) {
				if(!r1->segs[seg1_idx].is_problem) continue;
				int begin = r1->segs[seg1_idx].begin;
				int end = r1->segs[seg1_idx].end;
				if(end < begin) continue; // this is a closure segment (and so doesn't have a midpoint)
				if(end == begin+1) continue; // no midpoint

				// subdivide this segment
				int mid = (begin + end) / 2;
				r1->segs[seg1_idx].end = mid;
				r1->segs[r1->num_segs].begin = mid;
				r1->segs[r1->num_segs].end = end;
				r1->segs[r1->num_segs].is_problem = 1;
				r1->num_segs++;
				did_something = 1;
			} // seg loop
		} // ring loop

		have_problems = 0;
		// now test for resolved problems and new problems
		for(r1_idx=0; r1_idx < mpoly->num_rings; r1_idx++) {
			ring_t *c1 = &mpoly->rings[r1_idx];
			reduced_ring_t *r1 = &reduced_rings[r1_idx];
			for(seg1_idx=0; seg1_idx < r1->num_segs; seg1_idx++) {
				if(!r1->segs[seg1_idx].is_problem) continue;
				r1->segs[seg1_idx].is_problem = 0;
				for(r2_idx=0; r2_idx < mpoly->num_rings; r2_idx++) {
					ring_t *c2 = &mpoly->rings[r2_idx];
					reduced_ring_t *r2 = &reduced_rings[r2_idx];
					for(seg2_idx=0; seg2_idx < r2->num_segs; seg2_idx++) {
						char crosses = segs_cross(c1, &r1->segs[seg1_idx], c2, &r2->segs[seg2_idx]);
						if(crosses) {
							if(VERBOSE) {
								printf("found a crossing (still): %d,%d,%d,%d (%f,%f)-(%f,%f) (%f,%f)-(%f,%f)\n",
									r1_idx, seg1_idx, r2_idx, seg2_idx,
									c1->pts[r1->segs[seg1_idx].begin].x,
									c1->pts[r1->segs[seg1_idx].begin].y,
									c1->pts[r1->segs[seg1_idx].end].x,
									c1->pts[r1->segs[seg1_idx].end].y,
									c2->pts[r2->segs[seg2_idx].begin].x,
									c2->pts[r2->segs[seg2_idx].begin].y,
									c2->pts[r2->segs[seg2_idx].end].x,
									c2->pts[r2->segs[seg2_idx].end].y);
							}
							r1->segs[seg1_idx].is_problem = 1;
							r2->segs[seg2_idx].is_problem = 1;
							have_problems++;
						}
					} // seg loop
				} // ring loop
			} // seg loop
		} // ring loop
	} // while problems

	if(have_problems) {
		printf("WARNING: Could not fix all topology problems.\n  Please inspect output shapefile manually.\n");
	}
}

inline char segs_cross(ring_t *c1, segment_t *s1, ring_t *c2, segment_t *s2) {
	if(c1 == c2) {
		// don't test crossing if segments are identical or neighbors
		if(
			(s1->begin == s2->begin) ||
			(s1->begin == s2->end) ||
			(s1->end == s2->begin) ||
			(s1->end == s2->end)
		) return 0;
	}
	
	return line_intersects_line(
		c1->pts[s1->begin], c1->pts[s1->end  ],
		c2->pts[s2->begin], c2->pts[s2->end  ],
		1);
}
