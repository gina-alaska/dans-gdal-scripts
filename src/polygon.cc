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

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

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

typedef struct {
	int num_crossings;
	int array_size;
	int *crossings;
} row_crossings_t;

reduced_ring_t reduce_linestring_detail(ring_t *orig_string, double res);
double get_dist_to_seg(double seg_vec_x, double seg_vec_y, 
	vertex_t *seg_vert1, vertex_t *seg_vert2, vertex_t *test_vert);
ring_t make_ring_from_segs(ring_t *c_in, reduced_ring_t *r_in);
mpoly_t reduction_to_mpoly(mpoly_t *in_mpoly, reduced_ring_t *reduced_rings);
void fix_topology(mpoly_t *mpoly, reduced_ring_t *reduced_rings);
char segs_cross(ring_t *c1, segment_t *s1, ring_t *c2, segment_t *s2);

extern int VERBOSE;

void output_wkt_mpoly(char *wkt_fn, mpoly_t mpoly, int split_polys) {
	int num_rings = mpoly.num_rings;
	ring_t *rings = mpoly.rings;

	FILE *fout = fopen(wkt_fn, "w");
	if(!fout) fatal_error("cannot open output file for WKT");

	if(!split_polys) fprintf(fout, "MULTIPOLYGON(\n");

	int r_idx, h_idx, p_idx;
	int is_first_ring = 1;
	for(r_idx=0; r_idx<num_rings; r_idx++) {
		ring_t *ring = rings + r_idx;
		if(ring->is_hole) continue;

		if(!is_first_ring) fprintf(fout, split_polys ? "\n\n" : ", ");
		is_first_ring = 0;

		fprintf(fout, split_polys ? "POLYGON((\n" : "((\n");

		//fprintf(fout, "  ring:%d\n", r_idx);
		for(p_idx=0; p_idx<ring->npts+1; p_idx++) {
			if(!(p_idx%4)) {
				if(p_idx) fprintf(fout, "\n");
				fprintf(fout, "  ");
			}
			vertex_t v = ring->pts[p_idx % ring->npts];
			fprintf(fout, "%.15f %.15f", v.x, v.y);
			if(p_idx < ring->npts) fprintf(fout, ", ");
		}
		fprintf(fout, "\n)");

		for(h_idx=0; h_idx<num_rings; h_idx++) {
			ring_t *hole = rings + h_idx;
			if(hole->parent_id != r_idx) continue;

			fprintf(fout, ", (\n");
			//fprintf(fout, "  hole:%d\n", h_idx);
			for(p_idx=0; p_idx<hole->npts+1; p_idx++) {
				if(!(p_idx%4)) {
					if(p_idx) fprintf(fout, "\n");
					fprintf(fout, "  ");
				}
				vertex_t v = hole->pts[p_idx % hole->npts];
				fprintf(fout, "%.15f %.15f", v.x, v.y);
				if(p_idx < hole->npts) fprintf(fout, ", ");
			}
			fprintf(fout, "\n)");
		}
		fprintf(fout, ")");
	}

	if(!split_polys) fprintf(fout, ")\n");

	fclose(fout);
}

mpoly_t compute_reduced_pointset(mpoly_t *in_mpoly, double tolerance) {
	if(VERBOSE) fprintf(stderr, "reducing...\n");

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
// (and adapted from src/linework/dp.c in the sv_server module)

#define VECLEN(x,y) sqrt((x)*(x)+(y)*(y))

reduced_ring_t reduce_linestring_detail(ring_t *orig_string, double res) {
//fprintf(stderr, "enter dp\n");

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
//fprintf(stderr, "stack_ptr=%d / range=[%d,%d]\n", stack_ptr, seg_begin, seg_end);

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

//fprintf(stderr, "max=%.15f, toler=%.15f, idx=%i\n", max_dist, tolerance, idx_of_max);
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

//fprintf(stderr, "keeping %d of %d\n", num_keep_segs, num_in);
//fprintf(stderr, "exit dp\n");

	free(stack);

	reduced_ring_t ret;
	ret.segs = keep_segs;
	ret.num_segs = num_keep_segs;

	return ret;
}

double get_dist_to_seg(double seg_vec_x, double seg_vec_y, 
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

//fprintf(stderr, "keeping %d of %d\n", num_to_keep, num_in);

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

		if(VERBOSE) fprintf(stderr, "ring %d: %d => %d pts\n", c_idx, c_in->npts, r_in->num_segs);
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
			//fprintf(stderr, "map ring %d => %d\n", c_idx, out_idx);
			new_idx_map[c_idx] = out_idx++;
		} else {
			new_idx_map[c_idx] = -1;
		}
	}
	int num_out_rings = out_idx;

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
			//fprintf(stderr, "map parent_id %d => %d\n", old_parent_id, new_parent_id);
			new_string.parent_id = new_parent_id;
		}

		out_rings[out_idx] = new_string;
	}

	if(VERBOSE) fprintf(stderr, "reduced %d => %d rings, %d => %d pts\n",
		in_mpoly->num_rings, num_out_rings, total_npts_in, total_npts_out);

	return (mpoly_t){ num_out_rings, out_rings };
}

void fix_topology(mpoly_t *mpoly, reduced_ring_t *reduced_rings) {
	int c1_idx, c2_idx;
	int seg1_idx, seg2_idx;

	// clear problem flags
	for(c1_idx=0; c1_idx < mpoly->num_rings; c1_idx++) {
		reduced_ring_t *r1 = &reduced_rings[c1_idx];
		for(seg1_idx=0; seg1_idx < r1->num_segs; seg1_idx++) {
			r1->segs[seg1_idx].is_problem = 0;
		}
	}

	// flag segments that cross
	int have_problems = 0;
	for(c1_idx=0; c1_idx < mpoly->num_rings; c1_idx++) {
		ring_t *c1 = &mpoly->rings[c1_idx];
		reduced_ring_t *r1 = &reduced_rings[c1_idx];
		for(c2_idx=0; c2_idx < mpoly->num_rings; c2_idx++) {
			if(c2_idx > c1_idx) continue; // symmetry optimization

			ring_t *c2 = &mpoly->rings[c2_idx];
			reduced_ring_t *r2 = &reduced_rings[c2_idx];

			for(seg1_idx=0; seg1_idx < r1->num_segs; seg1_idx++) {
				for(seg2_idx=0; seg2_idx < r2->num_segs; seg2_idx++) {
					if(c2_idx == c1_idx && seg2_idx > seg1_idx) continue; // symmetry optimization

					char crosses = segs_cross(c1, &r1->segs[seg1_idx], c2, &r2->segs[seg2_idx]);
					if(crosses) {
						//fprintf(stderr, "found a crossing: %d,%d,%d,%d\n",
						//	c1_idx, seg1_idx, c2_idx, seg2_idx);
						r1->segs[seg1_idx].is_problem = 1;
						r2->segs[seg2_idx].is_problem = 1;
						have_problems += 2;
					}
				} // seg loop
			} // seg loop
		} // ring loop
	} // ring loop

	if(have_problems) {
		if(VERBOSE) fprintf(stderr, "fixing %d crossed segments from reduction\n", have_problems/2);
	}

	int did_something = 1;
	while(have_problems && did_something) {
		//fprintf(stderr, "%d crossings to fix\n", have_problems/2);
		did_something = 0;
		// subdivide problem segments
		for(c1_idx=0; c1_idx < mpoly->num_rings; c1_idx++) {
			reduced_ring_t *r1 = &reduced_rings[c1_idx];
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
		for(c1_idx=0; c1_idx < mpoly->num_rings; c1_idx++) {
			ring_t *c1 = &mpoly->rings[c1_idx];
			reduced_ring_t *r1 = &reduced_rings[c1_idx];
			for(seg1_idx=0; seg1_idx < r1->num_segs; seg1_idx++) {
				if(!r1->segs[seg1_idx].is_problem) continue;
				r1->segs[seg1_idx].is_problem = 0;
				for(c2_idx=0; c2_idx < mpoly->num_rings; c2_idx++) {
					ring_t *c2 = &mpoly->rings[c2_idx];
					reduced_ring_t *r2 = &reduced_rings[c2_idx];
					for(seg2_idx=0; seg2_idx < r2->num_segs; seg2_idx++) {
						char crosses = segs_cross(c1, &r1->segs[seg1_idx], c2, &r2->segs[seg2_idx]);
						if(crosses) {
							//fprintf(stderr, "found a crossing (still): %d,%d,%d,%d\n",
							//	c1_idx, seg1_idx, c2_idx, seg2_idx);
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
		fprintf(stderr, "WARNING: Could not fix all topology problems.\n  Please inspect output shapefile manually.\n");
	}
}

char segs_cross(ring_t *c1, segment_t *s1, ring_t *c2, segment_t *s2) {
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

int line_intersects_line(
	vertex_t p1, vertex_t p2,
	vertex_t p3, vertex_t p4,
	int fail_on_coincident
) {
	if(
		MAX(p1.x, p2.x) < MIN(p3.x, p4.x) ||
		MIN(p1.x, p2.x) > MAX(p3.x, p4.x) ||
		MAX(p1.y, p2.y) < MIN(p3.y, p4.y) ||
		MIN(p1.y, p2.y) > MAX(p3.y, p4.y)
	) return 0;
	double numer_a = (p4.x-p3.x)*(p1.y-p3.y) - (p4.y-p3.y)*(p1.x-p3.x);
	double numer_b = (p2.x-p1.x)*(p1.y-p3.y) - (p2.y-p1.y)*(p1.x-p3.x);
	double denom   = (p4.y-p3.y)*(p2.x-p1.x) - (p4.x-p3.x)*(p2.y-p1.y);
	if(denom == 0) {
		if(numer_a==0 && numer_b==0) { // coincident
			if(fail_on_coincident) {
				return 0;
			} else {
				// lines must touch because of min/max test above
				return 1;
			}
		} else { // parallel
			return 0;
		}
	}
	double ua = (double)numer_a / (double)denom;
	double ub = (double)numer_b / (double)denom;
	return ua>=0 && ua<=1 && ub>=0 && ub<=1;
}

void line_line_intersection(
	vertex_t p1, vertex_t p2,
	vertex_t p3, vertex_t p4,
	vertex_t *p_out
) {
	double numer_a = (p4.x-p3.x)*(p1.y-p3.y) - (p4.y-p3.y)*(p1.x-p3.x);
	//double numer_b = (p2.x-p1.x)*(p1.y-p3.y) - (p2.y-p1.y)*(p1.x-p3.x);
	double denom   = (p4.y-p3.y)*(p2.x-p1.x) - (p4.x-p3.x)*(p2.y-p1.y);
	if(!denom) fatal_error("lines are parallel");
	double ua = numer_a / denom;
	//double ub = numer_b / denom;
	(*p_out).x = p1.x + ua*(p2.x-p1.x);
	(*p_out).y = p1.y + ua*(p2.y-p1.y);
}

int polygon_contains(ring_t *c1, ring_t *c2) {
	// NOTE: by assumption, c1 and c2 don't cross

	int c2npts = c2->npts;
	// way faster, but can give false positive if the two polygons
	// share a common vertex
	//int c2npts = 1;

	int j;
	for(j=0; j<c2npts; j++) {
		double px = c2->pts[j].x;
		double py = c2->pts[j].y;
		//fprintf(stderr, "px=%.15f py=%.15f\n", px, py);
		int num_crossings = 0;
		int i;
		for(i=0; i<c1->npts; i++) {
			double x0 = c1->pts[i].x;
			double y0 = c1->pts[i].y;
			double x1 = c1->pts[(i+1) % c1->npts].x;
			double y1 = c1->pts[(i+1) % c1->npts].y;
			// we want to know whether a ray from (px,py) in the (1,0) direction
			// passes through this segment
			if(x0 < px && x1 < px) continue;

			int y0above = y0 >= py;
			int y1above = y1 >= py;
			if(y0above && y1above) continue;
			if(!y0above && !y1above) continue;

			double alpha = (py-y0)/(y1-y0);
			double cx = x0 + (x1-x0)*alpha;
			//fprintf(stderr, "x0=%.15f x1=%.15f\n", x0, x1);
			//fprintf(stderr, "y0=%.15f y1=%.15f\n", y0, y1);
			//fprintf(stderr, "alpha=%.15f cx=%.15f px=%.15f\n", alpha, cx, px);
			if(cx > px) num_crossings++;
		}
		//fprintf(stderr, "num_crossings=%d\n", num_crossings);
		// if there are an odd number of crossings, then (px,py) is within c1
		if(0 == (num_crossings & 1)) return 0;
	}
	return 1;
}

double polygon_area(ring_t *c) {
	double accum = 0;
	int i;
	for(i=0; i<c->npts; i++) {
		double x0 = c->pts[i].x;
		double y0 = c->pts[i].y;
		double x1 = c->pts[(i+1)%c->npts].x;
		double y1 = c->pts[(i+1)%c->npts].y;
		accum += x0*y1 - x1*y0;
	}
	return fabs(accum) / 2.0;
}

void compute_containments(mpoly_t *mp) {
	int i, j;

	unsigned char **containments = (unsigned char **)malloc_or_die(
		sizeof(unsigned char *) * mp->num_rings);
	for(i=0; i<mp->num_rings; i++) {
		containments[i] = (unsigned char *)malloc_or_die(mp->num_rings);
		for(j=0; j<mp->num_rings; j++) {
			if(i == j) {
				containments[i][j] = 0;
			} else {
				containments[i][j] = polygon_contains(&mp->rings[i], &mp->rings[j]);
				//fprintf(stderr, "containtments[%d][%d] = %d\n", i, j, containments[i][j]);
			}
		}
	}
	int *containment_levels = (int *)malloc_or_die(
		sizeof(int) * mp->num_rings);
	int max_level = 0;
	for(i=0; i<mp->num_rings; i++) {
		containment_levels[i] = 0;
		for(j=0; j<mp->num_rings; j++) {
			if(containments[i][j] && containments[j][i]) {
				fprintf(stderr, "topology error: %d and %d contain each other\n", i, j);
				fatal_error("topology error");
			}
			if(containments[j][i]) containment_levels[i]++;
		}
		if(VERBOSE) fprintf(stderr, "containment_levels[%d] = %d\n", i, containment_levels[i]);
		if(containment_levels[i] > max_level) max_level = containment_levels[i];
	}

	for(i=0; i<mp->num_rings; i++) {
		// only odd levels are holes
		mp->rings[i].is_hole = containment_levels[i] % 2;
	}

	for(i=0; i<mp->num_rings; i++) {
		mp->rings[i].parent_id = -1;

		for(j=0; j<mp->num_rings; j++) {
			if(
				containments[j][i] &&
				containment_levels[i] == containment_levels[j] + 1
			) {
				mp->rings[i].parent_id = j;
			}
		}
	}
}

void mask_from_mpoly(mpoly_t *mpoly, int w, int h, char *fn) {
	int i, j, y;

	fprintf(stderr, "mask draw: begin\n");

	row_crossings_t *rows = (row_crossings_t *)malloc_or_die(sizeof(row_crossings_t) * h);
	for(i=0; i<h; i++) {
		rows[i].num_crossings = 0;
		rows[i].array_size = 0;
		rows[i].crossings = NULL;
	}

	for(i=0; i<mpoly->num_rings; i++) {
		ring_t *c = mpoly->rings + i;
		for(j=0; j<c->npts; j++) {
			double x0 = c->pts[j].x;
			double y0 = c->pts[j].y;
			double x1 = c->pts[(j+1)%c->npts].x;
			double y1 = c->pts[(j+1)%c->npts].y;
			if(y0 == y1) continue;
			if(y0 > y1) {
				double tmp;
				tmp=x0; x0=x1; x1=tmp; 
				tmp=y0; y0=y1; y1=tmp; 
			}
			double alpha = (x1-x0) / (y1-y0);
			for(y=(int)y0; y<(int)y1; y++) {
				if(y<0 || y>h-1) continue;
				int x = (int)(x0 + ((double)y - y0)*alpha);
				row_crossings_t *r = rows+y;
				if(r->num_crossings == r->array_size) {
					r->array_size += 16;
					r->crossings = (int *)realloc_or_die(r->crossings,
						sizeof(int) * r->array_size);
				}
				r->crossings[r->num_crossings++] = x;
			}
		}
	}

	fprintf(stderr, "mask draw: write\n");

	FILE *fout = fopen(fn, "wb");
	if(!fout) fatal_error("cannot open mask output");
	fprintf(fout, "P4\n%d %d\n", w, h);
	unsigned char *buf = (unsigned char *)malloc_or_die((w+7)/8);
	for(y=0; y<h; y++) {
		memset(buf, 0, (w+7)/8);
		unsigned char *p = buf;
		unsigned char bitp = 128;
		row_crossings_t *r = rows+y;
		for(i=0; i<w; i++) {
			unsigned char v = 1;
			// not the fastest way...
			for(j=0; j<r->num_crossings; j++) {
				if(i >= r->crossings[j]) v = !v;
			}
			if(v) *p |= bitp;
			bitp >>= 1;
			if(!bitp) {
				p++;
				bitp = 128;
			}
		}
		fwrite(buf, (w+7)/8, 1, fout);
	}
	fclose(fout);

	fprintf(stderr, "mask draw: done\n");
}

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
					printf("touch at ring %d vert %d vs %d : xy=(%f,%f)\n", r_idx, v1out_idx, v2_idx, v1->x, v1->y);
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
