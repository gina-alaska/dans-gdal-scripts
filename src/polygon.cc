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
} reduced_contour_t;

typedef struct {
	int num_crossings;
	int array_size;
	int *crossings;
} row_crossings_t;

reduced_contour_t reduce_linestring_detail(contour_t *orig_string, double res);
double get_dist_to_seg(double seg_vec_x, double seg_vec_y, 
	vertex_t *seg_vert1, vertex_t *seg_vert2, vertex_t *test_vert);
contour_t make_contour_from_segs(contour_t *c_in, reduced_contour_t *r_in);
void fix_topology(mpoly_t *mpoly, reduced_contour_t *reduced_contours);
char segs_cross(contour_t *c1, segment_t *s1, contour_t *c2, segment_t *s2);
int line_intersects_line(
	double x1, double y1, double x2, double y2,
	double x3, double y3, double x4, double y4,
	int fail_on_coincident
);

extern int VERBOSE;

void output_wkt_mpoly(char *wkt_fn, mpoly_t mpoly) {
	int num_contours = mpoly.num_contours;
	contour_t *contours = mpoly.contours;

	FILE *fout = fopen(wkt_fn, "w");
	if(!fout) fatal_error("cannot open output file for WKT");

	fprintf(fout, "MULTIPOLYGON(\n");
	int r_idx, h_idx, p_idx;
	int is_first_ring = 1;
	for(r_idx=0; r_idx<num_contours; r_idx++) {
		contour_t *ring = contours + r_idx;
		if(ring->is_hole) continue;

		if(!is_first_ring) fprintf(fout, ", ");
		is_first_ring = 0;
		fprintf(fout, "((\n");
		//fprintf(fout, "  ring:%d\n", r_idx);
		for(p_idx=0; p_idx<ring->npts+1; p_idx++) {
			double x = ring->pts[p_idx % ring->npts].x;
			double y = ring->pts[p_idx % ring->npts].y;
			if(!(p_idx%4)) {
				if(p_idx) fprintf(fout, "\n");
				fprintf(fout, "  ");
			}
			fprintf(fout, "%.15f %.15f", x, y);
			if(p_idx < ring->npts) fprintf(fout, ", ");
		}
		fprintf(fout, "\n)");

		for(h_idx=0; h_idx<num_contours; h_idx++) {
			contour_t *hole = contours + h_idx;
			if(hole->parent_id != r_idx) continue;

			fprintf(fout, ", (\n");
			//fprintf(fout, "  hole:%d\n", h_idx);
			for(p_idx=0; p_idx<hole->npts+1; p_idx++) {
				double x = hole->pts[p_idx % hole->npts].x;
				double y = hole->pts[p_idx % hole->npts].y;
				if(!(p_idx%4)) {
					if(p_idx) fprintf(fout, "\n");
					fprintf(fout, "  ");
				}
				fprintf(fout, "%.15f %.15f", x, y);
				if(p_idx < hole->npts) fprintf(fout, ", ");
			}
			fprintf(fout, "\n)");
		}
		fprintf(fout, ")");
	}
	fprintf(fout, ")\n");

	fclose(fout);
}

mpoly_t compute_reduced_pointset(mpoly_t *in_mpoly, double tolerance) {
	if(VERBOSE) fprintf(stderr, "reducing...\n");

	reduced_contour_t *reduced_contours = (reduced_contour_t *)
		malloc_or_die(sizeof(reduced_contour_t) * in_mpoly->num_contours);
	int c_idx;
	for(c_idx=0; c_idx<in_mpoly->num_contours; c_idx++) {
		reduced_contours[c_idx] = reduce_linestring_detail(&in_mpoly->contours[c_idx], tolerance);
	}

	fix_topology(in_mpoly, reduced_contours);

	contour_t *out_contours = (contour_t *)malloc_or_die(sizeof(contour_t)*in_mpoly->num_contours);
	int num_out_contours = 0;
	int total_npts_in=0, total_npts_out=0;
	for(c_idx=0; c_idx<in_mpoly->num_contours; c_idx++) {
		contour_t *c_in = &in_mpoly->contours[c_idx];
		reduced_contour_t *r_in = &reduced_contours[c_idx];

		contour_t new_string = make_contour_from_segs(c_in, r_in);

		if(VERBOSE) fprintf(stderr, "contour %d: %d => %d pts\n", c_idx, c_in->npts, new_string.npts);
		total_npts_in += c_in->npts;
		// FIXME - if skipping an outer ring, skip its holes too
		if(new_string.npts > 2) {
			out_contours[num_out_contours++] = new_string;
			total_npts_out += new_string.npts;
		}
	}

	if(VERBOSE) fprintf(stderr, "reduced %d => %d contours, %d => %d pts\n",
		in_mpoly->num_contours, num_out_contours, total_npts_in, total_npts_out);

	return (mpoly_t){ num_out_contours, out_contours };
}

// Implementation of Douglas-Peucker polyline reduction algorithm
// rewrite of code from http://www.3dsoftware.com/Cartography/Programming/PolyLineReduction
// (and adapted from src/linework/dp.c in the sv_server module)

#define VECLEN(x,y) sqrt((x)*(x)+(y)*(y))

reduced_contour_t reduce_linestring_detail(contour_t *orig_string, double res) {
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

	reduced_contour_t ret;
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

contour_t make_contour_from_segs(contour_t *c_in, reduced_contour_t *r_in) {
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
//fprintf(stderr, "keeping %d of %d\n", num_to_keep, num_in);

	contour_t new_string = *c_in; // copy parent_id, is_hole, etc.
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

void fix_topology(mpoly_t *mpoly, reduced_contour_t *reduced_contours) {
	int c1_idx, c2_idx;
	int seg1_idx, seg2_idx;

	// clear problem flags
	for(c1_idx=0; c1_idx < mpoly->num_contours; c1_idx++) {
		reduced_contour_t *r1 = &reduced_contours[c1_idx];
		for(seg1_idx=0; seg1_idx < r1->num_segs; seg1_idx++) {
			r1->segs[seg1_idx].is_problem = 0;
		}
	}

	// flag segments that cross
	char have_problems = 0;
	for(c1_idx=0; c1_idx < mpoly->num_contours; c1_idx++) {
		contour_t *c1 = &mpoly->contours[c1_idx];
		reduced_contour_t *r1 = &reduced_contours[c1_idx];
		for(c2_idx=0; c2_idx < mpoly->num_contours; c2_idx++) {
			if(c2_idx > c1_idx) continue; // symmetry optimization

			contour_t *c2 = &mpoly->contours[c2_idx];
			reduced_contour_t *r2 = &reduced_contours[c2_idx];

			for(seg1_idx=0; seg1_idx < r1->num_segs; seg1_idx++) {
				for(seg2_idx=0; seg2_idx < r2->num_segs; seg2_idx++) {
					if(c2_idx == c1_idx && seg2_idx > seg1_idx) continue; // symmetry optimization

					char crosses = segs_cross(c1, &r1->segs[seg1_idx], c2, &r2->segs[seg2_idx]);
					if(crosses) {
						fprintf(stderr, "found a crossing: %d,%d,%d,%d\n",
							c1_idx, seg1_idx, c2_idx, seg2_idx);
						r1->segs[seg1_idx].is_problem = 1;
						r2->segs[seg2_idx].is_problem = 1;
						have_problems = 1;
					}
				} // seg loop
			} // seg loop
		} // contour loop
	} // contour loop

	int did_something = 1;
	while(have_problems && did_something) {
		did_something = 0;
		// subdivide problem segments
		for(c1_idx=0; c1_idx < mpoly->num_contours; c1_idx++) {
			reduced_contour_t *r1 = &reduced_contours[c1_idx];
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
		} // contour loop

		have_problems = 0;
		// now test for resolved problems and new problems
		for(c1_idx=0; c1_idx < mpoly->num_contours; c1_idx++) {
			contour_t *c1 = &mpoly->contours[c1_idx];
			reduced_contour_t *r1 = &reduced_contours[c1_idx];
			for(seg1_idx=0; seg1_idx < r1->num_segs; seg1_idx++) {
				if(!r1->segs[seg1_idx].is_problem) continue;
				r1->segs[seg1_idx].is_problem = 0;
				for(c2_idx=0; c2_idx < mpoly->num_contours; c2_idx++) {
					contour_t *c2 = &mpoly->contours[c2_idx];
					reduced_contour_t *r2 = &reduced_contours[c2_idx];
					for(seg2_idx=0; seg2_idx < r2->num_segs; seg2_idx++) {
						char crosses = segs_cross(c1, &r1->segs[seg1_idx], c2, &r2->segs[seg2_idx]);
						if(crosses) {
							fprintf(stderr, "found a crossing (still): %d,%d,%d,%d\n",
								c1_idx, seg1_idx, c2_idx, seg2_idx);
							r1->segs[seg1_idx].is_problem = 1;
							r2->segs[seg2_idx].is_problem = 1;
							have_problems = 1;
						}
					} // seg loop
				} // contour loop
			} // seg loop
		} // contour loop
	} // while problems

	if(have_problems) {
		fprintf(stderr, "WARNING: Could not fix all topology problems.\n  Please inspect output shapefile manually.\n");
	}
}

char segs_cross(contour_t *c1, segment_t *s1, contour_t *c2, segment_t *s2) {
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
		c1->pts[s1->begin].x, c1->pts[s1->begin].y,
		c1->pts[s1->end  ].x, c1->pts[s1->end  ].y,
		c2->pts[s2->begin].x, c2->pts[s2->begin].y,
		c2->pts[s2->end  ].x, c2->pts[s2->end  ].y,
		1
	);
}

int line_intersects_line(
	double x1, double y1, double x2, double y2,
	double x3, double y3, double x4, double y4,
	int fail_on_coincident
) {
	if(
		MAX(x1, x2) < MIN(x3, x4) ||
		MIN(x1, x2) > MAX(x3, x4) ||
		MAX(y1, y2) < MIN(y3, y4) ||
		MIN(y1, y2) > MAX(y3, y4)
	) return 0;
	double numer_a = (x4-x3)*(y1-y3) - (y4-y3)*(x1-x3);
	double numer_b = (x2-x1)*(y1-y3) - (y2-y1)*(x1-x3);
	double denom   = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1);
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

int polygon_contains(contour_t *c1, contour_t *c2) {
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

double polygon_area(contour_t *c) {
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

void mask_from_mpoly(mpoly_t *mpoly, int w, int h, char *fn) {
	int i, j, y;

	fprintf(stderr, "mask draw: begin\n");

	row_crossings_t *rows = (row_crossings_t *)malloc_or_die(sizeof(row_crossings_t) * h);
	for(i=0; i<h; i++) {
		rows[i].num_crossings = 0;
		rows[i].array_size = 0;
		rows[i].crossings = NULL;
	}

	for(i=0; i<mpoly->num_contours; i++) {
		contour_t *c = mpoly->contours + i;
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
