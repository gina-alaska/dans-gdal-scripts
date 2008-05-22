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
#include "georef.h"

// FIXME - memory leaks

typedef struct {
	int npts;
	int *point_ids;
} box_contents_t;

typedef struct {
	box_contents_t *cells;
	int num_x, num_y;
	double stepsize;
	double min_x, min_y;
} chopdata_t;

static int get_cell_id(chopdata_t *cd, vertex_t vert) {
	int x = (int)((vert.x - cd->min_x) / cd->stepsize);
	int y = (int)((vert.y - cd->min_y) / cd->stepsize);
	if(x<0 || y<0 || x>=cd->num_x || y>=cd->num_y) fatal_error("cell out of range");
	return (y*cd->num_x) + x;
}

static chopdata_t chop_ring_into_cells(ring_t *ring, double stepsize) {
	if(ring->npts < 3) fatal_error("ring must have at least three points");
	double min_x = ring->pts[0].x;
	double min_y = ring->pts[0].y;
	double max_x=min_x, max_y=min_y;
	for(int i=0; i<ring->npts; i++) {
		vertex_t v = ring->pts[i];
		if(v.x < min_x) min_x = v.x;
		if(v.x > max_x) max_x = v.x;
		if(v.y < min_y) min_y = v.y;
		if(v.y > max_y) max_y = v.y;
	}

	chopdata_t cd;
	cd.min_x = min_x;
	cd.min_y = min_y;
	cd.stepsize = stepsize;
	cd.num_x = (int)((max_x - cd.min_x) / cd.stepsize) + 1;
	cd.num_y = (int)((max_y - cd.min_y) / cd.stepsize) + 1;
	cd.cells = (box_contents_t *)malloc_or_die(sizeof(box_contents_t) * cd.num_x * cd.num_y);
	for(int i=0; i < cd.num_x * cd.num_y; i++) {
		cd.cells[i].npts = 0;
		cd.cells[i].point_ids = NULL;
	}

	for(int vid=0; vid<ring->npts; vid++) {
		int cid = get_cell_id(&cd, ring->pts[vid]);
		cd.cells[cid].point_ids = (int *)realloc_or_die(
			cd.cells[cid].point_ids,
			sizeof(int) * (cd.cells[cid].npts+1));
		cd.cells[cid].point_ids[ cd.cells[cid].npts++ ] = vid;
	}

	return cd;
}

static double vert_dist(vertex_t v1, vertex_t v2) {
	double dx = v1.x - v2.x;
	double dy = v1.y - v2.y;
	return sqrt(dx*dx + dy*dy);
}

//static int is_almost_colinear(ring_t *ring, int v1_idx, int v2_idx, int v3_idx, int v4_idx) {
//	vertex_t p1 = ring->pts[v1_idx];
//	vertex_t p2 = ring->pts[v2_idx];
//	vertex_t p3 = ring->pts[v3_idx];
//	vertex_t p4 = ring->pts[v4_idx];
//	double numer_a = (p4.x-p3.x)*(p1.y-p3.y) - (p4.y-p3.y)*(p1.x-p3.x);
//	double numer_b = (p2.x-p1.x)*(p1.y-p3.y) - (p2.y-p1.y)*(p1.x-p3.x);
//	double denom   = (p4.y-p3.y)*(p2.x-p1.x) - (p4.x-p3.x)*(p2.y-p1.y);
//	double scale   = (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) +
//	                 (p3.x-p4.x)*(p3.x-p4.x) + (p3.y-p4.y)*(p3.y-p4.y);
////	printf("is_colin: %d,%d,%d,%d :: %lf,%lf,%lf,%lf\n",
////		v1_idx, v2_idx, v3_idx, v4_idx, numer_a, numer_b, denom, scale);
//	double eps = 0.01;
//	if(
//		denom   < scale*eps && 
//		numer_a < scale*eps &&
//		numer_b < scale*eps
//	) {
//		return 1;
//	} else {
//		return 0;
//	}
//}

static int is_almost_colinear(ring_t *ring, int v1_idx, int v2_idx, int v3_idx, int v4_idx) {
	vertex_t p1 = ring->pts[v1_idx];
	vertex_t p2 = ring->pts[v2_idx];
	vertex_t p3 = ring->pts[v3_idx];
	vertex_t p4 = ring->pts[v4_idx];
	double d21x = p2.x - p1.x; double d21y = p2.y - p1.y;
	double d32x = p3.x - p2.x; double d32y = p3.y - p2.y;
	double d13x = p1.x - p3.x; double d13y = p1.y - p3.y;
	double d43x = p4.x - p3.x; double d43y = p4.y - p3.y;
	double len_a = sqrt(d21x*d21x + d21y*d21y);
	double len_b = sqrt(d43x*d43x + d43y*d43y);
	if(len_a < 10.0 || len_b < 10.0) return 0; // FIXME
	// from http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
	double dist_from_line = fabs(d21x*d13y - d13x*d21y) / sqrt(d21x*d21x + d21y*d21y);
	double dist_from_seg = sqrt(d32x*d32x + d32y*d32y);
	if(dist_from_line > dist_from_seg * 0.1) return 0; // FIXME - arbitrary
	double ang1 = atan2(d21y, d21x);
	double ang2 = atan2(d43y, d43x);
	double max_diff = PI / 180.0 * 5.0; // FIXME - arbitrary
	return 
		fabs(ang1 - ang2) < max_diff ||
		fabs(ang1 - ang2 + 2.0*PI) < max_diff ||
		fabs(ang1 - ang2 - 2.0*PI) < max_diff;
}

static int intcompare(const void *ap, const void *bp) {
	int a = *((int *)ap);
	int b = *((int *)bp);
	return (a<b) ? -1 : (a>b) ? 1 : 0;
}

static box_contents_t get_neighbors(ring_t *ring, chopdata_t *cd, int v1id, double radius) {
	box_contents_t out;
	out.npts = 0;
	out.point_ids = NULL;
	vertex_t v1 = ring->pts[v1id];
	int center_x = (int)((v1.x - cd->min_x) / cd->stepsize);
	int center_y = (int)((v1.y - cd->min_y) / cd->stepsize);
	int radius_cells = (int)(radius / cd->stepsize) + 1;
	for(int x=center_x-radius_cells; x<=center_x+radius_cells; x++) {
		if(x<0 || x>=cd->num_x) continue;
		for(int y=center_y-radius_cells; y<=center_y+radius_cells; y++) {
			if(y<0 || y>=cd->num_y) continue;
			box_contents_t *cont = cd->cells + (y*cd->num_x) + x;
			for(int cont_idx=0; cont_idx<cont->npts; cont_idx++) {
				int v2id = cont->point_ids[cont_idx];
				vertex_t v2 = ring->pts[v2id];
				double d = vert_dist(v1, v2);
				if(d > radius) continue;
				out.point_ids = (int *)realloc_or_die(out.point_ids,
					sizeof(int) * (out.npts + 1));
				out.point_ids[out.npts++] = v2id;
			}
		}
	}
	qsort(out.point_ids, out.npts, sizeof(int), intcompare);
	return out;
}

static int find_next_colinear(ring_t *ring, chopdata_t *cd, int v1id, double radius) {
	vertex_t v1 = ring->pts[v1id];
	int v1id_left = (v1id-1+ring->npts) % ring->npts;

	// FIXME - free list
	box_contents_t neig = get_neighbors(ring, cd, v1id, radius);
	//printf("neig =");
	//for(int neig_idx=0; neig_idx<neig.npts; neig_idx++) {
	//	int v2id = neig.point_ids[neig_idx];
	//	printf(" %d", v2id);
	//}
	//printf("\n");
	for(int neig_idx=0; neig_idx<neig.npts; neig_idx++) {
		int v2id = neig.point_ids[neig_idx];
		if(v2id <= v1id) continue;
		vertex_t v2 = ring->pts[v2id];
		int v2id_right = (v2id+1) % ring->npts;
		int is_colin = is_almost_colinear(ring, v1id_left, v1id, v2id, v2id_right);
		if(is_colin) return v2id;
	}
	return -1;
}

static double calc_ring_perimeter(ring_t *ring) {
	double perim = 0;
	for(int i=0; i<ring->npts; i++) {
		int i2 = (i+1) % ring->npts;
		perim += vert_dist(ring->pts[i], ring->pts[i2]);
	}
	return perim;
}

static double calc_string_perimeter(ring_t *ring, int from, int to) {
	double perim = 0;
	int i = from;
	while(i != to) {
		int i2 = (i+1) % ring->npts;
		perim += vert_dist(ring->pts[i], ring->pts[i2]);
		i = i2;
	}
	return perim;
}

static ring_t pinch_excursions_once(ring_t *ring) {
	double radius = 200; // FIXME
	chopdata_t cd = chop_ring_into_cells(ring, (int)(radius / 5.0));
	ring_t out = duplicate_ring(ring);
	out.npts = 0;

	double ring_perim = calc_ring_perimeter(ring);
	
	for(int vidx=0; vidx<ring->npts; vidx++) {
		out.pts[out.npts++] = ring->pts[vidx];
		int vnext = find_next_colinear(ring, &cd, vidx, radius);
		if(vnext < 0) continue;
		double string_perim = calc_string_perimeter(ring, vidx, vnext);
		double chord = vert_dist(ring->pts[vidx], ring->pts[vnext]);
		// don't allow excursions longer than a semicircle
		if(string_perim < ring_perim/2 && string_perim > chord * PI/2.0) {
			printf("pinch %d:(%lf,%lf) .. %d:(%lf,%lf) - (%lf and %lf)\n", 
				vidx, ring->pts[vidx].x, ring->pts[vidx].y,
				vnext, ring->pts[vnext].x, ring->pts[vnext].y,
				string_perim, chord);
			if(vnext > vidx) {
				vidx = vnext;
				out.pts[out.npts++] = ring->pts[vidx];
			} else break; // FIXME - handle excursions that cross idx==0
		}
	}

	out.pts = (vertex_t *)realloc_or_die(out.pts, sizeof(vertex_t) * out.npts);

	return out;
}

ring_t pinch_excursions(ring_t *ring) {
	ring_t last;
	ring_t next = *ring;
	do {
		last = next;
		next = pinch_excursions_once(&last);
		printf("pinched %d => %d pts\n", last.npts, next.npts);
	} while(last.npts != next.npts);
	return next;
}

/////////////////////////// version 2

/*
int polygon_orientation(ring_t *c) {
	double accum = 0;
	int i;
	for(i=0; i<c->npts; i++) {
		double x0 = c->pts[i].x;
		double y0 = c->pts[i].y;
		double x1 = c->pts[(i+1)%c->npts].x;
		double y1 = c->pts[(i+1)%c->npts].y;
		accum += x0*y1 - x1*y0;
	}
	return accum > 0;
}

int is_good_corner(vertex_t *v1, vertex_t *v2, vertex_t *v3) {
	double d21x = v2->x - v1->x; double d21y = v2->y - v1->y;
	double d32x = v3->x - v2->x; double d32y = v3->y - v2->y;
	double len1 = sqrt(d21x*d21x + d21y*d21y);
	double len2 = sqrt(d32x*d32x + d32y*d32y);
	double ca = (d21x*d32x + d21y*d32y) / len1 / len2;
	return ca > 0 || len1 > 100 || len2 > 100;
}

ring_t pinch_excursions(ring_t *ring) {
	ring_t in = duplicate_ring(ring);
	ring_t out = duplicate_ring(ring);
	int did_work = 1;
	while(did_work) {
		// swap buffers
		ring_t rtmp = out; out = in; in = rtmp;

		did_work = 0;
		out.npts = 0;
		int vidx;
		for(vidx=0; vidx<in.npts; vidx++) {
			int vleft  = (vidx-1+in.npts) % in.npts;
			int vright = (vidx+1+in.npts) % in.npts;
			int good = is_good_corner(in.pts+vleft, in.pts+vidx, in.pts+vright);
			if(good) {
				out.pts[out.npts++] = in.pts[vidx];
			} else {
				did_work = 1;
			}
		}
		printf("pinched %d => %d pts\n", in.npts, out.npts);
	}
	return out;
}
*/
