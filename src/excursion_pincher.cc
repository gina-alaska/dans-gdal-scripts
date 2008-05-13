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

static int is_almost_colinear(ring_t *ring, int v1_idx, int v2_idx, int v3_idx, int v4_idx) {
	double numer_a = (p4.x-p3.x)*(p1.y-p3.y) - (p4.y-p3.y)*(p1.x-p3.x);
	double numer_b = (p2.x-p1.x)*(p1.y-p3.y) - (p2.y-p1.y)*(p1.x-p3.x);
	double denom   = (p4.y-p3.y)*(p2.x-p1.x) - (p4.x-p3.x)*(p2.y-p1.y);
	if(denom == 0) {
		if(numer_a==0 && numer_b==0) { // coincident
		}
	}
}

static int find_next_colinear(ring_t *ring, chopdata_t *cd, int v1id, double radius) {
	vertex_t v1 = ring->pts[v1id];
	int v1id_left = (v1id-1+ring->npts) % ring->npts;
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
				if(v2id <= v1id) continue;
				vertex_t v2 = ring->pts[v2id];
				double d = vert_dist(v1, v2);
				if(d > radius) continue;
				int v2id_right = (v2id+1) % ring->npts;
				int is_colin = is_almost_colinear(ring, v1id_left, v1id, v2id, v2id_right);
				if(is_colin) return v2id;
			}
		}
	}
	return -1;
}
