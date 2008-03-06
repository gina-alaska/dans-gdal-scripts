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
#include "debugplot.h"
#include "georef.h"
#include "mask.h"

#define DIR_UP 0
#define DIR_RT 1
#define DIR_DN 2
#define DIR_LF 3

typedef struct {
	int x, y;
} intvert_t;

typedef struct {
	int npts;
	intvert_t *pts;
} intring_t;

typedef struct {
	int num_crossings;
	int array_size;
	int *crossings;
} row_crossings_t;

typedef int pixquad_t;

static intring_t *make_enclosing_ring(int w, int h) {
	intring_t *ring = (intring_t *)malloc_or_die(sizeof(intring_t));
	ring->npts = 4;
	ring->pts = (intvert_t *)malloc_or_die(sizeof(intvert_t *) * ring->npts);
	ring->pts[0].x = -1;
	ring->pts[0].y = -1;
	ring->pts[1].x = w;
	ring->pts[1].y = -1;
	ring->pts[2].x = w;
	ring->pts[2].y = h;
	ring->pts[3].x = -1;
	ring->pts[3].y = h;
	return ring;
}

static inline row_crossings_t *get_row_crossings(intring_t *ring, int min_y, int num_rows) {
	row_crossings_t *rows = (row_crossings_t *)malloc_or_die(sizeof(row_crossings_t) * num_rows);

	for(int i=0; i<num_rows; i++) {
		rows[i].num_crossings = 0;
		rows[i].array_size = 0;
		rows[i].crossings = NULL;
	}

	for(int v_idx=0; v_idx<ring->npts; v_idx++) {
		int x0 = ring->pts[v_idx].x;
		int y0 = ring->pts[v_idx].y;
		int x1 = ring->pts[(v_idx+1)%ring->npts].x;
		int y1 = ring->pts[(v_idx+1)%ring->npts].y;
		if(y0 == y1) continue;
		if(y0 > y1) {
			int tmp;
			tmp=x0; x0=x1; x1=tmp; 
			tmp=y0; y0=y1; y1=tmp; 
		}
		if(x0 != x1) fatal_error("segment was not horizontal or vertical");
		for(int y=y0; y<y1; y++) {
			int row = y - min_y;
			if(row<0 || row>num_rows-1) continue;
			row_crossings_t *r = rows+row;
			if(r->num_crossings == r->array_size) {
				r->array_size += 16;
				r->crossings = (int *)realloc_or_die(r->crossings,
					sizeof(int) * r->array_size);
			}
			r->crossings[r->num_crossings++] = x0;
		}
	}

	return rows;
}

static inline int is_inside_crossings(row_crossings_t *c, int x) {
	int inside = 0;
	for(int i=0; i<c->num_crossings; i++) {
		if(x >= c->crossings[i]) inside = !inside;
	}
	return inside;
}

static inline int get_pixel(unsigned char *mask, int w, int h, int x, int y) {
	if(x<0 || x>=w || y<0 || y>=h) return 0;
	int mask_rowlen = (w+7)/8;
	unsigned char mask_bitp = 1 << (x % 8);
	unsigned char *mask_bytep = mask + mask_rowlen*y + x/8;
	int val = *mask_bytep & mask_bitp;
	return val ? 1 : 0;
}

static inline void set_pixel(unsigned char *mask, int w, int h, int x, int y, int color) {
	if(x<0 || x>=w || y<0 || y>=h) return;
	int mask_rowlen = (w+7)/8;
	unsigned char mask_bitp = 1 << (x % 8);
	unsigned char *mask_bytep = mask + mask_rowlen*y + x/8;
	if(color) {
		*mask_bytep |= mask_bitp;
	} else {
		*mask_bytep &= ~mask_bitp;
	}
}

static inline pixquad_t get_quad(unsigned char *mask, int w, int h, int x, int y, int select_color) {
	// 0 1
	// 3 2
	int quad =
		(get_pixel(mask, w, h, x-1, y-1)     ) +
		(get_pixel(mask, w, h, x  , y-1) << 1) +
		(get_pixel(mask, w, h, x  , y  ) << 2) +
		(get_pixel(mask, w, h, x-1, y  ) << 3);
	if(!select_color) quad ^= 0xf;
	return quad;
}

static inline pixquad_t rotate_quad(pixquad_t q, int dir) {
	while(dir--) {
		q = (q>>1) + ((q&1)<<3);
	}
	return q;
}

int dbg_idx = 0;
static void debug_write_mask(unsigned char *mask, int w, int h) {
	char fn[1000];
	sprintf(fn, "zz-debug-%04d.pgm", dbg_idx++);
	FILE *fh = fopen(fn, "w");
	if(!fh) fatal_error("cannot open %s", fn);
	fprintf(fh, "P5\n%d %d\n255\n", w, h);
	for(int y=0; y<w; y++)
	for(int x=0; x<w; x++) {
		unsigned char pix = get_pixel(mask, w, h, x, y) ? 255 : 0;
		fwrite(&pix, 1, 1, fh);
	}
	fclose(fh);
}

static intring_t trace_single_mpoly(unsigned char *mask, int w, int h, int initial_x, int initial_y, int select_color) {
	//printf("trace_single_mpoly enter (%d,%d)\n", initial_x, initial_y);

	intring_t ring;
	ring.npts = 0;
	ring.pts = NULL;
	int ringbuf_size = 0;

	int x = initial_x;
	int y = initial_y;
	pixquad_t quad = get_quad(mask, w, h, x, y, select_color);
	int dir;
	for(dir=0; dir<4; dir++) {
		pixquad_t rq = rotate_quad(quad, dir);
		if((rq & 3) == 2) break;
	}
	if(dir == 4) fatal_error("couldn't choose a starting direction (q=%d)", quad);
	for(;;) {
		//printf("xy=(%d,%d)\n", x, y);

		if(ring.npts == ringbuf_size) {
			if(ringbuf_size) ringbuf_size *= 2;
			else ringbuf_size = 4;
			ring.pts = (intvert_t *)realloc_or_die(ring.pts, 
				sizeof(intvert_t *) * ringbuf_size);
		}
		ring.pts[ring.npts].x = x;
		ring.pts[ring.npts].y = y;
		ring.npts++;

		switch(dir) {
			case DIR_UP: y -= 1; break;
			case DIR_RT: x += 1; break;
			case DIR_DN: y += 1; break;
			case DIR_LF: x -= 1; break;
			default: fatal_error("bad direction");
		}
		if(x == initial_x && y == initial_y) break;
		pixquad_t quad = get_quad(mask, w, h, x, y, select_color);
		quad = rotate_quad(quad, dir);
		if((quad & 12) != 4) fatal_error("tracer was not on the right side of things");
		int rot;
		switch(quad & 3) {
			case 0: rot =  1; break; // N N
			case 1: rot =  1; break; // Y N
			case 2: rot =  0; break; // N Y
			case 3: rot = -1; break; // Y Y
			default: fatal_error("not possible");
		}
		dir = (dir + rot + 4) % 4;
	}

	if(ringbuf_size > ring.npts) {
		ring.pts = (intvert_t *)realloc_or_die(ring.pts, 
			sizeof(intvert_t *) * ring.npts);
	}

	return ring;
}

static ring_t ring_int2dbl(intring_t *in) {
	ring_t out;
	out.npts = in->npts;
	out.pts = (vertex_t *)malloc_or_die(sizeof(vertex_t) * in->npts);
	for(int i=0; i<in->npts; i++) {
		out.pts[i].x = in->pts[i].x;
		out.pts[i].y = in->pts[i].y;
	}
	return out;
}

static void recursive_trace(unsigned char *mask, int w, int h,
intring_t *bounds, int depth, mpoly_t *out_poly, int parent_id) {
	//printf("recursive_trace enter: depth=%d\n", depth);

	int select_color = (depth & 1) ? 0 : 1;

	if(!bounds) {
		bounds = make_enclosing_ring(w, h);
	}

	int bound_top, bound_bottom;
	bound_top = bound_bottom = bounds->pts[0].y;
	for(int v_idx=0; v_idx<bounds->npts; v_idx++) {
		intvert_t v = bounds->pts[v_idx];
		if(v.y < bound_top) bound_top = v.y;
		if(v.y > bound_bottom) bound_bottom = v.y;
	}

	row_crossings_t *crossings = get_row_crossings(bounds, bound_top, bound_bottom-bound_top);
	for(int y=bound_top+1; y<bound_bottom; y++) {
		if(!depth && !(y%1000)) printf("y=%d\n", y);

		row_crossings_t *c0 = crossings + (y-bound_top-1);
		row_crossings_t *c1 = crossings + (y-bound_top);
		for(int x=0; x<w; x++) {
			if(!is_inside_crossings(c0, x-1)) continue;
			if(!is_inside_crossings(c1, x-1)) continue;
			if(!is_inside_crossings(c0, x)) continue;
			if(!is_inside_crossings(c1, x)) continue;

			pixquad_t quad = get_quad(mask, w, h, x, y, select_color);
			if(quad == 0 || quad == 0xf) continue;

			intring_t outer_ring = trace_single_mpoly(mask, w, h, x, y, select_color);

			ring_t r = ring_int2dbl(&outer_ring);
			r.parent_id = parent_id;
			r.is_hole = depth % 2;
			out_poly->rings = (ring_t *)realloc_or_die(out_poly->rings,
				sizeof(ring_t) * (out_poly->num_rings + 1));
			int outer_ring_id = (out_poly->num_rings++);
			out_poly->rings[outer_ring_id] = r;

			recursive_trace(mask, w, h, &outer_ring, depth+1, out_poly, outer_ring_id);

			free(outer_ring.pts);
		}
	}

	for(int y=bound_top; y<bound_bottom; y++) {
		row_crossings_t *c = crossings + (y-bound_top);
		for(int x=0; x<w; x++) {
			if(!is_inside_crossings(c, x)) continue;
			set_pixel(mask, w, h, x, y, select_color);
		}
		free(c->crossings);
	}
	free(crossings);

	//debug_write_mask(mask, w, h);
}

mpoly_t *calc_ring_from_mask(unsigned char *mask, int w, int h,
report_image_t *dbuf, int major_ring_only, int no_donuts, 
double min_ring_area, double bevel_size) {
	//debug_write_mask(mask, w, h);
	mpoly_t *out_poly = (mpoly_t *)malloc_or_die(sizeof(mpoly_t));

	out_poly->num_rings = 0;
	out_poly->rings = NULL;

	recursive_trace(mask, w, h, NULL, 0, out_poly, -1);
	printf("nr=%d\n", out_poly->num_rings);
	return out_poly;
}
