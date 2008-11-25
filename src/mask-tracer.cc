/*
Copyright (c) 2008, Regents of the University of Alaska

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
#include "polygon-rasterizer.h"

#define DIR_UP 0
#define DIR_RT 1
#define DIR_DN 2
#define DIR_LF 3

typedef int pixquad_t;

int dbg_idx = 0;
static void debug_write_mask(uint8_t *mask, int w, int h) {
	char fn[1000];
	sprintf(fn, "zz-debug-%04d.pgm", dbg_idx++);
	FILE *fh = fopen(fn, "w");
	if(!fh) fatal_error("cannot open %s", fn);
	fprintf(fh, "P5\n%d %d\n255\n", w, h);
	for(int y=0; y<w; y++)
		fwrite(mask+(w+2)*(y+1)+1, w, 1, fh);
	fclose(fh);
}

static ring_t *make_enclosing_ring(int w, int h) {
	ring_t *ring = (ring_t *)malloc_or_die(sizeof(ring_t));
	ring->npts = 4;
	ring->pts = (vertex_t *)malloc_or_die(sizeof(vertex_t) * ring->npts);
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

static long compute_area(row_crossings_t *crossings, int num_rows) {
	long area = 0;
	for(int y=0; y<num_rows; y++) {
		int nc = crossings[y].num_crossings;
		int *rc = crossings[y].crossings;
		for(int cidx=0; cidx<nc/2; cidx++) {
			int from = rc[cidx*2  ];
			int to   = rc[cidx*2+1];
			area += to - from;
		}
	}
	return area;
}

/*
static int is_inside_crossings(row_crossings_t *c, int x) {
	int inside = 0;
	for(int i=0; i<c->num_crossings; i++) {
		if(x >= c->crossings[i]) inside = !inside;
	}
	return inside;
}
*/

static inline pixquad_t get_quad(uint8_t *mask, int w, int h, int x, int y, int select_color) {
	// 1 2
	// 8 4
	uint8_t *uprow = mask + (y  )*(w+2);
	uint8_t *dnrow = mask + (y+1)*(w+2);
	pixquad_t quad =
		(uprow[x  ]     ) + // y-1, x-1
		(uprow[x+1] << 1) + // y-1, x
		(dnrow[x+1] << 2) + // y  , x
		(dnrow[x  ] << 3);  // y  , x-1
	if(!select_color) quad ^= 0xf;
	return quad;
}

static inline pixquad_t rotate_quad(pixquad_t q, int dir) {
	return ((q + (q<<4)) >> dir) & 0xf;
}

static ring_t trace_single_mpoly(uint8_t *mask, int w, int h,
int initial_x, int initial_y, int select_color) {
	//printf("trace_single_mpoly enter (%d,%d)\n", initial_x, initial_y);

	ring_t ring;
	int ringbuf_size = 4;
	ring.pts = (vertex_t *)malloc_or_die(sizeof(vertex_t) * ringbuf_size);
	ring.pts[0].x = initial_x;
	ring.pts[0].y = initial_y;
	ring.npts = 1;

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

		switch(dir) {
			case DIR_UP: y -= 1; break;
			case DIR_RT: x += 1; break;
			case DIR_DN: y += 1; break;
			case DIR_LF: x -= 1; break;
			default: fatal_error("bad direction");
		}
		if(x == initial_x && y == initial_y) break;
		if(x<0 || y<0 || x>w || y>h) fatal_error("fell off edge (%d,%d)", x, y);
		pixquad_t quad = get_quad(mask, w, h, x, y, select_color);
		quad = rotate_quad(quad, dir);
		if((quad & 12) != 4) fatal_error("tracer was not on the right side of things (%d)", quad);
		int rot;
		switch(quad & 3) {
			case 0: rot =  1; break; // N N
			case 1: rot =  1; break; // Y N
			case 2: rot =  0; break; // N Y
			case 3: rot = -1; break; // Y Y
			default: fatal_error("not possible");
		}
		dir = (dir + rot + 4) % 4;

		if(rot) {
			if(ring.npts == ringbuf_size) {
				if(ringbuf_size) ringbuf_size *= 2;
				else ringbuf_size = 4;
				ring.pts = (vertex_t *)realloc_or_die(ring.pts, 
					sizeof(vertex_t) * ringbuf_size);
			}
			ring.pts[ring.npts].x = x;
			ring.pts[ring.npts].y = y;
			ring.npts++;
		}
	}

	if(ringbuf_size > ring.npts) {
		ring.pts = (vertex_t *)realloc_or_die(ring.pts, 
			sizeof(vertex_t) * ring.npts);
	}

	return ring;
}

static int recursive_trace(uint8_t *mask, int w, int h,
ring_t *bounds, int depth, mpoly_t *out_poly, int parent_id, 
long min_area, int no_donuts) {
	//printf("recursive_trace enter: depth=%d\n", depth);

	int select_color = (depth & 1) ? 0 : 1;

	int must_free_bounds = 0;
	if(!bounds) {
		bounds = make_enclosing_ring(w, h);
		must_free_bounds = 1;
	}

	int bound_left, bound_right;
	int bound_top, bound_bottom;
	bound_left = bound_right = (int)bounds->pts[0].x;
	bound_top = bound_bottom = (int)bounds->pts[0].y;
	for(int v_idx=0; v_idx<bounds->npts; v_idx++) {
		int x = (int)bounds->pts[v_idx].x;
		int y = (int)bounds->pts[v_idx].y;
		if(x < bound_left) bound_left = x;
		if(x > bound_right) bound_right = x;
		if(y < bound_top) bound_top = y;
		if(y > bound_bottom) bound_bottom = y;
	}

	mpoly_t bounds_mp;
	bounds_mp.num_rings = 1;
	bounds_mp.rings = bounds;

	row_crossings_t *crossings = get_row_crossings(&bounds_mp, bound_top, bound_bottom-bound_top);
	int skip_this = min_area && (compute_area(crossings, bound_bottom-bound_top) < min_area);
	int skip_child = skip_this || (depth && no_donuts);

	row_crossings_t cross_both;
	cross_both.array_size = 0;
	cross_both.crossings = NULL;

	if(!depth) {
		printf("Tracing: ");
		GDALTermProgress(0, NULL, NULL);
	}

	if(!skip_child) {
		for(int y=bound_top+1; y<bound_bottom; y++) {
			if(!depth) {
				GDALTermProgress((double)y/(double)(bound_bottom-bound_top-1), NULL, NULL);
			}

			uint8_t *uprow = mask + (y  )*(w+2);
			uint8_t *dnrow = mask + (y+1)*(w+2);

			// make sure the range (y-1,y)*(x-1,x) is in bounds
			crossings_intersection(&cross_both, 
				crossings + (y-bound_top-1),
				crossings + (y-bound_top));
			for(int cidx=0; cidx<cross_both.num_crossings/2; cidx++) {
				// make sure the range (y-1,y)*(x-1,x) is in bounds
				int from = 1+cross_both.crossings[cidx*2  ];
				int to   =   cross_both.crossings[cidx*2+1];

				// find the first possible quad that could match
				uint8_t *mc1 = (uint8_t *)memchr(uprow+from, select_color, to-from);
				uint8_t *mc2 = (uint8_t *)memchr(dnrow+from, select_color, to-from);
				int ic1 = mc1 ? (mc1-uprow)-1 : to;
				int ic2 = mc2 ? (mc2-dnrow)-1 : to;
				if(ic2 < ic1) ic1 = ic2;
				if(ic1 > from) from = ic1;

				for(int x=from; x<to; x++) {
					//pixquad_t quad = get_quad(mask, w, h, x, y, select_color);
					//int is_seed = (quad != 0 && quad != 0xf);
					int is_seed;
					if(select_color) {
						is_seed = 
							uprow[x  ] || // y-1, x-1
							uprow[x+1] || // y-1, x
							dnrow[x+1] || // y  , x
							dnrow[x  ];   // y  , x-1
					} else {
						is_seed = 
							(!uprow[x  ]) || // y-1, x-1
							(!uprow[x+1]) || // y-1, x
							(!dnrow[x+1]) || // y  , x
							(!dnrow[x  ]);   // y  , x-1
					}

					if(is_seed) {
						ring_t r = trace_single_mpoly(mask, w, h, x, y, select_color);

						r.parent_id = parent_id;
						r.is_hole = depth % 2;
						//r.parent_id = -1;
						//r.is_hole = 0;

						out_poly->rings = (ring_t *)realloc_or_die(out_poly->rings,
							sizeof(ring_t) * (out_poly->num_rings + 1));
						int outer_ring_id = (out_poly->num_rings++);
						out_poly->rings[outer_ring_id] = r;

						int was_skip = recursive_trace(
							mask, w, h, &r, depth+1, out_poly, outer_ring_id,
							min_area, no_donuts);
						if(was_skip) {
							free(r.pts);
							out_poly->num_rings--;
						}
					}
				} 
			}
		}
		if(cross_both.crossings != NULL) free(cross_both.crossings);
	}

	for(int y=bound_top; y<bound_bottom; y++) {
		row_crossings_t *r = crossings + (y-bound_top);
		if(depth>0) {
			uint8_t *maskrow = mask + (y+1)*(w+2);
			for(int cidx=0; cidx<r->num_crossings/2; cidx++) {
				int from = r->crossings[cidx*2  ];
				int to   = r->crossings[cidx*2+1];
				memset(maskrow+from+1, select_color, to-from);
			}
		}
	}
	free_row_crossings(crossings, bound_bottom-bound_top);

	if(must_free_bounds) {
		free(bounds->pts);
		free(bounds);
	}

	if(VERBOSE >= 4) debug_write_mask(mask, w, h);

	if(!depth) {
		GDALTermProgress(1, NULL, NULL);
	}

	return skip_this;
}

// this function has the side effect of erasing the mask
mpoly_t trace_mask(uint8_t *mask_8bit, int w, int h, long min_area, int no_donuts) {
	if(VERBOSE >= 4) debug_write_mask(mask_8bit, w, h);

	/*
	uint8_t *mask_8bit = (uint8_t *)malloc_or_die((w+2)*(h+2));
	// FIXME - only need to clear borders
	memset(mask_8bit, 0, (w+2)*(h+2));
	for(int y=0; y<h; y++) {
		int mask_rowlen = (w+7)/8;
		uint8_t mask_bitp = 1;
		uint8_t *mask_bytep = mask_1bit + mask_rowlen*y;
		uint8_t *outp = mask_8bit + (y+1)*(w+2) + 1;
		for(int x=0; x<w; x++) {
			*(outp++) = (*mask_bytep & mask_bitp) ? 1 : 0;
			mask_bitp <<= 1;
			if(!mask_bitp) {
				mask_bitp = 1;
				mask_bytep++;
			}
		}
	}
	*/

	mpoly_t out_poly;

	out_poly.num_rings = 0;
	out_poly.rings = NULL;

	recursive_trace(mask_8bit, w, h, NULL, 0, &out_poly, -1, min_area, no_donuts);
	printf("Trace found %d rings.\n", out_poly.num_rings);

	//free(mask_8bit);

	//fatal_error("OK");

	return out_poly;
}
