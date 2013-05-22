/*
Copyright (c) 2013, Regents of the University of Alaska

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

This code was developed by Dan Stahlke for the Geographic Information Network of Alaska.
*/



#include <vector>

#include "mask.h"
#include "mask-tracer.h"
#include "common.h"
#include "polygon.h"
#include "polygon-rasterizer.h"

namespace dangdal {

enum Direction {
	DIR_UP = 0,
	DIR_RT = 1,
	DIR_DN = 2,
	DIR_LF = 3
};

typedef int pixquad_t;

int dbg_idx = 0;
static void debug_write_mask(const BitGrid &mask, size_t w, size_t h) {
	char fn[1000];
	snprintf(fn, sizeof(fn), "zz-debug-%04d.pgm", dbg_idx++);

	uint8_t *row = new uint8_t[w];

	FILE *fh = fopen(fn, "w");
	if(!fh) fatal_error("cannot open %s", fn);
	fprintf(fh, "P5\n%zd %zd\n255\n", w, h);
	for(size_t y=0; y<h; y++) {
		for(size_t x=0; x<w; x++) {
			row[x] = mask(x, y) ? 255 : 0;
		}
		fwrite(row, w, 1, fh);
	}
	fclose(fh);

	delete[] row;
}

static Ring make_enclosing_ring(size_t w, size_t h) {
	Ring ring;
	ring.pts.reserve(4);
	ring.pts.push_back(Vertex(-1, -1));
	ring.pts.push_back(Vertex( w, -1));
	ring.pts.push_back(Vertex( w,  h));
	ring.pts.push_back(Vertex(-1,  h));
	return ring;
}

static int64_t compute_area(const std::vector<row_crossings_t> &crossings) {
	int64_t area = 0;
	for(size_t y=0; y<crossings.size(); y++) {
		const row_crossings_t &rc = crossings[y];
		size_t nc = rc.size();
		for(size_t cidx=0; cidx<nc/2; cidx++) {
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

static inline pixquad_t get_quad(const BitGrid &mask, int x, int y, bool select_color) {
	// 1 2
	// 8 4
	pixquad_t quad =
		(mask(x-1, y-1) ? 1 : 0) +
		(mask(x  , y-1) ? 2 : 0) +
		(mask(x  , y  ) ? 4 : 0) +
		(mask(x-1, y  ) ? 8 : 0);
	if(!select_color) quad ^= 0xf;
	return quad;
}

static inline pixquad_t rotate_quad(pixquad_t q, int dir) {
	return ((q + (q<<4)) >> dir) & 0xf;
}

static Ring trace_single_mpoly(const BitGrid &mask, size_t w, size_t h,
int initial_x, int initial_y, bool select_color) {
	//printf("trace_single_mpoly enter (%d,%d)\n", initial_x, initial_y);

	Ring ring;
	ring.pts.push_back(Vertex(initial_x, initial_y));

	int x = initial_x;
	int y = initial_y;
	pixquad_t quad = get_quad(mask, x, y, select_color);
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
		if(x<0 || y<0 || x>(int)w || y>(int)h) fatal_error("fell off edge (%d,%d)", x, y);
		pixquad_t quad = get_quad(mask, x, y, select_color);
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
			ring.pts.push_back(Vertex(x, y));
		}
	}

	return ring;
}

static int recursive_trace(BitGrid &mask, size_t w, size_t h,
const Ring &bounds, int depth, Mpoly &out_poly, int parent_id, 
int64_t min_area, bool no_donuts) {
	//printf("recursive_trace enter: depth=%d\n", depth);

	bool select_color = !(depth & 1);

	int bound_left, bound_right;
	int bound_top, bound_bottom;
	bound_left = bound_right = (int)bounds.pts[0].x;
	bound_top = bound_bottom = (int)bounds.pts[0].y;
	for(size_t v_idx=0; v_idx<bounds.pts.size(); v_idx++) {
		int x = (int)bounds.pts[v_idx].x;
		int y = (int)bounds.pts[v_idx].y;
		if(x < bound_left) bound_left = x;
		if(x > bound_right) bound_right = x;
		if(y < bound_top) bound_top = y;
		if(y > bound_bottom) bound_bottom = y;
	}

	Mpoly bounds_mp;
	bounds_mp.rings.push_back(bounds);

	std::vector<row_crossings_t> crossings = 
		get_row_crossings(bounds_mp, bound_top, bound_bottom-bound_top);
	assert(crossings.size() == size_t(bound_bottom-bound_top));
	int skip_this = min_area && (compute_area(crossings) < min_area);
	int skip_child = skip_this || (depth && no_donuts);

	if(!depth) {
		printf("Tracing: ");
		GDALTermProgress(0, NULL, NULL);
	}

	if(!skip_child) {
		for(int y=bound_top+1; y<bound_bottom; y++) {
			if(!depth) {
				GDALTermProgress((double)y/(double)(bound_bottom-bound_top-1), NULL, NULL);
			}

			// make sure the range (y-1,y)*(x-1,x) is in bounds
			row_crossings_t cross_both = crossings_intersection(
				crossings[y-bound_top-1], crossings[y-bound_top]);
			for(size_t cidx=0; cidx<cross_both.size()/2; cidx++) {
				// make sure the range (y-1,y)*(x-1,x) is in bounds
				int from = 1+cross_both[cidx*2  ];
				int to   =   cross_both[cidx*2+1];

				for(int x=from; x<to; x++) {
					//pixquad_t quad = get_quad(mask, x, y, select_color);
					//int is_seed = (quad != 0 && quad != 0xf);
					int is_seed;
					if(select_color) {
						is_seed = 
							mask(x-1, y-1) ||
							mask(x  , y-1) ||
							mask(x-1, y  ) ||
							mask(x  , y  );
					} else {
						is_seed = 
							(!mask(x-1, y-1)) ||
							(!mask(x  , y-1)) ||
							(!mask(x-1, y  )) ||
							(!mask(x  , y  ));
					}

					if(is_seed) {
						Ring r = trace_single_mpoly(mask, w, h, x, y, select_color);

						r.parent_id = parent_id;
						r.is_hole = depth % 2;
						//r.parent_id = -1;
						//r.is_hole = 0;

						size_t outer_ring_id = out_poly.rings.size();
						out_poly.rings.push_back(r);

						int was_skip = recursive_trace(
							mask, w, h, r, depth+1, out_poly, outer_ring_id,
							min_area, no_donuts);
						if(was_skip) {
							out_poly.rings.pop_back();
						}
					}
				} 
			}
		}
	}

	for(int y=bound_top; y<bound_bottom; y++) {
		const row_crossings_t &r = crossings[y-bound_top];
		if(depth>0) {
			for(size_t cidx=0; cidx<r.size()/2; cidx++) {
				int from = r[cidx*2  ];
				int to   = r[cidx*2+1];
				for(int x=from; x<=to; x++) {
					if(x>=0 && y>=0 && size_t(x)<w && size_t(y)<h) {
						mask.set(x, y, select_color);
					}
				}
			}
		}
	}

	if(VERBOSE >= 4) debug_write_mask(mask, w, h);

	if(!depth) {
		GDALTermProgress(1, NULL, NULL);
	}

	return skip_this;
}

// this function has the side effect of erasing the mask
Mpoly trace_mask(BitGrid &mask, size_t w, size_t h, int64_t min_area, bool no_donuts) {
	if(VERBOSE >= 4) debug_write_mask(mask, w, h);

	Mpoly out_poly;

	recursive_trace(mask, w, h, make_enclosing_ring(w, h), 0, out_poly, -1, min_area, no_donuts);
	printf("Trace found %zd rings.\n", out_poly.rings.size());

	//free(mask_8bit);

	//fatal_error("OK");

	return out_poly;
}

} // namespace dangdal
