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

#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a)<(b)?(a):(b))

#ifndef bool
#define uint8_t bool
#endif
#ifndef false
#define false 0
#endif
#ifndef true
#define true (!false)
#endif

/*
static inline double sqlen(vertex_t v0, vertex_t v1) {
	double dx = v1.x - v0.x;
	double dy = v1.y - v0.y;
	return dx*dx + dy*dy;
}

static int is_seg_good(ring_t *ring, int v0_idx) {
	double radius = 10.0; // FIXME
	double rad2 = radius*radius;

	vertex_t *pts = ring->pts;
	int npts = ring->npts;

	int vl_idx = v0_idx;
	int vr_idx = (v0_idx+1) % npts;
	vertex_t v0 = pts[vl_idx];
	vertex_t v1 = pts[vr_idx];
	vertex_t center;
	center.x = (v0.x + v1.x) / 2.0;
	center.y = (v0.y + v1.y) / 2.0;

	for(;;) {
		if(sqlen(center, pts[vr_idx]) > rad2) break;
		vr_idx = (vr_idx+1) % npts;
		if(vr_idx == vl_idx) return 0; // FIXME
	}
	for(;;) {
		if(sqlen(center, pts[vl_idx]) > rad2) break;
		vl_idx = (vl_idx+npts-1) % npts;
		if(vr_idx == vl_idx) return 0; // FIXME
	}

	double up=0, dn=0, lf=0, rt=0;
	for(int i=vl_idx; i!=vr_idx; i=(i+1)%npts) {
		vertex_t v0 = pts[i];
		vertex_t v1 = pts[(i+1)%npts];
		double dx = v1.x - v0.x;
		double dy = v1.y - v0.y;
		if(dx<0) lf += -dx;
		if(dx>0) rt +=  dx;
		if(dy<0) dn += -dy;
		if(dy>0) up +=  dy;
	}
	double loopback = MIN(up,dn) + MIN(lf,rt);

	return loopback == 0;
}

ring_t pinch_excursions(ring_t *ring, report_image_t *dbuf) {
	for(int i=0; i<ring->npts; i++) {
		vertex_t v0 = ring->pts[i];
		vertex_t v1 = ring->pts[(i+1)%ring->npts];
		int good = is_seg_good(ring, i);
		if(good) plot_line(dbuf, v0, v1, 255, 0, 0);
	}
	return *ring;
}
*/

static inline double seg_ang(vertex_t v0, vertex_t v1) {
	double dx = v1.x - v0.x;
	double dy = v1.y - v0.y;
	return atan2(dy, dx);
}

static int find_bottom_pt(ring_t *ring) {
	double min = 0;
	int min_idx = -1;
	for(int i=0; i<ring->npts; i++) {
		double y = ring->pts[i].y;
		if(i==0 || y<min) {
			min = y;
			min_idx = i;
		}
	}
	return min_idx;
}

static int find_next_convex(ring_t *ring, int start_idx, int limit_idx, double start_ang, double *ang_out) {
	int npts = ring->npts;
	vertex_t *pts = ring->pts;
	vertex_t v0 = pts[start_idx];
	double min_angdiff = 2.0 * M_PI;
	int best_vert = -1;
	double best_segang = 0;
	for(int i=(start_idx+1)%npts; i!=limit_idx; i=(i+1)%npts) {
		vertex_t v1 = pts[i];
		double segang = seg_ang(v0, v1);
		double angdiff = segang - start_ang;
		while(angdiff < 0) angdiff += 2.0 * M_PI;
		while(angdiff >= 2.0 * M_PI) angdiff -= 2.0 * M_PI;
		if(angdiff < min_angdiff) {
			min_angdiff = angdiff;
			best_vert = i;
			best_segang = segang;
		}
	}
	if(ang_out) *ang_out = best_segang;
	return best_vert;
}

static bool *find_chull(ring_t *ring) {
	bool *keep = (bool *)malloc_or_die(sizeof(bool) * ring->npts);
	for(int i=0; i<ring->npts; i++) keep[i] = false;

	int start_idx = find_bottom_pt(ring);
	keep[start_idx] = true;
	double ang = -M_PI;
	int idx = start_idx;
	for(;;) {
		idx = find_next_convex(ring, idx, start_idx, ang, &ang);
		if(idx < 0) break;
		keep[idx] = true;
	}

	return keep;
}

static double subring_area(ring_t *ring, int from, int to) {
	int npts = ring->npts;
	vertex_t *pts = ring->pts;

	double accum = 0;
	int i;
	for(i=from; ; i=(i+1)%npts) {
		int i2 = i==to ? from : (i+1)%npts;
		double x0 = pts[i].x;
		double y0 = pts[i].y;
		double x1 = pts[i2].x;
		double y1 = pts[i2].y;
		accum += x0*y1 - x1*y0;
		if(i == to) break;
	}
	return accum / 2.0;
}

static void refine(ring_t *ring, bool *keep) {
	int npts = ring->npts;
	vertex_t *pts = ring->pts;
	for(int i=0; i<npts; i++) {
		if(!keep[i]) continue;
		for(;;) {
			int j;
			for(j=(i+1)%npts; j!=i; j=(j+1)%npts) {
				if(keep[j]) break;
			}
			double area = subring_area(ring, i, j);
			area = fabs(area); // FIXME
			printf("area = %g\n", area);
			if(area > 1000) {
				if(j>i) keep[(i+j)/2] = true;
				else keep[((i+j+npts)/2)%npts] = true;
			} else {
				break;
			}
		}
	}
}

ring_t pinch_excursions(ring_t *ring, report_image_t *dbuf) {
	int npts = ring->npts;
	vertex_t *pts = ring->pts;

	if(!ring_is_ccw(ring)) {
		// reverse ring to make it CCW
		for(int i=0; i<npts/2; i++) {
			vertex_t tmp = pts[i];
			pts[i] = pts[npts-1-i];
			pts[npts-1-i] = tmp;
		}
	}

	bool *keep = find_chull(ring);

	refine(ring, keep);

	int nkeep = 0;
	for(int i=0; i<npts; i++) {
		if(keep[i]) nkeep++;
	}

	ring_t outring = *ring;
	outring.npts = nkeep;
	outring.pts = (vertex_t *)malloc_or_die(sizeof(vertex_t) * outring.npts);
	int oidx=0;
	for(int i=0; i<npts; i++) {
		if(keep[i]) outring.pts[oidx++] = pts[i];
	}

	for(int i=0; i<outring.npts; i++) {
		vertex_t v0 = outring.pts[i];
		vertex_t v1 = outring.pts[(i+1)%outring.npts];
		plot_line(dbuf, v0, v1, 255, 0, 0);
	}

	return outring;
}
