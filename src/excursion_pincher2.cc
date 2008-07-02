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
	double min_angdiff = M_PI;
	int best_vert = -1;
	double best_segang = 0;
	for(int i=(start_idx+1)%npts; ; i=(i+1)%npts) {
		if(i == start_idx) break;
		vertex_t v1 = pts[i];
		double segang = seg_ang(v0, v1);
//printf("start_ang=%g*PI, seg_ang=%g*PI\n", start_ang/M_PI, segang/M_PI);
		double angdiff = segang - start_ang;
		while(angdiff <= -M_PI) angdiff += 2.0 * M_PI;
		while(angdiff >   M_PI) angdiff -= 2.0 * M_PI;
		if(angdiff < min_angdiff) {
			min_angdiff = angdiff;
			best_vert = i;
			best_segang = segang;
		}
		if(i == limit_idx) break;
	}
	if(best_vert < 0) {
		return limit_idx;
	} else if(min_angdiff < 0) {
		//printf("point on wrong side of half-plane (ang=%g*PI) idx=%d\n", min_angdiff/M_PI, best_vert);
		return -1;
	} else {
		if(ang_out) *ang_out = best_segang;
		return best_vert;
	}
}

static bool *find_chull(ring_t *ring) {
	bool *keep = (bool *)malloc_or_die(sizeof(bool) * ring->npts);
	for(int i=0; i<ring->npts; i++) keep[i] = false;

	int start_idx = find_bottom_pt(ring);
	keep[start_idx] = true;
	double ang = 0;
	int idx = start_idx;
	for(;;) {
		idx = find_next_convex(ring, idx, start_idx, ang, &ang);
		if(idx < 0) fatal_error("could not get convex hull");
		if(idx == start_idx) break;
		keep[idx] = true;
	}

	return keep;
}

static inline double seg_len(vertex_t v0, vertex_t v1) {
	double dx = v1.x - v0.x;
	double dy = v1.y - v0.y;
	return sqrt(dx*dx + dy*dy);
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
		accum += x1*y0 - x0*y1;
		if(i == to) break;
	}
	//if(fabs(accum - round(accum)) > 1e-9) {
	//	fatal_error("accum not integer in subring_area (%g)", accum);
	//}
	if(accum < 0) {
		fatal_error("subring_area was negative");
	}
	return accum / 2.0;
}

static int next_keep(int npts, bool *keep, int i) {
	for(int j=(i+1)%npts; j!=i; j=(j+1)%npts) {
		if(keep[j]) return j;
	}
	return i;
}

static int prev_keep(int npts, bool *keep, int i) {
	for(int j=(i+npts-1)%npts; j!=i; j=(j+npts-1)%npts) {
		if(keep[j]) return j;
	}
	return i;
}

static int reach_point(ring_t *ring, bool *keep, int from, int to, double ang) {
	int npts = ring->npts;
	vertex_t *pts = ring->pts;

	int idx = from;
	for(;;) {
		idx = find_next_convex(ring, idx, to, ang, &ang);
		if(idx < 0) return 1;
		keep[idx] = true;
		if(idx == to) break;
	}
	for(int pk=from;;) {
		int nk = next_keep(npts, keep, pk);
		if(nk == to) return 0;

		double min_x = MIN(pts[pk].x, pts[nk].x);
		double max_x = MAX(pts[pk].x, pts[nk].x);
		double min_y = MIN(pts[pk].y, pts[nk].y);
		double max_y = MAX(pts[pk].y, pts[nk].y);

		// FIXME - this is the slowest part
		// FIXME - also test for crossings between keep-segs
		for(int i=0; i<npts; i++) {
			int i2 = i+1<npts ? i+1 : 0;
			vertex_t p1 = pts[i];
			vertex_t p2 = pts[i2];
			// this test is copied here from line_intersects_line to save
			// the time penalty of a function call for the common case of
			// the bbox not matching
			if(
				max_x < MIN(p1.x, p2.x) ||
				min_x > MAX(p1.x, p2.x) ||
				max_y < MIN(p1.y, p2.y) ||
				min_y > MAX(p1.y, p2.y)
			) continue;
			if(i==pk || i==nk || i2==pk || i2==nk) continue;
			if(line_intersects_line(pts[pk], pts[nk], p1, p2, 0)) return 1;
		}

		pk = nk;
	}
}

/*
static int chord_crosses_arc(ring_t *ring, int c0, int c1, int from, int to) {
	int npts = ring->npts;
	vertex_t *pts = ring->pts;

	for(int i=from; i!=to; i=(i+1)%npts) {
		int i2 = (i+1)%npts;
		if(i==c0 || i==c1 || i2==c0 || i2==c1) continue;
		if(line_intersects_line(pts[c0], pts[c1], pts[i], pts[i2], 0)) return 1;
	}
	return 0;
}
*/

static int refine_seg(ring_t *ring, bool *keep_orig, int from, int to) {
	int npts = ring->npts;
	vertex_t *pts = ring->pts;

	double start_area = subring_area(ring, from, to);
	double start_perim = seg_len(pts[from], pts[to]);

	bool *keep_new = (bool *)malloc_or_die(sizeof(bool) * npts);
	bool *keep_best = (bool *)malloc_or_die(sizeof(bool) * npts);

	//int mid; // FIXME
	//if(to>from) mid = (from+to)/2;
	//else mid = ((from+to+npts)/2)%npts;

	double best_improvement = 0;
	int best_touchpt = -1;

	for(int testpt=(from+1)%npts; testpt!=to; testpt=(testpt+1)%npts) {
		//if(testpt != mid) continue; // FIXME
		memcpy(keep_new, keep_orig, sizeof(bool) * npts);
		keep_new[testpt] = true;

		double ang = seg_ang(pts[from], pts[to]);
		int error = reach_point(ring, keep_new, from, testpt, ang);
		if(error) continue;

		int pk = prev_keep(npts, keep_new, testpt);
		if(pk == testpt) fatal_error("pk == testpt");
		ang = seg_ang(pts[testpt], pts[pk]);
		error = reach_point(ring, keep_new, testpt, to, ang);
		if(error) continue;

		double area = 0;
		double perim = 0;
		for(int pk=from;;) {
			int nk = next_keep(npts, keep_new, pk);
			area += subring_area(ring, pk, nk);
			perim += seg_len(pts[pk], pts[nk]);
			if(nk == to) break;
			pk = nk;
		}

		//double improvement = (start_area - area) - 50.0*(perim - start_perim);
		double improvement = ((start_area+2.0) / (area+2.0)) / pow(perim / start_perim, 2) - 2.0;
		//double improvement = (start_area - area) / pow(perim - start_perim, .5) - 100.0;
		//double improvement = (start_area - area) / pow(perim - start_perim, 2) - 100.0;

		if(VERBOSE >= 2) printf("improvement=%g, start_area=%g, area=%g, perim=%g, start_perim=%g\n",
			improvement, start_area, area, perim, start_perim);

		if(improvement > best_improvement) {
			best_improvement = improvement;
			memcpy(keep_best, keep_new, sizeof(bool) * npts);
			best_touchpt = testpt;
		}
	}
	if(best_improvement) {
		if(VERBOSE) {
			printf("best_improvement = %g\n", best_improvement);
			printf("tagged %d (%g,%g) as keep between %d and %d\n", best_touchpt,
				pts[best_touchpt].x, pts[best_touchpt].y, from, to);
			for(int i=0; i<npts; i++) {
				if(keep_best[i] && !keep_orig[i]) printf(
					"  rubberband touches %d (%g, %g)\n", i, pts[i].x, pts[i].y);
			}
		}

		memcpy(keep_orig, keep_best, sizeof(bool) * npts);
		return best_touchpt;
	} else {
		return -1;
	}
}


static void refine_ring(ring_t *ring, bool *keep, bool *touchpts) {
	int npts = ring->npts;
	for(int i=0; i<npts; i++) {
		if(!keep[i]) continue;
		for(;;) {
			int j = next_keep(npts, keep, i);
			double area = subring_area(ring, i, j);
			//printf("area = %g\n", area);
			if(area > 0) {
				int touchpt = refine_seg(ring, keep, i, j);
				if(touchpt < 0) break;
				touchpts[touchpt] = true;
			} else {
				break;
			}
		}
	}
}

static ring_t pinch_ring_excursions(ring_t *ring, report_image_t *dbuf) {
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
	bool *touchpts = (bool *)malloc_or_die(sizeof(bool) * ring->npts);
	for(int i=0; i<ring->npts; i++) touchpts[i] = false;

	refine_ring(ring, keep, touchpts);

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

	if(dbuf) {
		debug_plot_ring(dbuf, &outring, 128, 0, 0);
		for(int i=0; i<npts; i++) {
			if(touchpts[i]) {
				vertex_t p = pts[i];
				plot_point_big(dbuf, p.x, p.y, 255, 255, 255);
			}
		}
	}

	return outring;
}

mpoly_t pinch_excursions2(mpoly_t *mp_in, report_image_t *dbuf) {
	mpoly_t mp_out;
	mp_out.num_rings = mp_in->num_rings;
	mp_out.rings = (ring_t *)malloc_or_die(sizeof(ring_t) * mp_out.num_rings);
	for(int r_idx=0; r_idx<mp_in->num_rings; r_idx++) {
		mp_out.rings[r_idx] = pinch_ring_excursions(mp_in->rings+r_idx, dbuf);
	}
	// FIXME - fix topology using functions from dp.c
	return mp_out;
}
