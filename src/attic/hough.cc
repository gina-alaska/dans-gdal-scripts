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
#include "hough.h"

// See http://en.wikipedia.org/wiki/Hough_transform

typedef float accum_t;

typedef struct {
	int w_ang, w_rad;
	accum_t *accum;
	vertex_t center;
	double *sin_table, *cos_table;
} hough_accum_t;

static inline void hough_add_seg(hough_accum_t *ha, vertex_t p1, vertex_t p2) {
	int w_ang = ha->w_ang;
	int w_rad = ha->w_rad;

	p1.x -= ha->center.x;
	p1.y -= ha->center.y;
	p2.x -= ha->center.x;
	p2.y -= ha->center.y;

	double dx = p2.x - p1.x;
	double dy = p2.y - p1.y;
	double weight = sqrt(dx*dx + dy*dy);

	for(int theta_i=0; theta_i<w_ang; theta_i++) {
		accum_t *acc = ha->accum + theta_i * w_rad;
		double sin = ha->sin_table[theta_i];
		double cos = ha->cos_table[theta_i];
		double r1 = p1.x*cos + p1.y*sin + (double)(w_rad/2);
		double r2 = p2.x*cos + p2.y*sin + (double)(w_rad/2);
		if(r2 < r1) { double tmp = r1; r1 = r2; r2 = tmp; } // swap
		int r1_i = (int)floor(r1);
		int r2_i = (int)floor(r2);
		if(r1_i == r2_i) {
			acc[r1_i] += weight;
		} else {
			double w2 = weight / (r2-r1);
			acc[r1_i] += w2 * ((double)(r1_i+1) - r1);
			acc[r2_i] += w2 * (r2 - (double)r2_i);
			for(int i=r1_i+1; i<r2_i; i++) acc[i] += w2;
		}
	}
}

static inline accum_t hough_integrate_seg(hough_accum_t *ha, vertex_t p1, vertex_t p2) {
	int w_ang = ha->w_ang;
	int w_rad = ha->w_rad;

	p1.x -= ha->center.x;
	p1.y -= ha->center.y;
	p2.x -= ha->center.x;
	p2.y -= ha->center.y;

	double dx = p2.x - p1.x;
	double dy = p2.y - p1.y;
	double weight = 1.0; //sqrt(dx*dx + dy*dy);

	accum_t sum = 0;
	for(int theta_i=0; theta_i<w_ang; theta_i++) {
		accum_t *acc = ha->accum + theta_i * w_rad;
		double sin = ha->sin_table[theta_i];
		double cos = ha->cos_table[theta_i];
		double r1 = p1.x*cos + p1.y*sin + (double)(w_rad/2);
		double r2 = p2.x*cos + p2.y*sin + (double)(w_rad/2);
		if(r2 < r1) { double tmp = r1; r1 = r2; r2 = tmp; } // swap
		int r1_i = (int)floor(r1);
		int r2_i = (int)floor(r2);

		accum_t s2 = 0;
		if(r1_i == r2_i) {
			s2 += acc[r1_i] * weight;
		} else {
			double w2 = weight / (r2-r1);
			s2 += acc[r1_i] * w2 * ((double)(r1_i+1) - r1);
			s2 += acc[r2_i] * w2 * (r2 - (double)r2_i);
			for(int i=r1_i+1; i<r2_i; i++) {
				s2 += acc[i] * w2;
			}
		}
		sum += s2*s2;
	}
	return sum;
}

static void hough_add_poly(hough_accum_t *ha, mpoly_t *mp) {
	for(int r_idx=0; r_idx<mp->num_rings; r_idx++) {
		ring_t *ring = mp->rings + r_idx;
		for(int v_idx=0; v_idx<ring->npts; v_idx++) {
			vertex_t v1 = ring->pts[v_idx];
			vertex_t v2 = ring->pts[(v_idx+1) % ring->npts];
			hough_add_seg(ha, v1, v2);
		}
	}
}

static double color_step(double v) {
	double l = 1.0/4.0;
	v /= l;
	if(v > 0) v = 1.0 - v;
	else      v = 1.0 + v;
	return v<0 ? 0 : v;
}

static void hough_accum_dump(hough_accum_t *ha, const char *fn) {
	int w_ang = ha->w_ang;
	int w_rad = ha->w_rad;
	accum_t maxval = 0;
	for(int i=0; i<w_ang*w_rad; i++) {
		if(ha->accum[i] > maxval) maxval = ha->accum[i];
	}
	uint8_t *buf = (uint8_t *)malloc_or_die(w_ang*w_rad*3);
	float lmv = logf(maxval+1);
	for(int i=0; i<w_ang*w_rad; i++) {
		float v = logf(ha->accum[i]+1) / lmv;
		double b = color_step(v-1.0/4.0) + color_step(v-1.0);
		double g = color_step(v-2.0/4.0) + color_step(v-1.0);
		double r = color_step(v-3.0/4.0) + color_step(v-1.0);
		buf[i*3+0] = r>1 ? 255 : (int)(r*255.0);
		buf[i*3+1] = g>1 ? 255 : (int)(g*255.0);
		buf[i*3+2] = b>1 ? 255 : (int)(b*255.0);
	}
	FILE *fh = fopen(fn, "w");
	fprintf(fh, "P6\n%d %d\n255\n", w_ang, w_rad);
	fwrite(buf, 1, w_ang*w_rad*3, fh);
	fclose(fh);
}

static void hough_inverse(hough_accum_t *ha, mpoly_t *mp) {
	bbox_t bbox = get_polygon_bbox(mp);
	int w = (int)bbox.max_x + 1;
	int h = (int)bbox.max_y + 1;
	accum_t *invbuf = (accum_t *)malloc_or_die(sizeof(accum_t) * w*h);

	/*
	for(int r_idx=0; r_idx<mp->num_rings; r_idx++) {
		ring_t *ring = mp->rings + r_idx;
		for(int v_idx=0; v_idx<ring->npts; v_idx++) {
			vertex_t v1 = ring->pts[v_idx];
			vertex_t v2 = ring->pts[(v_idx+1) % ring->npts];
			double val = hough_integrate_seg(ha, v1, v2);

			double dx = v2.x - v1.x;
			double dy = v2.y - v1.y;
			double len = sqrt(dx*dx + dy*dy);
			int nsteps = (int)ceil(len);
			for(int i=0; i<nsteps; i++) {
				double alpha = (double)i / (double)nsteps;
				int x = (int)(v1.x + alpha*dx);
				int y = (int)(v1.y + alpha*dy);
				if(x>=0 && x<w && y>=0 && y<h) {
					invbuf[x+w*y] = val;
				}
			}
		}
	}
	*/

	for(int y=0; y<h; y++) for(int x=0; x<w; x++) {
		vertex_t v;
		v.x = x;
		v.y = y;
		double val = hough_integrate_seg(ha, v, v);
		invbuf[x+w*y] = val;
	}

	accum_t maxval = 0;
	for(int i=0; i<w*h; i++) {
		if(invbuf[i] > maxval) maxval = invbuf[i];
	}
	uint8_t *buf = (uint8_t *)malloc_or_die(w*h*3);
	for(int i=0; i<w*h; i++) {
		float v = invbuf[i] / maxval;
		double b = v;
		double g = v;
		double r = v;
		buf[i*3+0] = r>1 ? 255 : (int)(r*255.0);
		buf[i*3+1] = g>1 ? 255 : (int)(g*255.0);
		buf[i*3+2] = b>1 ? 255 : (int)(b*255.0);
	}
	FILE *fh = fopen("inv.ppm", "w");
	fprintf(fh, "P6\n%d %d\n255\n", w, h);
	fwrite(buf, 1, w*h*3, fh);
	fclose(fh);
}

void do_hough(mpoly_t *mp) {
	printf("hough init\n");

	hough_accum_t *ha = (hough_accum_t *)malloc_or_die(sizeof(hough_accum_t));

	bbox_t bbox = get_polygon_bbox(mp);
	ha->center.x = (bbox.max_x + bbox.min_x) / 2.0;
	ha->center.y = (bbox.max_y + bbox.min_y) / 2.0;

	double max_x = (bbox.max_x - bbox.min_x) / 2.0;
	double max_y = (bbox.max_y - bbox.min_y) / 2.0;
	max_x = max_y = MAX(max_x, max_y); // FIXME

	ha->w_ang = 2000; //(int)MAX(max_x, max_y); // FIXME
	ha->w_rad = 2000; //(int)MAX(max_x, max_y);
	ha->accum = (accum_t *)malloc_or_die(sizeof(accum_t) * ha->w_ang * ha->w_rad);
	for(int i=0; i<ha->w_ang*ha->w_rad; i++) ha->accum[i] = 0;

	double scale_x = (double)(ha->w_rad/2-1) / (max_x*sqrt(2.0));
	double scale_y = (double)(ha->w_rad/2-1) / (max_y*sqrt(2.0));
	ha->sin_table = (double *)malloc_or_die(sizeof(double) * ha->w_ang);
	ha->cos_table = (double *)malloc_or_die(sizeof(double) * ha->w_ang);
	for(int i=0; i<ha->w_ang; i++) {
		// FIXME - use full wave and take into account handedness of border
		double theta = (double)i / (double)ha->w_ang * M_PI;
		ha->cos_table[i] = cos(theta) * scale_x;
		ha->sin_table[i] = sin(theta) * scale_y;
	}

	//printf("bbox=[%g,%g,%g,%g] center=[%g,%g]\n",
	//	bbox.min_x, bbox.min_y, bbox.max_x, bbox.max_y, ha->center.x, ha->center.y);

	printf("hough accum\n");
	hough_add_poly(ha, mp);

	printf("hough dump\n");
	hough_accum_dump(ha, "hough.ppm");
	printf("hough inverse\n");
	hough_inverse(ha, mp);
}
