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
#include "hough.h"

// See http://en.wikipedia.org/wiki/Hough_transform

typedef struct {
	int w, h;
	int *accum;
	vertex_t center;
	double *sin_table, *cos_table;
} hough_accum_t;

static inline void hough_add_point(hough_accum_t *ha, double x0, double y0) {
	int w = ha->w;
	int h = ha->h;
	x0 -= ha->center.x;
	y0 -= ha->center.y;
	for(int theta_i=0; theta_i<w; theta_i++) {
		double sin = ha->sin_table[theta_i];
		double cos = ha->cos_table[theta_i];
		double r = x0*cos + y0*sin;
		int r_i = (int)round(r) + h/2;
		if(r_i >= 0 && r_i < h) {
			int p = theta_i + r_i * w;
			ha->accum[p]++;
		} else {
			printf("oob: %d\n", r_i);
		}
	}
}

static void hough_add_ring(hough_accum_t *ha, ring_t *ring) {
	for(int r_idx=0; r_idx<ring->npts; r_idx++) {
		vertex_t v1 = ring->pts[r_idx];
		vertex_t v2 = ring->pts[(r_idx+1) % ring->npts];
		double dx = v2.x - v1.x;
		double dy = v2.y - v1.y;
		double len = sqrt(dx*dx + dy*dy);
		int nsteps = (int)round(len);
		for(int s_idx=0; s_idx<nsteps; s_idx++) {
			double alpha = ((double)s_idx + 0.5) / (double)nsteps;
			double x = v1.x + dx * alpha;
			double y = v1.y + dy * alpha;
			//printf("alpha=%g x=%g y=%g\n", alpha, x, y);
			hough_add_point(ha, x, y);
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
	int w = ha->w;
	int h = ha->h;
	int maxval = 0;
	for(int i=0; i<w*h; i++) {
		if(ha->accum[i] > maxval) maxval = ha->accum[i];
	}
	uint8_t *buf = (uint8_t *)malloc_or_die(w*h*3);
	float lmv = logf(maxval+1);
	for(int i=0; i<w*h; i++) {
		float v = logf(ha->accum[i]+1) / lmv;
		double b = color_step(v-1.0/4.0) + color_step(v-1.0);
		double g = color_step(v-2.0/4.0) + color_step(v-1.0);
		double r = color_step(v-3.0/4.0) + color_step(v-1.0);
		buf[i*3+0] = r>1 ? 255 : (int)(r*255.0);
		buf[i*3+1] = g>1 ? 255 : (int)(g*255.0);
		buf[i*3+2] = b>1 ? 255 : (int)(b*255.0);
	}
	FILE *fh = fopen(fn, "w");
	fprintf(fh, "P6\n%d %d\n255\n", w, h);
	fwrite(buf, 1, w*h*3, fh);
	fclose(fh);
}

void do_hough(mpoly_t *mp) {
	printf("hough init\n");
	hough_accum_t *ha = (hough_accum_t *)malloc_or_die(sizeof(hough_accum_t));
	ha->w = 5000;
	ha->h = 5000;
	ha->accum = (int *)malloc_or_die(sizeof(int) * ha->w * ha->h);
	memset(ha->accum, 0, sizeof(int) * ha->w * ha->h);

	bbox_t bbox = get_polygon_bbox(mp);
	ha->center.x = (bbox.max_x + bbox.min_x) / 2.0;
	ha->center.y = (bbox.max_y + bbox.min_y) / 2.0;

	double max_x = (bbox.max_x - bbox.min_x) / 2.0;
	double max_y = (bbox.max_y - bbox.min_y) / 2.0;
	double scale_x = (double)(ha->h/2-1) / (max_x*sqrt(2.0));
	double scale_y = (double)(ha->h/2-1) / (max_y*sqrt(2.0));
	ha->sin_table = (double *)malloc_or_die(sizeof(double) * ha->w);
	ha->cos_table = (double *)malloc_or_die(sizeof(double) * ha->w);
	for(int i=0; i<ha->w; i++) {
		// FIXME - use full wave and take into account handedness of border
		double theta = (double)i / (double)ha->w * M_PI;
		ha->cos_table[i] = cos(theta) * scale_x;
		ha->sin_table[i] = sin(theta) * scale_y;
	}

	//printf("bbox=[%g,%g,%g,%g] center=[%g,%g]\n",
	//	bbox.min_x, bbox.min_y, bbox.max_x, bbox.max_y, ha->center.x, ha->center.y);

	printf("hough accum\n");
	for(int i=0; i<mp->num_rings; i++) {
		hough_add_ring(ha, mp->rings+i);
	}

	printf("hough dump\n");
	hough_accum_dump(ha, "hough.ppm");
}
