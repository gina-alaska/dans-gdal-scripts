/*
Copyright (c) 2009, Regents of the University of Alaska

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



#include "debugplot.h"

report_image_t *create_plot(double w, double h) {
	if(w<0 || h<0) fatal_error("negative size for debug plot (%g,%g)", w, h);

	report_image_t *dbuf = MYALLOC(report_image_t, 1);
	dbuf->canvas_w = w+2;
	dbuf->canvas_h = h+2;

	dbuf->img_w = size_t(w+1);
	dbuf->img_h = size_t(h+1);
	if(dbuf->img_w > 800) {
		dbuf->img_w = 800;
		dbuf->img_h = size_t(800.0 * (h+1) / (w+1));
	}
	if(dbuf->img_h > 800) {
		dbuf->img_w = size_t(800.0 * (w+1) / (h+1));
		dbuf->img_h = 800;
	}
	if(dbuf->img_w < 1) dbuf->img_w = 1;
	if(dbuf->img_h < 1) dbuf->img_h = 1;

	dbuf->stride_x = (int)floor(w / (double)dbuf->img_w);
	dbuf->stride_y = (int)floor(h / (double)dbuf->img_h);
	if(dbuf->stride_x < 1) dbuf->stride_x = 1;
	if(dbuf->stride_y < 1) dbuf->stride_y = 1;

	dbuf->img = MYALLOC(uint8_t, dbuf->img_w*dbuf->img_h*3);
	memset(dbuf->img, 0, dbuf->img_w*dbuf->img_h*3);
	
	return dbuf;
}

void write_plot(report_image_t *dbuf, const char *fn) {
	FILE *fout = fopen(fn, "wb");
	fprintf(fout, "P6\n%zd %zd\n255\n", dbuf->img_w, dbuf->img_h);
	fwrite(dbuf->img, dbuf->img_w*dbuf->img_h, 3, fout);
	fclose(fout);
}

void plot_point_big(report_image_t *dbuf, double x, double y, uint8_t r, uint8_t g, uint8_t b) {
	int center_x = (int)(x / dbuf->canvas_w * (double)(dbuf->img_w-1) + .5);
	int center_y = (int)(y / dbuf->canvas_h * (double)(dbuf->img_h-1) + .5);
	for(int dx=-1; dx<=1; dx++) for(int dy=-1; dy<=1; dy++) {
		int plot_x = center_x + dx;
		int plot_y = center_y + dy;
		if(plot_x>=0 && plot_y>=0 && size_t(plot_x)<dbuf->img_w && size_t(plot_y)<dbuf->img_h) {
			uint8_t *p = dbuf->img + (plot_x + dbuf->img_w*plot_y)*3;
			*(p++) = r; *(p++) = g; *(p++) = b;
		}
	}
}

void plot_point(report_image_t *dbuf, double x, double y, uint8_t r, uint8_t g, uint8_t b) {
	int plot_x = (int)round(x / dbuf->canvas_w * (double)(dbuf->img_w-1));
	int plot_y = (int)round(y / dbuf->canvas_h * (double)(dbuf->img_h-1));
	if(plot_x>=0 && plot_y>=0 && size_t(plot_x)<dbuf->img_w && size_t(plot_y)<dbuf->img_h) {
		uint8_t *p = dbuf->img + (plot_x + dbuf->img_w*plot_y)*3;
		*(p++) = r; *(p++) = g; *(p++) = b;
	}
}

void plot_line(report_image_t *dbuf, vertex_t p0, vertex_t p1, 
uint8_t r, uint8_t g, uint8_t b) {
	double dx = (p1.x-p0.x) / dbuf->canvas_w * (double)dbuf->img_w;
	double dy = (p1.y-p0.y) / dbuf->canvas_h * (double)dbuf->img_h;
	double len = sqrt(dx*dx + dy*dy) + 2.0;
	for(double alpha=0; alpha<=1; alpha+=1.0/len) {
		double x = p0.x+(p1.x-p0.x)*alpha;
		double y = p0.y+(p1.y-p0.y)*alpha;
		plot_point(dbuf, x, y, r, g, b);
	}
}

void debug_plot_ring(report_image_t *dbuf, ring_t *ring, uint8_t r, uint8_t g, uint8_t b) {
	for(int i=0; i<ring->npts; i++) {
		vertex_t p0 = ring->pts[i];
		vertex_t p1 = ring->pts[(i+1)%ring->npts];
		plot_line(dbuf, p0, p1, r, g, b);
	}
}


void debug_plot_mpoly(report_image_t *dbuf, mpoly_t *mpoly) {
	if(VERBOSE) printf("plotting...\n");

	for(int i=0; i<mpoly->num_rings; i++) {
		int v = (i%62)+1;
		uint8_t r = ((v&1) ? 150 : 0) + ((v&8) ? 100 : 0);
		uint8_t g = ((v&2) ? 150 : 0) + ((v&16) ? 100 : 0);
		uint8_t b = ((v&4) ? 150 : 0) + ((v&32) ? 100 : 0);
		ring_t *ring = mpoly->rings + i;
		if(ring->is_hole) { r=255; g=0; b=0; }
		else { r=255; g=255; b=0; }
		if(VERBOSE) printf("ring %d: %d pts color=%02x%02x%02x\n", i, ring->npts, r, g, b);
		debug_plot_ring(dbuf, ring, r, g, b);
		for(int j=0; j<ring->npts; j++) {
			vertex_t p = ring->pts[j];
			plot_point(dbuf, p.x, p.y, 255, 255, 255);
		}
	}
}

