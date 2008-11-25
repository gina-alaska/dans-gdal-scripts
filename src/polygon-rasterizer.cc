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

#define EPSILON 1e-9

typedef struct {
	int num_crossings;
	int array_size;
	double *crossings;
} row_crossings_dbl_t;

row_crossings_t crossings_dbl_to_int(row_crossings_dbl_t *in) {
	row_crossings_t out;
	out.num_crossings = 0;
	out.array_size = in->num_crossings;
	if(in->num_crossings == 0) {
		out.crossings = NULL;
	} else {
		out.crossings = (int *)malloc_or_die(sizeof(int) * out.array_size);
		for(int i=0; i<in->num_crossings; i+=2) {
			int from = (int)ceil(in->crossings[i] - EPSILON);
			int to = (int)floor(in->crossings[i+1] + EPSILON);
			if(to > from) {
				out.crossings[out.num_crossings++] = from;
				out.crossings[out.num_crossings++] = to;
			}
		}
	}
	return out;
}

static int dbl_compare(const void *ap, const void *bp) {
	double a = *((double *)ap);
	double b = *((double *)bp);
	return (a<b) ? -1 : (a>b) ? 1 : 0;
}

static void dbl_crossings_sort(row_crossings_dbl_t *r) {
	double *c = r->crossings;
	if(r->num_crossings == 2) {
		if(c[0] > c[1]) {
			double tmp = c[0];
			c[0] = c[1];
			c[1] = tmp;
		}
	} else {
		if(r->num_crossings % 2) {
			fatal_error("should not have an odd number of crossings");
		}
		qsort(c, r->num_crossings, sizeof(double), dbl_compare);
	}
}

// This function returns a list of pixel ranges for each row.  The ranges consist of pixels
// that are entirely contained within the polygon.  The results will be slightly wrong for
// polygons whose vertices are not integers.
row_crossings_t *get_row_crossings(mpoly_t *mpoly, int min_y, int num_rows) {
	row_crossings_dbl_t *rows_top = (row_crossings_dbl_t *)malloc_or_die(sizeof(row_crossings_dbl_t) * num_rows);
	row_crossings_dbl_t *rows_bot = (row_crossings_dbl_t *)malloc_or_die(sizeof(row_crossings_dbl_t) * num_rows);

	for(int i=0; i<num_rows; i++) {
		rows_top[i].num_crossings = 0;
		rows_top[i].array_size = 0;
		rows_top[i].crossings = NULL;
		rows_bot[i].num_crossings = 0;
		rows_bot[i].array_size = 0;
		rows_bot[i].crossings = NULL;
	}

	for(int i=0; i<mpoly->num_rings; i++) {
		ring_t *c = mpoly->rings + i;
		for(int j=0; j<c->npts; j++) {
			double x0 = c->pts[j].x;
			double y0 = c->pts[j].y;
			double x1 = c->pts[(j+1)%c->npts].x;
			double y1 = c->pts[(j+1)%c->npts].y;
			if(y0 == y1) continue;
			if(y0 > y1) {
				double tmp;
				tmp=x0; x0=x1; x1=tmp; 
				tmp=y0; y0=y1; y1=tmp; 
			}
			double alpha = (x1-x0) / (y1-y0);
			int y0i = (int)round(y0);
			int y1i = (int)round(y1);
			for(int y=y0i; y<=y1i; y++) {
				double x = x0 + ((double)y - y0)*alpha;

				int row = y - min_y - 1;
				if(y > y0i && row >= 0 && row < num_rows) {
					row_crossings_dbl_t *r = rows_bot + row;
					if(r->num_crossings == r->array_size) {
						r->array_size += 16;
						r->crossings = (double *)realloc_or_die(r->crossings,
							sizeof(double) * r->array_size);
					}
					r->crossings[r->num_crossings++] = x;
				}

				row = y - min_y;
				if(y < y1i && row >= 0 && row < num_rows) {
					row_crossings_dbl_t *r = rows_top + row;
					if(r->num_crossings == r->array_size) {
						r->array_size += 16;
						r->crossings = (double *)realloc_or_die(r->crossings,
							sizeof(double) * r->array_size);
					}
					r->crossings[r->num_crossings++] = x;
				}
			}
		}
	}

	row_crossings_t *rows_out = (row_crossings_t *)malloc_or_die(sizeof(row_crossings_t) * num_rows);

	for(int i=0; i<num_rows; i++) {
		rows_out[i].num_crossings = 0;
		rows_out[i].array_size = 0;
		rows_out[i].crossings = NULL;
	}

	for(int row=0; row<num_rows; row++) {
		dbl_crossings_sort(rows_top + row);
		dbl_crossings_sort(rows_bot + row);
		row_crossings_t top = crossings_dbl_to_int(rows_top + row);
		row_crossings_t bot = crossings_dbl_to_int(rows_bot + row);
		if(top.num_crossings && bot.num_crossings) {
			crossings_intersection(rows_out + row, &top, &bot);
			free(top.crossings);
			free(bot.crossings);
		} else if(top.num_crossings) {
			rows_out[row] = top;
		} else if(bot.num_crossings) {
			rows_out[row] = bot;
		} else {
			rows_out[row].num_crossings = 0;
			rows_out[row].array_size = 0;
			rows_out[row].crossings = NULL;
		}
		if(rows_top[row].crossings) free(rows_top[row].crossings);
		if(rows_bot[row].crossings) free(rows_bot[row].crossings);
	}
	free(rows_top);
	free(rows_bot);

	return rows_out;
}

void free_row_crossings(row_crossings_t *rc, int num_rows) {
	for(int i=0; i<num_rows; i++) {
		if(rc[i].crossings) free(rc[i].crossings);
	}
	free(rc);
}

void mask_from_mpoly(mpoly_t *mpoly, int w, int h, const char *fn) {
	int i, j, y;

	printf("mask draw: begin\n");

	row_crossings_t *rows = get_row_crossings(mpoly, 0, h);

	printf("mask draw: write\n");

	FILE *fout = fopen(fn, "wb");
	if(!fout) fatal_error("cannot open mask output");
	fprintf(fout, "P4\n%d %d\n", w, h);
	uint8_t *buf = (uint8_t *)malloc_or_die((w+7)/8);
	for(y=0; y<h; y++) {
		memset(buf, 0, (w+7)/8);
		uint8_t *p = buf;
		uint8_t bitp = 128;
		row_crossings_t *r = rows+y;
		for(i=0; i<w; i++) {
			uint8_t v = 1;
			// not the fastest way...
			for(j=0; j<r->num_crossings; j++) {
				if(i >= r->crossings[j]) v = !v;
			}
			if(v) *p |= bitp;
			bitp >>= 1;
			if(!bitp) {
				p++;
				bitp = 128;
			}
		}
		fwrite(buf, (w+7)/8, 1, fout);
	}
	fclose(fout);
	free_row_crossings(rows, h);
	printf("mask draw: done\n");
}

void crossings_intersection(row_crossings_t *out, row_crossings_t *in1, row_crossings_t *in2) {
	out->num_crossings = 0;
	int *c1 = in1->crossings;
	int  n1 = in1->num_crossings;
	int *c2 = in2->crossings;
	int  n2 = in2->num_crossings;
	int p1=0, p2=0;
	while(p1<n1 && p2<n2) {
		int open, close;
		if(c1[p1] > c2[p2]) {
			if(c1[p1] >= c2[p2+1]) {
				p2 += 2;
				continue;
			}
			open = c1[p1];
			if(c1[p1+1] < c2[p2+1]) {
				close = c1[p1+1];
				p1 += 2;
			} else {
				close = c2[p2+1];
				p2 += 2;
			}
		} else {
			if(c2[p2] >= c1[p1+1]) {
				p1 += 2;
				continue;
			}
			open = c2[p2];
			if(c2[p2+1] < c1[p1+1]) {
				close = c2[p2+1];
				p2 += 2;
			} else {
				close = c1[p1+1];
				p1 += 2;
			}
		}
		if(out->array_size < out->num_crossings+2) {
			out->array_size += 16;
			out->crossings = (int *)realloc_or_die(out->crossings,
				sizeof(int) * out->array_size);
		}
		out->crossings[out->num_crossings++] = open;
		out->crossings[out->num_crossings++] = close;
	}
}
