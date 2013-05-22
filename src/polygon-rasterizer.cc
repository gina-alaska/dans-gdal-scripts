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



#include <algorithm>
#include <string>
#include <vector>

#include "common.h"
#include "polygon.h"
#include "polygon-rasterizer.h"

static const double EPSILON = 1e-9;

namespace dangdal {

typedef std::vector<double> row_crossings_dbl_t;

row_crossings_t crossings_dbl_to_int(const row_crossings_dbl_t &in) {
	row_crossings_t out;
	for(size_t i=0; i<in.size(); i+=2) {
		int from = (int)ceil(in[i] - EPSILON);
		int to = (int)floor(in[i+1] + EPSILON);
		if(to > from) {
			out.push_back(from);
			out.push_back(to);
		}
	}
	return out;
}

// This function returns a list of pixel ranges for each row.  The ranges
// consist of pixels that are entirely contained within the polygon.  The
// results will be slightly wrong for polygons whose vertices are not integers.
std::vector<row_crossings_t> get_row_crossings(
	const Mpoly &mpoly, int min_y, int num_rows
) {
	std::vector<row_crossings_dbl_t> rows_top(num_rows);
	std::vector<row_crossings_dbl_t> rows_bot(num_rows);

	for(size_t i=0; i<mpoly.rings.size(); i++) {
		const Ring &c = mpoly.rings[i];
		size_t npts = c.pts.size();
		for(size_t j=0; j<npts; j++) {
			size_t j_plus1 = (j==npts-1) ? 0 : (j+1);
			double x0 = c.pts[j].x;
			double y0 = c.pts[j].y;
			double x1 = c.pts[j_plus1].x;
			double y1 = c.pts[j_plus1].y;
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
					row_crossings_dbl_t &r = rows_bot[row];
					r.push_back(x);
				}

				row = y - min_y;
				if(y < y1i && row >= 0 && row < num_rows) {
					row_crossings_dbl_t &r = rows_top[row];
					r.push_back(x);
				}
			}
		}
	}

	std::vector<row_crossings_t> rows_out(num_rows);

	for(int row=0; row<num_rows; row++) {
		std::sort(rows_top[row].begin(), rows_top[row].end());
		std::sort(rows_bot[row].begin(), rows_bot[row].end());
		row_crossings_t top = crossings_dbl_to_int(rows_top[row]);
		row_crossings_t bot = crossings_dbl_to_int(rows_bot[row]);
		if(top.size() && bot.size()) {
			row_crossings_t c = crossings_intersection(top, bot);
			std::swap(rows_out[row], c);
		} else if(!top.empty()) {
			std::swap(rows_out[row], top);
		} else if(!bot.empty()) {
			std::swap(rows_out[row], bot);
		} else {
			// no-op: leave rows_out[row] empty
		}
	}

	return rows_out;
}

void mask_from_mpoly(const Mpoly &mpoly, size_t w, size_t h, const std::string &fn) {
	printf("mask draw: begin\n");

	std::vector<row_crossings_t> rows = get_row_crossings(mpoly, 0, h);

	printf("mask draw: write\n");

	FILE *fout = fopen(fn.c_str(), "wb");
	if(!fout) fatal_error("cannot open mask output");
	fprintf(fout, "P4\n%zd %zd\n", w, h);
	std::vector<uint8_t> buf((w+7)/8);
	for(size_t y=0; y<h; y++) {
		buf.assign((w+7)/8, 0);
		uint8_t *p = &buf[0];
		uint8_t bitp = 128;
		const row_crossings_t &r = rows[y];
		for(size_t i=0; i<w; i++) {
			uint8_t v = 1;
			// not the fastest way...
			for(size_t j=0; j<r.size(); j++) {
				if(double(i) >= r[j]) v = !v;
			}
			if(v) *p |= bitp;
			bitp >>= 1;
			if(!bitp) {
				p++;
				bitp = 128;
			}
		}
		fwrite(&buf[0], (w+7)/8, 1, fout);
	}
	fclose(fout);
	printf("mask draw: done\n");
}

row_crossings_t crossings_intersection(
	const row_crossings_t &in1, const row_crossings_t &in2
) {
	row_crossings_t out;
	size_t n1 = in1.size();
	size_t n2 = in2.size();
	size_t p1=0, p2=0;
	while(p1<n1 && p2<n2) {
		int open, close;
		if(in1[p1] > in2[p2]) {
			if(in1[p1] >= in2[p2+1]) {
				p2 += 2;
				continue;
			}
			open = in1[p1];
			if(in1[p1+1] < in2[p2+1]) {
				close = in1[p1+1];
				p1 += 2;
			} else {
				close = in2[p2+1];
				p2 += 2;
			}
		} else {
			if(in2[p2] >= in1[p1+1]) {
				p1 += 2;
				continue;
			}
			open = in2[p2];
			if(in2[p2+1] < in1[p1+1]) {
				close = in2[p2+1];
				p2 += 2;
			} else {
				close = in1[p1+1];
				p1 += 2;
			}
		}
		out.push_back(open);
		out.push_back(close);
	}
	return out;
}

} // namespace dangdal
