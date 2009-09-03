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




#ifndef MASK_H
#define MASK_H

#include <boost/noncopyable.hpp>
#include <boost/smart_ptr.hpp>
#include "common.h"
#include "polygon.h"
#include "debugplot.h"
#include "ndv.h"

class BitGrid {
public:
	BitGrid(int _w, int _h) :
		w(_w), h(_h),
		arrlen((size_t(w)*h+7)/8),
		grid(new uint8_t[arrlen])
	{ }

// default dtor, copy, assign are OK

public:
	inline bool operator()(int x, int y) const { return get(x, y); }

	inline bool get(int x, int y) const {
		// out-of-bounds is OK and returns false
		if(x>=0 && y>=0 && x<w && y<h) {
			size_t p = size_t(y)*w + x;
			return grid[p/8] & (1 << (p&7));
		} else {
			return false;
		}
	}

	inline void set(int x, int y, bool val) {
		assert(x>=0 && y>=0 && x<w && y<h);

		size_t p = size_t(y)*w + x;
		if(val) {
			grid[p/8] |= (1 << (p&7));
		} else {
			grid[p/8] &= ~(1 << (p&7));
		}
	}

	void zero() {
		for(size_t i=0; i<arrlen; i++) {
			grid[i] = 0;
		}
	}

	void invert() {
		for(size_t i=0; i<arrlen; i++) {
			grid[i] = ~grid[i];
		}
	}

	void erode();

private:
	int w, h;
	size_t arrlen;
	boost::shared_array<uint8_t> grid;
};

uint8_t *get_mask_for_dataset(GDALDatasetH ds, int bandlist_size, int *bandlist, 
	ndv_def_t *ndv_def, report_image_t *dbuf);
uint8_t *read_dataset_8bit(GDALDatasetH ds, int band_idx, uint8_t *usage_array, report_image_t *dbuf);
uint8_t *get_mask_for_8bit_raster(size_t w, size_t h, const uint8_t *raster, uint8_t wanted);
void erode_mask(uint8_t *in_mask, size_t w, size_t h);
void invert_mask(uint8_t *in_mask, size_t w, size_t h);
vertex_t calc_centroid_from_mask(const uint8_t *mask, size_t w, size_t h);

BitGrid get_bitgrid_for_dataset(GDALDatasetH ds, int bandlist_size, int *bandlist, 
	ndv_def_t *ndv_def, report_image_t *dbuf);
BitGrid get_bitgrid_for_8bit_raster(size_t w, size_t h, const uint8_t *raster, uint8_t wanted);

#endif // ifndef MASK_H
