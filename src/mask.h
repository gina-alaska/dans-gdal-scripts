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




#ifndef DANGDAL_MASK_H
#define DANGDAL_MASK_H

#include <cassert>
#include <vector>

#include "common.h"
#include "polygon.h"
#include "debugplot.h"
#include "ndv.h"

namespace dangdal {

template <typename T>
class GridArray {
public:
	GridArray(int _w, int _h) :
		w(_w), h(_h),
		grid(w*h)
	{ }

// default dtor, copy, assign are OK

public:
	typename std::vector<T>::const_reference operator()(int x, int y) const {
		// out-of-bounds used to be okay and return 'false' but not now
		// FIXME - make sure that is okay
		assert(x>=0 && y>=0 && x<w && y<h);
		return grid[size_t(y)*w + x];
	}

	typename std::vector<T>::reference operator()(int x, int y) {
		// out-of-bounds used to be okay and return 'false' but not now
		// FIXME - make sure that is okay
		assert(x>=0 && y>=0 && x<w && y<h);
		return grid[size_t(y)*w + x];
	}

	T get(
		int x, int y, const T &default_val
	) const {
		if(x>=0 && y>=0 && x<w && y<h) {
			return (*this)(x, y);
		} else {
			return default_val;
		}
	}

	// FIXME - deprecate
	const typename std::vector<T>::const_reference get(int x, int y) const {
		return (*this)(x, y);
	}

	// FIXME - deprecate
	void set(int x, int y, const T &val) {
		(*this)(x, y) = val;
	}

	void zero() {
		for(size_t i=0; i<grid.size(); i++) {
			grid[i] = 0;
		}
	}

protected:
	int w, h;
	std::vector<T> grid;
};

class BitGrid : public GridArray<bool> {
public:
	BitGrid(int _w, int _h) : GridArray<bool>(_w, _h) { }

	void invert() {
		for(size_t i=0; i<grid.size(); i++) {
			grid[i] = !grid[i];
		}
	}

	void erode();

	Vertex centroid();
};

// Returns a BitGrid with 'true' values correspond to valid (not ndv) pixels.
BitGrid get_bitgrid_for_dataset(GDALDatasetH ds, const std::vector<size_t> &bandlist,
	const NdvDef &ndv_def, DebugPlot *dbuf);

} // namespace dangdal

#endif // ifndef DANGDAL_MASK_H
