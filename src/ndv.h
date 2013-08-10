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




#ifndef DANGDAL_NDV_H
#define DANGDAL_NDV_H

#include <gdal.h>

#include <string>
#include <utility>
#include <vector>
#include <complex>

#include "datatype_conversion.h"

namespace dangdal {

struct NdvInterval : std::pair<double, double> {
	NdvInterval() { }
	NdvInterval(double _min, double _max) :
		std::pair<double, double>(_min, _max) { }
	explicit NdvInterval(const std::string &s);

	template <typename T>
	bool contains(T v) const {
		return v >= first && v <= second;
	}

	template <typename T>
	bool contains(std::complex<T> v) const {
		return v.real() >= first && v.real() <= second;
	}

	bool contains(const void *, GDALDataType) const;
};

struct NdvSlab {
	NdvSlab() { }
	explicit NdvSlab(const std::string &s);

	// an interval for each band, or a single interval for all bands
	std::vector<NdvInterval> range_by_band;
};

class NdvDef {
public:
	static void printUsage();
	explicit NdvDef(std::vector<std::string> &arg_list);
	NdvDef(const GDALDatasetH ds, const std::vector<size_t> &bandlist);
	void debugPrint() const;
	bool empty() const { return slabs.empty(); }
	bool isInvert() const { return invert; }

	void getNdvMask(
		const void *band, GDALDataType dt,
		uint8_t *mask_out, size_t num_pixels
	) const;

	void getNdvMask(
		const std::vector<const void *> &bands,
		const std::vector<GDALDataType> &dt_list,
		uint8_t *mask_out, size_t num_pixels
	) const;

	template <typename T>
	void getNdvMask(
		const std::vector<std::vector<T> > &bands,
		uint8_t *mask_out, size_t num_pixels
	) const {
		GDALDataType dt = GetGDALDataTypeFor<T>::t;

		std::vector<const void *> band_p;
		std::vector<GDALDataType> dt_list;
		for(const std::vector<T> &v : bands) {
			band_p.push_back(reinterpret_cast<const void *>(&v[0]));
			dt_list.push_back(dt);
		}
		getNdvMask(band_p, dt_list, mask_out, num_pixels);
	}

	// For this one, it doesn't matter what T is, it is the GDALDataType that determines how
	// the data will be read from memory (i.e. we cast &bands[i][0] to (void *)).
	template <typename T>
	void getNdvMask(
		const std::vector<std::vector<T> > &bands,
		const std::vector<GDALDataType> &dt_list,
		uint8_t *mask_out, size_t num_pixels
	) const {
		std::vector<const void *> band_p;
		for(const std::vector<T> &v : bands) {
			band_p.push_back(reinterpret_cast<const void *>(&v[0]));
		}
		getNdvMask(band_p, dt_list, mask_out, num_pixels);
	}

	bool invert;
	std::vector<NdvSlab> slabs;
};

} // namespace dangdal

#endif // DANGDAL_NDV_H
