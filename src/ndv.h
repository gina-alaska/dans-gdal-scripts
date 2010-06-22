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




#ifndef DANGDAL_NDV_H
#define DANGDAL_NDV_H

#include <string>
#include <vector>

namespace dangdal {

struct NdvInterval : std::pair<double, double> {
	NdvInterval() { }
	NdvInterval(double _min, double _max) : 
		std::pair<double, double>(_min, _max) { }
	NdvInterval(const std::string &s);

	bool contains(double v) const {
		return v >= first && v <= second;
	}
};

struct NdvSlab {
	NdvSlab() { }
	NdvSlab(const std::string &s);

	// an interval for each band, or a single interval for all bands
	std::vector<NdvInterval> range_by_band;
};

class NdvDef {
public:
	static void printUsage();
	NdvDef(std::vector<std::string> &arg_list);
	NdvDef(const GDALDatasetH ds, const std::vector<size_t> &bandlist);
	void debugPrint() const;
	bool empty() const { return slabs.empty(); }
	bool isInvert() const { return invert; }

	template<class T>
	void arrayCheckNdv(
		size_t band, const T *in_data,
		uint8_t *mask_out, size_t nsamps
	) const;

	void aggregateMask(
		uint8_t *total_mask,
		const uint8_t *band_mask,
		size_t nsamps
	) const;

	bool invert;
	std::vector<NdvSlab> slabs;
};

} // namespace dangdal

#endif // DANGDAL_NDV_H
