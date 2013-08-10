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
#include <cstring>
#include <cassert>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include "common.h"
#include "ndv.h"
#include "datatype_conversion.h"

void usage(const std::string &cmdname); // externally defined

namespace dangdal {

void NdvDef::printUsage() {
	printf(
"No-data values:\n"
"  -ndv val                           Set a no-data value\n"
"  -ndv 'val val ...'                 Set a no-data value using all input bands\n"
"  -ndv 'min..max min..max ...'       Set a range of no-data values\n"
"                                     (-Inf and Inf are allowed; '*' == '-Inf..Inf')\n"
"  -valid-range 'min..max min..max ...'  Set a range of valid data values\n"
);
}

NdvInterval::NdvInterval(const std::string &s_in) {
	// copy the string in case we want to change it
	std::string s = s_in;

	if(VERBOSE >= 2) printf("minmax [%s]\n", s.c_str());

	if(s == "*") s = "-Inf..Inf";

	double min, max;
	size_t delim = s.find("..");
	try {
		if(delim == std::string::npos) {
			min = max = boost::lexical_cast<double>(s);
		} else {
			std::string s1 = s.substr(0, delim);
			std::string s2 = s.substr(delim+2);
			// This is capable of interpreting -Inf and Inf
			min = boost::lexical_cast<double>(s1);
			max = boost::lexical_cast<double>(s2);
		}
	} catch(boost::bad_lexical_cast &e) {
		fatal_error("NDV value was not a number");
	}

	first = min;
	second = max;
}

NdvSlab::NdvSlab(const std::string &s) {
	//printf("range [%s]\n", s.c_str());

	boost::char_separator<char> sep(" ");
	typedef boost::tokenizer<boost::char_separator<char> > toker;
	toker tok(s, sep);
	for(toker::iterator p=tok.begin(); p!=tok.end(); ++p) {
		range_by_band.push_back(NdvInterval(*p));
	}

	if(range_by_band.empty()) {
		fatal_error("could not parse given NDV term [%s]", s.c_str());
	}
}

NdvDef::NdvDef(std::vector<std::string> &arg_list) :
	invert(false)
{
	std::vector<std::string> args_out;
	const std::string cmdname = arg_list[0];
	args_out.push_back(cmdname);

	bool got_ndv=0, got_dv=0;

	size_t argp = 1;
	while(argp < arg_list.size()) {
		const std::string &arg = arg_list[argp++];
		if(arg[0] == '-') {
			if(arg == "-ndv") {
				if(argp == arg_list.size()) usage(cmdname);
				slabs.push_back(NdvSlab(arg_list[argp++]));
				got_ndv = 1;
			} else if(arg == "-valid-range") {
				if(argp == arg_list.size()) usage(cmdname);
				slabs.push_back(NdvSlab(arg_list[argp++]));
				got_dv = 1;
			} else {
				args_out.push_back(arg);
			}
		} else {
			args_out.push_back(arg);
		}
	}

	if(got_ndv && got_dv) {
		fatal_error("you cannot use both -ndv and -valid-range options");
	} else {
		invert = got_dv;
	}

	if(VERBOSE >= 2) debugPrint();

	arg_list = args_out;
}

NdvDef::NdvDef(const GDALDatasetH ds, const std::vector<size_t> &bandlist) :
	invert(false)
{
	bool got_error = 0;

	NdvSlab slab;

	size_t band_count = GDALGetRasterCount(ds);
	for(size_t bandlist_idx=0; bandlist_idx<bandlist.size(); bandlist_idx++) {
		size_t band_idx = bandlist[bandlist_idx];
		if(band_idx < 1 || band_idx > band_count) fatal_error("bandid out of range");

		GDALRasterBandH band = GDALGetRasterBand(ds, band_idx);

		int success;
		double val = GDALGetRasterNoDataValue(band, &success);
		if(success) {
			slab.range_by_band.push_back(NdvInterval(val, val));
		} else {
			got_error = true;
		}
	}

	if(!got_error) {
		slabs.push_back(slab);
	}
}

void NdvDef::debugPrint() const {
	printf("=== NDV\n");
	for(size_t i=0; i<slabs.size(); i++) {
		const NdvSlab &slab = slabs[i];
		for(size_t j=0; j<slab.range_by_band.size(); j++) {
			const NdvInterval &range = slab.range_by_band[j];
			printf("range %zd,%zd = [%g,%g]\n", i, j, range.first, range.second);
		}
	}
	printf("=== end NDV\n");
}

template <typename T>
static inline bool contains_templated(const NdvInterval *interval, const void *p) {
	T val = *(reinterpret_cast<const T *>(p));
	return interval->contains(val);
}

inline bool NdvInterval::contains(const void *p, GDALDataType dt) const {
	return DANGDAL_RUNTIME_TEMPLATE(dt, contains_templated, this, p);
}

void NdvDef::getNdvMask(
	const std::vector<const void *> &bands,
	const std::vector<GDALDataType> &dt_list,
	uint8_t *mask_out, size_t num_pixels
) const {
	assert(bands.size() == dt_list.size());
	// use uint8_t to allow incrementing the pointer by bytes
	std::vector<const uint8_t *> in_p;
	std::vector<size_t> dt_sizes;
	for(size_t i=0; i<bands.size(); i++) {
		in_p.push_back(reinterpret_cast<const uint8_t *>(bands[i]));
		dt_sizes.push_back(GDALGetDataTypeSize(dt_list[i]) / 8);
	}

	for(const NdvSlab &slab : slabs) {
		size_t num_intervals = slab.range_by_band.size();
		assert(num_intervals == bands.size() || num_intervals == 1);
	}

	memset(mask_out, 0, num_pixels);
	for(size_t pix_idx=0; pix_idx<num_pixels; pix_idx++) {
		mask_out[pix_idx] = 0;

		// A NaN value on any band makes this pixel NDV.
		for(size_t i=0; i<bands.size(); i++) {
			if(gdal_scalar_pointer_isnan(in_p[i], dt_list[i])) {
				mask_out[pix_idx] = 1;
			}
		}

		for(const NdvSlab &slab : slabs) {
			bool all_match = true;
			for(size_t i=0; i<bands.size(); i++) {
				size_t num_intervals = slab.range_by_band.size();
				// if only one interval is given, use it for all bands
				size_t j = num_intervals==1 ? 0 : i;
				all_match &= slab.range_by_band[j].contains(in_p[i], dt_list[i]);
			}
			mask_out[pix_idx] |= all_match;
		}

		if(invert) {
			mask_out[pix_idx] = mask_out[pix_idx] ? 0 : 1;
		}

		for(size_t band_idx=0; band_idx<bands.size(); band_idx++) {
			in_p[band_idx] += dt_sizes[band_idx];
		}
	}
}

void NdvDef::getNdvMask(
	const void *band, GDALDataType dt,
	uint8_t *mask_out, size_t num_pixels
) const {
	std::vector<const void *> bands;
	bands.push_back(band);
	std::vector<GDALDataType> dt_list;
	dt_list.push_back(dt);
	getNdvMask(bands, dt_list, mask_out, num_pixels);
}

} // namespace dangdal
