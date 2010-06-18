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



#include <boost/tokenizer.hpp>
#include "common.h"
#include "ndv.h"

void usage(const char *cmdname); // externally defined

namespace dangdal {

void NdvDef::printUsage() {
	// FIXME
	printf("\
No-data values:\n\
  -ndv val                           Set a no-data value\n\
  -ndv 'val val ...'                 Set a no-data value using all input bands\n\
  -ndv 'min..max min..max ...'       Set a range of no-data values\n\
                                     (-Inf and Inf are allowed)\n\
  -valid-range 'min..max min..max ...'  Set a range of valid data values\n\
");
}

static void add_arg_to_list(int *argc_ptr, char ***argv_ptr, char *new_arg) {
	*argv_ptr = REMYALLOC(char *, *argv_ptr, (*argc_ptr+1));
	(*argv_ptr)[*argc_ptr] = new_arg;
	(*argc_ptr)++;
}

NdvInterval::NdvInterval(const std::string &s) {
	printf("minmax [%s]\n", s.c_str());

	char *endptr;
	double min, max;
	size_t delim = s.find("..");
	if(delim == std::string::npos) {
		min = max = strtod(s.c_str(), &endptr);
		if(*endptr) fatal_error("NDV value was not a number");
	} else {
		std::string s1 = s.substr(0, delim);
		std::string s2 = s.substr(delim+2);
		// FIXME - allow -Inf, Inf
		min = strtod(s1.c_str(), &endptr);
		if(*endptr) fatal_error("NDV value was not a number");
		max = strtod(s2.c_str(), &endptr);
		if(*endptr) fatal_error("NDV value was not a number");
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
}

NdvDef::NdvDef(int *argc_ptr, char ***argv_ptr) : invert(false) {
	int argc = *argc_ptr;
	char **argv = *argv_ptr;

	int argc_out = 0;
	char **argv_out = NULL;
	add_arg_to_list(&argc_out, &argv_out, argv[0]);

	bool got_ndv=0, got_dv=0;

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		if(arg[0] == '-') {
			if(!strcmp(arg, "-ndv")) {
				if(argp == argc) usage(argv[0]);
				slabs.push_back(NdvSlab(argv[argp++]));
				got_ndv = 1;
			} else if(!strcmp(arg, "-valid-range")) {
				if(argp == argc) usage(argv[0]);
				slabs.push_back(NdvSlab(argv[argp++]));
				got_dv = 1;
			} else {
				add_arg_to_list(&argc_out, &argv_out, arg);
			}
		} else {
			add_arg_to_list(&argc_out, &argv_out, arg);
		}
	}

	if(got_ndv && got_dv) {
		fatal_error("you cannot use both -ndv and -valid-range options");
	} else {
		invert = got_dv;
	}

	debugPrint();

	*argc_ptr = argc_out;
	*argv_ptr = argv_out;
}

void NdvDef::debugPrint() const {
	for(size_t i=0; i<slabs.size(); i++) {
		const NdvSlab &slab = slabs[i];
		for(size_t j=0; j<slab.range_by_band.size(); j++) {
			const NdvInterval &range = slab.range_by_band[j];
			printf("range %zd,%zd = [%g,%g]\n", i, j, range.first, range.second);
		}
	}
}

NdvDef::NdvDef(const GDALDatasetH ds, const std::vector<size_t> &bandlist) {
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

template<class T>
void flagMatches(
	const NdvInterval range,
	const T *in_data,
	uint8_t *mask_out,
	size_t nsamps
) {
	for(size_t i=0; i<nsamps; i++) {
		T val = in_data[i];
		if(range.contains(val)) mask_out[i] = 1;
	}
}

template<>
void flagMatches<uint8_t>(
	const NdvInterval range,
	const uint8_t *in_data,
	uint8_t *mask_out,
	size_t nsamps
) {
	uint8_t min_byte = (uint8_t)std::max(ceil (range.first ), 0.0);
	uint8_t max_byte = (uint8_t)std::min(floor(range.second), 255.0);
	for(size_t i=0; i<nsamps; i++) {
		uint8_t v = in_data[i];
		uint8_t match = (v >= min_byte) && (v <= max_byte);
		if(match) mask_out[i] = 1;
	}
}

template<class T>
void flagNaN(
	const T *in_data,
	uint8_t *mask_out,
	size_t nsamps
) {
	for(size_t i=0; i<nsamps; i++) {
		if(isnan(in_data[i])) mask_out[i] = 1;
	}
}

template<>
void flagNaN<uint8_t>(
	const uint8_t *in_data __attribute__((unused)),
	uint8_t *mask_out __attribute__((unused)),
	size_t nsamps __attribute__((unused))
) { } // no-op

template<class T>
void NdvDef::arrayCheckNdv(
	size_t band, const T *in_data,
	uint8_t *mask_out, size_t nsamps
) const {
	for(size_t i=0; i<nsamps; i++) mask_out[i] = 0;
	for(size_t slab_idx=0; slab_idx<slabs.size(); slab_idx++) {
		const NdvSlab &slab = slabs[slab_idx];
		NdvInterval range;
		if(band > 0 && slab.range_by_band.size() == 1) {
			// if only a single range is defined, use it for all bands
			range = slab.range_by_band[0];
		} else if(band < slab.range_by_band.size()) {
			range = slab.range_by_band[band];
		} else {
			fatal_error("wrong number of bands in NDV def");
		}
		flagMatches(range, in_data, mask_out, nsamps);
	}
	if(invert) {
		for(size_t i=0; i<nsamps; i++) {
			mask_out[i] = mask_out[i] ? 0 : 1;
		}
	}
	flagNaN(in_data, mask_out, nsamps);
}

void NdvDef::aggregateMask(
	std::vector<uint8_t> &total_mask,
	const std::vector<uint8_t> &band_mask
) {
	assert(total_mask.size() == band_mask.size());
	size_t nsamps = total_mask.size();
	if(invert) {
		// pixel is valid only if all bands are within valid range
		for(size_t i=0; i<nsamps; i++) {
			if(band_mask[i]) total_mask[i] = 1;
		}
	} else {
		// pixel is NDV only if all bands are NDV
		for(size_t i=0; i<nsamps; i++) {
			if(!band_mask[i]) total_mask[i] = 0;
		}
	}
}

template void NdvDef::arrayCheckNdv<uint8_t>(
	size_t band, const uint8_t *in_data,
	uint8_t *mask_out, size_t nsamps
) const;

template void NdvDef::arrayCheckNdv<double>(
	size_t band, const double *in_data,
	uint8_t *mask_out, size_t nsamps
) const;

} // namespace dangdal
