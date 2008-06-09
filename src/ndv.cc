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
#include "ndv.h"

void usage(const char *cmdname); // externally defined

void print_ndv_usage() {
	// FIXME
	printf("\
No-data values:\n\
  -ndv\n\
  -valid\n\
");
}

static void add_arg_to_list(int *argc_ptr, char ***argv_ptr, char *new_arg) {
	*argv_ptr= (char **)realloc_or_die(*argv_ptr, sizeof(char *) * (*argc_ptr+1));
	(*argv_ptr)[*argc_ptr] = new_arg;
	(*argc_ptr)++;
}

// this modifies its input
static void add_minmax(ndv_range_t *r, char *minmax_string) {
	printf("minmax [%s]\n", minmax_string);

	char *endptr;
	double min, max;
	char *delim = strstr(minmax_string, "..");
	if(delim) {
		*delim = 0;
		char *s2 = delim + 2;

		// FIXME - allow -Inf, Inf
		min = strtod(minmax_string, &endptr);
		if(*endptr) usage("FIXME");
		max = strtod(s2, &endptr);
		if(*endptr) usage("FIXME");
	} else {
		min = max = strtod(minmax_string, &endptr);
		if(*endptr) usage("FIXME");
	}

	r->min = (double *)realloc_or_die(r->min, sizeof(double) * (r->nbands+1));
	r->max = (double *)realloc_or_die(r->max, sizeof(double) * (r->nbands+1));
	r->min[r->nbands] = min;
	r->max[r->nbands] = max;
	r->nbands++;
}

static void add_ndv_range(ndv_def_t *nd, const char *range_string) {
	printf("range [%s]\n", range_string);

	ndv_range_t r;
	r.min = r.max = NULL;
	r.nbands = 0;

	char *buf = strdup(range_string);
	char *p = buf;
	for(;;) {
		char *p2 = strstr(p, " ");
		if(p2) *p2 = 0;
		add_minmax(&r, p);
		if(p2) p = p2+1;
		else break;
	}
	free(buf);

	nd->ranges = (ndv_range_t *)realloc_or_die(nd->ranges, sizeof(ndv_range_t) * (nd->nranges+1));
	nd->ranges[nd->nranges++] = r;
}

ndv_def_t init_ndv_options(int *argc_ptr, char ***argv_ptr) {
	int argc = *argc_ptr;
	char **argv = *argv_ptr;

	int argc_out = 0;
	char **argv_out = NULL;
	add_arg_to_list(&argc_out, &argv_out, argv[0]);

	ndv_def_t nd;
	nd.ranges = NULL;
	nd.nranges = 0;
	int got_ndv=0, got_dv=0;

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(!strcmp(arg, "-ndv")) {
				if(argp == argc) usage(argv[0]);
				add_ndv_range(&nd, argv[argp++]);
				got_ndv = 1;
			} else if(!strcmp(arg, "-valid")) {
				if(argp == argc) usage(argv[0]);
				add_ndv_range(&nd, argv[argp++]);
				got_dv = 1;
			} else {
				add_arg_to_list(&argc_out, &argv_out, arg);
			}
		} else {
			add_arg_to_list(&argc_out, &argv_out, arg);
		}
	}

	if(got_ndv && got_dv) {
		fatal_error("you cannot use both -ndv and -valid options");
	} else {
		nd.invert = got_dv;
	}

	//for(int i=0; i<nd.nranges; i++) {
	//	ndv_range_t range = nd.ranges[i];
	//	printf("%d: %d bands\n", i, range.nbands);
	//	for(int j=0; j<range.nbands; j++) {
	//		printf("%d,%d = [%lf,%lf]\n", i, j, range.min[j], range.max[j]);
	//	}
	//}

	*argc_ptr = argc_out;
	*argv_ptr = argv_out;

	return nd;
}

void add_ndv_from_raster(ndv_def_t *nd, GDALDatasetH ds, int bandlist_size, int *bandlist) {
	ndv_range_t r;
	r.nbands = bandlist_size;
	r.min = (double *)malloc_or_die(sizeof(double) * r.nbands);
	r.max = (double *)malloc_or_die(sizeof(double) * r.nbands);
	int got_error = 0;

	int band_count = GDALGetRasterCount(ds);
	int bandlist_idx;
	for(bandlist_idx=0; bandlist_idx<bandlist_size; bandlist_idx++) {
		int band_idx = bandlist[bandlist_idx];
		if(band_idx < 1 || band_idx > band_count) fatal_error("bandid out of range");

		GDALRasterBandH band = GDALGetRasterBand(ds, band_idx);

		int success;
		double val = GDALGetRasterNoDataValue(band, &success);
		if(success) {
			r.min[bandlist_idx] = r.max[bandlist_idx] = val;
		} else {
			got_error = 1;
		}
	}

	if(got_error) {
		free(r.min);
		free(r.max);
		//fatal_error("could not determine NDV");
	} else {
		nd->ranges = (ndv_range_t *)realloc_or_die(nd->ranges, sizeof(ndv_range_t) * (nd->nranges+1));
		nd->ranges[nd->nranges++] = r;
	}
}

void array_check_ndv(
	ndv_def_t *nd, int band, double *in_dbl, uint8_t *in_byte,
	uint8_t *mask_out, int num_samples
) {
	memset(mask_out, 0, num_samples);
	for(int range_idx=0; range_idx<nd->nranges; range_idx++) {
		ndv_range_t *range = &nd->ranges[range_idx];
		double min, max;
		if(band > 0 && range->nbands == 1) {
			// if only a single range is defined, use it for all bands
			min = range->min[0];
			max = range->max[0];
		} else if(band < range->nbands) {
			min = range->min[band];
			max = range->max[band];
		} else {
			fatal_error("wrong number of bands in NDV def");
		}
		if(in_dbl) {
			for(int i=0; i<num_samples; i++) {
				uint8_t valid = (in_dbl[i] >= min) && (in_dbl[i] <= max);
				mask_out[i] |= (valid ? 1 : 0);
			}
		} else {
			uint8_t min_byte = (uint8_t)MAX(ceil (min), 0);
			uint8_t max_byte = (uint8_t)MIN(floor(max), 255);
			for(int i=0; i<num_samples; i++) {
				uint8_t valid = (in_byte[i] >= min_byte) && (in_byte[i] <= max_byte);
				mask_out[i] |= (valid ? 1 : 0);
			}
		}
	}
	if(nd->invert) {
		for(int i=0; i<num_samples; i++) {
			mask_out[i] = mask_out[i] ? 0 : 1;
		}
	}
}
