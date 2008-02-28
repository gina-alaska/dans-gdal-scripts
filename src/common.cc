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



#include <common.h>

int VERBOSE = 0;

void fatal_error(const char *fmt, ...) {
	va_list argp;
	
	fprintf(stderr, "\n\nerror:\n");
	
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);

	fprintf(stderr, "\n\n");

	exit(1);
}

void *malloc_or_die(size_t size) {
	if(size <= 0) fatal_error("size <= 0 in malloc_or_die (%ld)", (long)size);
	void *p = malloc(size);
	if(!p) fatal_error("out of memory");
	return p;
}

void *realloc_or_die(void *p, size_t size) {
	if(size <= 0) fatal_error("size <= 0 in realloc_or_die (%ld)", (long)size);
	p = realloc(p, size);
	if(!p) fatal_error("out of memory");
	return p;
}

int parse_list_of_doubles(char *input, int *num_out, double **list_out) {
	input = strdup(input);
	if(!input) fatal_error("out of memory");

	int num = 0;
	double *list = NULL;

	char *s1 = input;
	char *s2 = " \t\n\r";
	char *tok;
	while((tok = strtok(s1, s2))) {
		s1 = NULL;
		char *endptr;
		double val = strtod(tok, &endptr);
		if(*endptr) return 1;
		list = (double *)realloc_or_die(list, sizeof(double) * (num+1));
		list[num++] = val;
	}

	*num_out = num;
	*list_out = list;

	return 0;
}

void setup_ndv_list(GDALDatasetH ds, int bandlist_size, int *bandlist, int *num_ndv, double **ndv_list) {
	if(*num_ndv == 0) {
		*num_ndv = bandlist_size;
		*ndv_list = (double *)malloc_or_die(sizeof(double) * bandlist_size);

		int band_count = GDALGetRasterCount(ds);
		int bandlist_idx;
		for(bandlist_idx=0; bandlist_idx<bandlist_size; bandlist_idx++) {
			int band_idx = bandlist[bandlist_idx];
			if(band_idx < 1 || band_idx > band_count) fatal_error("bandid out of range");

			GDALRasterBandH band = GDALGetRasterBand(ds, band_idx);

			int success;
			(*ndv_list)[bandlist_idx] = GDALGetRasterNoDataValue(band, &success);
			if(!success) fatal_error("could not determine nodataval");
		}
	} else if(*num_ndv == 1) {
		double ndv = (*ndv_list)[0];
		*num_ndv = bandlist_size;
		*ndv_list = (double *)malloc_or_die(sizeof(double) * bandlist_size);
		int i;
		for(i=0; i<bandlist_size; i++) (*ndv_list)[i] = ndv;
	} else if(*num_ndv != bandlist_size) {
		fatal_error("number of vals passed to -nodataval must be one or equal to number of bands used");
	}
}
