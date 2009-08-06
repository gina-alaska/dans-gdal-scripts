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



#include "common.h"
#include "ndv.h"
#include "mask.h"

void usage(const char *cmdname) {
	printf("Usage:\n  %s [options] [image_name] [mask_name.pbm]\n", cmdname);
	printf("\n");
	
	print_ndv_usage();

	printf("\
\n\
Misc:\n\
  -b band_id -b band_id ...   Bands to inspect (default is all bands)\n\
  -invert              Make mask cover no-data pixels instead of data pixels\n\
  -erosion             Erode pixels that don't have two consecutive neighbors\n\
  -v                   Verbose\n\
\n\
");
	exit(1);
}

int main(int argc, char **argv) {
	char *input_raster_fn = NULL;

	char *mask_out_fn = NULL;
	int do_erosion = 0;
	int do_invert = 0;
	int inspect_numbands = 0;
	int *inspect_bandids = NULL;

	if(argc == 1) usage(argv[0]);

	ndv_def_t ndv_def = init_ndv_options(&argc, &argv);

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(!strcmp(arg, "-v")) {
				VERBOSE++;
			} else if(!strcmp(arg, "-b")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				int bandid = (int)strtol(argv[argp++], &endptr, 10);
				if(*endptr) usage(argv[0]);
				inspect_bandids = REMYALLOC(int, inspect_bandids, (inspect_numbands+1));
				inspect_bandids[inspect_numbands++] = bandid;
			} else if(!strcmp(arg, "-erosion")) {
				do_erosion++;
			} else if(!strcmp(arg, "-invert")) {
				do_invert = 1;
			} else if(!strcmp(arg, "-mask-out")) {
				if(argp == argc) usage(argv[0]);
				mask_out_fn = argv[argp++];
			} else usage(argv[0]);
		} else {
			if(!input_raster_fn) {
				input_raster_fn = arg;
			} else if(!mask_out_fn) {
				mask_out_fn = arg;
			} else {
				usage(argv[0]);
			}
		}
	}

	if(!input_raster_fn || !mask_out_fn) usage(argv[0]);

	GDALAllRegister();

	GDALDatasetH ds = GDALOpen(input_raster_fn, GA_ReadOnly);
	if(!ds) fatal_error("open failed");
	int w = GDALGetRasterXSize(ds);
	int h = GDALGetRasterYSize(ds);

	if(!inspect_numbands) {
		inspect_numbands = GDALGetRasterCount(ds);
		inspect_bandids = MYALLOC(int, inspect_numbands);
		for(int i=0; i<inspect_numbands; i++) inspect_bandids[i] = i+1;
	}

	CPLPushErrorHandler(CPLQuietErrorHandler);

	if(!ndv_def.nranges) {
		add_ndv_from_raster(&ndv_def, ds, inspect_numbands, inspect_bandids);
	}

	if(!ndv_def.nranges) {
		fatal_error("cannot determine no-data-value");
	}

	uint8_t *mask = get_mask_for_dataset(ds, inspect_numbands, inspect_bandids,
		&ndv_def, NULL);

	if(do_invert) {
		invert_mask(mask, w, h);
	}

	if(do_erosion) {
		erode_mask(mask, w, h);
	}

	FILE *fout = fopen(mask_out_fn, "wb");
	if(!fout) fatal_error("cannot open mask output");
	fprintf(fout, "P4\n%d %d\n", w, h);
	uint8_t *buf = MYALLOC(uint8_t, (w+7)/8);
	for(int y=0; y<h; y++) {
		memset(buf, 0, (w+7)/8);
		uint8_t *p = buf;
		uint8_t bitp = 128;
		uint8_t *mp = mask + (w+2)*(y+1) + 1;
		for(int x=0; x<w; x++) {
			uint8_t v = *(mp++);
			if(!v) *p |= bitp;
			bitp >>= 1;
			if(!bitp) {
				p++;
				bitp = 128;
			}
		}
		fwrite(buf, (w+7)/8, 1, fout);
	}
	fclose(fout);
}
