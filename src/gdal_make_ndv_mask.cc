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

using namespace dangdal;

void usage(const std::string &cmdname) {
	printf("Usage:\n  %s [options] [image_name] [mask_name.pbm]\n", cmdname.c_str());
	printf("\n");
	
	NdvDef::printUsage();

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
	const std::string cmdname = argv[0];
	if(argc == 1) usage(cmdname);
	std::vector<std::string> arg_list = argv_to_list(argc, argv);

	char *input_raster_fn = NULL;

	char *mask_out_fn = NULL;
	bool do_erosion = 0;
	bool do_invert = 0;
	std::vector<size_t> inspect_bandids;

	if(argc == 1) usage(argv[0]);

	NdvDef ndv_def = NdvDef(&argc, &argv);

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
				inspect_bandids.push_back(bandid);
			} else if(!strcmp(arg, "-erosion")) {
				do_erosion = 1;
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
	size_t w = GDALGetRasterXSize(ds);
	size_t h = GDALGetRasterYSize(ds);

	if(inspect_bandids.empty()) {
		size_t nbands = GDALGetRasterCount(ds);
		for(size_t i=0; i<nbands; i++) inspect_bandids.push_back(i+1);
	}

	CPLPushErrorHandler(CPLQuietErrorHandler);

	if(ndv_def.empty()) {
		ndv_def = NdvDef(ds, inspect_bandids);
	}

	if(ndv_def.empty()) {
		fatal_error("cannot determine no-data-value");
	}

	BitGrid mask = get_bitgrid_for_dataset(ds, inspect_bandids, ndv_def, NULL);

	if(do_invert) {
		mask.invert();
	}

	if(do_erosion) {
		mask.erode();
	}

	FILE *fout = fopen(mask_out_fn, "wb");
	if(!fout) fatal_error("cannot open mask output");
	fprintf(fout, "P4\n%zd %zd\n", w, h);
	std::vector<uint8_t> buf((w+7)/8);
	for(size_t y=0; y<h; y++) {
		buf.assign((w+7)/8, 0);
		uint8_t *p = &buf[0];
		uint8_t bitp = 128;
		for(size_t x=0; x<w; x++) {
			if(!mask(x, y)) *p |= bitp;
			bitp >>= 1;
			if(!bitp) {
				p++;
				bitp = 128;
			}
		}
		fwrite(&buf[0], (w+7)/8, 1, fout);
	}
	fclose(fout);
}
