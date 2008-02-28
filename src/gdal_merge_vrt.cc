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

void usage(char *cmdname) {
	fprintf(stderr, "Usage:\n");
	fprintf(stderr, "    %s -in <rgb.tif> -in <mask.tif> -out <out.vrt>\n", cmdname);
	fprintf(stderr, "\nMerges several images into one image with many bands.\n");
	exit(1);
}

int main(int argc, char *argv[]) {
	char *dst_fn = NULL;

	int src_ds_count = 0;
	char **src_fn = NULL;

	GDALAllRegister();

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(!strcmp(arg, "-out")) { if(argp == argc) usage(argv[0]); dst_fn = argv[argp++]; }
			else if(!strcmp(arg, "-in")) {
				if(argp == argc) usage(argv[0]);
				char *fn = argv[argp++];
				src_fn = (char**)realloc_or_die(src_fn, sizeof(char *) * (src_ds_count+1));
				src_fn[src_ds_count] = fn; 
				src_ds_count++;
			}
			else usage(argv[0]);
		} else {
			usage(argv[0]);
		}
	}

	if(src_ds_count < 1) usage(argv[0]);
	if(!dst_fn) usage(argv[0]);

	GDALDatasetH *src_ds = (GDALDatasetH *)malloc_or_die(sizeof(GDALDatasetH) * src_ds_count);

	int w=0, h=0;
	int ds_idx;
	for(ds_idx=0; ds_idx<src_ds_count; ds_idx++) {
		src_ds[ds_idx] = GDALOpen(src_fn[ds_idx], GA_ReadOnly);
		if(!src_ds[ds_idx]) fatal_error("open failed");

		int ds_w = GDALGetRasterXSize(src_ds[ds_idx]);
		int ds_h = GDALGetRasterYSize(src_ds[ds_idx]);
		if(!ds_w || !ds_h) fatal_error("missing width/height");

		if(ds_idx) {
			if(ds_w != w || ds_h != h) fatal_error("size mismatch for inputs");
		} else {
			w = ds_w;
			h = ds_h;
		}
	}

	GDALDriverH dst_driver = GDALGetDriverByName("VRT");
	if(!dst_driver) fatal_error("unrecognized output format");
	GDALDatasetH dst_ds = GDALCreateCopy(dst_driver, dst_fn, src_ds[0], 0, NULL, NULL, NULL);
	if(!dst_ds) fatal_error("could not create output");

	int band_idx = GDALGetRasterCount(src_ds[0]);
	for(ds_idx=1; ds_idx<src_ds_count; ds_idx++) {
		int nb = GDALGetRasterCount(src_ds[ds_idx]);
		int i;

		GDALDatasetH src_vrt_ds = GDALCreateCopy(dst_driver, "", src_ds[ds_idx], 0, NULL, NULL, NULL);
		if(!src_vrt_ds) fatal_error("could not create VRT copy");

		for(i=0; i<nb; i++) {
			GDALRasterBandH src_band = GDALGetRasterBand(src_vrt_ds, i+1);
			if(!src_band) fatal_error("could not get src_band");

			GDALAddBand(dst_ds, GDALGetRasterDataType(src_band), NULL);

			GDALRasterBandH dst_band = GDALGetRasterBand(dst_ds, band_idx+1);
			if(!dst_band) fatal_error("could not get dst_band");

			char **metadata = GDALGetMetadata(src_band, "vrt_sources");
			GDALSetMetadata(dst_band, metadata, "new_vrt_sources");

			band_idx++;
		}

		GDALClose(src_vrt_ds);
	}

	for(ds_idx=0; ds_idx<src_ds_count; ds_idx++) {
		GDALClose(src_ds[ds_idx]);
	}
	GDALClose(dst_ds);

	return 0;
}
