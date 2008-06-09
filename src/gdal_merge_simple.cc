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

void copyGeoCode(GDALDatasetH dst_ds, GDALDatasetH src_ds);

void usage(const char *cmdname) {
	printf("Usage:\n");
	printf("    %s -in <rgb.tif> -in <mask.tif> -out <out.tif>\n", cmdname);
	printf("\nMerges several images into one image with many bands.\n");
	printf("Currently only 8-bit data is supported.\n");
	exit(1);
}

int main(int argc, char *argv[]) {
	int i;

	const char *dst_fn = NULL;

	int src_ds_count = 0;
	int band_count = 0;
	GDALDatasetH *src_ds = NULL;

	const char *output_format = NULL;

	GDALAllRegister();

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(!strcmp(arg, "-out")) { if(argp == argc) usage(argv[0]); dst_fn = argv[argp++]; }
			else if(!strcmp(arg, "-of")) { if(argp == argc) usage(argv[0]); output_format = argv[argp++]; }
			else if(!strcmp(arg, "-in")) {
				if(argp == argc) usage(argv[0]);
				char *fn = argv[argp++];
				src_ds = (GDALDatasetH *)realloc_or_die(src_ds, sizeof(GDALDatasetH) * (src_ds_count+1));
				GDALDatasetH ds = GDALOpen(fn, GA_ReadOnly);
				if(!ds) fatal_error("open failed");
				src_ds[src_ds_count++] = ds; 
				band_count += GDALGetRasterCount(ds);
			}
			else usage(argv[0]);
		} else {
			usage(argv[0]);
		}
	}

	if(!band_count) usage(argv[0]);
	if(!dst_fn) usage(argv[0]);

	if(!output_format) output_format = "GTiff";

	//////// open source ////////

	GDALRasterBandH *src_bands = (GDALRasterBandH *)malloc_or_die(sizeof(GDALRasterBandH) * band_count);

	int w=0, h=0;

	int ds_idx;
	int band_idx;
	for(ds_idx=0, band_idx=0; ds_idx<src_ds_count; ds_idx++) {
		int ds_w = GDALGetRasterXSize(src_ds[ds_idx]);
		int ds_h = GDALGetRasterYSize(src_ds[ds_idx]);
		if(!ds_w || !ds_h) fatal_error("missing width/height");
		if(w) {
			if(ds_w != w || ds_h != h) fatal_error("size mismatch for inputs");
		} else {
			w = ds_w;
			h = ds_h;
		}

		int nb = GDALGetRasterCount(src_ds[ds_idx]);
		int i;
		for(i=0; i<nb; i++) {
			src_bands[band_idx++] = GDALGetRasterBand(src_ds[ds_idx], i+1);
		}
	}

	//////// open output ////////

	printf("Output size is %d x %d x %d\n", w, h, band_count);

	GDALDriverH dst_driver = GDALGetDriverByName(output_format);
	if(!dst_driver) fatal_error("unrecognized output format (%s)", output_format);
	GDALDatasetH dst_ds = GDALCreate(dst_driver, dst_fn, w, h, band_count, GDT_Byte, NULL);
	if(!dst_ds) fatal_error("could not create output");
	copyGeoCode(dst_ds, src_ds[0]);

	GDALRasterBandH *dst_bands = (GDALRasterBandH *)malloc_or_die(sizeof(GDALRasterBandH) * band_count);
	for(i=0; i<band_count; i++) {
		dst_bands[i] = GDALGetRasterBand(dst_ds, i+1);
	}

	//////// process data ////////

	int chunk_size = 200;

	// FIXME - handle other datatypes
	uint8_t *buf = (uint8_t *)malloc_or_die(w * chunk_size);

	int row;
	for(row=0; row<h; row+=chunk_size) {
		GDALTermProgress((double)row/(double)h, NULL, NULL);

		int num_lines = chunk_size;
		if(row+num_lines > h) num_lines = h - row;
		int band_idx;
		for(band_idx=0; band_idx<band_count; band_idx++) {
			if(GDALRasterIO(
				src_bands[band_idx], GF_Read,
				0, row, w, num_lines,
				buf, w, num_lines,
				GDT_Byte, 0, 0
			) != CE_None) fatal_error("read error");
			if(GDALRasterIO(
				dst_bands[band_idx], GF_Write,
				0, row, w, num_lines,
				buf, w, num_lines,
				GDT_Byte, 0, 0
			) != CE_None) fatal_error("write error");
		}
	}

	free(buf);

	//////// shutdown ////////

	for(ds_idx=0; ds_idx<src_ds_count; ds_idx++) {
		GDALClose(src_ds[ds_idx]);
	}
	GDALClose(dst_ds);

	GDALTermProgress(1, NULL, NULL);

	return 0;
}

void copyGeoCode(GDALDatasetH dst_ds, GDALDatasetH src_ds) {
	double affine[6];
	if(GDALGetGeoTransform(src_ds, affine) == CE_None) {
		GDALSetGeoTransform(dst_ds, affine);
	}
	GDALSetProjection(dst_ds, GDALGetProjectionRef(src_ds));
}
