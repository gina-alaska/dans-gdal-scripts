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



#include <vector>

#include "common.h"

using namespace dangdal;

void copyGeoCode(GDALDatasetH dst_ds, GDALDatasetH src_ds);

void usage(const std::string &cmdname) {
	printf("Usage:\n");
	printf("    %s -in <rgb.tif> -in <mask.tif> -out <out.tif>\n", cmdname.c_str());
	printf("\nMerges several images into one image with many bands.\n");
	printf("Currently only 8-bit data is supported.\n");
	exit(1);
}

int main(int argc, char *argv[]) {
	const std::string cmdname = argv[0];
	if(argc == 1) usage(cmdname);
	std::vector<std::string> arg_list = argv_to_list(argc, argv);

	std::string dst_fn;
	std::vector<GDALDatasetH> src_ds;
	std::string output_format;

	GDALAllRegister();

	size_t argp = 1;
	while(argp < arg_list.size()) {
		const std::string &arg = arg_list[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(arg == "-out") {
				if(argp == arg_list.size()) usage(cmdname); 
				dst_fn = arg_list[argp++];
			} else if(arg == "-of") { 
				if(argp == arg_list.size()) usage(cmdname);
				output_format = arg_list[argp++];
			} else if(arg == "-in") {
				if(argp == arg_list.size()) usage(cmdname);
				std::string fn = arg_list[argp++];
				GDALDatasetH ds = GDALOpen(fn.c_str(), GA_ReadOnly);
				if(!ds) fatal_error("open failed");
				src_ds.push_back(ds); 
			} else {
				usage(cmdname);
			}
		} else {
			usage(cmdname);
		}
	}

	if(dst_fn.empty()) usage(cmdname);

	if(output_format.empty()) output_format = "GTiff";

	//////// open source ////////

	std::vector<GDALRasterBandH> src_bands;

	size_t w=0, h=0;

	for(size_t ds_idx=0; ds_idx<src_ds.size(); ds_idx++) {
		size_t ds_w = GDALGetRasterXSize(src_ds[ds_idx]);
		size_t ds_h = GDALGetRasterYSize(src_ds[ds_idx]);
		if(!ds_w || !ds_h) fatal_error("missing width/height");
		if(w) {
			if(ds_w != w || ds_h != h) fatal_error("size mismatch for inputs");
		} else {
			w = ds_w;
			h = ds_h;
		}

		int nb = GDALGetRasterCount(src_ds[ds_idx]);
		for(int i=0; i<nb; i++) {
			src_bands.push_back(GDALGetRasterBand(src_ds[ds_idx], i+1));
		}
	}

	size_t band_count = src_bands.size();
	if(!band_count) usage(cmdname);

	//////// open output ////////

	printf("Output size is %zd x %zd x %zd\n", w, h, band_count);

	GDALDriverH dst_driver = GDALGetDriverByName(output_format.c_str());
	if(!dst_driver) fatal_error("unrecognized output format (%s)", output_format.c_str());
	GDALDatasetH dst_ds = GDALCreate(dst_driver, dst_fn.c_str(), w, h, band_count, GDT_Byte, NULL);
	if(!dst_ds) fatal_error("could not create output");
	copyGeoCode(dst_ds, src_ds[0]);

	std::vector<GDALRasterBandH> dst_bands;
	for(size_t i=0; i<band_count; i++) {
		dst_bands.push_back(GDALGetRasterBand(dst_ds, i+1));
	}

	//////// process data ////////

	int chunk_size = 200;

	// FIXME - handle other datatypes
	std::vector<uint8_t> buf(w * chunk_size);

	for(size_t row=0; row<h; row+=chunk_size) {
		GDALTermProgress((double)row/(double)h, NULL, NULL);

		int num_lines = chunk_size;
		if(row+num_lines > h) num_lines = h - row;
		for(size_t band_idx=0; band_idx<band_count; band_idx++) {
			if(GDALRasterIO(
				src_bands[band_idx], GF_Read,
				0, row, w, num_lines,
				&buf[0], w, num_lines,
				GDT_Byte, 0, 0
			) != CE_None) fatal_error("read error");
			if(GDALRasterIO(
				dst_bands[band_idx], GF_Write,
				0, row, w, num_lines,
				&buf[0], w, num_lines,
				GDT_Byte, 0, 0
			) != CE_None) fatal_error("write error");
		}
	}

	//////// shutdown ////////

	for(size_t ds_idx=0; ds_idx<src_ds.size(); ds_idx++) {
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
