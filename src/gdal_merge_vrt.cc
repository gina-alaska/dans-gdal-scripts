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



#include <vector>

#include "common.h"

void usage(const std::string &cmdname) {
	printf("Usage:\n");
	printf("    %s -in <rgb.tif> -in <mask.tif> -out <out.vrt>\n", cmdname.c_str());
	printf("\nMerges several images into one image with many bands.\n");
	printf("This program is obsoleted by \"gdalbuildvrt -separate\" from GDAL 1.7.\n");
	exit(1);
}

int main(int argc, char *argv[]) {
	const std::string cmdname = argv[0];
	if(argc == 1) usage(cmdname);
	std::vector<std::string> arg_list = argv_to_list(argc, argv);

	std::string dst_fn;
	std::vector<std::string> src_fn;

	GDALAllRegister();

	size_t argp = 1;
	while(argp < arg_list.size()) {
		const std::string &arg = arg_list[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(arg == "-out") {
				if(argp == arg_list.size()) usage(cmdname);
				dst_fn = arg_list[argp++];
			} else if(arg == "-in") {
				if(argp == arg_list.size()) usage(cmdname);
				src_fn.push_back(arg_list[argp++]);
			}
			else usage(cmdname);
		} else {
			usage(cmdname);
		}
	}

	if(src_fn.empty()) usage(cmdname);
	if(dst_fn.empty()) usage(cmdname);

	std::vector<GDALDatasetH> src_ds;

	size_t w=0, h=0;
	for(size_t ds_idx=0; ds_idx<src_fn.size(); ds_idx++) {
		GDALDatasetH ds = GDALOpen(src_fn[ds_idx].c_str(), GA_ReadOnly);
		if(!ds) fatal_error("open failed");
		src_ds.push_back(ds);

		size_t ds_w = GDALGetRasterXSize(ds);
		size_t ds_h = GDALGetRasterYSize(ds);
		if(!ds_w || !ds_h) fatal_error("missing width/height");

		if(ds_idx) {
			if(ds_w != w || ds_h != h) fatal_error("size mismatch for inputs");
		} else {
			w = ds_w;
			h = ds_h;
		}
	}

	GDALDriverH dst_driver = GDALGetDriverByName("VRT");
	if(!dst_driver) fatal_error("unrecognized output format (VRT)");
	GDALDatasetH dst_ds = GDALCreateCopy(dst_driver, dst_fn.c_str(), src_ds[0], 0, NULL, NULL, NULL);
	if(!dst_ds) fatal_error("could not create output");

	int band_idx = GDALGetRasterCount(src_ds[0]);
	for(size_t ds_idx=1; ds_idx<src_fn.size(); ds_idx++) {
		int nb = GDALGetRasterCount(src_ds[ds_idx]);

		GDALDatasetH src_vrt_ds = GDALCreateCopy(dst_driver, "", src_ds[ds_idx], 0, NULL, NULL, NULL);
		if(!src_vrt_ds) fatal_error("could not create VRT copy");

		for(int i=0; i<nb; i++) {
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

	for(size_t ds_idx=0; ds_idx<src_fn.size(); ds_idx++) {
		GDALClose(src_ds[ds_idx]);
	}
	GDALClose(dst_ds);

	return 0;
}
