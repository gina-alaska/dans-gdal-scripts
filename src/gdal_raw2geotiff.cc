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
	printf("Usage: %s\n", cmdname.c_str());
	printf("\t-wh <width> <height>\n");
	printf("\t-origin <left easting> <top northing>\n");
	printf("\t-res <pixel_size>\n");
	printf("\t-srs <proj4>\n");
	printf("\t[-datatype { UINT8 | UINT16 | INT16 | UINT32 | INT32 | FLOAT32 | FLOAT64 }]\n");
	printf("\t[-lsb | -msb]\n");
	printf("\t<input.bil> <output.tif>\n");
	exit(1);
}

int main(int argc, char *argv[]) {
	const std::string cmdname = argv[0];
	if(argc == 1) usage(cmdname);
	std::vector<std::string> arg_list = argv_to_list(argc, argv);

	std::string src_fn;
	std::string dst_fn;
	size_t w=0, h=0;
	double res=0;
	double origin_e=0, origin_n=0;
	double got_en=0;
	std::string srs;
	double ndv=0;
	bool got_ndv=0;
	std::string datatype = "UINT8";
	double affine[6];
	bool got_affine=0;
	char endian=0;

	size_t argp = 1;
	while(argp < arg_list.size()) {
		const std::string &arg = arg_list[argp++];
		// FIXME - check for duplicate values
		if(arg[0] == '-' && arg[1]) {
			if(arg == "-wh") {
				char *endptr;

				if(argp == arg_list.size()) usage(cmdname);
				w = (int)strtol(arg_list[argp++].c_str(), &endptr, 10);
				if(*endptr) usage(cmdname);

				if(argp == arg_list.size()) usage(cmdname);
				h = (int)strtol(arg_list[argp++].c_str(), &endptr, 10);
				if(*endptr) usage(cmdname);
			} else if(arg == "-affine") {
				for(int i=0; i<6; i++) {
					if(argp == arg_list.size()) usage(cmdname);
					char *endptr;
					affine[i] = strtod(arg_list[argp++].c_str(), &endptr);
					if(*endptr) usage(cmdname);
				}
				got_affine = 1;
			} else if(arg == "-res") {
				if(argp == arg_list.size()) usage(cmdname);
				char *endptr;
				res = strtod(arg_list[argp++].c_str(), &endptr);
				if(*endptr) usage(cmdname);
			} else if(arg == "-origin") {
				char *endptr;

				if(argp == arg_list.size()) usage(cmdname);
				origin_e = strtod(arg_list[argp++].c_str(), &endptr);
				if(*endptr) usage(cmdname);

				if(argp == arg_list.size()) usage(cmdname);
				origin_n = strtod(arg_list[argp++].c_str(), &endptr);
				if(*endptr) usage(cmdname);

				got_en++;
			} else if(arg == "-srs") {
				if(argp == arg_list.size()) usage(cmdname);
				srs = arg_list[argp++];
			} else if(arg == "-ndv") {
				if(argp == arg_list.size()) usage(cmdname);
				char *endptr;
				ndv = strtod(arg_list[argp++].c_str(), &endptr);
				if(*endptr) usage(cmdname);

				got_ndv = 1;
			} else if(arg == "-datatype") {
				if(argp == arg_list.size()) usage(cmdname);
				datatype = arg_list[argp++];
			} else if(arg == "-lsb") {
				endian = 'L';
			} else if(arg == "-msb") {
				endian = 'M';
			} else {
				usage(cmdname);
			}
		} else {
			if(src_fn.size() && dst_fn.size()) {
				usage(cmdname);
			} else if(src_fn.size()) {
				dst_fn = arg;
			} else {
				src_fn = arg;
			}
		}
	}

	if(!(
		src_fn.size() && dst_fn.size() &&
		w && h && srs.size()
	)) usage(cmdname);

	if(!got_affine) {
		if(!(res && got_en)) usage(cmdname);
		affine[0] = origin_e;
		affine[1] = res;
		affine[2] = 0;
		affine[3] = origin_n;
		affine[4] = 0;
		affine[5] = -res;
	}

	GDALDataType gdal_dt = GDT_Unknown;
	int bytes_per_pixel = 0;
	if(datatype == "UINT8") {
		gdal_dt = GDT_Byte;
		bytes_per_pixel = 1;
	} else if(datatype == "UINT16") {
		gdal_dt = GDT_UInt16;
		bytes_per_pixel = 2;
	} else if(datatype == "INT16") {
		gdal_dt = GDT_Int16;
		bytes_per_pixel = 2;
	} else if(datatype == "UINT32") {
		gdal_dt = GDT_UInt32;
		bytes_per_pixel = 4;
	} else if(datatype == "INT32") {
		gdal_dt = GDT_Int32;
		bytes_per_pixel = 4;
	} else if(datatype == "FLOAT32") {
		gdal_dt = GDT_Float32;
		bytes_per_pixel = 4;
	} else if(datatype == "FLOAT64") {
		gdal_dt = GDT_Float64;
		bytes_per_pixel = 8;
	}

	bool endian_mismatch;
	if(bytes_per_pixel == 1) {
		endian_mismatch = 0;
	} else {
		if(!endian) fatal_error("must specify endian");
		if(CPL_IS_LSB) {
			endian_mismatch = (endian == 'M');
		} else {
			endian_mismatch = (endian == 'L');
		}
	}

	//////////// open input

	FILE *fin;
	if(src_fn == "-") {
		fin = stdin;
	} else {
		fin = fopen(src_fn.c_str(), "r");
	}
	if(!fin) fatal_error("could not open input");

	//////////// open output

	GDALAllRegister();

	GDALDriverH dst_driver = GDALGetDriverByName("GTiff");
	if(!dst_driver) fatal_error("unrecognized output format (GTiff)");
	GDALDatasetH dst_ds = GDALCreate(dst_driver, dst_fn.c_str(), w, h, 1, gdal_dt, NULL);
	if(!dst_ds) fatal_error("couldn't create dst_dataset");

	GDALSetGeoTransform(dst_ds, affine);

	OGRSpatialReference output_srs;
	if(output_srs.SetFromUserInput(srs.c_str()) != OGRERR_NONE) {
		fatal_error("could not parse SRS");
	}
	char *wkt = NULL;
	output_srs.exportToWkt(&wkt);
	if(!wkt) fatal_error("could not convert SRS to WKT");
	GDALSetProjection(dst_ds, wkt);

	GDALRasterBandH dst_band = GDALGetRasterBand(dst_ds, 1);

	if(got_ndv) {
		GDALSetRasterNoDataValue(dst_band, ndv);
	}

	//////////// transfer data

	std::vector<uint8_t> linebuf(w * bytes_per_pixel);
	for(size_t row=0; row<h; row++) {
		GDALTermProgress((double)row / h, NULL, NULL);
		if(w != fread(&linebuf[0], bytes_per_pixel, w, fin)) {
			fatal_error("input was short");
		}
		if(endian_mismatch) {
			GDALSwapWords(&linebuf[0], bytes_per_pixel, w, 0);
		}
		GDALRasterIO(dst_band, GF_Write, 0, row, w, 1, &linebuf[0], w, 1, gdal_dt, 0, 0);
	}

	//////////// shutdown

	GDALClose(dst_ds);

	GDALTermProgress(1, NULL, NULL);

	// This error is checked after the output is closed.
	// The script exits with error but the output is
	// still saved to disk.
	if(fread(&linebuf[0], 1, 1, fin)) {
		fatal_error("warning: input had extra data at end\n");
	}

	return 0;
}
