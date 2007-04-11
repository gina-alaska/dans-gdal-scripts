/*
    Copyright 2006 University of Alaska,
    Geographic Information Network of Alaska
    All Rights Reserved
    UAF Confidential Information
*/


#include <common.h>

void usage(char *cmdname) {
	fprintf(stderr, "Usage: %s\n", cmdname);
	fprintf(stderr, "\t-wh <width> <height>\n");
	fprintf(stderr, "\t-origin <left easting> <top northing>\n");
	fprintf(stderr, "\t-res <pixel_size>\n");
	fprintf(stderr, "\t-srs <proj4>\n");
	fprintf(stderr, "\t[-lsb | -msb]\n");
	fprintf(stderr, "\t<input.bil> <output.tif>\n");
	exit(1);
}

int main(int argc, char *argv[]) {
	int i;

	char *src_fn = NULL;
	char *dst_fn = NULL;
	int w=0, h=0;
	double res=0;
	double origin_e=0, origin_n=0;
	double got_en=0;
	char *srs = NULL;
	double ndv=0;
	int got_ndv=0;
	char *datatype = "UINT8";
	double affine[6];
	int got_affine=0;
	char endian=0;

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check for duplicate values
		if(arg[0] == '-' && arg[1]) {
			if(!strcmp(arg, "-wh")) {
				char *endptr;

				if(argp == argc) usage(argv[0]);
				w = strtol(argv[argp++], &endptr, 10);
				if(*endptr) usage(argv[0]);

				if(argp == argc) usage(argv[0]);
				h = strtol(argv[argp++], &endptr, 10);
				if(*endptr) usage(argv[0]);
			} else if(!strcmp(arg, "-affine")) {
				for(i=0; i<6; i++) {
					if(argp == argc) usage(argv[0]);
					char *endptr;
					affine[i] = strtod(argv[argp++], &endptr);
					if(*endptr) usage(argv[0]);
				}
				got_affine++;
			} else if(!strcmp(arg, "-res")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				res = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
			} else if(!strcmp(arg, "-origin")) {
				char *endptr;

				if(argp == argc) usage(argv[0]);
				origin_e = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);

				if(argp == argc) usage(argv[0]);
				origin_n = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);

				got_en++;
			} else if(!strcmp(arg, "-srs")) {
				if(argp == argc) usage(argv[0]);
				srs = argv[argp++];
			} else if(!strcmp(arg, "-ndv")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				ndv = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);

				got_ndv++;
			} else if(!strcmp(arg, "-datatype")) {
				if(argp == argc) usage(argv[0]);
				datatype = argv[argp++];
			} else if(!strcmp(arg, "-lsb")) {
				endian = 'L';
			} else if(!strcmp(arg, "-msb")) {
				endian = 'M';
			} else usage(argv[0]);
		} else {
			if(src_fn && dst_fn) {
				usage(argv[0]);
			} else if(src_fn) {
				dst_fn = arg;
			} else {
				src_fn = arg;
			}
		}
	}

	if(!(
		src_fn && dst_fn &&
		w && h && srs
	)) usage(argv[0]);

	if(!got_affine) {
		if(!(res && got_en)) usage(argv[0]);
		affine[0] = origin_e;
		affine[1] = res;
		affine[2] = 0;
		affine[3] = origin_n;
		affine[4] = 0;
		affine[5] = -res;
	}

	GDALDataType gdal_dt = GDT_Unknown;
	int bytes_per_pixel = 0;
	if(!strcmp(datatype, "UINT8")) {
		gdal_dt = GDT_Byte;
		bytes_per_pixel = 1;
	} else if(!strcmp(datatype, "UINT16")) {
		gdal_dt = GDT_UInt16;
		bytes_per_pixel = 2;
	} else if(!strcmp(datatype, "INT16")) {
		gdal_dt = GDT_Int16;
		bytes_per_pixel = 2;
	} else if(!strcmp(datatype, "UINT32")) {
		gdal_dt = GDT_UInt32;
		bytes_per_pixel = 4;
	} else if(!strcmp(datatype, "INT32")) {
		gdal_dt = GDT_Int32;
		bytes_per_pixel = 4;
	} else if(!strcmp(datatype, "FLOAT32")) {
		gdal_dt = GDT_Float32;
		bytes_per_pixel = 4;
	} else if(!strcmp(datatype, "FLOAT64")) {
		gdal_dt = GDT_Float64;
		bytes_per_pixel = 8;
	}

	int endian_mismatch;
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
	if(!strcmp(src_fn, "-")) {
		fin = stdin;
	} else {
		fin = fopen(src_fn, "r");
	}
	if(!fin) fatal_error("could not open input");

	//////////// open output

	GDALAllRegister();

	GDALDriverH dst_driver = GDALGetDriverByName("GTiff");
	if(!dst_driver) fatal_error("unrecognized output format");
	GDALDatasetH dst_ds = GDALCreate(dst_driver, dst_fn, w, h, 1, gdal_dt, NULL);
	if(!dst_ds) fatal_error("could create dst_dataset");

	GDALSetGeoTransform(dst_ds, affine);

	OGRSpatialReference output_srs;
	if(output_srs.SetFromUserInput(srs) != OGRERR_NONE) {
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

	void *linebuf = malloc_or_die(w * bytes_per_pixel);
	int row;
	for(row=0; row<h; row++) {
		GDALTermProgress((double)row / (double)h, NULL, NULL);
		if((size_t)w != fread(linebuf, bytes_per_pixel, w, fin)) {
			fatal_error("input was short");
		}
		if(endian_mismatch) {
			GDALSwapWords(linebuf, bytes_per_pixel, w, 0);
		}
		GDALRasterIO(dst_band, GF_Write, 0, row, w, 1, linebuf, w, 1, gdal_dt, 0, 0);
	}

	//////////// shutdown

	GDALClose(dst_ds);

	GDALTermProgress(1, NULL, NULL);

	// This error is checked after the output is closed.
	// The script exits with error but the output is
	// still saved to disk.
	if(fread(linebuf, 1, 1, fin)) {
		fatal_error("warning: input had extra data at end\n");
	}

	return 0;
}
