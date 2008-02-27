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
#include "georef.h"

void usage(char *cmdname); // externally defined

void print_georef_usage(FILE *fh) {
	fprintf(fh, "\
Geocoding: \n\
  -s_srs 's_srs'                  Set or override source SRS \n\
  -ll_en left_east lower_north    Set or override lower-left coordinate \n\
  -ul_en left_east lower_north    Set or override upper-left coordinate (don't use both ll_en and ul_en)\n\
  -wh width height                Set or override image size \n\
  -res res_x res_y                Set or override resolution \n\
");
}

void _add_arg_to_list(int *argc_ptr, char ***argv_ptr, char *new_arg) {
	*argv_ptr= (char **)realloc_or_die(*argv_ptr, sizeof(char *) * (*argc_ptr+1));
	(*argv_ptr)[*argc_ptr] = new_arg;
	(*argc_ptr)++;
}

geo_opts_t init_geo_options(int *argc_ptr, char ***argv_ptr) {
	int argc = *argc_ptr;
	char **argv = *argv_ptr;

	int argc_out = 0;
	char **argv_out = NULL;
	_add_arg_to_list(&argc_out, &argv_out, argv[0]);

	geo_opts_t opt;
	opt.s_srs = NULL;
	opt.w=0, opt.h=0;
	opt.got_ll_en = 0;
	opt.got_ul_en = 0;
	opt.given_left_e=0, opt.given_lower_n=0, opt.given_upper_n=0;
	opt.res_x=0, opt.res_y=0;

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(!strcmp(arg, "-s_srs")) {
				if(argp == argc) usage(argv[0]);
				opt.s_srs = argv[argp++];
			} else if(!strcmp(arg, "-ll_en")) {
				if(argp+2 > argc) usage(argv[0]);
				char *endptr;
				opt.given_left_e = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
				opt.given_lower_n = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
				opt.got_ll_en++;
			} else if(!strcmp(arg, "-ul_en")) {
				if(argp+2 > argc) usage(argv[0]);
				char *endptr;
				opt.given_left_e = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
				opt.given_upper_n = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
				opt.got_ul_en++;
			} else if(!strcmp(arg, "-wh")) {
				if(argp+2 > argc) usage(argv[0]);
				char *endptr;
				opt.w = strtol(argv[argp++], &endptr, 10);
				if(*endptr) usage(argv[0]);
				opt.h = strtol(argv[argp++], &endptr, 10);
				if(*endptr) usage(argv[0]);
			} else if(!strcmp(arg, "-res")) {
				char *endptr;
				if(argp == argc) usage(argv[0]);
				opt.res_x = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);

				if(argp == argc) usage(argv[0]);
				opt.res_y = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
			} else {
				_add_arg_to_list(&argc_out, &argv_out, arg);
			}
		} else {
			_add_arg_to_list(&argc_out, &argv_out, arg);
		}
	}

	if(opt.got_ll_en && opt.got_ul_en) fatal_error("don't specify both -ll_en and -ul_en");

	*argc_ptr = argc_out;
	*argv_ptr = argv_out;

	return opt;
}

georef_t init_georef(geo_opts_t *opt, GDALDatasetH ds) {
	if(!ds && !(opt->s_srs && (opt->got_ll_en || opt->got_ul_en) && 
		opt->w && opt->h && opt->res_x && opt->res_y)) fatal_error("not enough information to determine geolocation");

	georef_t georef;
	georef.fwd_xform = NULL;
	georef.fwd_affine = NULL;

	OGRSpatialReferenceH p1 = NULL;
	if(opt->s_srs) {
		p1 = OSRNewSpatialReference(NULL);
		if(OSRImportFromProj4(p1, opt->s_srs) != OGRERR_NONE) fatal_error("cannot parse proj4 definition");
	} else if(ds) {
		const char *wkt = GDALGetProjectionRef(ds);
		if(wkt && strlen(wkt)) {
			//if(VERBOSE) fprintf(stderr, "%s\n", wkt);
			p1 = OSRNewSpatialReference(wkt);
		}
	}

	if(p1) {
		opt->s_srs = NULL;
		OSRExportToProj4(p1, &opt->s_srs);

		OGRSpatialReferenceH p2 = OSRNewSpatialReference(NULL);
		OSRImportFromProj4(p2, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs");

		georef.fwd_xform = OCTNewCoordinateTransformation(p1, p2);
		georef.inv_xform = OCTNewCoordinateTransformation(p2, p1);
	}

	if(ds) {
		if(!opt->w) opt->w = GDALGetRasterXSize(ds);
		if(!opt->h) opt->h = GDALGetRasterYSize(ds);
	}
	if(!opt->w || !opt->h) fatal_error("missing width/height");

	opt->res_x = fabs(opt->res_x);
	opt->res_y = fabs(opt->res_y);

	if((opt->got_ul_en || opt->got_ll_en) && opt->res_x && opt->res_y) {
		if(opt->got_ll_en) {
			opt->given_upper_n = opt->given_lower_n + (double)opt->h*opt->res_y;
			opt->got_ul_en = 1;
		}

		if(!opt->got_ul_en) fatal_error("impossibility");
		georef.fwd_affine = (double *)malloc_or_die(sizeof(double) * 6);
		georef.fwd_affine[0] = opt->given_left_e;  
		georef.fwd_affine[3] = opt->given_upper_n; 
		georef.fwd_affine[1] = opt->res_x; georef.fwd_affine[2] =      0;
		georef.fwd_affine[4] =     0; georef.fwd_affine[5] = -opt->res_y;
	} else if(ds) {
		int has_rotation = 0;
		georef.fwd_affine = (double *)malloc_or_die(sizeof(double) * 6);
		if(GDALGetGeoTransform(ds, georef.fwd_affine) == CE_None) {
			has_rotation = georef.fwd_affine[2] || georef.fwd_affine[4];
		} else {
			georef.fwd_affine = NULL;
		}

		if(opt->res_x && opt->res_y) {
			// if corner coordinate were specified, the first branch of the outer if
			// statement would have been taken
			if(!georef.fwd_affine) fatal_error("missing ll_en/ul_en parameter");

			if(has_rotation) fatal_error("cannot override resolution if source image has rotation");

			georef.fwd_affine[1] =  opt->res_x;
			georef.fwd_affine[5] = -opt->res_y;
		} else {
			if(has_rotation || !georef.fwd_affine) {
				opt->res_x = opt->res_y = 0;
			} else {
				opt->res_x = fabs(georef.fwd_affine[1]);
				opt->res_y = fabs(georef.fwd_affine[5]);
			}
		}

		if(opt->got_ul_en || opt->got_ll_en) {
			if(has_rotation) fatal_error("cannot override ll_en/ul_en if source image has rotation");

			if(!georef.fwd_affine || !opt->res_x || !opt->res_y) fatal_error("missing -res parameter");

			if(opt->got_ll_en) {
				opt->given_upper_n = opt->given_lower_n + (double)opt->h*opt->res_y;
				opt->got_ul_en = 1;
			}

			if(!opt->got_ul_en) fatal_error("impossibility");
			georef.fwd_affine[0] = opt->given_left_e;
			georef.fwd_affine[3] = opt->given_upper_n;
		}
	}

	georef.s_srs = opt->s_srs;
	georef.res_x = opt->res_x;
	georef.res_y = opt->res_y;
	georef.w = opt->w;
	georef.h = opt->h;

	if(georef.fwd_affine) {
		georef.inv_affine = (double *)malloc_or_die(sizeof(double) * 6);
		if(!GDALInvGeoTransform(georef.fwd_affine, georef.inv_affine)) {
			fatal_error("affine is not invertible");
		}
	}

	return georef;
}

void xy2en(
	georef_t *georef,
	double xpos, double ypos,
	double *e_out, double *n_out
) {
	double *affine = georef->fwd_affine;
	if(!affine) fatal_error("missing affine");
	*e_out = affine[0] + affine[1] * xpos + affine[2] * ypos;
	*n_out = affine[3] + affine[4] * xpos + affine[5] * ypos;
}

void en2xy(
	georef_t *georef,
	double east, double north,
	double *x_out, double *y_out
) {
	double *affine = georef->inv_affine;
	if(!affine) fatal_error("missing affine");
	*x_out = affine[0] + affine[1] * east + affine[2] * north;
	*y_out = affine[3] + affine[4] * east + affine[5] * north;
}

void en2ll(
	georef_t *georef,
	double east, double north,
	double *lon_out, double *lat_out
) {
	if(!georef->fwd_xform) fatal_error("missing xform");

	if(!OCTTransform(georef->fwd_xform, 1, &east, &north, NULL)) {
		fatal_error("OCTTransform failed");
	}

	if(north < -90.0 || north > 90.0) fatal_error("latitude out of range");
	// images in latlong projection that cross the dateline can
	// have numbers outside of this range...
	//if(east < -180.0 || east > 180.0) fatal_error("longitude out of range");
	// but it shouldn't be outside of *this* range no matter what!
	if(east < -360.0 || east > 540.0) fatal_error("longitude out of range");

	*lon_out = east;
	*lat_out = north;
}

void ll2en(
	georef_t *georef,
	double lon, double lat,
	double *e_out, double *n_out
) {
	if(!georef->inv_xform) fatal_error("missing xform");

	if(lat < -90.0 || lat > 90.0) fatal_error("latitude out of range");
	// images in latlong projection that cross the dateline can
	// have numbers outside of this range...
	//if(lon < -180.0 || lon > 180.0) fatal_error("longitude out of range");
	// but it shouldn't be outside of *this* range no matter what!
	if(lon < -360.0 || lon > 540.0) fatal_error("longitude out of range");

	if(!OCTTransform(georef->inv_xform, 1, &lon, &lat, NULL)) {
		fatal_error("OCTTransform failed");
	}

	*e_out = lon;
	*n_out = lat;
}

void xy2ll(
	georef_t *georef,
	double x, double y,
	double *lon_out, double *lat_out
) {
	double east, north;
	xy2en(georef, x, y, &east, &north);
	en2ll(georef, east, north, lon_out, lat_out);
}

void ll2xy(
	georef_t *georef,
	double lon, double lat,
	double *x_out, double *y_out
) {
	double east, north;
	ll2en(georef, lon, lat, &east, &north);
	en2xy(georef, east, north, x_out, y_out);
}
