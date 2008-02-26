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
#include "polygon.h"
#include "debugplot.h"
#include "georef.h"

#ifndef PI
#define PI 3.141592653
#endif

int VERBOSE = 0;

void usage(char *cmdname) {
	fprintf(stderr, "Usage:\n  %s image_name [options]\n", cmdname);
	fprintf(stderr, "\
\n\
Geocoding: \n\
  -s_srs 's_srs'                  Set or override source SRS \n\
  -ll_en left_east lower_north    Set or override lower-left coordinate \n\
  -ul_en left_east lower_north    Set or override upper-left coordinate (don't use both ll_en and ul_en)\n\
  -wh width height                Set or override image size \n\
  -res res_x res_y                Set or override resolution \n\
\n\
Inspection: \n\
  -inspect-rect4                  Attempt to find 4-sided bounding polygon  \n\
  -inspect-contour                Trace out entire polygon of image \n\
  -nodataval 'val [val ...]'      Specify value of no-data pixels \n\
  -ndv-toler val                  Tolerance for deciding if a pixel matches nodataval \n\
  -b band_id -b band_id ...       Bands to inspect (default is all bands) \n\
  -skip-erosion                   Don't use erosion filter \n\
  -major-ring                     Output only the biggest outer ring \n\
  -no-donuts                      Output only top-level rings \n\
  -min-ring-area val              Drop rings with less than this area (in square pixels) \n\
  -dp-toler val                   Tolerance for point reduction (in pixel units) \n\
                                      default is 2 pixels\n\
  -report fn.ppm                  Output graphical report of bounds found \n\
  -wkt-xy fn.wkt                  Output bounds in WKT format (x/y coords) \n\
  -wkt-en fn.wkt                  Output bounds in WKT format (east/north coords) \n\
  -wkt-ll fn.wkt                  Output bounds in WKT format (lat/lon coords) \n\
  -wkt-split-polys                Output several polygons rather than one multipolygon\n\
  -mask-out fn.pbm                Output mask of bounding polygon in PBM format \n\
\n\
Misc: \n\
  -v                              Verbose\n\
\n\
Examples:\n\
  Output basic geocoding info:\n\
    gdal_list_corners raster.tif > geocode.yaml \n\
  Inspect image to find corners of actual data (arbitrary four-sided region):\n\
    gdal_list_corners raster.tif -inspect-rect4 -nodataval 0 > geocode.yaml\n\
  Inspect image and output contour of data region:\n\
    gdal_list_corners raster.tif -inspect-contour -nodataval 0 -wkt-ll outline.wkt > geocode.yaml \n\
  Same as above but polygon actually follows border pixel-by-pixel:\n\
    gdal_list_corners raster.tif -inspect-contour -nodataval 0 -dp-toler 0 -wkt-ll outline.wkt > geocode.yaml \n\
\n\
");
	exit(1);
}

int parse_list_of_doubles(char *input, int *num_out, double **list_out);
void setup_ndv_list(GDALDatasetH ds, int bandlist_size, int *bandlist, int *num_ndv, double **ndv_list);
unsigned char *get_mask_for_dataset(GDALDatasetH ds, int bandlist_size, int *bandlist, 
	int num_ndv, double *ndv_list, double ndv_tolerance, report_image_t *dbuf);
vertex_t calc_centroid_from_mask(unsigned char *mask, int w, int h);
ring_t calc_rect4_from_mask(unsigned char *mask, int w, int h, report_image_t *dbuf);
mpoly_t calc_ring_from_mask(unsigned char *mask, int w, int h,
	report_image_t *dbuf, int major_ring_only, int no_donuts, double min_ring_area);
unsigned char *erode_mask(unsigned char *in_mask, int w, int h);

int main(int argc, char **argv) {
	char *input_raster_fn = NULL;
	char *s_srs = NULL;

	int w=0, h=0;
	int got_ll_en = 0;
	int got_ul_en = 0;
	double given_left_e=0, given_lower_n=0, given_upper_n=0;
	double res_x=0, res_y=0;
	int inspect_rect4 = 0;
	int instpect_contour = 0;
	int num_ndv = 0;
	double *ndv_list = NULL;
	double ndv_tolerance = 0;
	char *debug_report = NULL;
	int inspect_numbands = 0;
	int *inspect_bandids = NULL;
	int split_polys = 0;
	char *wkt_xy_fn = NULL;
	char *wkt_en_fn = NULL;
	char *wkt_ll_fn = NULL;
	char *mask_out_fn = NULL;
	int major_ring_only = 0;
	int no_donuts = 0;
	double min_ring_area = 0;
	double reduction_tolerance = 2;
	int skip_erosion = 0;

	int i, j;

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(!strcmp(arg, "-v")) {
				VERBOSE++;
			} else if(!strcmp(arg, "-s_srs")) {
				if(argp == argc) usage(argv[0]);
				s_srs = argv[argp++];
			} else if(!strcmp(arg, "-ll_en")) {
				if(argp+2 > argc) usage(argv[0]);
				char *endptr;
				given_left_e = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
				given_lower_n = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
				got_ll_en++;
			} else if(!strcmp(arg, "-ul_en")) {
				if(argp+2 > argc) usage(argv[0]);
				char *endptr;
				given_left_e = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
				given_upper_n = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
				got_ul_en++;
			} else if(!strcmp(arg, "-wh")) {
				if(argp+2 > argc) usage(argv[0]);
				char *endptr;
				w = strtol(argv[argp++], &endptr, 10);
				if(*endptr) usage(argv[0]);
				h = strtol(argv[argp++], &endptr, 10);
				if(*endptr) usage(argv[0]);
			} else if(!strcmp(arg, "-res")) {
				char *endptr;
				if(argp == argc) usage(argv[0]);
				res_x = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);

				if(argp == argc) usage(argv[0]);
				res_y = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
			} else if(!strcmp(arg, "-inspect-rect4")) {
				inspect_rect4++;
			} else if(!strcmp(arg, "-inspect-contour")) {
				instpect_contour++;
			} else if(!strcmp(arg, "-nodataval")) {
				if(argp == argc) usage(argv[0]);
				int result = parse_list_of_doubles(argv[argp++], &num_ndv, &ndv_list);
				if(result) fatal_error("input to -nodataval must be space-separated list of numbers");
			} else if(!strcmp(arg, "-report")) {
				if(argp == argc) usage(argv[0]);
				debug_report = argv[argp++];
			} else if(!strcmp(arg, "-b")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				int bandid = strtol(argv[argp++], &endptr, 10);
				if(*endptr) usage(argv[0]);
				inspect_bandids = (int *)realloc_or_die(inspect_bandids,
					sizeof(int)*(inspect_numbands+1));
				inspect_bandids[inspect_numbands++] = bandid;
			} else if(!strcmp(arg, "-skip-erosion")) {
				skip_erosion++;
			} else if(!strcmp(arg, "-wkt-split-polys")) {
				split_polys++;
			} else if(!strcmp(arg, "-wkt-xy")) {
				if(argp == argc) usage(argv[0]);
				wkt_xy_fn = argv[argp++];
			} else if(!strcmp(arg, "-wkt-en")) {
				if(argp == argc) usage(argv[0]);
				wkt_en_fn = argv[argp++];
			} else if(!strcmp(arg, "-wkt-ll")) {
				if(argp == argc) usage(argv[0]);
				wkt_ll_fn = argv[argp++];
			} else if(!strcmp(arg, "-mask-out")) {
				if(argp == argc) usage(argv[0]);
				mask_out_fn = argv[argp++];
			} else if(!strcmp(arg, "-major-ring")) {
				major_ring_only = 1;
			} else if(!strcmp(arg, "-no-donuts")) {
				no_donuts = 1;
			} else if(!strcmp(arg, "-min-ring-area")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				min_ring_area = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
			} else if(!strcmp(arg, "-dp-toler")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				reduction_tolerance = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
			} else if(!strcmp(arg, "-ndv-toler")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				ndv_tolerance = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
			} else usage(argv[0]);
		} else {
			if(input_raster_fn) usage(argv[0]);
			input_raster_fn = arg;
		}
	}

	if(!input_raster_fn && !(s_srs && (got_ll_en || got_ul_en) && w && h && res_x && res_y)) usage(argv[0]);
	if(got_ll_en && got_ul_en) usage(argv[0]);

	int do_wkt_output = wkt_xy_fn || wkt_en_fn || wkt_ll_fn;
	int do_inspect = inspect_rect4 || instpect_contour;
	if(do_inspect && !input_raster_fn) fatal_error("must specify filename of image");
	if((do_wkt_output || mask_out_fn) && !do_inspect) fatal_error(
		"must specify -inspect-rect4 or -inspect-contour");

	if(major_ring_only && min_ring_area) fatal_error(
		"-major-ring and -min-ring-area options cannot both be used at the same time");
	if(major_ring_only && no_donuts) fatal_error(
		"-major-ring and -no-donuts options cannot both be used at the same time");

	GDALAllRegister();

	GDALDatasetH ds = NULL;
	if(input_raster_fn) {
		ds = GDALOpen(input_raster_fn, GA_ReadOnly);
		if(!ds) fatal_error("open failed");
	}

	if(do_inspect && !inspect_numbands) {
		inspect_numbands = GDALGetRasterCount(ds);
		inspect_bandids = (int *)malloc_or_die(sizeof(int)*inspect_numbands);
		for(i=0; i<inspect_numbands; i++) inspect_bandids[i] = i+1;
	}

	CPLPushErrorHandler(CPLQuietErrorHandler);

	georef_t georef;
	georef.fwd_xform = NULL;
	georef.fwd_affine = NULL;

	OGRSpatialReferenceH p1 = NULL;
	if(s_srs) {
		p1 = OSRNewSpatialReference(NULL);
		if(OSRImportFromProj4(p1, s_srs) != OGRERR_NONE) fatal_error("cannot parse proj4 definition");
	} else if(ds) {
		const char *wkt = GDALGetProjectionRef(ds);
		if(wkt && strlen(wkt)) {
			//if(VERBOSE) fprintf(stderr, "%s\n", wkt);
			p1 = OSRNewSpatialReference(wkt);
		}
	}

	if(p1) {
		s_srs = NULL;
		OSRExportToProj4(p1, &s_srs);

		OGRSpatialReferenceH p2 = OSRNewSpatialReference(NULL);
		OSRImportFromProj4(p2, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs");

		georef.fwd_xform = OCTNewCoordinateTransformation(p1, p2);
		georef.inv_xform = OCTNewCoordinateTransformation(p2, p1);
	}

	if(ds) {
		if(!w) w = GDALGetRasterXSize(ds);
		if(!h) h = GDALGetRasterYSize(ds);
	}
	if(!w || !h) fatal_error("missing width/height");

	res_x = fabs(res_x);
	res_y = fabs(res_y);

	if((got_ul_en || got_ll_en) && res_x && res_y) {
		if(got_ll_en) {
			given_upper_n = given_lower_n + (double)h*res_y;
			got_ul_en = 1;
		}

		if(!got_ul_en) fatal_error("impossibility");
		georef.fwd_affine = (double *)malloc_or_die(sizeof(double) * 6);
		georef.fwd_affine[0] = given_left_e;  
		georef.fwd_affine[3] = given_upper_n; 
		georef.fwd_affine[1] = res_x; georef.fwd_affine[2] =      0;
		georef.fwd_affine[4] =     0; georef.fwd_affine[5] = -res_y;
	} else if(ds) {
		int has_rotation = 0;
		georef.fwd_affine = (double *)malloc_or_die(sizeof(double) * 6);
		if(GDALGetGeoTransform(ds, georef.fwd_affine) == CE_None) {
			has_rotation = georef.fwd_affine[2] || georef.fwd_affine[4];
		} else {
			georef.fwd_affine = NULL;
		}

		if(res_x && res_y) {
			// if corner coordinate were specified, the first branch of the outer if
			// statement would have been taken
			if(!georef.fwd_affine) fatal_error("missing ll_en/ul_en parameter");

			if(has_rotation) fatal_error("cannot override resolution if source image has rotation");

			georef.fwd_affine[1] =  res_x;
			georef.fwd_affine[5] = -res_y;
		} else {
			if(has_rotation || !georef.fwd_affine) {
				res_x = res_y = 0;
			} else {
				res_x = fabs(georef.fwd_affine[1]);
				res_y = fabs(georef.fwd_affine[5]);
			}
		}

		if(got_ul_en || got_ll_en) {
			if(has_rotation) fatal_error("cannot override ll_en/ul_en if source image has rotation");

			if(!georef.fwd_affine || !res_x || !res_y) fatal_error("missing -res parameter");

			if(got_ll_en) {
				given_upper_n = given_lower_n + (double)h*res_y;
				got_ul_en = 1;
			}

			if(!got_ul_en) fatal_error("impossibility");
			georef.fwd_affine[0] = given_left_e;
			georef.fwd_affine[3] = given_upper_n;
		}
	}

	if(georef.fwd_affine) {
		georef.inv_affine = (double *)malloc_or_die(sizeof(double) * 6);
		if(!GDALInvGeoTransform(georef.fwd_affine, georef.inv_affine)) {
			fatal_error("affine is not invertible");
		}
	}

	if((wkt_en_fn || wkt_ll_fn) && !georef.fwd_affine) fatal_error("missing affine transform");
	if(wkt_ll_fn && !georef.fwd_xform) fatal_error("missing coordinate transform");

	report_image_t *dbuf = NULL;
	unsigned char *mask = NULL;
	if(do_inspect) {
		setup_ndv_list(ds, inspect_numbands, inspect_bandids, &num_ndv, &ndv_list);

		if(debug_report) {
			dbuf = create_plot(w, h);
			if(inspect_rect4) dbuf->mode = PLOT_RECT4;
			else dbuf->mode = PLOT_CONTOURS;
		}

		mask = get_mask_for_dataset(ds, inspect_numbands, inspect_bandids,
			num_ndv, ndv_list, ndv_tolerance, dbuf);
		if(!skip_erosion) {
			unsigned char *eroded_mask = erode_mask(mask, w, h);
			free(mask);
			mask = eroded_mask;
		}
	}

	// output phase

	printf("width: %d\nheight: %d\n", w, h);

	if(ds) {
		int band_count = GDALGetRasterCount(ds);
		const char *datatypes = "";
		for(i=0; i<band_count; i++) {
			GDALRasterBandH band = GDALGetRasterBand(ds, i+1);
			GDALDataType gdt = GDALGetRasterDataType(band);
			const char *dt = GDALGetDataTypeName(gdt);
			if(i) {
				char *join_str = (char *)malloc_or_die(
					strlen(datatypes) + 1 + strlen(dt) + 1);
				sprintf(join_str, "%s,%s", datatypes, dt);
				datatypes = join_str;
			} else {
				datatypes = dt;
			}
		}
		printf("num_bands: %d\n", band_count);
		printf("datatype: %s\n", datatypes);
	}

	if(s_srs && strlen(s_srs)) {
		printf("s_srs: '%s'\n", s_srs);
	}
	if(res_x && res_y) printf("res: %.15f %.15f\n", res_x, res_y);
	if(georef.fwd_affine) {
		printf("affine:\n");
		for(i=0; i<6; i++) printf("  - %.15f\n", georef.fwd_affine[i]);
	}

	double lon, lat;
	double east, north;

	vertex_t center;
	printf("center:\n");
	center = (vertex_t){ (double)w/2.0, (double)h/2.0 };
	if(georef.fwd_xform && georef.fwd_affine) {
		xy2ll(&georef, center.x, center.y, &lon, &lat);
		printf("  lon: %.15f\n", lon);
		printf("  lat: %.15f\n", lat);
	}
	if(georef.fwd_affine) {
		xy2en(&georef, center.x, center.y, &east, &north);
		printf("  east: %.15f\n", east);
		printf("  north: %.15f\n", north);
	}
	printf("  x: %.15f\n", center.x);
	printf("  y: %.15f\n", center.y);

	if(do_inspect) {
		vertex_t centroid;
		printf("centroid:\n");
		centroid = calc_centroid_from_mask(mask, w, h);
		if(georef.fwd_xform && georef.fwd_affine) {
			xy2ll(&georef, centroid.x, centroid.y, &lon, &lat);
			printf("  lon: %.15f\n", lon);
			printf("  lat: %.15f\n", lat);
		}
		if(georef.fwd_affine) {
			xy2en(&georef, centroid.x, centroid.y, &east, &north);
			printf("  east: %.15f\n", east);
			printf("  north: %.15f\n", north);
		}
		printf("  x: %.15f\n", centroid.x);
		printf("  y: %.15f\n", centroid.y);
	}

	mpoly_t *bpoly = NULL;

	if(inspect_rect4) {
		ring_t rect4 = calc_rect4_from_mask(mask, w, h, dbuf);

		bpoly = (mpoly_t *)malloc_or_die(sizeof(mpoly_t));
		bpoly->num_rings = 1;
		bpoly->rings = (ring_t *)malloc_or_die(sizeof(ring_t));
		bpoly->rings[0] = rect4;

		if(rect4.npts == 4) {
			char *labels[] = { "upper_left", "upper_right", "lower_right", "lower_left" };
			if(georef.fwd_xform && georef.fwd_affine) {
				printf("geometry_ll:\n  type: rectangle4\n");
				for(i=0; i<4; i++) {
					xy2ll(&georef, rect4.pts[i].x, rect4.pts[i].y, &lon, &lat);
					printf("  %s_lon: %.15f\n", labels[i], lon);
					printf("  %s_lat: %.15f\n", labels[i], lat);
				}
			}
			if(georef.fwd_affine) {
				printf("geometry_en:\n  type: rectangle4\n");
				for(i=0; i<4; i++) {
					xy2en(&georef, rect4.pts[i].x, rect4.pts[i].y, &east, &north);
					printf("  %s_east: %.15f\n", labels[i], east);
					printf("  %s_north: %.15f\n", labels[i], north);
				}
			}
			printf("geometry_xy:\n  type: rectangle4\n");
			for(i=0; i<4; i++) {
				printf("  %s_x: %.15f\n", labels[i], rect4.pts[i].x);
				printf("  %s_y: %.15f\n", labels[i], rect4.pts[i].y);
			}
		}
	} else {
		char *e_labels[] = { "left", "mid", "right" };
		double e_pos[] = { 0, (double)w/2.0, w };
		char *n_labels[] = { "upper", "mid", "lower" };
		double n_pos[] = { 0, (double)h/2.0, h };
		if(georef.fwd_xform && georef.fwd_affine) {
			printf("geometry_ll:\n  type: rectangle8\n");
			for(i=0; i<3; i++) for(j=0; j<3; j++) {
				if(!strcmp(e_labels[i], "mid") && !strcmp(n_labels[j], "mid")) continue;
				xy2ll(&georef, e_pos[i], n_pos[j], &lon, &lat);
				printf("  %s_%s_lon: %.15f\n", n_labels[j], e_labels[i], lon);
				printf("  %s_%s_lat: %.15f\n", n_labels[j], e_labels[i], lat);
			}
		}
		if(georef.fwd_affine) {
			printf("geometry_en:\n  type: rectangle8\n");
			for(i=0; i<3; i++) for(j=0; j<3; j++) {
				if(!strcmp(e_labels[i], "mid") && !strcmp(n_labels[j], "mid")) continue;
				xy2en(&georef, e_pos[i], n_pos[j], &east, &north);
				printf("  %s_%s_east: %.15f\n", n_labels[j], e_labels[i], east);
				printf("  %s_%s_north: %.15f\n", n_labels[j], e_labels[i], north);
			}
		}
		printf("geometry_xy:\n  type: rectangle8\n");
		for(i=0; i<3; i++) for(j=0; j<3; j++) {
			if(!strcmp(e_labels[i], "mid") && !strcmp(n_labels[j], "mid")) continue;
			printf("  %s_%s_x: %.15f\n", n_labels[j], e_labels[i], e_pos[i]);
			printf("  %s_%s_y: %.15f\n", n_labels[j], e_labels[i], n_pos[j]);
		}
	}

	if(instpect_contour) {
		bpoly = (mpoly_t *)malloc_or_die(sizeof(mpoly_t));
		*bpoly = calc_ring_from_mask(mask, w, h, dbuf, major_ring_only, no_donuts, min_ring_area);
		if(bpoly->num_rings == 0) bpoly = 0;
	}

	if(mask_out_fn) {
		if(!bpoly) fatal_error("missing bound polygon");
		mask_from_mpoly(bpoly, w, h, mask_out_fn);
	}

	if(bpoly && reduction_tolerance > 0) {
		*bpoly = compute_reduced_pointset(bpoly, reduction_tolerance);
	}

	if(dbuf && dbuf->mode == PLOT_CONTOURS) {
		if(!bpoly) fatal_error("missing bound polygon");
		debug_plot_rings(bpoly, dbuf);
	}

	if(do_wkt_output) {
		if(!bpoly) fatal_error("missing bound polygon");

		int num_outer=0, num_inner=0, total_pts=0;
		for(i=0; i<bpoly->num_rings; i++) {
			if(bpoly->rings[i].is_hole) num_inner++;
			else num_outer++;
			total_pts += bpoly->rings[i].npts;
		}
		fprintf(stderr, "Found %d outer rings and %d holes with a total of %d vertices.\n", num_outer, num_inner, total_pts);

		if(wkt_xy_fn) output_wkt_mpoly(wkt_xy_fn, *bpoly, split_polys);

		if(wkt_en_fn) {
			if(!georef.fwd_affine) fatal_error("missing affine transform");
			mpoly_t *en_poly = mpoly_xy2en(&georef, bpoly);
			output_wkt_mpoly(wkt_en_fn, *en_poly, split_polys);
		}

		if(wkt_ll_fn) {
			if(!georef.fwd_affine) fatal_error("missing affine transform");
			if(!georef.fwd_xform) fatal_error("missing coordinate transform");
			mpoly_t *ll_poly = mpoly_xy2ll_with_interp(&georef, bpoly, .2);
			output_wkt_mpoly(wkt_ll_fn, *ll_poly, split_polys);
		}
	}

	if(dbuf) write_plot(dbuf, debug_report);

	CPLPopErrorHandler();

	return 0;
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

vertex_t calc_centroid_from_mask(unsigned char *mask, int w, int h) {
	int mask_rowlen = (w+7)/8;

	long weight_x=0, weight_y=0, num_datavals=0;
	int i, j;
	for(j=0; j<h; j++) {
		unsigned char mask_bitp = 1;
		unsigned char *mask_bytep = mask + mask_rowlen*j;
		for(i=0; i<w; i++) {
			if(*mask_bytep & mask_bitp) {
				weight_x += i;
				weight_y += j;
				num_datavals++;
			}
			mask_bitp <<= 1;
			if(!mask_bitp) {
				mask_bitp = 1;
				mask_bytep++;
			}
		}
	}

	return (vertex_t){
		(double)weight_x / (double)num_datavals,
		(double)weight_y / (double)num_datavals
	};
}

typedef struct {
	vertex_t p0, p1;
	double angle;
	double seg_len;
	int group;
} edge_t;

typedef struct {
	double arc_len;
	double wx, wy;
	double avg_ang;
	int use;
	edge_t best_edge;
	double sort_key;
} edge_group_t;

double ang_diff(double a1, double a2) {
	double d = fabs(a1 - a2);
	return d<=180.0 ? d : 360.0-d;
}

ring_t calc_rect4_from_mask(unsigned char *mask, int w, int h, report_image_t *dbuf) {
	int i, j;

	int mask_rowlen = (w+7)/8;
	int *chrows_left = (int *)malloc_or_die(sizeof(int) * h);
	int *chrows_right = (int *)malloc_or_die(sizeof(int) * h);
	for(i=0; i<h; i++) {
		chrows_left[i] = w;
		chrows_right[i] = -1;
	}
	for(j=0; j<h; j++) {
		int left = w;
		int right = -1;
		unsigned char mask_bitp = 1;
		unsigned char *mask_bytep = mask + mask_rowlen*j;
		for(i=0; i<w; i++) {
			if(*mask_bytep & mask_bitp) {
				if(left > i) left = i;
				if(right < i) right = i;
			}
			mask_bitp <<= 1;
			if(!mask_bitp) {
				mask_bitp = 1;
				mask_bytep++;
			}
		}
		if(chrows_left[j] > left) chrows_left[j] = left;
		if(chrows_right[j] < right) chrows_right[j] = right;
	}

	int fulcrum_x=-1, fulcrum_y=-1;
	for(j=0; j<h; j++) {
		if(chrows_left[j] <= chrows_right[j]) {
			fulcrum_x = chrows_right[j];
			fulcrum_y = j;
			break;
		}
	}
	if(fulcrum_x<0) fatal_error("image was empty");
	//if(VERBOSE) fprintf(stderr, "start point: %d,%d\n", fulcrum_x, fulcrum_y);

	edge_t *all_edges = NULL;
	int num_edges = 0;

	int chop_dx = 1, chop_dy = 0;
	for(;;) {
		if(dbuf && dbuf->mode == PLOT_RECT4) plot_point_big(dbuf, fulcrum_x, fulcrum_y, 0, 255, 0);

		int best_dx = -chop_dx, best_dy = -chop_dy;
		int best_x=-1, best_y=-1;
		for(j=0; j<h; j++) {
			int l = chrows_left[j];
			int r = chrows_right[j];
			if(l > r) continue;
			int pix_dy = j-fulcrum_y;

			int lr;
			for(lr=0; lr<2; lr++) {
				i = lr ? r : l;
				int pix_dx = i-fulcrum_x;
				if(pix_dx*chop_dy >= pix_dy*chop_dx) continue;
				if(pix_dx*best_dy <  pix_dy*best_dx) continue;
				if(pix_dx*best_dy == pix_dy*best_dx) {
					int pdist = pix_dx*pix_dx + pix_dy*pix_dy;
					int bdist = best_dx*best_dx + best_dy*best_dy;
					if(pdist < bdist) continue;
				}
				//if(VERBOSE) fprintf(stderr, "%d %d   %.15f\n", i, j, atan2((double)pix_dy, (double)pix_dx)*180.0/PI);
				best_dx = pix_dx; best_dy = pix_dy;
				best_x = i; best_y = j;
			}
		}
		//if(VERBOSE) fprintf(stderr, "  f=[%3d,%3d] ", best_x, best_y);
		//if(VERBOSE) fprintf(stderr, "  a=[% 3.1f]\n", angle);

		all_edges = (edge_t *)realloc_or_die(all_edges, (num_edges+1)*sizeof(edge_t));
		all_edges[num_edges].p0.x = fulcrum_x;
		all_edges[num_edges].p0.y = fulcrum_y;
		all_edges[num_edges].p1.x = best_x;
		all_edges[num_edges].p1.y = best_y;
		num_edges++;

		if(best_x<0) fatal_error("could not find new fulcrum");
		fulcrum_x = best_x; fulcrum_y = best_y;
		if(chop_dy < 0 && best_dy >= 0) break;
		chop_dx = best_dx; chop_dy = best_dy;
	}

	for(i=0; i<num_edges; i++) {
		edge_t e = all_edges[i];
		double dx = e.p1.x - e.p0.x;
		double dy = e.p1.y - e.p0.y;
		e.angle = atan2(dy, dx)*180.0/PI;
		e.seg_len = sqrt(dx*dx + dy*dy);
		e.group = -1;
		all_edges[i] = e;
	}

	// input consists of a single point or line in this case
	// (in other words it is zero or one dimensional)
	if(num_edges < 3) fatal_error("convex hull has less than three sides");

	int num_groups = 0;
	all_edges[0].group = (num_groups++);
	for(i=0; i<num_edges; i++) {
		edge_t l = all_edges[i];
		edge_t r = all_edges[(i+1) % num_edges];
		double len = l.seg_len + r.seg_len;
		double adiff = ang_diff(l.angle, r.angle);

		// check whether l and r should be in the same group
		if(len > adiff * 5.0 && adiff < 15.0) { // FIXME - arbitrary
			if(i<num_edges-1) {
				// set r.group = l.group
				all_edges[i+1].group = l.group;
			} else {
				// wrapped around... set the tail group to be
				// equal to the head group
				for(j=num_edges-1; j; j--) {
					if(all_edges[j].group != l.group) {
						j++; break;
					}
				}
				for(; j<num_edges; j++) {
					all_edges[j].group = all_edges[0].group;
				}
			}
		} else {
			// create new group for next edge (if we are at the
			// end of the loop then no action is necessary)
			if(i<num_edges-1) {
				all_edges[i+1].group = (num_groups++);
			}
		}

		if(VERBOSE) {
			fprintf(stderr, "a=%.15f  l=%.15f  ", l.angle, l.seg_len);
			fprintf(stderr, "  l2=%.15f  ad=%.15f  ratio=%.15f", len, adiff, len/adiff);
			l = all_edges[i];
			r = all_edges[(i+1) % num_edges];
			fprintf(stderr, "  lg=%d rg=%d\n", l.group, r.group);
		}
	}
	if(VERBOSE) fprintf(stderr, "num groups: %d\n", num_groups);
	for(i=0; i<num_edges; i++) {
		if(all_edges[i].group < 0) fatal_error("edge not assigned to a group");
		//all_edges[i].group = (num_groups++);
	}
	//if(VERBOSE) fprintf(stderr, "num groups: %d\n", num_groups);

	if(VERBOSE) for(i=0; i<num_edges; i++) {
		edge_t l = all_edges[i];
		edge_t r = all_edges[(i+1) % num_edges];
		double len = l.seg_len + r.seg_len;
		double adiff = ang_diff(l.angle, r.angle);
		fprintf(stderr, "a=%.15f  l=%.15f  ", l.angle, l.seg_len);
		fprintf(stderr, "  l2=%.15f  ad=%.15f  ratio=%.15f", len, adiff, len/adiff);
		fprintf(stderr, "  group=%d\n", l.group);
	}

	edge_group_t *groups = (edge_group_t *)malloc_or_die(sizeof(edge_group_t) * num_groups);
	for(i=0; i<num_groups; i++) {
		groups[i].arc_len = 0;
		groups[i].wx = 0;
		groups[i].wy = 0;
	}
	for(i=0; i<num_edges; i++) {
		edge_t e = all_edges[i];
		int eg = e.group;
		if(eg < 0 || eg >= num_groups) {
			fprintf(stderr, "i=%d, g=%d, num_groups=%d\n", i, eg, num_groups);
			fatal_error("group out of range");
		}
		groups[eg].arc_len += e.seg_len;
		groups[eg].wx += e.seg_len * cos(e.angle / 180.0 * PI);
		groups[eg].wy += e.seg_len * sin(e.angle / 180.0 * PI);
	}
	for(i=0; i<num_groups; i++) {
		if(VERBOSE) fprintf(stderr, "group %d: l=%.15f\n", i, groups[i].arc_len);
		if(groups[i].arc_len > (w+h)/10) { // FIXME - arbitrary
			groups[i].use = 1;
			groups[i].avg_ang = atan2(groups[i].wy, groups[i].wx) * 180.0 / PI;
		} else {
			groups[i].use = 0;
		}
	}
	j=0;
	for(i=0; i<num_groups; i++) {
		if(groups[i].use) {
			groups[j++] = groups[i];
		}
	}
	num_groups = j;
	double top_edge_angle = 0;
	if(VERBOSE) fprintf(stderr, "num groups: %d\n", num_groups);
	for(i=0; i<num_groups; i++) {
		// FIXME - instead of choosing an existing edge close to avg_ang, it would
		// be better to create a new edge with the desired angle and with proper
		// p0.x,p0.y value
		groups[i].best_edge = all_edges[0];
		for(j=0; j<num_edges; j++) {
			double d1 = ang_diff(groups[i].avg_ang, all_edges[j].angle);
			double d2 = ang_diff(groups[i].avg_ang, groups[i].best_edge.angle);
			if(d1 < d2) groups[i].best_edge = all_edges[j];
		}
		if(VERBOSE) fprintf(stderr, "group %d: l=%.15f  a=%.15f  b=%.15f\n", i, groups[i].arc_len, groups[i].avg_ang, groups[i].best_edge.angle);
		double ang = groups[i].best_edge.angle;
		if(i==0 || (fabs(ang) < fabs(top_edge_angle))) top_edge_angle = ang;
	}
	for(i=0; i<num_groups; i++) {
		double a = groups[i].best_edge.angle - top_edge_angle;
		if(a < 0.0) a += 360.0;
		groups[i].sort_key = a;
	}
	// bubble sort - start at top edge and go clockwise
	for(i=0; i<num_groups; i++) {
		for(j=num_groups-1; j>i; j--) {
			edge_group_t e1 = groups[j-1];
			edge_group_t e2 = groups[j];
			if(e1.sort_key > e2.sort_key) {
				groups[j-1] = e2;
				groups[j] = e1;
			}
		}
	}
	if(VERBOSE) fprintf(stderr, "sorted:\n");
	if(VERBOSE) for(i=0; i<num_groups; i++) {
		fprintf(stderr, "group %d: l=%.15f  a=%.15f  s=%.15f\n", i, groups[i].arc_len, groups[i].best_edge.angle, groups[i].sort_key);
	}

	//if(VERBOSE) fprintf(stderr, "%d edges\n", num_groups);
	vertex_t *verts = (vertex_t *)malloc_or_die(sizeof(vertex_t) * num_groups);
	for(i=0; i<num_groups; i++) {
		j = i ? i-1 : num_groups-1;
		edge_t e1 = groups[i].best_edge;
		edge_t e2 = groups[j].best_edge;
		line_line_intersection(
			e1.p0, e1.p1,
			e2.p0, e2.p1,
			&verts[i]);
		if(VERBOSE) fprintf(stderr, "vert[%d] = %.15f, %.15f\n", i, verts[i].x, verts[i].y);
	}

	if(dbuf && dbuf->mode == PLOT_RECT4) {
		for(i=0; i<num_groups; i++) {
			j = i<num_groups-1 ? i+1 : 0;
			plot_line(dbuf, verts[i], verts[j], 255, 0, 0);
			plot_point_big(dbuf, verts[i].x, verts[i].y, 255, 255, 0);
			plot_point_big(dbuf, verts[j].x, verts[j].y, 255, 255, 0);
		}
	}

	ring_t rect4;
	if(num_groups == 4) {
		rect4.npts = 4;
		rect4.pts = verts;
	} else {
		rect4.npts = 0;
		rect4.pts = NULL;
		fprintf(stderr, "could not find a 4-sided bounding polygon\n");
	}
	return rect4;
}

typedef struct {
	int top_y;
	int bottom_y;
	int *pts;
	int top_linkage;
	int bottom_linkage;
} descender_t;

typedef struct {
	int num_transitions;
	int *openings;
	int *closings;
	int *descender_ids;
} rowstat_t;

int create_descender_pair(int *num_descenders, descender_t **descenders, int y, int max_pts) {
	int n = *num_descenders;
	*descenders = (descender_t *)realloc_or_die(*descenders, sizeof(descender_t)*(n+2));
	descender_t *d1 = (*descenders) + n;
	descender_t *d2 = (*descenders) + n+1;
	d1->top_y = d2->top_y = y;
	d1->top_linkage = n+1;
	d2->top_linkage = n;
	d1->bottom_y = d2->bottom_y = -1;
	d1->bottom_linkage = -1;
	d2->bottom_linkage = -1;
	d1->pts = (int *)malloc_or_die(sizeof(int) * max_pts);
	d2->pts = (int *)malloc_or_die(sizeof(int) * max_pts);
	(*num_descenders) += 2;
	if(VERBOSE) fprintf(stderr, "num_descenders = %d (y=%d)\n", *num_descenders, y);
	return n;
}

mpoly_t calc_ring_from_mask(unsigned char *mask, int w, int h,
report_image_t *dbuf, int major_ring_only, int no_donuts, double min_ring_area) {
	int x, y;

	fprintf(stderr, "finding rings: begin\n");

	rowstat_t up_row, down_row;
	down_row.num_transitions = 0; // prevent compiler warning;

	up_row.openings        = (int *)malloc_or_die(sizeof(int) * (w+1)/2);
	up_row.closings        = (int *)malloc_or_die(sizeof(int) * (w+1)/2);
	up_row.descender_ids   = (int *)malloc_or_die(sizeof(int) * (w+1));
	down_row.openings      = (int *)malloc_or_die(sizeof(int) * (w+1)/2);
	down_row.closings      = (int *)malloc_or_die(sizeof(int) * (w+1)/2);
	down_row.descender_ids = (int *)malloc_or_die(sizeof(int) * (w+1));

	int num_descenders = 0;
	descender_t *descenders = NULL;

	for(y=0; y<=h; y++) {
		if(y) {
			up_row.num_transitions = down_row.num_transitions;
			memcpy(up_row.openings, down_row.openings, sizeof(int)*down_row.num_transitions);
			memcpy(up_row.closings, down_row.closings, sizeof(int)*down_row.num_transitions);
			memcpy(up_row.descender_ids, down_row.descender_ids, sizeof(int)*down_row.num_transitions*2);
		} else {
			up_row.num_transitions = 0;
		}

		if(y==h) {
			down_row.num_transitions = 0;
		} else {
			down_row.num_transitions = 0;

			int mask_rowlen = (w+7)/8;
			unsigned char mask_bitp = 1;
			unsigned char *mask_bytep = mask + mask_rowlen*y;

			int mask_left = 0;
			for(x=0; x<=w; x++) {
				int mask_right;
				if(x == w) {
					mask_right = 0;
				} else {
					mask_right = (*mask_bytep & mask_bitp) ? 1 : 0;
					mask_bitp <<= 1;
					if(!mask_bitp) {
						mask_bitp = 1;
						mask_bytep++;
					}
				}
				if(mask_right && !mask_left) {
					down_row.openings[down_row.num_transitions] = x;
					//if(VERBOSE) fprintf(stderr, "[%d,", x);
				}
				if(mask_left && !mask_right) {
					down_row.closings[down_row.num_transitions++] = x;
					//if(VERBOSE) fprintf(stderr, "%d] ", x);
				}
				mask_left = mask_right;
			}
			//if(VERBOSE) fprintf(stderr, "\n");
		}

		int up_tid=0, down_tid=0;
		while(up_tid < up_row.num_transitions || down_tid < down_row.num_transitions) {
			//if(VERBOSE) fprintf(stderr, "%d/%d:[%d,%d]  %d/%d:[%d,%d]\n",
			//	up_tid, up_row.num_transitions, up_row.openings[up_tid], up_row.closings[up_tid],
			//	down_tid, down_row.num_transitions, down_row.openings[down_tid], down_row.closings[down_tid]);
			if(
				(down_tid == down_row.num_transitions) ||
				(up_tid < up_row.num_transitions && up_row.closings[up_tid] <= down_row.openings[down_tid])
			) {
				//      \----/ 
				if(dbuf && dbuf->mode == PLOT_DESCENDERS) {
					for(x=up_row.openings[up_tid]; x<=up_row.closings[up_tid]; x++)
						plot_point(dbuf, x, y, 255, 0, 0);
				}
				int d1 = up_row.descender_ids[up_tid*2];
				int d2 = up_row.descender_ids[up_tid*2+1];
				descenders[d1].bottom_linkage = d2;
				descenders[d1].bottom_y = y-1;
				descenders[d1].pts = (int *)realloc_or_die(descenders[d1].pts,
					sizeof(int) * (descenders[d1].bottom_y - descenders[d1].top_y + 1));
				descenders[d2].bottom_linkage = d1;
				descenders[d2].bottom_y = y-1;
				descenders[d2].pts = (int *)realloc_or_die(descenders[d2].pts,
					sizeof(int) * (descenders[d2].bottom_y - descenders[d2].top_y + 1));
				up_tid++;
			} else if(
				(up_tid == up_row.num_transitions) ||
				(down_tid < down_row.num_transitions && down_row.closings[down_tid] <= up_row.openings[up_tid])
			) {
				//      /----\  .
				if(dbuf && dbuf->mode == PLOT_DESCENDERS) {
					for(x=down_row.openings[down_tid]; x<=down_row.closings[down_tid]; x++)
						plot_point(dbuf, x, y, 0, 255, 0);
				}
				int d = create_descender_pair(&num_descenders, &descenders, y, h-y+1);
				descenders[d  ].pts[0] = down_row.openings[down_tid];
				descenders[d+1].pts[0] = down_row.closings[down_tid];
				down_row.descender_ids[down_tid*2  ] = d;
				down_row.descender_ids[down_tid*2+1] = d+1;
				down_tid++;
			} else if(up_tid < up_row.num_transitions && down_tid < down_row.num_transitions) {
				if(dbuf && dbuf->mode == PLOT_DESCENDERS) {
					plot_point(dbuf, up_row.openings[up_tid], y, 255, 255, 0);
					plot_point(dbuf, down_row.openings[down_tid], y, 255, 255, 0);
				}
				int dl = up_row.descender_ids[up_tid*2];
//if(VERBOSE) fprintf(stderr, "dl=%d\n", dl);
				down_row.descender_ids[down_tid*2] = dl;
				descenders[dl].pts[y-descenders[dl].top_y] = down_row.openings[down_tid];
				for(;;) {
					if(
						(up_tid < up_row.num_transitions-1) &&
						(up_row.openings[up_tid+1] < down_row.closings[down_tid])
					) {
						if(dbuf && dbuf->mode == PLOT_DESCENDERS) {
							for(x=up_row.closings[up_tid]; x<=up_row.openings[up_tid+1]; x++)
								plot_point(dbuf, x, y, 255, 0, 255);
						}
						int d1 = up_row.descender_ids[up_tid*2+1];
						int d2 = up_row.descender_ids[up_tid*2+2];
						descenders[d1].bottom_linkage = d2;
						descenders[d1].bottom_y = y-1;
						descenders[d1].pts = (int *)realloc_or_die(descenders[d1].pts,
							sizeof(int) * (descenders[d1].bottom_y - descenders[d1].top_y + 1));
						descenders[d2].bottom_linkage = d1;
						descenders[d2].bottom_y = y-1;
						descenders[d2].pts = (int *)realloc_or_die(descenders[d2].pts,
							sizeof(int) * (descenders[d2].bottom_y - descenders[d2].top_y + 1));
						up_tid++;
					} else if(
						(down_tid < down_row.num_transitions-1) &&
						(down_row.openings[down_tid+1] < up_row.closings[up_tid])
					) {
						if(dbuf && dbuf->mode == PLOT_DESCENDERS) {
							for(x=down_row.closings[down_tid]; x<=down_row.openings[down_tid+1]; x++)
								plot_point(dbuf, x, y, 0, 255, 255);
						}
						int d = create_descender_pair(&num_descenders, &descenders, y, h-y+1);
						descenders[d  ].pts[0] = down_row.closings[down_tid];
						descenders[d+1].pts[0] = down_row.openings[down_tid+1];
						down_row.descender_ids[down_tid*2+1] = d;
						down_row.descender_ids[down_tid*2+2] = d+1;
						down_tid++;
					} else break;
				}
				if(dbuf && dbuf->mode == PLOT_DESCENDERS) {
					plot_point(dbuf, up_row.closings[up_tid], y, 255, 255, 0);
					plot_point(dbuf, down_row.closings[down_tid], y, 255, 255, 0);
				}
				int dr = up_row.descender_ids[up_tid*2+1];
				down_row.descender_ids[down_tid*2+1] = dr;
				descenders[dr].pts[y-descenders[dr].top_y] = down_row.closings[down_tid];
				up_tid++;
				down_tid++;
			}
		}
	}

	int i;

	if(VERBOSE) {
		for(i=0; i<num_descenders; i++) fprintf(stderr, "%d: top_link=%d bottom_link=%d\n",
			i, descenders[i].top_linkage, descenders[i].bottom_linkage);
	}

	if(!num_descenders) {
		fprintf(stderr, "image was completely blank - therefore there is no bounding polygon\n");
		return (mpoly_t){ 0, NULL };
	}

	unsigned char *used_desc = (unsigned char *)malloc_or_die(num_descenders);
	for(i=0; i<num_descenders; i++) used_desc[i] = 0;

	mpoly_t mp = (mpoly_t){ 0, NULL };

	for(;;) {
		int start_d = -1;
		for(i=0; i<num_descenders; i++) {
			if(!used_desc[i]) {
				start_d = i;
				break;
			}
		}
		if(start_d < 0) break;
		ring_t ring;
		ring.npts = 0;
		ring.pts = NULL;
		// prevent compiler warning - these fields are filled in by compute_containments
		ring.parent_id = ring.is_hole = 0;

		int cur_d = start_d;
		do {
			if(VERBOSE) fprintf(stderr, "d:%d ", cur_d);
			descender_t *d = descenders + cur_d;
			if(used_desc[cur_d]) fatal_error("descender used twice");
			used_desc[cur_d]++;
			int n = d->bottom_y - d->top_y + 1;
			if(n <= 0) fatal_error("npts <= 0 in ring segment");
			ring.pts = (vertex_t *)realloc_or_die(ring.pts, 
				sizeof(vertex_t)*(ring.npts+n*2));
			ring.pts[ring.npts++] = (vertex_t){ d->pts[0], d->top_y };
			for(i=1; i<n; i++) {
				if(d->pts[i] != d->pts[i-1]) {
					ring.pts[ring.npts++] = 
						(vertex_t){ d->pts[i-1], d->top_y + i };
					ring.pts[ring.npts++] =
						(vertex_t){ d->pts[i  ], d->top_y + i };
				}
			}
			ring.pts[ring.npts++] = (vertex_t){ d->pts[n-1], d->bottom_y+1 };
			cur_d = d->bottom_linkage;
			if(cur_d < 0) fatal_error("uninitialized val in descender linkage");

			if(VERBOSE) fprintf(stderr, "u:%d ", cur_d);
			d = descenders + cur_d;
			if(used_desc[cur_d]) fatal_error("descender used twice");
			used_desc[cur_d]++;
			n = d->bottom_y - d->top_y + 1;
			if(n <= 0) fatal_error("npts <= 0 in ring segment");
			ring.pts = (vertex_t *)realloc_or_die(ring.pts, 
				sizeof(vertex_t)*(ring.npts+n*2));
			ring.pts[ring.npts++] = (vertex_t){ d->pts[n-1], d->bottom_y + 1 };
			for(i=n-2; i>=0; i--) {
				if(d->pts[i] != d->pts[i+1]) {
					ring.pts[ring.npts++] = 
						(vertex_t){ d->pts[i+1], d->top_y + i+1 };
					ring.pts[ring.npts++] = 
						(vertex_t){ d->pts[i  ], d->top_y + i+1 };
				}
			}
			ring.pts[ring.npts++] = (vertex_t){ d->pts[0], d->top_y };
			cur_d = d->top_linkage;
			if(cur_d < 0) fatal_error("uninitialized val in descender linkage");
		} while(cur_d != start_d);
		if(VERBOSE) fprintf(stderr, "\n");

		if(VERBOSE) fprintf(stderr, "ring %d: %d pts\n", mp.num_rings, ring.npts);

		mp.rings = (ring_t *)realloc_or_die(mp.rings, sizeof(ring_t)*(mp.num_rings+1));
		mp.rings[mp.num_rings++] = ring;
	}
	fprintf(stderr, "finding rings: end\n");

	// the topology cannot be resolved by us or by geos/jump/postgis if
	// there are self-intersections
	bevel_self_intersections(&mp);

	if(min_ring_area > 0) {
		if(VERBOSE) fprintf(stderr, "removing small rings...\n");

		ring_t *filtered_rings = (ring_t *)malloc_or_die(sizeof(ring_t)*mp.num_rings);
		int num_filtered_rings = 0;
		for(i=0; i<mp.num_rings; i++) {
			double area = polygon_area(mp.rings+i);
			if(VERBOSE) if(area > 10) fprintf(stderr, "ring %d has area %.15f\n", i, area);
			if(area >= min_ring_area) {
				filtered_rings[num_filtered_rings++] = mp.rings[i];
			}
		}
		fprintf(stderr, "filtered by area %d => %d rings\n",
			mp.num_rings, num_filtered_rings);

		mp.rings = filtered_rings;
		mp.num_rings = num_filtered_rings;
	}

	if(major_ring_only) {
		double biggest_area = 0;
		int best_idx = 0;
		for(i=0; i<mp.num_rings; i++) {
			double area = polygon_area(mp.rings+i);
			if(area > biggest_area) {
				biggest_area = area;
				best_idx = i;
			}
		}
		if(VERBOSE) fprintf(stderr, "major ring was %d with %d pts, %.1f area\n",
			best_idx, mp.rings[best_idx].npts, biggest_area);
		mp.rings = mp.rings+best_idx;
		mp.num_rings = 1;
	}

	//fprintf(stderr, "computing containments: begin\n");
	compute_containments(&mp);
	//fprintf(stderr, "computing containments: end\n");

	if(no_donuts) {
		int out_idx = 0;
		for(i=0; i<mp.num_rings; i++) {
			if(mp.rings[i].parent_id < 0) {
				mp.rings[out_idx++] = mp.rings[i];
			}
		}
		mp.num_rings = out_idx;
		if(VERBOSE) fprintf(stderr, "number of non-donut rings is %d", mp.num_rings);
	}

	return mp;
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

unsigned char *get_mask_for_dataset(GDALDatasetH ds, int bandlist_size, int *bandlist, 
int num_ndv, double *ndv_list, double ndv_tolerance, report_image_t *dbuf) {
	int i, j;

	int w = GDALGetRasterXSize(ds);
	int h = GDALGetRasterYSize(ds);
	int band_count = GDALGetRasterCount(ds);
	if(VERBOSE) fprintf(stderr, "input is %d x %d x %d\n", w, h, band_count);

	int mask_rowlen = (w+7)/8;
	if(VERBOSE) fprintf(stderr, "mask array is %.1f megabytes\n", (double)mask_rowlen*h/1024.0/1024.0);
	unsigned char *mask = (unsigned char *)malloc_or_die(mask_rowlen*h);
	for(i=0; i<mask_rowlen*h; i++) mask[i] = 0;

	int bandlist_idx;
	int last_progress = 0;
	for(bandlist_idx=0; bandlist_idx<bandlist_size; bandlist_idx++) {
		int band_idx = bandlist[bandlist_idx];
		if(band_idx < 1 || band_idx > band_count) fatal_error("bandid out of range");

		GDALRasterBandH band = GDALGetRasterBand(ds, band_idx);

		int blocksize_x, blocksize_y;
		GDALGetBlockSize(band, &blocksize_x, &blocksize_y);

		double nodataval = ndv_list[bandlist_idx];

		if(VERBOSE) fprintf(stderr, "band %d: block size = %d,%d  nodataval = %.15f\n",
			band_idx, blocksize_x, blocksize_y, nodataval);

		double *buf = (double *)malloc_or_die(blocksize_x*blocksize_y*sizeof(double));
		int boff_x, boff_y;
		for(boff_y=0; boff_y<h; boff_y+=blocksize_y) {
			int bsize_y = blocksize_y;
			if(bsize_y + boff_y > h) bsize_y = h - boff_y;
			for(boff_x=0; boff_x<w; boff_x+=blocksize_x) {
				int bsize_x = blocksize_x;
				if(bsize_x + boff_x > w) bsize_x = w - boff_x;

				GDALRasterIO(band, GF_Read, boff_x, boff_y, bsize_x, bsize_y, 
					buf, bsize_x, bsize_y, GDT_Float64, 0, 0);

				/*
				if(!boff_x && !boff_y) {
					if(VERBOSE) fprintf(stderr, "band %d: pixel[0] = %.15f\n", band_idx, buf[0]);
					nodataval = buf[0];
				}
				*/

				double *p = buf;
				for(j=0; j<bsize_y; j++) {
					int y = j + boff_y;
					int is_dbuf_stride_y = dbuf && ((y % dbuf->stride_y) == 0);
					unsigned char mask_bitp = 1 << (boff_x % 8);
					unsigned char *mask_bytep = mask + mask_rowlen*y + boff_x/8;
					for(i=0; i<bsize_x; i++) {
						int is_dbuf_stride = is_dbuf_stride_y && ((i % dbuf->stride_x) == 0);
						double val = *(p++);
						if(fabs(val - nodataval) > ndv_tolerance) {
							*mask_bytep |= mask_bitp;

							if(is_dbuf_stride) {
								int x = i + boff_x;
								unsigned char db_v = 100 + (unsigned char)(val/2);
								if(db_v < 100) db_v = 100;
								if(db_v > 254) db_v = 254;
								unsigned char r = (unsigned char)(db_v*.75);
								plot_point(dbuf, x, y, r, db_v, db_v);
							}
						}
						mask_bitp <<= 1;
						if(!mask_bitp) {
							mask_bitp = 1;
							mask_bytep++;
						}
					}
				}

				int progress = (int)(100L * (
					(long)(bandlist_idx) * (long)w * (long)h +
					(long)boff_y * (long)w +
					(long)(boff_x+bsize_x) * (long)bsize_y) /
					((long)bandlist_size * (long)w * (long)h));
				if(progress != last_progress) {
					fprintf(stderr, "reading: %d%%\r", progress);
					fflush(stderr);
					last_progress = progress;
				}
			}
		}

		free(buf);
	}
	fprintf(stderr, "\n");

	return mask;
}

unsigned char *erode_mask(unsigned char *in_mask, int w, int h) {
	int i, j;
	int mask_rowlen = (w+7)/8;

	unsigned char *out_mask = (unsigned char *)malloc_or_die(mask_rowlen*h);
	for(i=0; i<mask_rowlen*h; i++) out_mask[i] = 0;

	for(j=0; j<h; j++) {
		unsigned char *in_bytep = in_mask + mask_rowlen*j;
		unsigned char *out_bytep = out_mask + mask_rowlen*j;
		unsigned char ul = 0, um = j ? *(in_bytep-mask_rowlen) & 1 : 0;
		unsigned char ml = 0, mm = *in_bytep & 1;
		unsigned char ll = 0, lm = (j<h-1) ? *(in_bytep+mask_rowlen) & 1 : 0;
		unsigned char in_bitp = 2;
		unsigned char out_bitp = 1;
		for(i=0; i<w; i++) {
			unsigned char ur = (j && i<w-1) ? *(in_bytep-mask_rowlen) & in_bitp : 0;
			unsigned char mr = (i<w-1) ? *in_bytep & in_bitp : 0;
			unsigned char lr = (j<h-1 && i<w-1) ? *(in_bytep+mask_rowlen) & in_bitp : 0;

			// remove pixels that don't have two consecutive filled neighbors
			if(mm && (
				ul&&um || um&&ur || ur&&mr || mr&&lr ||
				lr&&lm || lm&&ll || ll&&ml || ml&&ul
			)) {
				*out_bytep |= out_bitp;
			}
			// fill pixels that don't have two consecutive empty neighbors
			/*
			if(!mm && (
				(ul||um) && (um||ur) && (ur||mr) && (mr||lr) &&
				(lr||lm) && (lm||ll) && (ll||ml) && (ml||ul)
			)) {
				*out_bytep |= out_bitp;
			}
			*/
			
			ul=um; ml=mm; ll=lm;
			um=ur; mm=mr; lm=lr;

			in_bitp <<= 1;
			if(!in_bitp) {
				in_bitp = 1;
				in_bytep++;
			}

			out_bitp <<= 1;
			if(!out_bitp) {
				out_bitp = 1;
				out_bytep++;
			}
		}
	}

	return out_mask;
}
