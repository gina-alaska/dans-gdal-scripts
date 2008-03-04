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
#include "mask.h"

#include <ogrsf_frmts.h>

#define CS_UNKNOWN 0
#define CS_XY 1
#define CS_EN 2
#define CS_LL 3

void usage(char *cmdname) {
	printf("Usage:\n  %s [options] [image_name]\n", cmdname);
	printf("\n");
	
	print_georef_usage();

	printf("\
\n\
Behavior:\n\
  -classify                    Output a polygon for each value of an 8-bit band\n\
                               (default is to generate a single polygon that\n\
                               surrounds all pixels that don't match\n\
                               the no-data-value)\n\
  -nodataval 'val [val ...]'   Specify value of no-data pixels\n\
  -ndv-toler val               Tolerance for deciding if a pixel\n\
                               matches nodataval\n\
  -b band_id -b band_id ...    Bands to inspect (default is all bands)\n\
  -invert                      Trace no-data pixels rather than data pixels\n\
  -erosion                     Erode pixels that don't have two consecutive\n\
                               neighbors\n\
  -major-ring                  Take only the biggest outer ring\n\
  -no-donuts                   Take only top-level rings\n\
  -min-ring-area val           Drop rings with less than this area\n\
                               (in square pixels)\n\
  -dp-toler val                Tolerance for point reduction\n\
                               (in pixels, default is 2.0)\n\
  -bevel-size                  How much to shave off corners at\n\
                               self-intersection points\n\
                               (in pixels, default it 0.1)\n\
                               (this is done to make geometries that\n\
                               PostGIS/GEOS/Jump can handle)\n\
\n\
Output:\n\
  -report fn.ppm               Output graphical report of bounds found\n\
  -mask-out fn.pbm             Output mask of bounding polygon in PBM format\n\
  -out-cs [xy | en | ll]       Coordinate system for output\n\
                               (pixel coords, easting/northing, or lon/lat)\n\
  -llproj-toler val            Error tolerance for curved lines when\n\
                               using '-out-cs ll'\n\
                               (in pixels, default is 0.2)\n\
  -wkt-out fn.wkt              Output bounds in WKT format\n\
  -ogr-out fn.shp              Output bounds using an OGR format\n\
  -ogr-fmt                     OGR format to use (default is 'ESRI Shapefile')\n\
  -split-polys                 Output several polygons rather than one\n\
                               multipolygon\n\
\n\
Misc:\n\
  -v                           Verbose\n\
\n\
Examples:\n\
\n\
Inspect image and output contour of data region:\n\
gdal_trace_outline raster.tif -nodataval 0 -erosion -out-cs ll -wkt-out outline.wkt\n\
\n\
Same as above but polygon actually follows border pixel-by-pixel:\n\
gdal_trace_outline raster.tif -nodataval 0 -dp-toler 0 -out-cs ll -wkt-out outline.wkt\n\
\n\
Output ESRI Shapefile in projection of input image:\n\
gdal_trace_outline raster.tif -nodataval 0 -erosion -out-cs en -ogr-out outline.shp\n\
\n\
Generate one shape for each value in input image:\n\
gdal_trace_outline raster.tif -classify -out-cs en -ogr-out outline.shp\n\
\n\
");
	exit(1);
}

mpoly_t calc_ring_from_mask(unsigned char *mask, int w, int h,
	report_image_t *dbuf, int major_ring_only, int no_donuts,
	double min_ring_area, double bevel_size);

int main(int argc, char **argv) {
	char *input_raster_fn = NULL;
	int classify = 0;
	int num_ndv = 0;
	double *ndv_list = NULL;
	double ndv_tolerance = 0;
	char *debug_report = NULL;
	int inspect_numbands = 0;
	int *inspect_bandids = NULL;
	int split_polys = 0;
	int out_cs = CS_UNKNOWN;
	char *wkt_fn = NULL;
	char *ogr_fn = NULL;
	char *ogr_fmt = NULL;
	char *mask_out_fn = NULL;
	int major_ring_only = 0;
	int no_donuts = 0;
	double min_ring_area = 0;
	double reduction_tolerance = 2;
	int do_erosion = 0;
	int do_invert = 0;
	double llproj_toler = .2;
	double bevel_size = .1;

	if(argc == 1) usage(argv[0]);

	geo_opts_t geo_opts = init_geo_options(&argc, &argv);

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(!strcmp(arg, "-v")) {
				VERBOSE++;
			} else if(!strcmp(arg, "-classify")) {
				classify = 1;
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
			} else if(!strcmp(arg, "-erosion")) {
				do_erosion = 1;
			} else if(!strcmp(arg, "-invert")) {
				do_invert = 1;
			} else if(!strcmp(arg, "-split-polys")) {
				split_polys++;
			} else if(!strcmp(arg, "-wkt-out")) {
				if(argp == argc) usage(argv[0]);
				wkt_fn = argv[argp++];
			} else if(!strcmp(arg, "-ogr-out")) {
				if(argp == argc) usage(argv[0]);
				ogr_fn = argv[argp++];
			} else if(!strcmp(arg, "-ogr-fmt")) {
				if(argp == argc) usage(argv[0]);
				ogr_fmt = argv[argp++];
			} else if(!strcmp(arg, "-out-cs")) {
				if(argp == argc) usage(argv[0]);
				char *cs = argv[argp++];
				if(!strcmp(cs, "xy")) out_cs = CS_XY;
				else if(!strcmp(cs, "en")) out_cs = CS_EN;
				else if(!strcmp(cs, "ll")) out_cs = CS_LL;
				else fatal_error("unrecognized value for -out-cs option (%s)", cs);
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
			} else if(!strcmp(arg, "-bevel-size")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				bevel_size = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
				if(bevel_size < 0 || bevel_size >= 1) fatal_error(
					"-bevel-size must be in the range 0 <= bevel < 1");
			} else if(!strcmp(arg, "-llproj-toler")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				llproj_toler = strtod(argv[argp++], &endptr);
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

	if(!input_raster_fn) fatal_error("must specify filename of image");

	int do_geom_output = wkt_fn || ogr_fn;
	if(do_geom_output && (out_cs == CS_UNKNOWN)) 
		fatal_error("must specify output coordinate system with -out-cs option");

	if(major_ring_only && min_ring_area) fatal_error(
		"-major-ring and -min-ring-area options cannot both be used at the same time");
	if(major_ring_only && no_donuts) fatal_error(
		"-major-ring and -no-donuts options cannot both be used at the same time");

	if(classify) {
		if(num_ndv) fatal_error("-classify option is not compatible with -nodataval option");
		if(do_invert) fatal_error("-classify option is not compatible with -invert option");
		if(mask_out_fn) fatal_error("-classify option is not compatible with -mask-out option");
	}

	GDALAllRegister();

	GDALDatasetH ds = GDALOpen(input_raster_fn, GA_ReadOnly);
	if(!ds) fatal_error("open failed");

	if(!inspect_numbands) {
		inspect_numbands = classify ? 1 : GDALGetRasterCount(ds);
		inspect_bandids = (int *)malloc_or_die(sizeof(int)*inspect_numbands);
		int i;
		for(i=0; i<inspect_numbands; i++) inspect_bandids[i] = i+1;
	}

	// FIXME - optional NDV for classify
	if(!classify) {
		setup_ndv_list(ds, inspect_numbands, inspect_bandids, &num_ndv, &ndv_list);
	}

	CPLPushErrorHandler(CPLQuietErrorHandler);

	georef_t georef = init_georef(&geo_opts, ds);

	if((out_cs == CS_EN || out_cs == CS_LL) && !georef.fwd_affine) fatal_error("missing affine transform");
	if((out_cs == CS_LL) && !georef.fwd_xform) fatal_error("missing coordinate transform");

	report_image_t *dbuf = NULL;
	if(debug_report) {
		dbuf = create_plot(georef.w, georef.h);
		dbuf->mode = PLOT_CONTOURS;
	}

	unsigned char *raster = NULL;
	unsigned char *mask = NULL;
	unsigned char usage_array[256];
	GDALColorTableH color_table = NULL;
	if(classify) {
		if(inspect_numbands != 1) {
			fatal_error("only one band may be used in classify mode");
		}

		raster = read_dataset_8bit(ds, inspect_bandids[0], usage_array, dbuf);

		GDALRasterBandH band = GDALGetRasterBand(ds, inspect_bandids[0]);
		if(GDALGetRasterColorInterpretation(band) == GCI_PaletteIndex) {
			color_table = GDALGetRasterColorTable(band);
		}
	} else {
		mask = get_mask_for_dataset(ds, inspect_numbands, inspect_bandids,
			num_ndv, ndv_list, ndv_tolerance, dbuf);
	}

	FILE *wkt_fh = NULL;
	OGRDataSourceH ogr_ds = NULL;
	OGRLayerH ogr_layer = NULL;
	int class_fld_idx = -1;
	int color_fld_idx[4] = { -1, -1, -1, -1 };
	if(do_geom_output) {
		if(wkt_fn) {
			wkt_fh = fopen(wkt_fn, "w");
			if(!wkt_fh) fatal_error("cannot open output file for WKT");
		}

		if(ogr_fn) {
			OGRRegisterAll();
			if(!ogr_fmt) ogr_fmt = "ESRI Shapefile";
			OGRSFDriverH ogr_driver = OGRGetDriverByName(ogr_fmt);
			if(!ogr_driver) fatal_error("cannot get OGR driver (%s)", ogr_fmt);
			ogr_ds = OGR_Dr_CreateDataSource(ogr_driver, ogr_fn, NULL);
			if(!ogr_ds) fatal_error("cannot create OGR data source");

			char *layer_name = ogr_fn;

			OGRSpatialReferenceH sref = NULL;
			if(out_cs == CS_EN) {
				sref = georef.spatial_ref;
			} else if(out_cs == CS_LL) {
				sref = OSRNewSpatialReference(NULL);
				OSRImportFromProj4(sref, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs");
			}

			ogr_layer = OGR_DS_CreateLayer(ogr_ds, layer_name, sref, 
				(split_polys ? wkbPolygon : wkbMultiPolygon), NULL);
			if(!ogr_layer) fatal_error("cannot create OGR layer");

			if(classify) {
				OGRFieldDefnH fld = OGR_Fld_Create("value", OFTInteger);
				OGR_Fld_SetWidth(fld, 4);
				OGR_L_CreateField(ogr_layer, fld, TRUE);
				class_fld_idx = OGR_FD_GetFieldIndex(OGR_L_GetLayerDefn(ogr_layer), "value");

				if(color_table) {
					char *names[4] = { "c1", "c2", "c3", "c4" };
					int i;
					for(i=0; i<4; i++) {
						fld = OGR_Fld_Create(names[i], OFTInteger);
						OGR_Fld_SetWidth(fld, 4);
						OGR_L_CreateField(ogr_layer, fld, TRUE);
						color_fld_idx[i] = OGR_FD_GetFieldIndex(OGR_L_GetLayerDefn(ogr_layer), names[i]);
					}
				}
			}
		}
	}

	int num_shapes_written = 0;

	int class_id;
	for(class_id=0; class_id<256; class_id++) {
		const GDALColorEntry *color = NULL;
		if(classify) {
			if(!usage_array[class_id]) continue;
			printf("\nTracing feature class %d\n", class_id);

			if(color_table) {
				color = GDALGetColorEntry(color_table, class_id);
				if(color) printf("  Color=%d,%d,%d,%d\n",
					color->c1, color->c2, color->c3, color->c4);
			}

			mask = get_mask_for_8bit_raster(georef.w, georef.h,
				raster, (unsigned char)class_id);
		} else {
			if(class_id != 0) continue;
			printf("\nTracing data pixels\n");
		}

		if(do_invert) {
			invert_mask(mask, georef.w, georef.h);
		}

		if(do_erosion) {
			unsigned char *eroded_mask = erode_mask(mask, georef.w, georef.h);
			free(mask);
			mask = eroded_mask;
		}

		mpoly_t *bounds_poly = (mpoly_t *)malloc_or_die(sizeof(mpoly_t));
		*bounds_poly = calc_ring_from_mask(mask, georef.w, georef.h, dbuf, 
			major_ring_only, no_donuts, min_ring_area, bevel_size);
		free(mask);

		if(mask_out_fn) {
			mask_from_mpoly(bounds_poly, georef.w, georef.h, mask_out_fn);
		}

		if(reduction_tolerance > 0) {
			mpoly_t reduced_poly = compute_reduced_pointset(bounds_poly, reduction_tolerance);
			free_mpoly(bounds_poly);
			*bounds_poly = reduced_poly;
		}

		int num_outer=0, num_inner=0, total_pts=0;
		int r_idx;
		for(r_idx=0; r_idx<bounds_poly->num_rings; r_idx++) {
			if(bounds_poly->rings[r_idx].is_hole) num_inner++;
			else num_outer++;
			total_pts += bounds_poly->rings[r_idx].npts;
		}
		printf("Found %d outer rings and %d holes with a total of %d vertices.\n",
			num_outer, num_inner, total_pts);

		if(dbuf && dbuf->mode == PLOT_CONTOURS) {
			debug_plot_rings(bounds_poly, dbuf);
		}

		if(do_geom_output && bounds_poly->num_rings) {
			int num_shapes;
			mpoly_t *shapes;
			int shape_is_copy;
			if(split_polys) {
				split_mpoly_to_polys(bounds_poly, &num_shapes, &shapes);
				shape_is_copy = 0;
			} else {
				num_shapes = 1;
				shapes = (mpoly_t *)malloc_or_die(sizeof(mpoly_t));
				shapes[0] = *bounds_poly;
				shape_is_copy = 1;
			}

			int shape_idx;
			for(shape_idx=0; shape_idx<num_shapes; shape_idx++) {
				mpoly_t *poly_in = shapes+shape_idx;

				mpoly_t *proj_poly;
				int proj_is_copy;
				if(out_cs == CS_XY) {
					proj_poly = poly_in;
					proj_is_copy = 1;
				} else if(out_cs == CS_EN) {
					proj_poly = mpoly_xy2en(&georef, poly_in);
					proj_is_copy = 0;
				} else if(out_cs == CS_LL) {
					proj_poly = mpoly_xy2ll_with_interp(&georef, poly_in, llproj_toler);
					proj_is_copy = 0;
				} else {
					fatal_error("bad val for out_cs");
				}

				OGRGeometryH ogr_geom = mpoly_to_ogr(proj_poly);

				if(wkt_fh) {
					char *wkt_out;
					OGR_G_ExportToWkt(ogr_geom, &wkt_out);
					fprintf(wkt_fh, "%s\n", wkt_out);
				}

				if(ogr_ds) {
					OGRFeatureH ogr_feat = OGR_F_Create(OGR_L_GetLayerDefn(ogr_layer));
					if(class_fld_idx >= 0) OGR_F_SetFieldInteger(ogr_feat, class_fld_idx, class_id);
					if(color) {
						if(color_fld_idx[0] >= 0) OGR_F_SetFieldInteger(ogr_feat, color_fld_idx[0], color->c1);
						if(color_fld_idx[1] >= 0) OGR_F_SetFieldInteger(ogr_feat, color_fld_idx[1], color->c2);
						if(color_fld_idx[2] >= 0) OGR_F_SetFieldInteger(ogr_feat, color_fld_idx[2], color->c3);
						if(color_fld_idx[3] >= 0) OGR_F_SetFieldInteger(ogr_feat, color_fld_idx[3], color->c4);
					}
					OGR_F_SetGeometryDirectly(ogr_feat, ogr_geom); // assumes ownership of geom
					OGR_L_CreateFeature(ogr_layer, ogr_feat);
					OGR_F_Destroy(ogr_feat);
				} else {
					OGR_G_DestroyGeometry(ogr_geom);
				}

				if(!proj_is_copy) free_mpoly(proj_poly);
				if(!shape_is_copy) free_mpoly(poly_in);

				num_shapes_written++;
			}
		}

		free_mpoly(bounds_poly);
	}

	printf("\n");

	if(wkt_fh) fclose(wkt_fh);
	if(ogr_ds) OGR_DS_Destroy(ogr_ds);
	if(dbuf) write_plot(dbuf, debug_report);

	if(do_geom_output) {
		if(num_shapes_written) printf("Wrote %d shapes.\n", num_shapes_written);
		else printf("Wrote empty shapefile.\n");
	}

	CPLPopErrorHandler();

	return 0;
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
	if(VERBOSE >= 2) printf("num_descenders = %d (y=%d)\n", *num_descenders, y);
	return n;
}

mpoly_t calc_ring_from_mask(unsigned char *mask, int w, int h,
report_image_t *dbuf, int major_ring_only, int no_donuts, 
double min_ring_area, double bevel_size) {
	int x, y;

	if(VERBOSE) printf("finding rings: begin\n");

	// up_row is previous row, down_row is current row
	rowstat_t up_row, down_row;
	down_row.num_transitions = 0; // prevent compiler warning;

	// Openings and closings come in pairs having the same index.
	// An opening/closing pair represents entering and then leaving
	// a region of marked pixels in a given row.
	up_row.openings        = (int *)malloc_or_die(sizeof(int) * (w+1)/2);
	up_row.closings        = (int *)malloc_or_die(sizeof(int) * (w+1)/2);
	// Descenders trace the border between marked/unmarked pixels.
	up_row.descender_ids   = (int *)malloc_or_die(sizeof(int) * (w+1));
	down_row.openings      = (int *)malloc_or_die(sizeof(int) * (w+1)/2);
	down_row.closings      = (int *)malloc_or_die(sizeof(int) * (w+1)/2);
	down_row.descender_ids = (int *)malloc_or_die(sizeof(int) * (w+1));

	int num_descenders = 0;
	descender_t *descenders = NULL;

	int show_progress = 1;
	if(show_progress) printf("Tracing: ");

	for(y=0; y<=h; y++) {
		if(show_progress) GDALTermProgress((double)y/(double)(h+1), NULL, NULL);
		if(y) {
			// the previous down_row becomes the new up_row
			up_row.num_transitions = down_row.num_transitions;
			memcpy(up_row.openings, down_row.openings, sizeof(int)*down_row.num_transitions);
			memcpy(up_row.closings, down_row.closings, sizeof(int)*down_row.num_transitions);
			memcpy(up_row.descender_ids, down_row.descender_ids, sizeof(int)*down_row.num_transitions*2);
		} else {
			up_row.num_transitions = 0;
		}

		if(y==h) {
			// This is below the last row of the image.
			// This area is considered to be all unmarked pixels.
			down_row.num_transitions = 0;
		} else {
			// Tabulate the opening/closing pairs for this row.
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
					//if(VERBOSE) printf("[%d,", x);
				}
				if(mask_left && !mask_right) {
					down_row.closings[down_row.num_transitions++] = x;
					//if(VERBOSE) printf("%d] ", x);
				}
				mask_left = mask_right;
			}
			//if(VERBOSE) printf("\n");
		}

		int up_tid=0, down_tid=0;
		while(up_tid < up_row.num_transitions || down_tid < down_row.num_transitions) {
			//if(VERBOSE) printf("%d/%d:[%d,%d]  %d/%d:[%d,%d]\n",
			//	up_tid, up_row.num_transitions, up_row.openings[up_tid], up_row.closings[up_tid],
			//	down_tid, down_row.num_transitions, down_row.openings[down_tid], down_row.closings[down_tid]);
			if(
				(down_tid == down_row.num_transitions) ||
				(up_tid < up_row.num_transitions && up_row.closings[up_tid] <= down_row.openings[down_tid])
			) {
				// Close off a pair of descenders:     \____/ 
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
				// Create a new pair of descenders:     /^^^^\  .
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
				// Move descender to new location of boundary for this row
				if(dbuf && dbuf->mode == PLOT_DESCENDERS) {
					plot_point(dbuf, up_row.openings[up_tid], y, 255, 255, 0);
					plot_point(dbuf, down_row.openings[down_tid], y, 255, 255, 0);
				}
				int dl = up_row.descender_ids[up_tid*2];
//if(VERBOSE) printf("dl=%d\n", dl);
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
	if(show_progress) GDALTermProgress(1, NULL, NULL);

	int i;

	if(VERBOSE >= 2) {
		for(i=0; i<num_descenders; i++) printf("%d: top_link=%d bottom_link=%d\n",
			i, descenders[i].top_linkage, descenders[i].bottom_linkage);
	}

	if(!num_descenders) {
		printf("Mask was completely blank - therefore there is no bounding polygon\n");
		return empty_polygon();
	}

	unsigned char *used_desc = (unsigned char *)malloc_or_die(num_descenders);
	for(i=0; i<num_descenders; i++) used_desc[i] = 0;
	int num_desc_used = 0;

	mpoly_t mp = empty_polygon();

	if(show_progress) printf("Joining segments: ");
	int start_d = 0;
	for(;;) {
		if(show_progress) GDALTermProgress((double)num_desc_used/(double)num_descenders, NULL, NULL);
		for(; start_d<num_descenders; start_d++) {
			if(!used_desc[start_d]) {
				break;
			}
		}
		if(start_d == num_descenders) break;
		ring_t ring;
		ring.npts = 0;
		ring.pts = NULL;
		// prevent compiler warning - these fields are filled in by compute_containments
		ring.parent_id = ring.is_hole = 0;

		int cur_d = start_d;
		do {
			if(VERBOSE >= 2) printf("d:%d ", cur_d);
			descender_t *d = descenders + cur_d;
			if(used_desc[cur_d]) fatal_error("descender used twice");
			used_desc[cur_d]++; num_desc_used++;
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

			if(VERBOSE >= 2) printf("u:%d ", cur_d);
			d = descenders + cur_d;
			if(used_desc[cur_d]) fatal_error("descender used twice");
			used_desc[cur_d]++; num_desc_used++;
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
		if(VERBOSE >= 2) printf("\n");

		if(VERBOSE >= 2) printf("ring %d: %d pts\n", mp.num_rings, ring.npts);

		mp.rings = (ring_t *)realloc_or_die(mp.rings, sizeof(ring_t)*(mp.num_rings+1));
		mp.rings[mp.num_rings++] = ring;
	}
	if(show_progress) GDALTermProgress(1, NULL, NULL);

	free(used_desc);
	for(i=0; i<num_descenders; i++) {
		free(descenders[i].pts);
	}
	free(descenders);

	if(VERBOSE) printf("finding rings: end\n");

	if(VERBOSE) {
		int total_pts = 0;
		int r_idx;
		for(r_idx=0; r_idx<mp.num_rings; r_idx++) {
			total_pts += mp.rings[r_idx].npts;
		}
		printf("traced produced %d rings with a total of %d points\n",
			mp.num_rings, total_pts);
	}

	if(bevel_size > 0) {
		// the topology cannot be resolved by us or by geos/jump/postgis if
		// there are self-intersections
		bevel_self_intersections(&mp, bevel_size);
	}

	if(min_ring_area > 0) {
		if(VERBOSE) printf("removing small rings...\n");

		ring_t *filtered_rings = (ring_t *)malloc_or_die(sizeof(ring_t)*mp.num_rings);
		int num_filtered_rings = 0;
		for(i=0; i<mp.num_rings; i++) {
			double area = polygon_area(mp.rings+i);
			if(VERBOSE) if(area > 10) printf("ring %d has area %.15f\n", i, area);
			if(area >= min_ring_area) {
				filtered_rings[num_filtered_rings++] = mp.rings[i];
			} else {
				free_ring(mp.rings + i);
			}
		}
		printf("filtered by area %d => %d rings\n",
			mp.num_rings, num_filtered_rings);

		free(mp.rings);
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
		if(VERBOSE) printf("major ring was %d with %d pts, %.1f area\n",
			best_idx, mp.rings[best_idx].npts, biggest_area);
		mp.rings = mp.rings+best_idx;
		mp.num_rings = 1;
	}

	//printf("computing containments: begin\n");
	compute_containments(&mp);
	//printf("computing containments: end\n");

	if(no_donuts) {
		int out_idx = 0;
		for(i=0; i<mp.num_rings; i++) {
			if(mp.rings[i].parent_id < 0) {
				mp.rings[out_idx++] = mp.rings[i];
			}
		}
		mp.num_rings = out_idx;
		if(VERBOSE) printf("number of non-donut rings is %d", mp.num_rings);
	}

	return mp;
}
