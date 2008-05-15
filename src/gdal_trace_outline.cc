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
#include "mask-tracer.h"
#include "dp.h"
#include "excursion_pincher.h"

#include <ogrsf_frmts.h>
#include <cpl_string.h>
#include <cpl_conv.h>
#include <cpl_port.h>

#ifdef CPL_MSB 
#define WKB_BYTE_ORDER wkbNDR
#else
#define WKB_BYTE_ORDER wkbXDR
#endif

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
  -wkb-out fn.wkb              Output bounds in WKB format\n\
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
	long min_ring_area, double bevel_size);

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
	char *wkb_fn = NULL;
	char *ogr_fn = NULL;
	char *ogr_fmt = NULL;
	char *mask_out_fn = NULL;
	int major_ring_only = 0;
	int no_donuts = 0;
	long min_ring_area = 0;
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
			} else if(!strcmp(arg, "-wkb-out")) {
				if(argp == argc) usage(argv[0]);
				wkb_fn = argv[argp++];
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
				min_ring_area = strtol(argv[argp++], &endptr, 10);
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

	int do_geom_output = wkb_fn || wkt_fn || ogr_fn;
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
	FILE *wkb_fh = NULL;
	OGRDataSourceH ogr_ds = NULL;
	OGRLayerH ogr_layer = NULL;
	int class_fld_idx = -1;
	int color_fld_idx[4] = { -1, -1, -1, -1 };
	if(do_geom_output) {
		if(wkt_fn) {
			wkt_fh = fopen(wkt_fn, "w");
			if(!wkt_fh) fatal_error("cannot open output file for WKT");
		}
		if(wkb_fn) {
			wkb_fh = fopen(wkb_fn, "w");
			if(!wkb_fh) fatal_error("cannot open output file for WKB");
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
			printf("\nFeature class %d\n", class_id);

			if(color_table) {
				color = GDALGetColorEntry(color_table, class_id);
				if(color) printf("  Color=%d,%d,%d,%d\n",
					color->c1, color->c2, color->c3, color->c4);
			}

			mask = get_mask_for_8bit_raster(georef.w, georef.h,
				raster, (unsigned char)class_id);
		} else {
			if(class_id != 0) continue;
		}

		if(do_invert) {
			invert_mask(mask, georef.w, georef.h);
		}

		if(do_erosion) {
			erode_mask(mask, georef.w, georef.h);
		}

		mpoly_t feature_poly = calc_ring_from_mask(mask, georef.w, georef.h, dbuf, 
			major_ring_only, no_donuts, min_ring_area, bevel_size);
		free(mask);

		if(mask_out_fn) {
			mask_from_mpoly(&feature_poly, georef.w, georef.h, mask_out_fn);
		}

		if(feature_poly.num_rings && reduction_tolerance > 0) {
			mpoly_t reduced_poly = compute_reduced_pointset(&feature_poly, reduction_tolerance);
			free_mpoly(&feature_poly);
			feature_poly = reduced_poly;
		}

//		for(int ridx=0; ridx<feature_poly.num_rings; ridx++) {
//			feature_poly.rings[ridx] = pinch_excursions(feature_poly.rings + ridx);
//		}

		if(feature_poly.num_rings) {
			int num_outer=0, num_inner=0, total_pts=0;
			int r_idx;
			for(r_idx=0; r_idx<feature_poly.num_rings; r_idx++) {
				if(feature_poly.rings[r_idx].is_hole) num_inner++;
				else num_outer++;
				total_pts += feature_poly.rings[r_idx].npts;
			}
			printf("Found %d outer rings and %d holes with a total of %d vertices.\n",
				num_outer, num_inner, total_pts);

			if(dbuf && dbuf->mode == PLOT_CONTOURS) {
				debug_plot_rings(&feature_poly, dbuf);
			}

			if(do_geom_output && feature_poly.num_rings) {
				printf("Writing output\n");

				int num_shapes;
				mpoly_t *shapes;
				int shape_is_copy;
				if(split_polys) {
					split_mpoly_to_polys(&feature_poly, &num_shapes, &shapes);
					shape_is_copy = 0;
				} else {
					num_shapes = 1;
					shapes = (mpoly_t *)malloc_or_die(sizeof(mpoly_t));
					shapes[0] = feature_poly;
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
					if(wkb_fh) {
						int wkb_size = OGR_G_WkbSize(ogr_geom);
						printf("WKB size = %d\n", wkb_size);
						unsigned char *wkb_out = (unsigned char *)malloc_or_die(wkb_size);
						OGR_G_ExportToWkb(ogr_geom, WKB_BYTE_ORDER, wkb_out);
						fwrite(wkb_out, wkb_size, 1, wkb_fh);
						free(wkb_out);
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
		}

		free_mpoly(&feature_poly);
	}

	printf("\n");

	if(wkt_fh) fclose(wkt_fh);
	if(wkb_fh) fclose(wkb_fh);
	if(ogr_ds) OGR_DS_Destroy(ogr_ds);
	if(dbuf) write_plot(dbuf, debug_report);

	if(do_geom_output) {
		if(num_shapes_written) printf("Wrote %d shapes.\n", num_shapes_written);
		else printf("Wrote empty shapefile.\n");
	}

	CPLPopErrorHandler();

	return 0;
}

mpoly_t calc_ring_from_mask(unsigned char *mask, int w, int h,
report_image_t *dbuf, int major_ring_only, int no_donuts, 
long min_ring_area, double bevel_size) {
	if(major_ring_only) no_donuts = 1;

	mpoly_t mp = trace_mask(mask, w, h, min_ring_area, no_donuts);

	if(VERBOSE) {
		int total_pts = 0;
		for(int r_idx=0; r_idx<mp.num_rings; r_idx++) {
			total_pts += mp.rings[r_idx].npts;
		}
		printf("tracer produced %d rings with a total of %d points\n",
			mp.num_rings, total_pts);
	}

	// this is now done directly by tracer
	/*
	if(min_ring_area > 0) {
		if(VERBOSE) printf("removing small rings...\n");

		ring_t *filtered_rings = (ring_t *)malloc_or_die(sizeof(ring_t) * mp.num_rings);
		int *parent_map = (int *)malloc_or_die(sizeof(int) * mp.num_rings);
		for(int i=0; i<mp.num_rings; i++) {
			parent_map[i] = -1;
		}
		int num_filtered_rings = 0;
		for(int i=0; i<mp.num_rings; i++) {
			double area = polygon_area(mp.rings+i);
			if(VERBOSE) if(area > 10) printf("ring %d has area %.15f\n", i, area);
			if(area >= min_ring_area) {
				parent_map[i] = num_filtered_rings;
				filtered_rings[num_filtered_rings++] = mp.rings[i];
			} else {
				free_ring(mp.rings + i);
			}
		}
		for(int i=0; i<num_filtered_rings; i++) {
			int old_parent = filtered_rings[i].parent_id;
			if(old_parent >= 0) {
				int new_parent = parent_map[old_parent];
				if(new_parent < 0) fatal_error("took a ring but not its parent");
				filtered_rings[i].parent_id = new_parent;
			}
		}
		printf("filtered by area %d => %d rings\n",
			mp.num_rings, num_filtered_rings);

		free(mp.rings);
		mp.rings = filtered_rings;
		mp.num_rings = num_filtered_rings;
	}
	*/

	if(major_ring_only && mp.num_rings > 1) {
		double biggest_area = 0;
		int best_idx = 0;
		for(int i=0; i<mp.num_rings; i++) {
			double area = polygon_area(mp.rings+i);
			if(area > biggest_area) {
				biggest_area = area;
				best_idx = i;
			}
		}
		if(VERBOSE) printf("major ring was %d with %d pts, %.1f area\n",
			best_idx, mp.rings[best_idx].npts, biggest_area);
		if(mp.rings[best_idx].parent_id >= 0) fatal_error("largest ring should not have a parent");

		mpoly_t new_mp;
		new_mp.num_rings = 1;
		new_mp.rings = (ring_t *)malloc_or_die(sizeof(ring_t));
		new_mp.rings[0] = mp.rings[best_idx];

		for(int i=0; i<mp.num_rings; i++) {
			if(i != best_idx) free_ring(mp.rings+i);
		}
		free(mp.rings);

		mp = new_mp;
	}

	// this is now done directly by tracer
	/*
	if(no_donuts) {
		// we don't have to worry about remapping parent_id in this
		// case because we only take rings with no parent
		int out_idx = 0;
		for(int i=0; i<mp.num_rings; i++) {
			if(mp.rings[i].parent_id < 0) {
				mp.rings[out_idx++] = mp.rings[i];
			}
		}
		mp.num_rings = out_idx;
		if(VERBOSE) printf("number of non-donut rings is %d", mp.num_rings);
	}
	*/

	if(mp.num_rings && bevel_size > 0) {
		// the topology cannot be resolved by us or by geos/jump/postgis if
		// there are self-intersections
		bevel_self_intersections(&mp, bevel_size);
	}

	return mp;
}
