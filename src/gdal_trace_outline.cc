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



#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <map>
#include <vector>

#include "common.h"
#include "polygon.h"
#include "polygon-rasterizer.h"
#include "debugplot.h"
#include "georef.h"
#include "ndv.h"
#include "mask.h"
#include "mask-tracer.h"
#include "dp.h"
#include "excursion_pincher.h"
#include "beveler.h"
#include "raster_features.h"

#include <ogrsf_frmts.h>
#include <cpl_string.h>
#include <cpl_conv.h>
#include <cpl_port.h>

#ifdef CPL_MSB 
#define WKB_BYTE_ORDER wkbNDR
#else
#define WKB_BYTE_ORDER wkbXDR
#endif

using namespace dangdal;

void usage(const std::string &cmdname) {
	printf("Usage:\n  %s [options] [image_name]\n", cmdname.c_str());
	printf("\n");
	
	GeoOpts::printUsage();
	printf("\n");
	NdvDef::printUsage();

	printf(
"\n"
"Behavior:\n"
"  -classify                    Output a polygon for each value of an 8-bit band\n"
"                               (default is to generate a single polygon that\n"
"                               surrounds all pixels that don't match\n"
"                               the no-data-value)\n"
"  -b band_id -b band_id ...    Bands to inspect (default is all bands)\n"
"  -invert                      Trace no-data pixels rather than data pixels\n"
"  -erosion                     Erode pixels that don't have two consecutive\n"
"                               neighbors\n"
"  -major-ring                  Take only the biggest outer ring\n"
"  -no-donuts                   Take only top-level rings\n"
"  -min-ring-area val           Drop rings with less than this area\n"
"                               (in square pixels)\n"
"  -containing [xy | en | ll | percent] xval yval\n"
"                               Take only rings that contain the given point.\n"
"                               This option can be specified more than once.\n"
"                               'percent' means interpret coordinate as percent\n"
"                               of image width/height.\n"
"  -not-containing [xy | en | ll | percent] xval yval\n"
"                               Take rings that don't contain the given point.\n"
"                               This option can be specified more than once.\n"
"                               'percent' means interpret coordinate as percent\n"
"                               of image width/height.\n"
"  -dp-toler val                Tolerance for polygon simplification\n"
"                               (in pixels, default is 2.0)\n"
"  -bevel-size                  How much to shave off corners at\n"
"                               self-intersection points\n"
"                               (in pixels, default is 0.1)\n"
"                               (this is done to make geometries that\n"
"                               PostGIS/GEOS/Jump can handle)\n"
"  -pinch-excursions            Remove all the complicated 'mouse bites' that\n"
"                               occur in the outline when lossy compression\n"
"                               has been used (experimental)\n"
"\n"
"Output:\n"
"  -report fn.ppm               Output graphical report of polygons found\n"
"  -mask-out fn.pbm             Output mask of bounding polygon in PBM format\n"
"  -out-cs [xy | en | ll]       Set coordinate system for following outputs\n"
"                               (pixel coords, easting/northing, or lon/lat)\n"
"                               Must be specified before -{wkt,wkb,ogr}-out options\n"
"  -llproj-toler val            Error tolerance for curved lines when\n"
"                               using '-out-cs ll' (in pixels, default is 1.0)\n"
"  -wkt-out fn.wkt              Output polygons in WKT format\n"
"  -wkb-out fn.wkb              Output polygons in WKB format\n"
"  -ogr-out fn.shp              Output polygons using an OGR format\n"
"  -ogr-fmt                     OGR format to use (default is 'ESRI Shapefile')\n"
"                               Must be specified before -ogr-out option\n"
"  -split-polys                 Output several polygons rather than one\n"
"                               multipolygon\n"
"\n"
"Misc:\n"
"  -v                           Verbose\n"
"\n"
"Examples:\n"
"\n"
"Inspect image and output contour of data region:\n"
"gdal_trace_outline raster.tif -ndv 0 -erosion -out-cs ll -wkt-out outline.wkt\n"
"\n"
"Same as above but polygon actually follows border pixel-by-pixel:\n"
"gdal_trace_outline raster.tif -ndv 0 -dp-toler 0 -out-cs ll -wkt-out outline.wkt\n"
"\n"
"Output ESRI Shapefile in projection of input image:\n"
"gdal_trace_outline raster.tif -ndv 0 -erosion -out-cs en -ogr-out outline.shp\n"
"\n"
"Generate one shape for each value in input image:\n"
"gdal_trace_outline raster.tif -classify -out-cs en -ogr-out outline.shp\n"
"\n"
	);
	exit(1);
}

enum CoordSystem {
	CS_UNKNOWN,
	CS_XY,
	CS_EN,
	CS_LL,
	CS_PERCENT
};

struct GeomOutput {
	explicit GeomOutput(CoordSystem _out_cs=CS_UNKNOWN) :
		out_cs(_out_cs),
		wkt_fh(NULL),
		wkb_fh(NULL),
		ogr_ds(NULL),
		ogr_layer(NULL)
	{ }

	CoordSystem out_cs;

	std::string wkt_fn;
	FILE *wkt_fh;

	std::string wkb_fn;
	FILE *wkb_fh;

	std::string ogr_fn;
	std::string ogr_fmt;
	OGRDataSourceH ogr_ds;
	OGRLayerH ogr_layer;
};

struct ContainingOption {
	ContainingOption() : x(0), y(0), cs(CS_UNKNOWN), wanted_point(false) { }
	double x, y;
	CoordSystem cs;
	bool wanted_point;
};

Mpoly take_largest_ring(const Mpoly &mp_in);

Mpoly remove_holes(const Mpoly &mp_in);

Mpoly containment_filters(
	const Mpoly &mp_in,
	const std::vector<ContainingOption> &containing_options,
	const GeoRef &georef,
	DebugPlot *dbuf
);

int main(int argc, char **argv) {
	const std::string cmdname = argv[0];
	if(argc == 1) usage(cmdname);
	std::vector<std::string> arg_list = argv_to_list(argc, argv);

	std::string input_raster_fn;
	bool classify = 0;
	std::string debug_report;
	std::vector<size_t> inspect_bandids;
	bool split_polys = 0;
	CoordSystem cur_out_cs = CS_UNKNOWN;
	std::string cur_ogr_fmt = "ESRI Shapefile";
	std::vector<GeomOutput> geom_outputs;
	std::string mask_out_fn;
	bool major_ring_only = 0;
	bool trace_no_donuts = 0;
	bool output_no_donuts = 0;
	int64_t min_ring_area = 0;
	double reduction_tolerance = 2;
	bool do_erosion = 0;
	bool do_invert = 0;
	double llproj_toler = 1;
	double bevel_size = .1;
	bool do_pinch_excursions = 0;
	std::vector<ContainingOption> containing_options;

	GeoOpts geo_opts = GeoOpts(arg_list);
	NdvDef ndv_def = NdvDef(arg_list);

	size_t argp = 1;
	while(argp < arg_list.size()) {
		const std::string &arg = arg_list[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			try {
				if(arg == "-v") {
					VERBOSE++;
				} else if(arg == "-classify") {
					classify = 1;
				} else if(arg == "-report") {
					if(argp == arg_list.size()) usage(cmdname);
					debug_report = arg_list[argp++];
				} else if(arg == "-b") {
					if(argp == arg_list.size()) usage(cmdname);
					int bandid = boost::lexical_cast<int>(arg_list[argp++]);
					inspect_bandids.push_back(bandid);
				} else if(arg == "-erosion") {
					do_erosion = 1;
				} else if(arg == "-invert") {
					do_invert = 1;
				} else if(arg == "-split-polys") {
					split_polys = 1;
				} else if(arg == "-wkt-out") {
					if(argp == arg_list.size()) usage(cmdname);
					GeomOutput go(cur_out_cs);
					go.wkt_fn = arg_list[argp++];
					geom_outputs.push_back(go);
				} else if(arg == "-wkb-out") {
					if(argp == arg_list.size()) usage(cmdname);
					GeomOutput go(cur_out_cs);
					go.wkb_fn = arg_list[argp++];
					geom_outputs.push_back(go);
				} else if(arg == "-ogr-out") {
					if(argp == arg_list.size()) usage(cmdname);
					GeomOutput go(cur_out_cs);
					go.ogr_fmt = cur_ogr_fmt;
					go.ogr_fn = arg_list[argp++];
					geom_outputs.push_back(go);
				} else if(arg == "-ogr-fmt") {
					if(argp == arg_list.size()) usage(cmdname);
					cur_ogr_fmt = arg_list[argp++];
				} else if(arg == "-out-cs") {
					if(argp == arg_list.size()) usage(cmdname);
					std::string cs = arg_list[argp++];
					if(cs == "xy") cur_out_cs = CS_XY;
					else if(cs == "en") cur_out_cs = CS_EN;
					else if(cs == "ll") cur_out_cs = CS_LL;
					else fatal_error("unrecognized value for -out-cs option (%s)", cs.c_str());
				} else if(arg == "-mask-out") {
					if(argp == arg_list.size()) usage(cmdname);
					mask_out_fn = arg_list[argp++];
				} else if(arg == "-major-ring") {
					major_ring_only = 1;
				} else if(arg == "-no-donuts") {
					output_no_donuts = 1;
				} else if(arg == "-min-ring-area") {
					if(argp == arg_list.size()) usage(cmdname);
					min_ring_area = boost::lexical_cast<int64_t>(arg_list[argp++]);
				} else if(arg == "-dp-toler") {
					if(argp == arg_list.size()) usage(cmdname);
					reduction_tolerance = boost::lexical_cast<double>(arg_list[argp++]);
				} else if(arg == "-bevel-size") {
					if(argp == arg_list.size()) usage(cmdname);
					bevel_size = boost::lexical_cast<double>(arg_list[argp++]);
					if(bevel_size < 0 || bevel_size >= 1) fatal_error(
						"-bevel-size must be in the range 0 <= bevel < 1");
				} else if(arg == "-pinch-excursions") {
					do_pinch_excursions = 1;
				} else if(arg == "-llproj-toler") {
					if(argp == arg_list.size()) usage(cmdname);
					llproj_toler = boost::lexical_cast<double>(arg_list[argp++]);
				} else if(arg == "-containing" || arg == "-not-containing") {
					if(argp+3 > arg_list.size()) usage(cmdname);
					ContainingOption opt;
					opt.wanted_point = (arg == "-containing");

					std::string cs = arg_list[argp++];
					if     (cs == "xy") opt.cs = CS_XY;
					else if(cs == "en") opt.cs = CS_EN;
					else if(cs == "ll") opt.cs = CS_LL;
					else if(cs == "percent") opt.cs = CS_PERCENT;
					else fatal_error("unrecognized coordinate system for -containing option (%s)", cs.c_str());

					opt.x = boost::lexical_cast<double>(arg_list[argp++]);
					opt.y = boost::lexical_cast<double>(arg_list[argp++]);

					containing_options.push_back(opt);
				} else if(arg == "-h" || arg == "--help") {
					usage(cmdname);
				} else {
					fatal_error("unrecognized option: %s", arg.c_str());
				}
			} catch(boost::bad_lexical_cast &e) {
				fatal_error("cannot parse number given on command line");
			}
		} else {
			if(input_raster_fn.size()) usage(cmdname);
			input_raster_fn = arg;
		}
	}

	if(input_raster_fn.empty()) fatal_error("must specify filename of image");

	bool do_geom_output = geom_outputs.size();

	if(major_ring_only && output_no_donuts) fatal_error(
		"-major-ring and -no-donuts options cannot both be used at the same time");

	if(do_pinch_excursions && !(output_no_donuts || major_ring_only)) {
		// some extra logic would be needed in pinch_excursions2 in order to
		// support holes
		fatal_error("the -pinch-excursions option requires the -no-donuts or -major-ring options");
	}

	if(classify) {
		if(do_invert) fatal_error("-classify option is not compatible with -invert option");
		if(mask_out_fn.size()) fatal_error("-classify option is not compatible with -mask-out option");
	}

	GDALAllRegister();

	GDALDatasetH ds = GDALOpen(input_raster_fn.c_str(), GA_ReadOnly);
	if(!ds) fatal_error("open failed");

	if(inspect_bandids.empty()) {
		size_t nbands = GDALGetRasterCount(ds);
		for(size_t i=0; i<nbands; i++) inspect_bandids.push_back(i+1);
	}

	if(ndv_def.empty()) {
		ndv_def = NdvDef(ds, inspect_bandids);
	}

	CPLPushErrorHandler(CPLQuietErrorHandler);

	GeoRef georef = GeoRef(geo_opts, ds);

	for(size_t i=0; i<geom_outputs.size(); i++) {
		CoordSystem out_cs = geom_outputs[i].out_cs;
		if(out_cs == CS_UNKNOWN) fatal_error(
			"must specify output coordinate system with -out-cs option before specifying output");
		if((out_cs == CS_EN || out_cs == CS_LL) && !georef.hasAffine()) 
			fatal_error("missing affine transform");
		if((out_cs == CS_LL) && !georef.fwd_xform) 
			fatal_error("missing coordinate transform");
	}

	DebugPlot *dbuf = NULL;
	if(debug_report.size()) {
		dbuf = new DebugPlot(georef.w, georef.h,
			do_pinch_excursions ? PLOT_PINCH : PLOT_CONTOURS);
	}

	// only used if classify==true, but doesn't hurt otherwise
	const FeatureInterpreter feature_interp(ds, inspect_bandids);

	FeatureBitmap *features_bitmap = NULL;
	BitGrid mask(0, 0);
	if(classify) {
		features_bitmap = read_raster_features(ds, inspect_bandids, ndv_def, dbuf);
	} else {
		mask = get_bitgrid_for_dataset(ds, inspect_bandids, ndv_def, dbuf);
	}

	for(size_t go_idx=0; go_idx<geom_outputs.size(); go_idx++) {
		GeomOutput &go = geom_outputs[go_idx];

		if(go.wkt_fn.size()) {
			go.wkt_fh = fopen(go.wkt_fn.c_str(), "w");
			if(!go.wkt_fh) fatal_error("cannot open output file for WKT");
		}
		if(go.wkb_fn.size()) {
			go.wkb_fh = fopen(go.wkb_fn.c_str(), "w");
			if(!go.wkb_fh) fatal_error("cannot open output file for WKB");
		}

		if(go.ogr_fn.size()) {
			OGRRegisterAll();
			if(go.ogr_fmt.empty()) fatal_error("no OGR format was specified");
			OGRSFDriverH ogr_driver = OGRGetDriverByName(go.ogr_fmt.c_str());
			if(!ogr_driver) fatal_error("cannot get OGR driver (%s)", go.ogr_fmt.c_str());
			go.ogr_ds = OGR_Dr_CreateDataSource(ogr_driver, go.ogr_fn.c_str(), NULL);
			if(!go.ogr_ds) fatal_error(
				"cannot create OGR data source (does output file already exist?)");

			std::string layer_name = go.ogr_fn;

			OGRSpatialReferenceH sref = NULL;
			if(go.out_cs == CS_EN) {
				sref = georef.spatial_ref;
			} else if(go.out_cs == CS_LL) {
				sref = georef.geo_sref;
			}

			go.ogr_layer = OGR_DS_CreateLayer(go.ogr_ds, layer_name.c_str(), sref, 
				(split_polys ? wkbPolygon : wkbMultiPolygon), NULL);
			if(!go.ogr_layer) fatal_error("cannot create OGR layer");

			if(classify) {
				feature_interp.create_ogr_fields(go.ogr_layer);
			}
		}
	}

	int num_shapes_written = 0;

	std::map<FeatureRawVal, FeatureBitmap::Index> features_list;
	if(classify) {
		features_list = features_bitmap->feature_table();
	} else {
		// If we are not running in feature classify mode, then just add an arbitrary element
		// to this list so that the loop below runs once.
		features_list[FeatureRawVal()];
	}

	for(const std::pair<FeatureRawVal, FeatureBitmap::Index> &feature : features_list) {
		if(classify) {
			printf("\nTracing feature %s\n", feature_interp.pixel_to_string(feature.first).c_str());
			mask = features_bitmap->get_mask_for_feature(feature.second);
		}

		if(do_invert) {
			mask.invert();
		}

		if(do_erosion) {
			mask.erode();
		}

		if(!containing_options.empty()) {
			// We need to trace donuts even if not outputting them, in order to
			// see if the polygons satisfy the containment options.  Ideally
			// the user should be able to specify which happens first, hole
			// removal or containment options.  Maybe there needs to be a
			// rudimentary scripting language?  Or maybe just process the
			// options in the order they are specified on the command line.

			// Note: in this case, donuts must be removed later on!
			trace_no_donuts = 0;
		} else {
			trace_no_donuts = output_no_donuts;
			// If taking only the major ring, no holes are needed.
			trace_no_donuts |= major_ring_only;
		}
		// If we are only taking the largest ring, and don't need to compute
		// containments, then skip donuts for speed.
		if(major_ring_only && containing_options.empty()) {
			trace_no_donuts = 1;
		}

		Mpoly feature_poly = trace_mask(mask, georef.w, georef.h, min_ring_area, trace_no_donuts);
		mask = BitGrid(0, 0); // free some memory

		if(VERBOSE) {
			size_t num_inner = 0, num_outer = 0, total_pts = 0;
			for(size_t r_idx=0; r_idx<feature_poly.rings.size(); r_idx++) {
				if(feature_poly.rings[r_idx].is_hole) num_inner++;
				else num_outer++;
				total_pts += feature_poly.rings[r_idx].pts.size();
			}
			printf("tracer produced %zd rings (%zd outer, %zd holes) with a total of %zd points\n",
				feature_poly.rings.size(), num_outer, num_inner, total_pts);
		}

		if(!feature_poly.rings.empty() && !containing_options.empty()) {
			feature_poly = containment_filters(feature_poly, containing_options, georef, dbuf);
		}

		if(major_ring_only && feature_poly.rings.size() > 1) {
			printf("Taking largest ring.\n");
			feature_poly = take_largest_ring(feature_poly);
		}

		if(output_no_donuts && !trace_no_donuts) {
			// Hole removal was deferred until now.
			printf("Removing donut holes.\n");
			feature_poly = remove_holes(feature_poly);
		}

		if(!feature_poly.rings.empty() && bevel_size > 0) {
			// the topology cannot be resolved by us or by geos/jump/postgis if
			// there are self-intersections
			bevel_self_intersections(feature_poly, bevel_size);
		}

		if(feature_poly.rings.size() && do_pinch_excursions) {
			printf("Pinching excursions...\n");
			feature_poly = pinch_excursions2(feature_poly, dbuf);
			printf("Done pinching excursions.\n");
		}

		if(mask_out_fn.size()) {
			mask_from_mpoly(feature_poly, georef.w, georef.h, mask_out_fn);
		}

		if(feature_poly.rings.size() && reduction_tolerance > 0) {
			Mpoly reduced_poly = compute_reduced_pointset(feature_poly, reduction_tolerance);
			feature_poly = reduced_poly;
		}

		if(feature_poly.rings.empty()) {
			printf("WARNING! No rings found!\n");
		} else {
			size_t num_outer=0, num_inner=0, total_pts=0;
			for(size_t r_idx=0; r_idx<feature_poly.rings.size(); r_idx++) {
				if(feature_poly.rings[r_idx].is_hole) num_inner++;
				else num_outer++;
				total_pts += feature_poly.rings[r_idx].pts.size();
			}
			printf("Found %zd outer rings and %zd holes with a total of %zd vertices.\n",
				num_outer, num_inner, total_pts);

			if(dbuf && dbuf->mode == PLOT_CONTOURS) {
				dbuf->debugPlotMpoly(feature_poly);
			}

			if(do_geom_output && feature_poly.rings.size()) {
				printf("Writing output\n");

				std::vector<Mpoly> shapes;
				if(split_polys) {
					shapes = split_mpoly_to_polys(feature_poly);
				} else {
					shapes.push_back(feature_poly);
				}

				for(size_t shape_idx=0; shape_idx<shapes.size(); shape_idx++) {
					const Mpoly &poly_in = shapes[shape_idx];

					for(size_t go_idx=0; go_idx<geom_outputs.size(); go_idx++) {
						GeomOutput &go = geom_outputs[go_idx];

						Mpoly proj_poly = poly_in;
						if(go.out_cs == CS_XY) {
							// no-op
						} else if(go.out_cs == CS_EN) {
							proj_poly.xy2en(georef);
						} else if(go.out_cs == CS_LL) {
							proj_poly.xy2ll_with_interp(georef, llproj_toler);
						} else {
							fatal_error("bad val for out_cs");
						}

						OGRGeometryH ogr_geom = mpoly_to_ogr(proj_poly);

						if(go.wkt_fh) {
							char *wkt_out;
							OGR_G_ExportToWkt(ogr_geom, &wkt_out);
							fprintf(go.wkt_fh, "%s\n", wkt_out);
						}
						if(go.wkb_fh) {
							size_t wkb_size = OGR_G_WkbSize(ogr_geom);
							printf("WKB size = %zd\n", wkb_size);
							std::vector<unsigned char> wkb_out(wkb_size);
							OGR_G_ExportToWkb(ogr_geom, WKB_BYTE_ORDER, &wkb_out[0]);
							fwrite(&wkb_out[0], wkb_size, 1, go.wkb_fh);
						}

						if(go.ogr_ds) {
							OGRFeatureH ogr_feat = OGR_F_Create(OGR_L_GetLayerDefn(go.ogr_layer));

							if(classify) {
								feature_interp.set_ogr_fields(go.ogr_layer, ogr_feat, feature.first);
							}

							OGR_F_SetGeometryDirectly(ogr_feat, ogr_geom); // assumes ownership of geom
							OGR_L_CreateFeature(go.ogr_layer, ogr_feat);
							OGR_F_Destroy(ogr_feat);
						} else {
							OGR_G_DestroyGeometry(ogr_geom);
						}
					}

					num_shapes_written++;
				}
			}
		}
	}

	printf("\n");

	for(size_t go_idx=0; go_idx<geom_outputs.size(); go_idx++) {
		GeomOutput &go = geom_outputs[go_idx];
		if(go.wkt_fh) fclose(go.wkt_fh);
		if(go.wkb_fh) fclose(go.wkb_fh);
		if(go.ogr_ds) OGR_DS_Destroy(go.ogr_ds);
	}

	if(dbuf) dbuf->writePlot(debug_report);

	if(do_geom_output) {
		if(num_shapes_written) printf("Wrote %d shapes.\n", num_shapes_written);
		else printf("Wrote empty shapefile.\n");
	}

	GDALClose(ds);

	CPLPopErrorHandler();

	return 0;
}

Mpoly take_largest_ring(const Mpoly &mp_in) {
	double biggest_area = 0;
	size_t best_idx = 0;
	for(size_t i=0; i<mp_in.rings.size(); i++) {
		double area = mp_in.rings[i].area();
		if(area > biggest_area) {
			biggest_area = area;
			best_idx = i;
		}
	}
	if(VERBOSE) printf("major ring was %zd with %zd pts, %.1f area\n",
		best_idx, mp_in.rings[best_idx].pts.size(), biggest_area);
	if(mp_in.rings[best_idx].parent_id >= 0) fatal_error("largest ring should not have a parent");

	Mpoly new_mp;
	new_mp.rings.push_back(mp_in.rings[best_idx]);
	return new_mp;
}

Mpoly remove_holes(const Mpoly &mp_in) {
	Mpoly new_mp;

	for(size_t i=0; i<mp_in.rings.size(); i++) {
		const Ring &ring = mp_in.rings[i];
		// Take only top-level rings.  Since we are filling holes, it doesn't
		// make sense to keep an island within a hole.
		if(ring.parent_id < 0) {
			new_mp.rings.push_back(ring);
		}
	}

	return new_mp;
}

Mpoly containment_filters(
	const Mpoly &mp_in,
	const std::vector<ContainingOption> &containing_options,
	const GeoRef &georef,
	DebugPlot *dbuf
) {
	std::vector<Vertex> wanted_pts;
	std::vector<Vertex> unwanted_pts;
	BOOST_FOREACH(const ContainingOption &opt, containing_options) {
		Vertex v;
		switch(opt.cs) {
			case CS_XY:
				v.x = opt.x;
				v.y = opt.y;
				break;
			case CS_PERCENT:
				v.x = opt.x / 100.0 * georef.w;
				v.y = opt.y / 100.0 * georef.h;
				break;
			case CS_EN:
				georef.en2xy(opt.x, opt.y, &v.x, &v.y);
				break;
			case CS_LL:
				georef.ll2xy(opt.x, opt.y, &v.x, &v.y);
				break;
			default:
				fatal_error("coord system not implemented");
		};
		if(opt.wanted_point) {
			wanted_pts.push_back(v);
			dbuf->plotPointBig(v.x, v.y, 0, 255, 0);
		} else {
			unwanted_pts.push_back(v);
			dbuf->plotPointBig(v.x, v.y, 255, 0, 0);
		}
	}

	if(wanted_pts.empty() && unwanted_pts.empty()) {
		fatal_error("no wanted/unwanted pts given");
	}

	printf("Looking for polygons");
	if(!wanted_pts.empty()) {
		printf(" containing:");
		BOOST_FOREACH(const Vertex &v, wanted_pts) {
			printf(" (%.1f,%.1f)", v.x, v.y);
		}
	}
	if(!unwanted_pts.empty()) {
		printf(" not containing:");
		BOOST_FOREACH(const Vertex &v, unwanted_pts) {
			printf(" (%.1f,%.1f)", v.x, v.y);
		}
	}
	printf("\n");

	Mpoly new_mp;
	std::map<int, int> relabeling;

	int num_outer=0, num_holes=0;

	for(size_t outer_idx=0; outer_idx<mp_in.rings.size(); outer_idx++) {
		const Ring &outer = mp_in.rings[outer_idx];
		if(outer.is_hole) continue;

		bool contains_wanted_pt = false;
		BOOST_FOREACH(const Vertex &v, wanted_pts) {
			if(mp_in.component_contains(v, outer_idx)) {
				if(VERBOSE) printf("ring %zd contains wanted point %g,%g\n",
					outer_idx, v.x, v.y);
				contains_wanted_pt = true;
				break;
			}
		}
		if(!wanted_pts.empty() && !contains_wanted_pt) {
			// not interested in this ring
			continue;
		}

		bool contains_unwanted_pt = false;
		BOOST_FOREACH(const Vertex &v, unwanted_pts) {
			if(mp_in.component_contains(v, outer_idx)) {
				if(VERBOSE) printf("ring %zd contains unwanted point %g,%g\n",
					outer_idx, v.x, v.y);
				contains_unwanted_pt = true;
				break;
			}
		}
		if(contains_unwanted_pt) {
			// not interested in this ring
			continue;
		}

		relabeling[outer_idx] = int(new_mp.rings.size());
		new_mp.rings.push_back(outer);
		num_outer++;

		for(size_t j=0; j<mp_in.rings.size(); j++) {
			Ring inner = mp_in.rings[j];
			// take children of outer ring
			if(inner.parent_id != int(outer_idx)) continue;

			relabeling[j] = int(new_mp.rings.size());
			new_mp.rings.push_back(inner);
			num_holes++;
		}
	}

	for(size_t i=0; i<new_mp.rings.size(); i++) {
		Ring &ring = new_mp.rings[i];
		// Compute new parent.  Iterate outwards through the ring hierarchy
		// until a ring is found that has been kept.
		while(ring.parent_id >= 0) {
			if(relabeling.count(ring.parent_id)) {
				ring.parent_id = relabeling[ring.parent_id];
				break;
			} else {
				ring.parent_id = mp_in.rings[ring.parent_id].parent_id;
			}
		}
	}

	if(new_mp.rings.empty()) {
		printf("   None found!\n");
	} else {
		printf("   Found %d connected components (with %d holes).\n", num_outer, num_holes);
	}

	return new_mp;
}
