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



#include "common.h"
#include "polygon.h"
#include "polygon-rasterizer.h"
#include "debugplot.h"
#include "georef.h"
#include "ndv.h"
#include "mask.h"
#include "rectangle_finder.h"

using namespace dangdal;

void usage(const char *cmdname) {
	printf("Usage:\n  %s [options] [image_name]\n", cmdname);
	printf("\n");
	
	GeoOpts::printUsage();
	printf("\n");
	NdvDef::printUsage();

	printf("\
\n\
Inspection:\n\
  -inspect-rect4              Attempt to find 4-sided bounding polygon\n\
  -fuzzy-match                Try to exclude logos and other extraneous\n\
                              pixels from bounding polygon\n\
  -b band_id -b band_id ...   Bands to inspect (default is all bands)\n\
  -erosion                    Erode pixels that don't have two consecutive\n\
                              neighbors\n\
  -report fn.ppm              Output graphical report of bounds found\n\
  -mask-out fn.pbm            Output mask of bounding polygon in PBM format\n\
\n\
Misc:\n\
  -v                          Verbose\n\
\n\
Examples:\n\
  Output basic geocoding info:\n\
    gdal_list_corners raster.tif > geocode.yaml\n\
  Inspect image to find corners of actual data (arbitrary four-sided region):\n\
    gdal_list_corners raster.tif -inspect-rect4 -nodataval 0 > geocode.yaml\n\
\n\
");
	exit(1);
}

int main(int argc, char **argv) {
	char *input_raster_fn = NULL;

	bool inspect_rect4 = 0;
	bool fuzzy_match = 0;
	char *debug_report = NULL;
	char *mask_out_fn = NULL;
	std::vector<size_t> inspect_bandids;
	bool do_erosion = 0;

	if(argc == 1) usage(argv[0]);

	// We will be sending YAML to stdout, so stuff that would normally
	// go to stdout (such as debug messages or progress bars) should
	// go to stderr.
	// See http://forums.devshed.com/c-programming-42/redirect-standard-error-and-assert-how-to-52650.html
	FILE *yaml_fh = fdopen(dup(1), "w");
	close(1);
	dup2(2, 1);

	GeoOpts geo_opts = GeoOpts(&argc, &argv);
	NdvDef ndv_def = NdvDef(&argc, &argv);

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(!strcmp(arg, "-v")) {
				VERBOSE++;
			} else if(!strcmp(arg, "-inspect-rect4")) {
				inspect_rect4 = 1;
			} else if(!strcmp(arg, "-fuzzy-match")) {
				fuzzy_match = 1;
			} else if(!strcmp(arg, "-b")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				int bandid = strtol(argv[argp++], &endptr, 10);
				if(*endptr) usage(argv[0]);
				inspect_bandids.push_back(bandid);
			} else if(!strcmp(arg, "-erosion")) {
				do_erosion = 1;
			} else if(!strcmp(arg, "-report")) {
				if(argp == argc) usage(argv[0]);
				debug_report = argv[argp++];
			} else if(!strcmp(arg, "-mask-out")) {
				if(argp == argc) usage(argv[0]);
				mask_out_fn = argv[argp++];
			} else usage(argv[0]);
		} else {
			if(input_raster_fn) usage(argv[0]);
			input_raster_fn = arg;
		}
	}

	bool do_inspect = inspect_rect4;
	if(do_inspect && !input_raster_fn) fatal_error("must specify filename of image");

	GDALAllRegister();

	GDALDatasetH ds = NULL;
	if(input_raster_fn) {
		ds = GDALOpen(input_raster_fn, GA_ReadOnly);
		if(!ds) fatal_error("open failed");
	}

	if(do_inspect && inspect_bandids.empty()) {
		size_t nbands = GDALGetRasterCount(ds);
		for(size_t i=0; i<nbands; i++) inspect_bandids.push_back(i+1);
	}

	if(!do_inspect) {
		if(fuzzy_match)      fatal_error("-fuzzy-match option can only be used with -inspect-rect4 option");
		if(!ndv_def.empty()) fatal_error("NDV options can only be used with -inspect-rect4 option");
		if(debug_report)     fatal_error("-report option can only be used with -inspect-rect4 option");
		if(mask_out_fn)      fatal_error("-mask-out option can only be used with -inspect-rect4 option");
		if(inspect_bandids.size()) fatal_error("-b option can only be used with -inspect-rect4 option");
		if(do_erosion)       fatal_error("-erosion option can only be used with -inspect-rect4 option");
	}

	CPLPushErrorHandler(CPLQuietErrorHandler);

	GeoRef georef = GeoRef(geo_opts, ds);

	DebugPlot *dbuf = NULL;
	BitGrid mask(0, 0);
	if(do_inspect) {
		if(ndv_def.empty()) {
			ndv_def = NdvDef(ds, inspect_bandids);
		}

		if(debug_report) {
			dbuf = new DebugPlot(georef.w, georef.h);
			dbuf->mode = PLOT_RECT4;
		}

		mask = get_bitgrid_for_dataset(ds, inspect_bandids, ndv_def, dbuf);

		if(do_erosion) {
			mask.erode();
		}
	}

	// output phase

	fprintf(yaml_fh, "---\n"); // begin YAML document

	fprintf(yaml_fh, "width: %zd\nheight: %zd\n", georef.w, georef.h);

	if(ds) {
		int band_count = GDALGetRasterCount(ds);
		const char *datatypes = "";
		for(int i=0; i<band_count; i++) {
			GDALRasterBandH band = GDALGetRasterBand(ds, i+1);
			GDALDataType gdt = GDALGetRasterDataType(band);
			const char *dt = GDALGetDataTypeName(gdt);
			if(i) {
				char *join_str = MYALLOC(char,
					strlen(datatypes) + 1 + strlen(dt) + 1);
				sprintf(join_str, "%s,%s", datatypes, dt);
				datatypes = join_str;
			} else {
				datatypes = dt;
			}
		}
		fprintf(yaml_fh, "num_bands: %d\n", band_count);
		fprintf(yaml_fh, "datatype: %s\n", datatypes);

		char **metadata = GDALGetMetadata(ds, "");
		if(metadata) {
			fprintf(yaml_fh, "metadata:\n");
			for(char **p=metadata; *p; p++) {
				fprintf(yaml_fh, "  - '%s'\n", *p);
			}
		}
	}

	if(georef.s_srs && strlen(georef.s_srs)) {
		fprintf(yaml_fh, "s_srs: '%s'\n", georef.s_srs);
	}
	if(georef.units_name) {
		fprintf(yaml_fh, "units_name: '%s'\n", georef.units_name);
	}
	if(georef.units_val) {
		fprintf(yaml_fh, "units_val: %lf\n", georef.units_val);
	}
	if(georef.res_x && georef.res_y) {
		fprintf(yaml_fh, "res: %.15f %.15f\n", georef.res_x, georef.res_y);
	}
	if(georef.res_meters_x && georef.res_meters_y) {
		fprintf(yaml_fh, "res_meters: %.15f %.15f\n", georef.res_meters_x, georef.res_meters_y);
	}
	if(georef.hasAffine()) {
		fprintf(yaml_fh, "affine:\n");
		for(int i=0; i<6; i++) fprintf(yaml_fh, "  - %.15f\n", georef.fwd_affine[i]);
	}

	double lon, lat;
	double east, north;

	Vertex center;
	fprintf(yaml_fh, "center:\n");
	center = Vertex((double)georef.w/2.0, (double)georef.h/2.0);
	if(georef.fwd_xform && georef.hasAffine()) {
		georef.xy2ll_or_die(center.x, center.y, &lon, &lat);
		fprintf(yaml_fh, "  lon: %.15f\n", lon);
		fprintf(yaml_fh, "  lat: %.15f\n", lat);
	}
	if(georef.hasAffine()) {
		georef.xy2en(center.x, center.y, &east, &north);
		fprintf(yaml_fh, "  east: %.15f\n", east);
		fprintf(yaml_fh, "  north: %.15f\n", north);
	}
	fprintf(yaml_fh, "  x: %.15f\n", center.x);
	fprintf(yaml_fh, "  y: %.15f\n", center.y);

	if(do_inspect) {
		Vertex centroid;
		fprintf(yaml_fh, "centroid:\n");
		centroid = mask.centroid();
		if(georef.fwd_xform && georef.hasAffine()) {
			georef.xy2ll_or_die(centroid.x, centroid.y, &lon, &lat);
			fprintf(yaml_fh, "  lon: %.15f\n", lon);
			fprintf(yaml_fh, "  lat: %.15f\n", lat);
		}
		if(georef.hasAffine()) {
			georef.xy2en(centroid.x, centroid.y, &east, &north);
			fprintf(yaml_fh, "  east: %.15f\n", east);
			fprintf(yaml_fh, "  north: %.15f\n", north);
		}
		fprintf(yaml_fh, "  x: %.15f\n", centroid.x);
		fprintf(yaml_fh, "  y: %.15f\n", centroid.y);
	}

	if(inspect_rect4) {
		Ring rect4 = calc_rect4_from_mask(mask, georef.w, georef.h, dbuf, fuzzy_match);

		if(rect4.pts.size() != 4) {
			fatal_error("could not find four-sided region");
		}

		if(mask_out_fn) {
			Mpoly bpoly;
			bpoly.rings.push_back(rect4);

			mask_from_mpoly(bpoly, georef.w, georef.h, mask_out_fn);
		}

		const char *labels[] = { "upper_left", "upper_right", "lower_right", "lower_left" };
		if(georef.fwd_xform && georef.hasAffine()) {
			fprintf(yaml_fh, "geometry_ll:\n  type: rectangle4\n");
			for(int i=0; i<4; i++) {
				georef.xy2ll_or_die(rect4.pts[i].x, rect4.pts[i].y, &lon, &lat);
				fprintf(yaml_fh, "  %s_lon: %.15f\n", labels[i], lon);
				fprintf(yaml_fh, "  %s_lat: %.15f\n", labels[i], lat);
			}
		}
		if(georef.hasAffine()) {
			fprintf(yaml_fh, "geometry_en:\n  type: rectangle4\n");
			for(int i=0; i<4; i++) {
				georef.xy2en(rect4.pts[i].x, rect4.pts[i].y, &east, &north);
				fprintf(yaml_fh, "  %s_east: %.15f\n", labels[i], east);
				fprintf(yaml_fh, "  %s_north: %.15f\n", labels[i], north);
			}
		}
		fprintf(yaml_fh, "geometry_xy:\n  type: rectangle4\n");
		for(int i=0; i<4; i++) {
			fprintf(yaml_fh, "  %s_x: %.15f\n", labels[i], rect4.pts[i].x);
			fprintf(yaml_fh, "  %s_y: %.15f\n", labels[i], rect4.pts[i].y);
		}
	} else {
		const char *e_labels[] = { "left", "mid", "right" };
		double e_pos[] = { 0, (double)georef.w/2.0, georef.w };
		const char *n_labels[] = { "upper", "mid", "lower" };
		double n_pos[] = { 0, (double)georef.h/2.0, georef.h };
		if(georef.fwd_xform && georef.hasAffine()) {
			fprintf(yaml_fh, "geometry_ll:\n  type: rectangle8\n");
			for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
				if(!strcmp(e_labels[i], "mid") && !strcmp(n_labels[j], "mid")) continue;
				georef.xy2ll_or_die(e_pos[i], n_pos[j], &lon, &lat);
				fprintf(yaml_fh, "  %s_%s_lon: %.15f\n", n_labels[j], e_labels[i], lon);
				fprintf(yaml_fh, "  %s_%s_lat: %.15f\n", n_labels[j], e_labels[i], lat);
			}
		}
		if(georef.hasAffine()) {
			fprintf(yaml_fh, "geometry_en:\n  type: rectangle8\n");
			for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
				if(!strcmp(e_labels[i], "mid") && !strcmp(n_labels[j], "mid")) continue;
				georef.xy2en(e_pos[i], n_pos[j], &east, &north);
				fprintf(yaml_fh, "  %s_%s_east: %.15f\n", n_labels[j], e_labels[i], east);
				fprintf(yaml_fh, "  %s_%s_north: %.15f\n", n_labels[j], e_labels[i], north);
			}
		}
		fprintf(yaml_fh, "geometry_xy:\n  type: rectangle8\n");
		for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
			if(!strcmp(e_labels[i], "mid") && !strcmp(n_labels[j], "mid")) continue;
			fprintf(yaml_fh, "  %s_%s_x: %.15f\n", n_labels[j], e_labels[i], e_pos[i]);
			fprintf(yaml_fh, "  %s_%s_y: %.15f\n", n_labels[j], e_labels[i], n_pos[j]);
		}
	}

	if(dbuf) dbuf->writePlot(debug_report);

	CPLPopErrorHandler();

	return 0;
}
