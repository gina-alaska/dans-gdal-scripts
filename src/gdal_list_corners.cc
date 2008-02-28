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

void usage(char *cmdname) {
	printf("Usage:\n  %s [options] [image_name]\n", cmdname);
	printf("\n");
	
	print_georef_usage();

	printf("\
\n\
Inspection:\n\
  -inspect-rect4                  Attempt to find 4-sided bounding polygon\n\
  -nodataval 'val [val ...]'      Specify value of no-data pixels\n\
  -ndv-toler val                  Tolerance for deciding if a pixel\n\
                                  matches nodataval\n\
  -b band_id -b band_id ...       Bands to inspect (default is all bands)\n\
  -skip-erosion                   Don't use erosion filter\n\
  -report fn.ppm                  Output graphical report of bounds found\n\
\n\
Misc:\n\
  -v                              Verbose\n\
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

ring_t calc_rect4_from_mask(unsigned char *mask, int w, int h, report_image_t *dbuf);

int main(int argc, char **argv) {
	char *input_raster_fn = NULL;

	int inspect_rect4 = 0;
	int num_ndv = 0;
	double *ndv_list = NULL;
	double ndv_tolerance = 0;
	char *debug_report = NULL;
	int inspect_numbands = 0;
	int *inspect_bandids = NULL;
	int skip_erosion = 0;

	int i, j;

	if(argc == 1) usage(argv[0]);

	// We will be sending YAML to stdout, so stuff that would normally
	// go to stdout (such as debug messages or progress bars) should
	// go to stderr.
	// See http://forums.devshed.com/c-programming-42/redirect-standard-error-and-assert-how-to-52650.html
	FILE *yaml_fh = fdopen(dup(1), "w");
	close(1);
	dup2(2, 1);

	geo_opts_t geo_opts = init_geo_options(&argc, &argv);

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(!strcmp(arg, "-v")) {
				VERBOSE++;
			} else if(!strcmp(arg, "-inspect-rect4")) {
				inspect_rect4++;
			} else if(!strcmp(arg, "-nodataval")) {
				if(argp == argc) usage(argv[0]);
				int result = parse_list_of_doubles(argv[argp++], &num_ndv, &ndv_list);
				if(result) fatal_error("input to -nodataval must be space-separated list of numbers");
			} else if(!strcmp(arg, "-ndv-toler")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				ndv_tolerance = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);
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
			} else if(!strcmp(arg, "-report")) {
				if(argp == argc) usage(argv[0]);
				debug_report = argv[argp++];
			} else usage(argv[0]);
		} else {
			if(input_raster_fn) usage(argv[0]);
			input_raster_fn = arg;
		}
	}

	int do_inspect = inspect_rect4;
	if(do_inspect && !input_raster_fn) fatal_error("must specify filename of image");

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

	georef_t georef = init_georef(&geo_opts, ds);

	report_image_t *dbuf = NULL;
	unsigned char *mask = NULL;
	if(do_inspect) {
		setup_ndv_list(ds, inspect_numbands, inspect_bandids, &num_ndv, &ndv_list);

		if(debug_report) {
			dbuf = create_plot(georef.w, georef.h);
			dbuf->mode = PLOT_RECT4;
		}

		mask = get_mask_for_dataset(ds, inspect_numbands, inspect_bandids,
			num_ndv, ndv_list, ndv_tolerance, dbuf);
		if(!skip_erosion) {
			unsigned char *eroded_mask = erode_mask(mask, georef.w, georef.h);
			free(mask);
			mask = eroded_mask;
		}
	}

	// output phase

	fprintf(yaml_fh, "width: %d\nheight: %d\n", georef.w, georef.h);

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
		fprintf(yaml_fh, "num_bands: %d\n", band_count);
		fprintf(yaml_fh, "datatype: %s\n", datatypes);
	}

	if(georef.s_srs && strlen(georef.s_srs)) {
		fprintf(yaml_fh, "s_srs: '%s'\n", georef.s_srs);
	}
	if(georef.res_x && georef.res_y) fprintf(yaml_fh, "res: %.15f %.15f\n", georef.res_x, georef.res_y);
	if(georef.fwd_affine) {
		fprintf(yaml_fh, "affine:\n");
		for(i=0; i<6; i++) fprintf(yaml_fh, "  - %.15f\n", georef.fwd_affine[i]);
	}

	double lon, lat;
	double east, north;

	vertex_t center;
	fprintf(yaml_fh, "center:\n");
	center = (vertex_t){ (double)georef.w/2.0, (double)georef.h/2.0 };
	if(georef.fwd_xform && georef.fwd_affine) {
		xy2ll(&georef, center.x, center.y, &lon, &lat);
		fprintf(yaml_fh, "  lon: %.15f\n", lon);
		fprintf(yaml_fh, "  lat: %.15f\n", lat);
	}
	if(georef.fwd_affine) {
		xy2en(&georef, center.x, center.y, &east, &north);
		fprintf(yaml_fh, "  east: %.15f\n", east);
		fprintf(yaml_fh, "  north: %.15f\n", north);
	}
	fprintf(yaml_fh, "  x: %.15f\n", center.x);
	fprintf(yaml_fh, "  y: %.15f\n", center.y);

	if(do_inspect) {
		vertex_t centroid;
		fprintf(yaml_fh, "centroid:\n");
		centroid = calc_centroid_from_mask(mask, georef.w, georef.h);
		if(georef.fwd_xform && georef.fwd_affine) {
			xy2ll(&georef, centroid.x, centroid.y, &lon, &lat);
			fprintf(yaml_fh, "  lon: %.15f\n", lon);
			fprintf(yaml_fh, "  lat: %.15f\n", lat);
		}
		if(georef.fwd_affine) {
			xy2en(&georef, centroid.x, centroid.y, &east, &north);
			fprintf(yaml_fh, "  east: %.15f\n", east);
			fprintf(yaml_fh, "  north: %.15f\n", north);
		}
		fprintf(yaml_fh, "  x: %.15f\n", centroid.x);
		fprintf(yaml_fh, "  y: %.15f\n", centroid.y);
	}

	if(inspect_rect4) {
		ring_t rect4 = calc_rect4_from_mask(mask, georef.w, georef.h, dbuf);

		if(rect4.npts == 4) {
			char *labels[] = { "upper_left", "upper_right", "lower_right", "lower_left" };
			if(georef.fwd_xform && georef.fwd_affine) {
				fprintf(yaml_fh, "geometry_ll:\n  type: rectangle4\n");
				for(i=0; i<4; i++) {
					xy2ll(&georef, rect4.pts[i].x, rect4.pts[i].y, &lon, &lat);
					fprintf(yaml_fh, "  %s_lon: %.15f\n", labels[i], lon);
					fprintf(yaml_fh, "  %s_lat: %.15f\n", labels[i], lat);
				}
			}
			if(georef.fwd_affine) {
				fprintf(yaml_fh, "geometry_en:\n  type: rectangle4\n");
				for(i=0; i<4; i++) {
					xy2en(&georef, rect4.pts[i].x, rect4.pts[i].y, &east, &north);
					fprintf(yaml_fh, "  %s_east: %.15f\n", labels[i], east);
					fprintf(yaml_fh, "  %s_north: %.15f\n", labels[i], north);
				}
			}
			fprintf(yaml_fh, "geometry_xy:\n  type: rectangle4\n");
			for(i=0; i<4; i++) {
				fprintf(yaml_fh, "  %s_x: %.15f\n", labels[i], rect4.pts[i].x);
				fprintf(yaml_fh, "  %s_y: %.15f\n", labels[i], rect4.pts[i].y);
			}
		}
	} else {
		char *e_labels[] = { "left", "mid", "right" };
		double e_pos[] = { 0, (double)georef.w/2.0, georef.w };
		char *n_labels[] = { "upper", "mid", "lower" };
		double n_pos[] = { 0, (double)georef.h/2.0, georef.h };
		if(georef.fwd_xform && georef.fwd_affine) {
			fprintf(yaml_fh, "geometry_ll:\n  type: rectangle8\n");
			for(i=0; i<3; i++) for(j=0; j<3; j++) {
				if(!strcmp(e_labels[i], "mid") && !strcmp(n_labels[j], "mid")) continue;
				xy2ll(&georef, e_pos[i], n_pos[j], &lon, &lat);
				fprintf(yaml_fh, "  %s_%s_lon: %.15f\n", n_labels[j], e_labels[i], lon);
				fprintf(yaml_fh, "  %s_%s_lat: %.15f\n", n_labels[j], e_labels[i], lat);
			}
		}
		if(georef.fwd_affine) {
			fprintf(yaml_fh, "geometry_en:\n  type: rectangle8\n");
			for(i=0; i<3; i++) for(j=0; j<3; j++) {
				if(!strcmp(e_labels[i], "mid") && !strcmp(n_labels[j], "mid")) continue;
				xy2en(&georef, e_pos[i], n_pos[j], &east, &north);
				fprintf(yaml_fh, "  %s_%s_east: %.15f\n", n_labels[j], e_labels[i], east);
				fprintf(yaml_fh, "  %s_%s_north: %.15f\n", n_labels[j], e_labels[i], north);
			}
		}
		fprintf(yaml_fh, "geometry_xy:\n  type: rectangle8\n");
		for(i=0; i<3; i++) for(j=0; j<3; j++) {
			if(!strcmp(e_labels[i], "mid") && !strcmp(n_labels[j], "mid")) continue;
			fprintf(yaml_fh, "  %s_%s_x: %.15f\n", n_labels[j], e_labels[i], e_pos[i]);
			fprintf(yaml_fh, "  %s_%s_y: %.15f\n", n_labels[j], e_labels[i], n_pos[j]);
		}
	}

	if(dbuf) write_plot(dbuf, debug_report);

	CPLPopErrorHandler();

	return 0;
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
	//if(VERBOSE) printf("start point: %d,%d\n", fulcrum_x, fulcrum_y);

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
				//if(VERBOSE) printf("%d %d   %.15f\n", i, j, atan2((double)pix_dy, (double)pix_dx)*180.0/PI);
				best_dx = pix_dx; best_dy = pix_dy;
				best_x = i; best_y = j;
			}
		}
		//if(VERBOSE) printf("  f=[%3d,%3d] ", best_x, best_y);
		//if(VERBOSE) printf("  a=[% 3.1f]\n", angle);

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
			printf("a=%.15f  l=%.15f  ", l.angle, l.seg_len);
			printf("  l2=%.15f  ad=%.15f  ratio=%.15f", len, adiff, len/adiff);
			l = all_edges[i];
			r = all_edges[(i+1) % num_edges];
			printf("  lg=%d rg=%d\n", l.group, r.group);
		}
	}
	if(VERBOSE) printf("num groups: %d\n", num_groups);
	for(i=0; i<num_edges; i++) {
		if(all_edges[i].group < 0) fatal_error("edge not assigned to a group");
		//all_edges[i].group = (num_groups++);
	}
	//if(VERBOSE) printf("num groups: %d\n", num_groups);

	if(VERBOSE) for(i=0; i<num_edges; i++) {
		edge_t l = all_edges[i];
		edge_t r = all_edges[(i+1) % num_edges];
		double len = l.seg_len + r.seg_len;
		double adiff = ang_diff(l.angle, r.angle);
		printf("a=%.15f  l=%.15f  ", l.angle, l.seg_len);
		printf("  l2=%.15f  ad=%.15f  ratio=%.15f", len, adiff, len/adiff);
		printf("  group=%d\n", l.group);
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
			fatal_error("group out of range (i=%d, g=%d, num_groups=%d)", i, eg, num_groups);
		}
		groups[eg].arc_len += e.seg_len;
		groups[eg].wx += e.seg_len * cos(e.angle / 180.0 * PI);
		groups[eg].wy += e.seg_len * sin(e.angle / 180.0 * PI);
	}
	for(i=0; i<num_groups; i++) {
		if(VERBOSE) printf("group %d: l=%.15f\n", i, groups[i].arc_len);
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
	if(VERBOSE) printf("num groups: %d\n", num_groups);
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
		if(VERBOSE) printf("group %d: l=%.15f  a=%.15f  b=%.15f\n", i, groups[i].arc_len, groups[i].avg_ang, groups[i].best_edge.angle);
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
	if(VERBOSE) printf("sorted:\n");
	if(VERBOSE) for(i=0; i<num_groups; i++) {
		printf("group %d: l=%.15f  a=%.15f  s=%.15f\n", i, groups[i].arc_len, groups[i].best_edge.angle, groups[i].sort_key);
	}

	//if(VERBOSE) printf("%d edges\n", num_groups);
	vertex_t *verts = (vertex_t *)malloc_or_die(sizeof(vertex_t) * num_groups);
	for(i=0; i<num_groups; i++) {
		j = i ? i-1 : num_groups-1;
		edge_t e1 = groups[i].best_edge;
		edge_t e2 = groups[j].best_edge;
		line_line_intersection(
			e1.p0, e1.p1,
			e2.p0, e2.p1,
			&verts[i]);
		if(VERBOSE) printf("vert[%d] = %.15f, %.15f\n", i, verts[i].x, verts[i].y);
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
		printf("could not find a 4-sided bounding polygon\n");
	}
	return rect4;
}
