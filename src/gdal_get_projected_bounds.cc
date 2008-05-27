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

void plot_points(ring_t *pl, const char *fn);

void usage(const char *cmdname) {
	printf("Usage: %s [options] \n", cmdname);
	printf("  -s_wkt <fn>           File containing WKT of source region\n");
	printf("  -t_bounds_wkt <fn>    File containing WKT for valid region of target SRS (optional)\n");
	printf("  -s_srs <srs_def>      Source SRS\n");
	printf("  -t_srs <srs_def>      Target SRS\n");
	printf("  -report <out.ppm>     Ouput a graphical report (optional)\n");
	printf("\nOutput is the envelope of the source region projected into the target SRS.\n");
	printf("If the -t_bounds_wkt option is given it will be used as a clip mask in the\n");
	printf("projected space.\n");
	printf("\n");
	
	exit(1);
}

int main(int argc, char **argv) {
	const char *src_wkt_fn = NULL;
	const char *t_bounds_wkt_fn = NULL;
	const char *s_srs = NULL;
	const char *t_srs = NULL;
	const char *report_fn = NULL;

	if(argc == 1) usage(argv[0]);

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(!strcmp(arg, "-v")) {
				VERBOSE++;
			} else if(!strcmp(arg, "-s_wkt")) {
				if(argp == argc) usage(argv[0]);
				src_wkt_fn = argv[argp++];
			} else if(!strcmp(arg, "-t_bounds_wkt")) {
				if(argp == argc) usage(argv[0]);
				t_bounds_wkt_fn = argv[argp++];
			} else if(!strcmp(arg, "-s_srs")) {
				if(argp == argc) usage(argv[0]);
				s_srs = argv[argp++];
			} else if(!strcmp(arg, "-t_srs")) {
				if(argp == argc) usage(argv[0]);
				t_srs = argv[argp++];
			} else if(!strcmp(arg, "-report")) {
				if(argp == argc) usage(argv[0]);
				report_fn = argv[argp++];
			} else usage(argv[0]);
		} else {
			usage(argv[0]);
		}
	}

	if(!src_wkt_fn || !s_srs || !t_srs) usage(argv[0]);

	GDALAllRegister();

	CPLPushErrorHandler(CPLQuietErrorHandler);

	///////////////////////
	
	OGRSpatialReferenceH s_sref = OSRNewSpatialReference(NULL);
	if(OSRImportFromProj4(s_sref, s_srs) != OGRERR_NONE)
		fatal_error("cannot parse proj4 definition for -s_srs");

	OGRSpatialReferenceH t_sref = OSRNewSpatialReference(NULL);
	if(OSRImportFromProj4(t_sref, t_srs) != OGRERR_NONE)
		fatal_error("cannot parse proj4 definition for -t_srs");

	OGRCoordinateTransformationH fwd_xform = 
		OCTNewCoordinateTransformation(s_sref, t_sref);
	OGRCoordinateTransformationH inv_xform = 
		OCTNewCoordinateTransformation(t_sref, s_sref);

	mpoly_t src_mp = mpoly_from_wktfile(src_wkt_fn);
	bbox_t src_bbox = get_polygon_bbox(&src_mp);

	mpoly_t t_bounds_mp;
	int use_t_bounds;
	bbox_t t_bounds_bbox;
	if(t_bounds_wkt_fn) {
		use_t_bounds = 1;
		t_bounds_mp = mpoly_from_wktfile(t_bounds_wkt_fn);
		t_bounds_bbox = get_polygon_bbox(&t_bounds_mp);
	} else {
		use_t_bounds = 0;
	}

	ring_t pl;
	pl.pts = NULL; 
	pl.npts = 0;

	int num_grid_steps = 100;
	for(int grid_xi=0; grid_xi<=num_grid_steps; grid_xi++) {
		vertex_t src_pt;
		double alpha_x = (double)grid_xi / (double)num_grid_steps;
		src_pt.x = src_bbox.min_x + (src_bbox.max_x - src_bbox.min_x) * alpha_x;
		for(int grid_yi=0; grid_yi<=num_grid_steps; grid_yi++) {
			double alpha_y = (double)grid_yi / (double)num_grid_steps;
			src_pt.y = src_bbox.min_y + (src_bbox.max_y - src_bbox.min_y) * alpha_y;
			if(!polygon_contains_point(&src_mp, src_pt.x, src_pt.y)) continue;

			vertex_t tgt_pt = src_pt;
			if(!OCTTransform(fwd_xform, 1, &tgt_pt.x, &tgt_pt.y, NULL)) {
				continue;
			}

			if(!use_t_bounds || polygon_contains_point(&t_bounds_mp, tgt_pt.x, tgt_pt.y)) {
				add_point_to_ring(&pl, tgt_pt);
			}
		}
	}

	double max_step_len = MAX(
		src_bbox.max_x - src_bbox.min_x,
		src_bbox.max_y - src_bbox.min_y) / 1000.0;
	for(int r_idx=0; r_idx<src_mp.num_rings; r_idx++) {
		ring_t *ring = src_mp.rings + r_idx;
		for(int v_idx=0; v_idx<ring->npts; v_idx++) {
			vertex_t v1 = ring->pts[v_idx];
			vertex_t v2 = ring->pts[(v_idx+1) % ring->npts];
			double dx = v2.x - v1.x;
			double dy = v2.y - v1.y;
			double len = sqrt(dx*dx + dy*dy);
			int num_steps = 1 + (int)(len / max_step_len);
			for(int step=0; step<=num_steps; step++) {
				double alpha = (double)step / (double)num_steps;
				vertex_t src_pt;
				src_pt.x = v1.x + dx * alpha;
				src_pt.y = v1.y + dy * alpha;

				vertex_t tgt_pt = src_pt;
				if(!OCTTransform(fwd_xform, 1, &tgt_pt.x, &tgt_pt.y, NULL)) {
					continue;
				}

				if(!use_t_bounds || polygon_contains_point(&t_bounds_mp, tgt_pt.x, tgt_pt.y)) {
					add_point_to_ring(&pl, tgt_pt);
				}
			}
		}
	}

	if(use_t_bounds) {
		double max_step_len = MAX(
			t_bounds_bbox.max_x - t_bounds_bbox.min_x,
			t_bounds_bbox.max_y - t_bounds_bbox.min_y) / 1000.0;
		for(int r_idx=0; r_idx<t_bounds_mp.num_rings; r_idx++) {
			ring_t *ring = t_bounds_mp.rings + r_idx;
			for(int v_idx=0; v_idx<ring->npts; v_idx++) {
				vertex_t v1 = ring->pts[v_idx];
				vertex_t v2 = ring->pts[(v_idx+1) % ring->npts];
				double dx = v2.x - v1.x;
				double dy = v2.y - v1.y;
				double len = sqrt(dx*dx + dy*dy);
				int num_steps = 1 + (int)(len / max_step_len);
				for(int step=0; step<=num_steps; step++) {
					double alpha = (double)step / (double)num_steps;
					vertex_t tgt_pt;
					tgt_pt.x = v1.x + dx * alpha;
					tgt_pt.y = v1.y + dy * alpha;

					vertex_t src_pt = tgt_pt;
					if(!OCTTransform(inv_xform, 1, &src_pt.x, &src_pt.y, NULL)) {
						continue;
					}

					if(polygon_contains_point(&src_mp, src_pt.x, src_pt.y)) {
						add_point_to_ring(&pl, tgt_pt);
					}
				}
			}
		}
	}

	//fprintf(stderr, "got %d points\n", pl.npts);
	bbox_t bbox = get_ring_bbox(&pl);
	printf("bounds:\n");
	printf("  min_e: %.15f\n", bbox.min_x);
	printf("  min_n: %.15f\n", bbox.min_y);
	printf("  max_e: %.15f\n", bbox.max_x);
	printf("  max_n: %.15f\n", bbox.max_y);

	if(report_fn) plot_points(&pl, report_fn);

	return 0;
}

void plot_points(ring_t *pl, const char *fn) {
	bbox_t bbox = get_ring_bbox(pl);
	bbox.min_x -= (bbox.max_x - bbox.min_x) * .05;
	bbox.max_x += (bbox.max_x - bbox.min_x) * .05;
	bbox.min_y -= (bbox.max_y - bbox.min_y) * .05;
	bbox.max_y += (bbox.max_y - bbox.min_y) * .05;
	double W = bbox.max_x - bbox.min_x;
	double H = bbox.max_y - bbox.min_y;
	report_image_t *dbuf = create_plot(W, H);
	for(int i=0; i<pl->npts; i++) {
		vertex_t v = pl->pts[i];
		double x = v.x - bbox.min_x;
		double y = bbox.max_y - v.y;
		plot_point(dbuf, x, y, 255, 255, 255);
	}
	write_plot(dbuf, fn);
}
