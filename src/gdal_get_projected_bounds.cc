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
#include "debugplot.h"
#include "polygon.h"

using namespace dangdal;

void plot_points(const Ring &pl, const char *fn);

struct PointStats {
	PointStats() : 
		total(0),
		proj_ok(0),
		contained(0)
	{ }

	void printYaml(const char *label) {
		printf("%s:\n", label);
		printf("  total: %zd\n", total);
		printf("  proj_ok: %zd\n", proj_ok);
		printf("  contained: %zd\n", contained);
	}

	size_t total;
	size_t proj_ok;
	size_t contained;
};

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

// This function transforms a point, and then as a check transforms it back to
// see if it comes back to the same place.  This allows us to detect cases
// where OCTTransform reports success when really it just returned some
// meaningless result.
bool picky_transform(
	OGRCoordinateTransformationH fwd_xform,
	OGRCoordinateTransformationH inv_xform,
	Vertex *v_in
) {
	// tolerance in meters, could probably be much smaller
	const double toler = 1.0;

	Vertex v_out = *v_in;
	if(!OCTTransform(fwd_xform, 1, &v_out.x, &v_out.y, NULL)) {
		return 0;
	}

	Vertex v_back = v_out;
	if(!OCTTransform(inv_xform, 1, &v_back.x, &v_back.y, NULL)) {
		return 0;
	}

	double err = hypot(v_in->x - v_back.x, v_in->y - v_back.y);
	//fprintf(stderr, "err=%g\n", err);
	if(err > toler) {
		return 0;
	}

	*v_in = v_out;
	return 1;
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

	Mpoly src_mp = mpoly_from_wktfile(src_wkt_fn);
	Bbox src_bbox = src_mp.getBbox();

	Mpoly t_bounds_mp;
	bool use_t_bounds;
	Bbox t_bounds_bbox;
	if(t_bounds_wkt_fn) {
		use_t_bounds = 1;
		t_bounds_mp = mpoly_from_wktfile(t_bounds_wkt_fn);
		t_bounds_bbox = t_bounds_mp.getBbox();
	} else {
		use_t_bounds = 0;
	}

	Ring pl;

	PointStats ps_border;
	PointStats ps_interior;
	PointStats ps_bounds;

	// Sample a regular grid of points, take the ones within the source region,
	// and project them to the target projection.  This is done to handle the
	// cases where the projected border does not necessarily encircle the
	// source region (such as would be the case for a source region that
	// encircles the pole with a target lonlat projection).
	int num_grid_steps = 100;
	for(int grid_xi=0; grid_xi<=num_grid_steps; grid_xi++) {
		Vertex src_pt;
		double alpha_x = (double)grid_xi / (double)num_grid_steps;
		src_pt.x = src_bbox.min_x + (src_bbox.max_x - src_bbox.min_x) * alpha_x;
		for(int grid_yi=0; grid_yi<=num_grid_steps; grid_yi++) {
			double alpha_y = (double)grid_yi / (double)num_grid_steps;
			src_pt.y = src_bbox.min_y + (src_bbox.max_y - src_bbox.min_y) * alpha_y;
			if(!src_mp.contains(src_pt)) continue;

			ps_interior.total++;

			Vertex tgt_pt = src_pt;
			if(!picky_transform(fwd_xform, inv_xform, &tgt_pt)) {
				continue;
			}

			ps_interior.proj_ok++;

			if(!use_t_bounds || t_bounds_mp.contains(tgt_pt)) {
				ps_interior.contained++;
				pl.pts.push_back(tgt_pt);
			}
		}
	}

	// Project points along the source region border to the target projection.
	double max_step_len = MAX(
		src_bbox.max_x - src_bbox.min_x,
		src_bbox.max_y - src_bbox.min_y) / 1000.0;
	for(size_t r_idx=0; r_idx<src_mp.rings.size(); r_idx++) {
		const Ring &ring = src_mp.rings[r_idx];
		for(size_t v_idx=0; v_idx<ring.pts.size(); v_idx++) {
			Vertex v1 = ring.pts[v_idx];
			Vertex v2 = ring.pts[(v_idx+1) % ring.pts.size()];
			double dx = v2.x - v1.x;
			double dy = v2.y - v1.y;
			double len = sqrt(dx*dx + dy*dy);
			int num_steps = 1 + (int)(len / max_step_len);
			for(int step=0; step<=num_steps; step++) {
				double alpha = (double)step / (double)num_steps;
				Vertex src_pt;
				src_pt.x = v1.x + dx * alpha;
				src_pt.y = v1.y + dy * alpha;

				ps_border.total++;

				Vertex tgt_pt = src_pt;
				if(!picky_transform(fwd_xform, inv_xform, &tgt_pt)) {
					continue;
				}

				ps_border.proj_ok++;

				if(!use_t_bounds || t_bounds_mp.contains(tgt_pt)) {
					ps_border.contained++;
					pl.pts.push_back(tgt_pt);
				}
			}
		}
	}

	// Take points along the border of the t_bounds clip shape that lie within the
	// source region.
	if(use_t_bounds) {
		double max_step_len = MAX(
			t_bounds_bbox.max_x - t_bounds_bbox.min_x,
			t_bounds_bbox.max_y - t_bounds_bbox.min_y) / 1000.0;
		for(size_t r_idx=0; r_idx<t_bounds_mp.rings.size(); r_idx++) {
			const Ring &ring = t_bounds_mp.rings[r_idx];
			for(size_t v_idx=0; v_idx<ring.pts.size(); v_idx++) {
				Vertex v1 = ring.pts[v_idx];
				Vertex v2 = ring.pts[(v_idx+1) % ring.pts.size()];
				double dx = v2.x - v1.x;
				double dy = v2.y - v1.y;
				double len = sqrt(dx*dx + dy*dy);
				int num_steps = 1 + (int)(len / max_step_len);
				for(int step=0; step<=num_steps; step++) {
					double alpha = (double)step / (double)num_steps;
					Vertex tgt_pt;
					tgt_pt.x = v1.x + dx * alpha;
					tgt_pt.y = v1.y + dy * alpha;

					ps_bounds.total++;

					Vertex src_pt = tgt_pt;
					if(!picky_transform(inv_xform, fwd_xform, &src_pt)) {
						continue;
					}

					ps_bounds.proj_ok++;

					if(src_mp.contains(src_pt)) {
						ps_bounds.contained++;
						pl.pts.push_back(tgt_pt);
					}
				}
			}
		}
	}

	//bool debug = 1;
	//if(debug) {
	//	ps_border.printYaml("stats_border");
	//	ps_interior.printYaml("stats_interior");
	//	ps_bounds.printYaml("stats_bounds");
	//}

	//fprintf(stderr, "got %zd points\n", pl.npts);
	Bbox bbox = pl.getBbox();
	printf("bounds:\n");
	printf("  min_e: %.15f\n", bbox.min_x);
	printf("  min_n: %.15f\n", bbox.min_y);
	printf("  max_e: %.15f\n", bbox.max_x);
	printf("  max_n: %.15f\n", bbox.max_y);

	if(report_fn) plot_points(pl, report_fn);

	return 0;
}

void plot_points(const Ring &pl, const char *fn) {
	Bbox bbox = pl.getBbox();
	bbox.min_x -= (bbox.max_x - bbox.min_x) * .05;
	bbox.max_x += (bbox.max_x - bbox.min_x) * .05;
	bbox.min_y -= (bbox.max_y - bbox.min_y) * .05;
	bbox.max_y += (bbox.max_y - bbox.min_y) * .05;
	double W = bbox.max_x - bbox.min_x;
	double H = bbox.max_y - bbox.min_y;
	report_image_t *dbuf = create_plot(W, H);
	for(size_t i=0; i<pl.pts.size(); i++) {
		Vertex v = pl.pts[i];
		double x = v.x - bbox.min_x;
		double y = bbox.max_y - v.y;
		plot_point(dbuf, x, y, 255, 255, 255);
	}
	write_plot(dbuf, fn);
}
