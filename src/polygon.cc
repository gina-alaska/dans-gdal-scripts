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



#include <string>

#include "common.h"
#include "polygon.h"
#include "georef.h"

namespace dangdal {

Bbox Ring::getBbox() const {
	Bbox bbox;
	for(size_t i=0; i<pts.size(); i++) {
		bbox.expand(pts[i]);
	}
	return bbox;
}

Bbox Mpoly::getBbox() const {
	Bbox bbox;
	for(size_t i=0; i<rings.size(); i++) {
		bbox.expand(rings[i].getBbox());
	}
	return bbox;
}

std::vector<Bbox> Mpoly::getRingBboxes() const {
	const size_t nrings = rings.size();
	std::vector<Bbox> bboxes(nrings);
	for(size_t i=0; i<nrings; i++) {
		bboxes[i] = rings[i].getBbox();
	}
	return bboxes;
}

void Bbox::expand(const Bbox bb) {
	if(bb.empty) {
		// no-op
	} else if(empty) {
		*this = bb;
	} else {
		expand(Vertex(min_x, min_y));
		expand(Vertex(max_x, max_y));
	}
}

Bbox box_union(const Bbox bb1, const Bbox bb2) {
	Bbox u = bb1;
	u.expand(bb2);
	return u;
}

OGRGeometryH ring_to_ogr(const Ring &ring) {
	OGRGeometryH ogr = OGR_G_CreateGeometry(wkbLinearRing);
	const size_t npts = ring.pts.size();;
	for(size_t i=0; i<npts+1; i++) {
		size_t j_in = (i==npts) ? 0 : i;
		OGR_G_AddPoint_2D(ogr, ring.pts[j_in].x, ring.pts[j_in].y);
	}
	return ogr;
}

Ring ogr_to_ring(OGRGeometryH ogr) {
	//OGRwkbGeometryType type = OGR_G_GetGeometryType(ogr);
	//if(type != wkbLinearRing) {
	//	fatal_error("type != wkbLinearRing in ogr_to_ring (it was %s)",
	//		OGR_G_GetGeometryName(ogr));
	//}
	Ring ring;
	size_t npts = OGR_G_GetPointCount(ogr);
	if(!npts) fatal_error("ring has no points");
	ring.pts.resize(npts);
	for(size_t i=0; i<npts; i++) {
		ring.pts[i].x = OGR_G_GetX(ogr, i);
		ring.pts[i].y = OGR_G_GetY(ogr, i);
	}
	return ring;
}

OGRGeometryH mpoly_to_ogr(const Mpoly mpoly_in) {
	size_t num_rings_in = mpoly_in.rings.size();

	std::vector<std::vector<int> > holes;
	holes.resize(num_rings_in);

	size_t num_geom_out = 0;
	for(size_t outer_idx=0; outer_idx<num_rings_in; outer_idx++) {
		const Ring &ring = mpoly_in.rings[outer_idx];
		if(ring.is_hole) {
			int parent = ring.parent_id;
			holes[parent].push_back(outer_idx);
		} else {
			num_geom_out++;
		}
	}

	bool use_multi = num_geom_out > 1;

	OGRGeometryH geom_out = OGR_G_CreateGeometry(
		use_multi ? wkbMultiPolygon : wkbPolygon);

	for(size_t outer_idx=0; outer_idx<num_rings_in; outer_idx++) {
		const Ring &ring = mpoly_in.rings[outer_idx];
		if(ring.is_hole) continue;

		OGRGeometryH poly_out = use_multi ?
			OGR_G_CreateGeometry(wkbPolygon) : geom_out;
		OGR_G_AddGeometry(poly_out, ring_to_ogr(ring));

		for(size_t hole_idx=0; hole_idx<holes[outer_idx].size(); hole_idx++) {
			const Ring &hole = mpoly_in.rings[holes[outer_idx][hole_idx]];
			if(hole.parent_id != int(outer_idx)) fatal_error("could not sort out holes");
			OGR_G_AddGeometry(poly_out, ring_to_ogr(hole));
		}

		if(size_t(OGR_G_GetGeometryCount(poly_out)) != holes[outer_idx].size()+1) {
			fatal_error("GeometryCount != num_holes+1 (%d vs. %zd)", 
				OGR_G_GetGeometryCount(poly_out), holes[outer_idx].size()+1);
		}

		if(use_multi) OGR_G_AddGeometry(geom_out, poly_out);
	}

	if(use_multi && size_t(OGR_G_GetGeometryCount(geom_out)) != num_geom_out) {
		fatal_error("GeometryCount != num_geom_out (%d vs. %zd)",
			OGR_G_GetGeometryCount(geom_out), num_geom_out);
	}

	return geom_out;
}

Mpoly ogr_to_mpoly(OGRGeometryH geom_in) {
	OGRwkbGeometryType type = OGR_G_GetGeometryType(geom_in);
	if(type == wkbPolygon) {
		Mpoly mpoly_out;
		const size_t num_rings = OGR_G_GetGeometryCount(geom_in);
		if(num_rings < 1) fatal_error("num_rings<1 in ogr_to_mpoly");

		mpoly_out.rings.resize(num_rings);

		for(size_t i=0; i<num_rings; i++) {
			mpoly_out.rings[i] = ogr_to_ring(OGR_G_GetGeometryRef(geom_in, i));
			mpoly_out.rings[i].is_hole = (i>0);
			mpoly_out.rings[i].parent_id = (i>0) ? 0 : -1;
		}

		return mpoly_out;
	} else if(type == wkbMultiPolygon || type == wkbGeometryCollection) {
		size_t num_geom = OGR_G_GetGeometryCount(geom_in);
		std::vector<Mpoly> polys(num_geom);
		size_t total_rings = 0;
		for(size_t i=0; i<num_geom; i++) {
			OGRGeometryH g = OGR_G_GetGeometryRef(geom_in, i);
			polys[i] = ogr_to_mpoly(g);
			total_rings += polys[i].rings.size();
		}
		if(total_rings < 1) fatal_error("num_rings<1 in ogr_to_mpoly");

		Mpoly mpoly_out;
		mpoly_out.rings.resize(total_rings);

		int o = 0;
		for(size_t i=0; i<num_geom; i++) {
			for(size_t j=0; j<polys[i].rings.size(); j++) {
				Ring &ring = polys[i].rings[j];
				if(ring.is_hole) ring.parent_id += o;
				mpoly_out.rings[o+j] = ring;
			}
			o += polys[i].rings.size();
		}
		return mpoly_out;
	} else {
		fatal_error("not a polygon type: %d\n", type);
	}
}

std::vector<Mpoly> split_mpoly_to_polys(const Mpoly &mpoly_in) {
	size_t num_rings_in = mpoly_in.rings.size();

	std::vector<std::vector<int> > holes;
	holes.resize(num_rings_in);

	size_t num_geom_out = 0;
	for(size_t outer_idx=0; outer_idx<num_rings_in; outer_idx++) {
		const Ring &ring = mpoly_in.rings[outer_idx];
		if(ring.is_hole) {
			int parent = ring.parent_id;
			holes[parent].push_back(outer_idx);
		} else {
			num_geom_out++;
		}
	}

	std::vector<Mpoly> polys(num_geom_out);
	size_t poly_out_idx = 0;

	for(size_t outer_idx=0; outer_idx<num_rings_in; outer_idx++) {
		const Ring &ring = mpoly_in.rings[outer_idx];
		if(ring.is_hole) continue;

		Mpoly &out_poly = polys[poly_out_idx];
		out_poly.rings.resize(holes[outer_idx].size()+1);
		size_t ring_out_idx = 0;

		Ring dup_ring(ring);
		dup_ring.parent_id = -1;
		out_poly.rings[ring_out_idx++] = dup_ring;

		for(size_t hole_idx=0; hole_idx<holes[outer_idx].size(); hole_idx++) {
			const Ring &hole = mpoly_in.rings[holes[outer_idx][hole_idx]];
			if(hole.parent_id != int(outer_idx)) fatal_error("could not sort out holes");

			Ring dup_ring(hole);
			dup_ring.parent_id = 0;
			out_poly.rings[ring_out_idx++] = dup_ring;
		}

		poly_out_idx++;
	}

	return polys;
}

/*
void write_mpoly_wkt(char *wkt_fn, mpoly_t mpoly, int split_polys) {
	int num_rings = mpoly.num_rings;
	ring_t *rings = mpoly.rings;

	FILE *fout = fopen(wkt_fn, "w");
	if(!fout) fatal_error("cannot open output file for WKT");

	if(!split_polys) fprintf(fout, "MULTIPOLYGON(\n");

	int r_idx, h_idx, p_idx;
	int is_first_ring = 1;
	for(r_idx=0; r_idx<num_rings; r_idx++) {
		ring_t *ring = rings + r_idx;
		if(ring->is_hole) continue;

		if(!is_first_ring) fprintf(fout, split_polys ? "\n\n" : ", ");
		is_first_ring = 0;

		fprintf(fout, split_polys ? "POLYGON((\n" : "((\n");

		//fprintf(fout, "  ring:%d\n", r_idx);
		for(p_idx=0; p_idx<ring->npts+1; p_idx++) {
			if(!(p_idx%4)) {
				if(p_idx) fprintf(fout, "\n");
				fprintf(fout, "  ");
			}
			vertex_t v = ring->pts[p_idx % ring->npts];
			fprintf(fout, "%.15f %.15f", v.x, v.y);
			if(p_idx < ring->npts) fprintf(fout, ", ");
		}
		fprintf(fout, "\n)");

		for(h_idx=0; h_idx<num_rings; h_idx++) {
			ring_t *hole = rings + h_idx;
			if(hole->parent_id != r_idx) continue;

			fprintf(fout, ", (\n");
			//fprintf(fout, "  hole:%d\n", h_idx);
			for(p_idx=0; p_idx<hole->npts+1; p_idx++) {
				if(!(p_idx%4)) {
					if(p_idx) fprintf(fout, "\n");
					fprintf(fout, "  ");
				}
				vertex_t v = hole->pts[p_idx % hole->npts];
				fprintf(fout, "%.15f %.15f", v.x, v.y);
				if(p_idx < hole->npts) fprintf(fout, ", ");
			}
			fprintf(fout, "\n)");
		}
		fprintf(fout, ")");
	}

	if(!split_polys) fprintf(fout, ")\n");

	fclose(fout);
}
*/

bool line_intersects_line(
	Vertex p1, Vertex p2,
	Vertex p3, Vertex p4,
	bool fail_on_coincident
) {
	if(
		std::max(p1.x, p2.x) < std::min(p3.x, p4.x) ||
		std::min(p1.x, p2.x) > std::max(p3.x, p4.x) ||
		std::max(p1.y, p2.y) < std::min(p3.y, p4.y) ||
		std::min(p1.y, p2.y) > std::max(p3.y, p4.y)
	) return 0;
	double numer_a = (p4.x-p3.x)*(p1.y-p3.y) - (p4.y-p3.y)*(p1.x-p3.x);
	double numer_b = (p2.x-p1.x)*(p1.y-p3.y) - (p2.y-p1.y)*(p1.x-p3.x);
	double denom   = (p4.y-p3.y)*(p2.x-p1.x) - (p4.x-p3.x)*(p2.y-p1.y);
	if(denom == 0) {
		if(numer_a==0 && numer_b==0) { // coincident
			if(fail_on_coincident) {
				return 0;
			} else {
				// lines must touch because of min/max test above
				return 1;
			}
		} else { // parallel
			return 0;
		}
	}
	double ua = (double)numer_a / (double)denom;
	double ub = (double)numer_b / (double)denom;
	return ua>=0 && ua<=1 && ub>=0 && ub<=1;
}

Vertex line_line_intersection(
	Vertex p1, Vertex p2,
	Vertex p3, Vertex p4
) {
	double numer_a = (p4.x-p3.x)*(p1.y-p3.y) - (p4.y-p3.y)*(p1.x-p3.x);
	//double numer_b = (p2.x-p1.x)*(p1.y-p3.y) - (p2.y-p1.y)*(p1.x-p3.x);
	double denom   = (p4.y-p3.y)*(p2.x-p1.x) - (p4.x-p3.x)*(p2.y-p1.y);
	if(!denom) fatal_error("lines are parallel");
	double ua = numer_a / denom;
	//double ub = numer_b / denom;
	double x = p1.x + ua*(p2.x-p1.x);
	double y = p1.y + ua*(p2.y-p1.y);
	return Vertex(x, y);
}

double Ring::orientedArea() const {
	double accum = 0;
	const size_t npts = pts.size();
	for(size_t i=0; i<npts; i++) {
		size_t j = (i==npts-1) ? 0 : (i+1);
		double x0 = pts[i].x;
		double y0 = pts[i].y;
		double x1 = pts[j].x;
		double y1 = pts[j].y;
		accum += x0*y1 - x1*y0;
	}
	return accum / 2.0;
}

bool Ring::isCCW() const {
	return orientedArea() > 0;
}

double Ring::area() const {
	return fabs(orientedArea());
}

bool Ring::contains(Vertex p) const {
	const double px = p.x;
	const double py = p.y;

	const size_t npts = pts.size();
	int num_crossings = 0;
	for(size_t i=0; i<npts; i++) {
		size_t i2 = (i==npts-1) ? 0 : (i+1);
		double x0 = pts[i].x;
		double y0 = pts[i].y;
		double x1 = pts[i2].x;
		double y1 = pts[i2].y;
		// we want to know whether a ray from (px,py) in the (1,0) direction
		// passes through this segment
		if(x0 < px && x1 < px) continue;

		int y0above = y0 >= py;
		int y1above = y1 >= py;
		if(y0above && y1above) continue;
		if(!y0above && !y1above) continue;

		double alpha = (py-y0)/(y1-y0);
		double cx = x0 + (x1-x0)*alpha;
		//printf("x0=%.15f x1=%.15f\n", x0, x1);
		//printf("y0=%.15f y1=%.15f\n", y0, y1);
		//printf("alpha=%.15f cx=%.15f px=%.15f\n", alpha, cx, px);
		if(cx > px) num_crossings++;
	}
	//printf("num_crossings=%d\n", num_crossings);
	// if there are an odd number of crossings, then (px,py) is within c1
	return num_crossings & 1;
}

bool Mpoly::contains(Vertex p) const {
	int num_crossings = 0;
	for(size_t r_idx=0; r_idx<rings.size(); r_idx++) {
		if(rings[r_idx].contains(p)) num_crossings++;
	}
	// if it is within an odd number of rings it is not in a hole
	return num_crossings & 1;
}

void Mpoly::deleteRing(size_t idx) {
	rings.erase(rings.begin() + idx);
}

RingRelation ring_ring_relation(const Ring &r1, const Ring &r2) {
	Bbox bb1 = r1.getBbox();
	Bbox bb2 = r2.getBbox();
	if(is_disjoint(bb1, bb2)) return RINGREL_DISJOINT;

	size_t n1 = r1.pts.size();
	size_t n2 = r2.pts.size();
	if(n1==0 || n2==0) return RINGREL_DISJOINT;

	// test for crossings
	for(size_t i1=0; i1<n1; i1++) {
		size_t i1_plus1 = (i1==n1-1) ? 0 : i1+1;
		Vertex p1a = r1.pts[i1];
		Vertex p1b = r1.pts[i1_plus1];
		if(
			(p1a.x < bb2.min_x && p1b.x < bb2.min_x) ||
			(p1a.y < bb2.min_y && p1b.y < bb2.min_y) ||
			(p1a.x > bb2.max_x && p1b.x > bb2.max_x) ||
			(p1a.y > bb2.max_y && p1b.y > bb2.max_y)
		) continue;
		for(size_t i2=0; i2<n2; i2++) {
			size_t i2_plus1 = (i2==n2-1) ? 0 : i2+1;
			Vertex p2a = r2.pts[i2];
			Vertex p2b = r2.pts[i2_plus1];
			if(line_intersects_line(p1a, p1b, p2a, p2b, false)) {
				return RINGREL_CROSSES;
			}
		}
	}

	if(r1.contains(r2.pts[0])) return RINGREL_CONTAINS;
	else if(r2.contains(r1.pts[0])) return RINGREL_CONTAINED_BY;
	else return RINGREL_DISJOINT;
}

/*
void compute_containments(mpoly_t *mp) {
	int i;
	int nrings = mp->num_rings;

	bbox_t *bboxes = make_bboxes(mp);

	int **ancestors = MYALLOC(int *, nrings);
	int *num_ancestors = MYALLOC(int, nrings);

	long num_hits = 0;

	printf("Computing containments: ");
	for(i=0; i<nrings; i++) {
		GDALTermProgress((double)i/(double)nrings, NULL, NULL);

		bbox_t bbox1 = bboxes[i];
		if(bbox1.empty) continue;

		num_ancestors[i] = 0;
		ancestors[i] = NULL;

		int j;
		for(j=0; j<nrings; j++) {
			bbox_t bbox2 = bboxes[j];
			if(bboxes_disjoint(&bbox1, &bbox2)) continue;

			int contains = (i != j) && ring_contains_ring(
				mp->rings + j, mp->rings + i);
			if(contains) {
				if(VERBOSE) printf("%d contains %d\n", j, i);
				ancestors[i] = REMYALLOC(int, ancestors[i], (num_ancestors[i]+1));
				ancestors[i][ num_ancestors[i]++ ] = j;
				num_hits++;
			}
		}
	}
	GDALTermProgress(1, NULL, NULL);

	free(bboxes);

	if(VERBOSE) {
		printf("containment hits = %ld/%ld\n", num_hits, (long)nrings*(long)nrings);
		if(VERBOSE >= 2) {
			for(i=0; i<nrings; i++) {
				printf("num_ancestors[%d] = %d\n", i, num_ancestors[i]);
			}
		}
	}

	for(i=0; i<nrings; i++) {
		// only odd levels are holes
		mp->rings[i].is_hole = num_ancestors[i] % 2;
	}

	for(i=0; i<nrings; i++) {
		mp->rings[i].parent_id = -1;

		int a;
		for(a=0; a<num_ancestors[i]; a++) {
			int j = ancestors[i][a];
			if(num_ancestors[i] == num_ancestors[j] + 1) {
				mp->rings[i].parent_id = j;
			}
		}
	}

	for(i=0; i<nrings; i++) {
		if(ancestors[i]) free(ancestors[i]);
	}
	free(ancestors);
	free(num_ancestors);
}
*/

void Mpoly::xy2en(const GeoRef &georef) {
	for(size_t r_idx=0; r_idx<rings.size(); r_idx++) {
		Ring &ring = rings[r_idx];
		for(size_t v_idx=0; v_idx<ring.pts.size(); v_idx++) {
			double x = ring.pts[v_idx].x;
			double y = ring.pts[v_idx].y;
			double east, north;
			georef.xy2en(x, y, &east, &north);
			ring.pts[v_idx].x = east;
			ring.pts[v_idx].y = north;
		}
	}
}

void Mpoly::en2xy(const GeoRef &georef) {
	for(size_t r_idx=0; r_idx<rings.size(); r_idx++) {
		Ring &ring = rings[r_idx];
		for(size_t v_idx=0; v_idx<ring.pts.size(); v_idx++) {
			double east  = ring.pts[v_idx].x;
			double north = ring.pts[v_idx].y;
			double x, y;
			georef.en2xy(east, north, &x, &y);
			ring.pts[v_idx].x = x;
			ring.pts[v_idx].y = y;
		}
	}
}

// This way is very simple and robust, but doesn't know anything
// about the actual distortion of the projection.
/*
mpoly_t *mpoly_xy2ll_with_interp(const GeoRef &georef, mpoly_t *xy_poly, double toler_pixels) {
	mpoly_t *ll_poly = MYALLOC(mpoly_t, 1);
	ll_poly->num_rings = xy_poly->num_rings;
	ll_poly->rings = MYALLOC(ring_t, ll_poly->num_rings);

	OGRErr err = OGRERR_NONE;
	double earth_radius = OSRGetSemiMajor(georef.spatial_ref, &err);
	if(err != OGRERR_NONE) fatal_error("could not determine globe radius");

	double toler_radians = toler_pixels * 
		std::min(georef.res_meters_x, georef.res_meters_y) / earth_radius;
	// error is (approximately) proportional to segment length squared
	// FIXME - need to multiply this by some constant
	double max_seg_len = sqrt(toler_radians);

	int total_midpoints = 0;
	
	int r_idx;
	for(r_idx=0; r_idx<xy_poly->num_rings; r_idx++) {
		ring_t *xy_ring = xy_poly->rings + r_idx;
		int npts_in = xy_ring->npts;
		// this will be the output
		ring_t ll_ring = *xy_ring;
		ll_ring.npts = 0;
		ll_ring.pts = NULL;

		int v_idx;
		for(v_idx=0; v_idx<npts_in; v_idx++) {
			vertex_t xy1 = xy_ring->pts[v_idx];
			vertex_t xy2 = xy_ring->pts[(v_idx+1) % npts_in];

			// compute segment length in radians
			double dx = (xy1.x - xy2.x) * georef.res_meters_x / earth_radius;
			double dy = (xy1.y - xy2.y) * georef.res_meters_y / earth_radius;
			double seg_len = sqrt(dx*dx + dy*dy);

			int nmid = (int)floor(seg_len / max_seg_len);
			total_midpoints += nmid;
			for(int i=0; i<=nmid; i++) {
				double alpha = (double)i / (double)(nmid+1);
				double x = xy1.x + alpha * (xy2.x - xy1.x);
				double y = xy1.y + alpha * (xy2.y - xy1.y);
				vertex_t ll;
				xy2ll_or_die(georef, x, y, &ll.x, &ll.y);
				add_point_to_ring(&ll_ring, ll);
			}
		}

		ll_poly->rings[r_idx] = ll_ring;
	}

	if(VERBOSE) {
		printf("Added %d interpolation midpoints\n", total_midpoints);
	}

	return ll_poly;
}
*/

// returns length of longest raster side in meters, squared
static double estimate_canvas_size_sq(const GeoRef &georef, double semi_major) {
	double e1, n1, e2, n2;
	georef.xy2en(0, 0, &e1, &n1);
	georef.xy2en(georef.w, 0, &e2, &n2);
	double dx = (e1 - e2);
	double dy = (n1 - n2);
	double size1 = dx*dx + dy*dy;
	georef.xy2en(0, georef.h, &e2, &n2);
	dx = (e1 - e2);
	dy = (n1 - n2);
	double size2 = dx*dx + dy*dy;

	double size = std::max(size1, size2);
	size *= georef.units_val * georef.units_val;
	if(OSRIsGeographic(georef.spatial_ref)) {
		size *= semi_major * semi_major;
	}
	return size;
}

void Mpoly::xy2ll_with_interp(const GeoRef &georef, double toler) {
	size_t nrings = rings.size();
	Mpoly ll_poly;
	ll_poly.rings.resize(nrings);
	
	double semi_major;
	if(georef.have_semi_major) {
		semi_major = georef.semi_major;
	} else {
		semi_major = 6370997.0;
		fprintf(stderr, "Warning: could not get globe size, assuming %lf\n", semi_major);
	}
	double canvas_size_sq = estimate_canvas_size_sq(georef, semi_major);
	//printf("canvas_size_sq = %lf\n", canvas_size_sq);

	// FIXME - now that we don't use ll2xy this is probably not needed:
	// This is a kludge that shrinks the canvas by a millionth of a pixel
	// to avoid problems with images that span an entire 360 degrees of
	// longitude.  Without this the map xy -> ll -> xy is not single-valued.
	double epsilon = 5e-7;
	double shrink = ((double)georef.w - 2.0*epsilon) / (double)georef.w;

	for(size_t r_idx=0; r_idx<nrings; r_idx++) {
		// note that this is *not* a reference, since it will potentially
		// be modified
		Ring xy_ring = rings[r_idx];
		// this will be the output
		Ring &ll_ring = ll_poly.rings[r_idx];
		ll_ring = xy_ring.copyMetadata();
		ll_ring.pts.resize(xy_ring.pts.size());

		for(size_t v_idx=0; v_idx<ll_ring.pts.size(); v_idx++) {
			double x = xy_ring.pts[v_idx].x;
			double y = xy_ring.pts[v_idx].y;
			double lon, lat;
			georef.xy2ll_or_die(x*shrink+epsilon, y, &lon, &lat);
			ll_ring.pts[v_idx].x = lon;
			ll_ring.pts[v_idx].y = lat;
		}

		int num_consec = 0;

		for(size_t v_idx=0; v_idx<ll_ring.pts.size(); ) {
			if(xy_ring.pts.size() != ll_ring.pts.size()) {
				fatal_error("xy_ring.npts != ll_ring.npts");
			}
			size_t npts = xy_ring.pts.size();

			const Vertex xy1 = xy_ring.pts[v_idx];
			const Vertex xy2 = xy_ring.pts[(v_idx + 1) % npts];
			const Vertex xy_m(
				(xy1.x + xy2.x)/2.0,
				(xy1.y + xy2.y)/2.0);

			Vertex ll_m_proj;
			georef.xy2ll_or_die(
				xy_m.x*shrink+epsilon, xy_m.y,
				&ll_m_proj.x, &ll_m_proj.y);

			Vertex &ll1 = ll_ring.pts[v_idx];
			Vertex &ll2 = ll_ring.pts[(v_idx + 1) % npts];

			while(ll1.x - ll2.x >  180) ll1.x -= 360;
			while(ll1.x - ll2.x < -180) ll1.x += 360;

			bool need_midpt = 0;

			// avoid topological errors by not allowing segments
			// to be longer than 90 degrees
			if(!need_midpt) {
				double dx = fabs(ll1.x - ll2.x);
				double dy = ll1.y - ll2.y;
				double sqr_error = dx*dx + dy*dy;
				double max_error = 90; // degrees

				need_midpt = sqr_error > max_error*max_error;

				if(VERBOSE && need_midpt) {
					printf("  inserting midpoint at %zd,%zd (lon/lat difference %g > %g)\n",
						r_idx, v_idx, sqrt(sqr_error), max_error);
				}
			}

			if(!need_midpt) {
				Vertex ll_m_interp(
					(ll1.x + ll2.x)/2.0,
					(ll1.y + ll2.y)/2.0);

				while(ll_m_interp.x - ll_m_proj.x >  180) ll_m_interp.x -= 360;
				while(ll_m_interp.x - ll_m_proj.x < -180) ll_m_interp.x += 360;
				double lonscale = cos(D2R * 
					std::max(fabs(ll_m_interp.y), fabs(ll_m_proj.y)));
				double dx = (ll_m_interp.x - ll_m_proj.x) * lonscale;
				double dy = ll_m_interp.y - ll_m_proj.y;
				dx *= D2R * semi_major;
				dy *= D2R * semi_major;
				double sqr_error = dx*dx + dy*dy;

				// if the midpoint is this far off then something is seriously wrong
				if(sqr_error > canvas_size_sq) {
					fprintf(stderr, "\nInfo on bad point:\n");
					fprintf(stderr, "xy1 = %lf,%lf\n", xy1.x, xy1.y);
					fprintf(stderr, "xy2 = %lf,%lf\n", xy2.x, xy2.y);
					fprintf(stderr, "xy_m = %lf,%lf\n", xy_m.x, xy_m.y);
					fprintf(stderr, "ll1 = %lf,%lf\n", ll1.x, ll1.y);
					fprintf(stderr, "ll2 = %lf,%lf\n", ll2.x, ll2.y);
					fprintf(stderr, "ll_m_proj = %lf,%lf\n", ll_m_proj.x, ll_m_proj.y);
					fprintf(stderr, "ll_m_interp = %lf,%lf\n", ll_m_interp.x, ll_m_interp.y);
					fprintf(stderr, "error = %lf > %lf\n", sqr_error, canvas_size_sq);
					fatal_error("projection error in mpoly_xy2ll_with_interp");
				}

				need_midpt = toler && sqr_error > toler*toler;

				if(VERBOSE && need_midpt) {
					printf("  inserting midpoint at %zd,%zd (delta=%lf,%lf > %lf)\n",
						r_idx, v_idx, dx, dy, toler);
					printf("    testll=[%lf,%lf] projll=[%lf,%lf]\n", 
						ll_m_interp.x, ll_m_interp.y, ll_m_proj.x, ll_m_proj.y);
				}
			}

			if(need_midpt) {
				if(num_consec++ > 10) fatal_error("convergence error in mpoly_xy2ll_with_interp");

				if(VERBOSE) {
					printf("    xy=[%lf,%lf]:[%lf,%lf]\n", xy1.x, xy1.y, xy2.x, xy2.y);
					printf("    ll=[%lf,%lf]:[%lf,%lf]\n", ll1.x, ll1.y, ll2.x, ll2.y);
					printf("    midxy=[%lf,%lf] midll=[%lf,%lf]\n", 
						xy_m.x, xy_m.y, ll_m_proj.x, ll_m_proj.y);
				}
				xy_ring.pts.insert(xy_ring.pts.begin()+v_idx+1, xy_m);
				ll_ring.pts.insert(ll_ring.pts.begin()+v_idx+1, ll_m_proj);
			} else {
				v_idx++;
				num_consec = 0;
			}
		}

		// Now reproject everything again, without using the shrink
		// kludge.  We no longer care whether the projection is
		// single-valued and we don't want the loss of accuracy
		// that comes from multiplying by shrink.
		for(size_t v_idx=0; v_idx<ll_ring.pts.size(); v_idx++) {
			double x = xy_ring.pts[v_idx].x;
			double y = xy_ring.pts[v_idx].y;

			double lon, lat;
			georef.xy2ll_or_die(x, y, &lon, &lat);

			ll_ring.pts[v_idx].x = lon;
			ll_ring.pts[v_idx].y = lat;
		}
	}

	*this = ll_poly;
}

static std::string read_whole_file(FILE *fin) {
	size_t chunk_size = 65536;
	std::string accum;
	std::string chunk(chunk_size, 0);
	size_t num_read;
	while(0 < (num_read = fread(&chunk[0], 1, chunk_size, fin))) {
		accum.append(chunk, 0, num_read);
	}
	return accum;
}

Mpoly mpoly_from_wktfile(const char *fn) {
	FILE *fh = fopen(fn, "r");
	if(!fh) fatal_error("cannot read file [%s]", fn);
	std::string wkt = read_whole_file(fh);
	for(size_t i=0; i<wkt.size(); i++) {
		if(wkt[i] == '\r') wkt[i] = ' ';
		if(wkt[i] == '\n') wkt[i] = ' ';
		if(wkt[i] == '\t') wkt[i] = ' ';
	}

	OGRGeometryH geom;
	char *wkt_str = strdup(wkt.c_str());
	char *wkt_str2 = wkt_str;
	OGRErr err = OGR_G_CreateFromWkt(&wkt_str2, NULL, &geom);
	free(wkt_str);
	if(OGRERR_NONE != err) {
		fatal_error("OGR_G_CreateFromWkt failed: %d", err);
	}
	
	return ogr_to_mpoly(geom);
}

} // namespace dangdal
