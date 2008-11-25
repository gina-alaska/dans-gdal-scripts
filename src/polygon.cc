/*
Copyright (c) 2008, Regents of the University of Alaska

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
#include "georef.h"

mpoly_t empty_polygon() {
	mpoly_t empty;
	empty.num_rings = 0;
	empty.rings = NULL;
	return empty;
}

ring_t duplicate_ring(ring_t *in_ring) {
	ring_t out_ring = *in_ring; // copy metadata
	out_ring.pts = (vertex_t *)malloc_or_die(sizeof(vertex_t) * out_ring.npts);
	memcpy(out_ring.pts, in_ring->pts, sizeof(vertex_t) * out_ring.npts);
	return out_ring;
}

void free_ring(ring_t *ring) {
	if(ring->pts) free(ring->pts);
}

void free_mpoly(mpoly_t *mpoly) {
	int i;
	for(i=0; i<mpoly->num_rings; i++) {
		free_ring(mpoly->rings + i);
	}
	if(mpoly->rings) free(mpoly->rings);
}

void insert_point_into_ring(ring_t *ring, int idx) {
	if(idx<0 || idx>ring->npts) fatal_error("idx out of range in insert_point_into_ring");
	ring->pts = (vertex_t *)realloc_or_die(ring->pts, sizeof(vertex_t) * (ring->npts+1));
	memmove(ring->pts + idx + 1, ring->pts + idx, sizeof(vertex_t) * (ring->npts - idx));
	ring->npts++;
}

void add_point_to_ring(ring_t *ring, vertex_t v) {
	ring->pts = (vertex_t *)realloc_or_die(ring->pts, sizeof(vertex_t) * (ring->npts+1));
	ring->pts[ring->npts] = v;
	ring->npts++;
}

void delete_ring_from_mpoly(mpoly_t *mp, int idx) {
	free_ring(&mp->rings[idx]);
	memmove(&mp->rings[idx], &mp->rings[idx+1], 
		sizeof(ring_t) * (mp->num_rings-1-idx));
	mp->num_rings--;
}

bbox_t get_ring_bbox(ring_t *ring) {
	bbox_t bbox;

	if(ring->npts == 0) {
		bbox.min_x = bbox.max_x = 0;
		bbox.min_y = bbox.max_y = 0;
		bbox.empty = 1;
	} else {
		bbox.empty = 0;

		bbox.min_x = bbox.max_x = ring->pts[0].x;
		bbox.min_y = bbox.max_y = ring->pts[0].y;
		int i;
		for(i=0; i<ring->npts; i++) {
			vertex_t v = ring->pts[i];
			if(v.x < bbox.min_x) bbox.min_x = v.x;
			if(v.y < bbox.min_y) bbox.min_y = v.y;
			if(v.x > bbox.max_x) bbox.max_x = v.x;
			if(v.y > bbox.max_y) bbox.max_y = v.y;
		}
	}
	return bbox;
}

bbox_t get_polygon_bbox(mpoly_t *mp) {
	bbox_t bbox;
	bbox.min_x = bbox.max_x = 0;
	bbox.min_y = bbox.max_y = 0;
	bbox.empty = 1;
	for(int i=0; i<mp->num_rings; i++) {
		bbox_t rbb = get_ring_bbox(mp->rings + i);
		bbox = union_bbox(bbox, rbb);
	}
	return bbox;
}

bbox_t *make_bboxes(mpoly_t *mp) {
	int nrings = mp->num_rings;
	bbox_t *bboxes = (bbox_t *)malloc_or_die(sizeof(bbox_t) * nrings);
	int i;
	for(i=0; i<nrings; i++) {
		bboxes[i] = get_ring_bbox(mp->rings + i);
	}
	return bboxes;
}

bbox_t union_bbox(bbox_t bb1, bbox_t bb2) {
	bbox_t bb;
	if(bb1.empty) {
		// bb2 may be empty also...
		bb = bb2;
	} else if(bb2.empty) {
		bb = bb1;
	} else {
		bb.empty = 0;
		bb.min_x = MIN(bb1.min_x, bb2.min_x);
		bb.min_y = MIN(bb1.min_y, bb2.min_y);
		bb.max_x = MAX(bb1.max_x, bb2.max_x);
		bb.max_y = MAX(bb1.max_y, bb2.max_y);
	}
	return bb;
}

int bboxes_disjoint(bbox_t *bbox1, bbox_t *bbox2) {
	return
		bbox1->empty || bbox2->empty ||
		bbox1->min_x > bbox2->max_x ||
		bbox1->min_y > bbox2->max_y ||
		bbox2->min_x > bbox1->max_x ||
		bbox2->min_y > bbox1->max_y;
}

OGRGeometryH ring_to_ogr(ring_t *ring) {
	OGRGeometryH ogr = OGR_G_CreateGeometry(wkbLinearRing);
	int i;
	for(i=0; i<ring->npts+1; i++) {
		int j_in = i==ring->npts ? 0 : i;
		OGR_G_AddPoint_2D(ogr, ring->pts[j_in].x, ring->pts[j_in].y);
	}
	return ogr;
}

ring_t ogr_to_ring(OGRGeometryH ogr) {
	//OGRwkbGeometryType type = OGR_G_GetGeometryType(ogr);
	//if(type != wkbLinearRing) {
	//	fatal_error("type != wkbLinearRing in ogr_to_ring (it was %s)",
	//		OGR_G_GetGeometryName(ogr));
	//}
	ring_t ring;
	ring.npts = OGR_G_GetPointCount(ogr);
	if(!ring.npts) fatal_error("ring has no points");
	ring.pts = (vertex_t *)malloc_or_die(sizeof(vertex_t) * ring.npts);
	for(int i=0; i<ring.npts; i++) {
		ring.pts[i].x = OGR_G_GetX(ogr, i);
		ring.pts[i].y = OGR_G_GetY(ogr, i);
	}
	return ring;
}

OGRGeometryH mpoly_to_ogr(mpoly_t *mpoly_in) {
	int num_rings_in = mpoly_in->num_rings;

	int *num_holes = (int *)malloc_or_die(sizeof(int) * num_rings_in);
	int **holes = (int **)malloc_or_die(sizeof(int *) * num_rings_in);
	for(int outer_idx=0; outer_idx<num_rings_in; outer_idx++) {
		num_holes[outer_idx] = 0;
		holes[outer_idx] = NULL;
	}

	int num_geom_out = 0;
	for(int outer_idx=0; outer_idx<num_rings_in; outer_idx++) {
		ring_t *ring = mpoly_in->rings + outer_idx;
		if(ring->is_hole) {
			int parent = ring->parent_id;
			holes[parent] = (int *)realloc_or_die(holes[parent],
				sizeof(int) * (num_holes[parent]+1));
			holes[parent][num_holes[parent]++] = outer_idx;
		} else {
			num_geom_out++;
		}
	}

	int use_multi = num_geom_out > 1;

	OGRGeometryH geom_out = OGR_G_CreateGeometry(
		use_multi ? wkbMultiPolygon : wkbPolygon);

	for(int outer_idx=0; outer_idx<num_rings_in; outer_idx++) {
		ring_t *ring = mpoly_in->rings + outer_idx;
		if(ring->is_hole) continue;

		OGRGeometryH poly_out = use_multi ? OGR_G_CreateGeometry(wkbPolygon) : geom_out;
		OGR_G_AddGeometry(poly_out, ring_to_ogr(ring));

		for(int hole_idx=0; hole_idx<num_holes[outer_idx]; hole_idx++) {
			ring_t *hole = mpoly_in->rings + holes[outer_idx][hole_idx];
			if(hole->parent_id != outer_idx) fatal_error("could not sort out holes");
			OGR_G_AddGeometry(poly_out, ring_to_ogr(hole));
		}

		if(OGR_G_GetGeometryCount(poly_out) != num_holes[outer_idx]+1) {
			fatal_error("GeometryCount != num_holes+1 (%d vs. %d)", 
				OGR_G_GetGeometryCount(poly_out), num_holes[outer_idx]+1);
		}

		if(use_multi) OGR_G_AddGeometry(geom_out, poly_out);
	}

	if(use_multi && OGR_G_GetGeometryCount(geom_out) != num_geom_out) {
		fatal_error("GeometryCount != num_geom_out (%d vs. %d)",
			OGR_G_GetGeometryCount(geom_out), num_geom_out);
	}

	for(int i=0; i<num_rings_in; i++) {
		if(holes[i]) free(holes[i]);
	}
	free(holes);
	free(num_holes);

	return geom_out;
}

mpoly_t ogr_to_mpoly(OGRGeometryH geom_in) {
	OGRwkbGeometryType type = OGR_G_GetGeometryType(geom_in);
	if(type == wkbPolygon) {
		mpoly_t mpoly_out;
		mpoly_out.num_rings = OGR_G_GetGeometryCount(geom_in);
		if(mpoly_out.num_rings < 1) fatal_error("num_rings<1 in ogr_to_mpoly");

		mpoly_out.rings = (ring_t *)malloc_or_die(sizeof(ring_t) * mpoly_out.num_rings);

		for(int i=0; i<mpoly_out.num_rings; i++) {
			mpoly_out.rings[i] = ogr_to_ring(OGR_G_GetGeometryRef(geom_in, i));
			mpoly_out.rings[i].is_hole = (i>0);
			mpoly_out.rings[i].parent_id = (i>0) ? 0 : -1;
		}

		return mpoly_out;
	} else if(type == wkbMultiPolygon || type == wkbGeometryCollection) {
		int num_geom = OGR_G_GetGeometryCount(geom_in);
		mpoly_t *polys = (mpoly_t *)malloc_or_die(sizeof(mpoly_t) * num_geom);
		int total_rings = 0;
		for(int i=0; i<num_geom; i++) {
			OGRGeometryH g = OGR_G_GetGeometryRef(geom_in, i);
			polys[i] = ogr_to_mpoly(g);
			total_rings += polys[i].num_rings;
		}

		mpoly_t mpoly_out;
		mpoly_out.num_rings = total_rings;
		if(mpoly_out.num_rings < 1) fatal_error("num_rings<1 in ogr_to_mpoly");
		mpoly_out.rings = (ring_t *)malloc_or_die(sizeof(ring_t) * mpoly_out.num_rings);

		int o = 0;
		for(int i=0; i<num_geom; i++) {
			for(int j=0; j<polys[i].num_rings; j++) {
				ring_t ring = polys[i].rings[j];
				if(ring.is_hole) ring.parent_id += o;
				mpoly_out.rings[o+j] = ring;
			}
			o += polys[i].num_rings;
			free(polys[i].rings);
		}
		free(polys);
		return mpoly_out;
	} else {
		fatal_error("not a polygon type: %d\n", type);
	}
}

void split_mpoly_to_polys(mpoly_t *mpoly_in, int *num_polys, mpoly_t **polys) {
	int num_rings_in = mpoly_in->num_rings;

	int *num_holes = (int *)malloc_or_die(sizeof(int) * num_rings_in);
	int **holes = (int **)malloc_or_die(sizeof(int *) * num_rings_in);
	for(int outer_idx=0; outer_idx<num_rings_in; outer_idx++) {
		num_holes[outer_idx] = 0;
		holes[outer_idx] = NULL;
	}

	int num_geom_out = 0;
	for(int outer_idx=0; outer_idx<num_rings_in; outer_idx++) {
		ring_t *ring = mpoly_in->rings + outer_idx;
		if(ring->is_hole) {
			int parent = ring->parent_id;
			holes[parent] = (int *)realloc_or_die(holes[parent],
				sizeof(int) * (num_holes[parent]+1));
			holes[parent][num_holes[parent]++] = outer_idx;
		} else {
			num_geom_out++;
		}
	}

	*num_polys = num_geom_out;
	*polys = (mpoly_t *)malloc_or_die(sizeof(mpoly_t) * num_geom_out);
	int poly_out_idx = 0;

	for(int outer_idx=0; outer_idx<num_rings_in; outer_idx++) {
		ring_t *ring = mpoly_in->rings + outer_idx;
		if(ring->is_hole) continue;

		mpoly_t out_poly;
		out_poly.num_rings = num_holes[outer_idx]+1;
		out_poly.rings = (ring_t *)malloc_or_die(sizeof(ring_t) * out_poly.num_rings);
		int ring_out_idx = 0;

		ring_t dup_ring = duplicate_ring(ring);
		dup_ring.parent_id = -1;
		out_poly.rings[ring_out_idx++] = dup_ring;

		for(int hole_idx=0; hole_idx<num_holes[outer_idx]; hole_idx++) {
			ring_t *hole = mpoly_in->rings + holes[outer_idx][hole_idx];
			if(hole->parent_id != outer_idx) fatal_error("could not sort out holes");

			ring_t dup_ring = duplicate_ring(hole);
			dup_ring.parent_id = 0;
			out_poly.rings[ring_out_idx++] = dup_ring;
		}

		(*polys)[poly_out_idx++] = out_poly;
	}

	for(int i=0; i<num_rings_in; i++) {
		if(holes[i]) free(holes[i]);
	}
	free(holes);
	free(num_holes);
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

int line_intersects_line(
	vertex_t p1, vertex_t p2,
	vertex_t p3, vertex_t p4,
	int fail_on_coincident
) {
	if(
		MAX(p1.x, p2.x) < MIN(p3.x, p4.x) ||
		MIN(p1.x, p2.x) > MAX(p3.x, p4.x) ||
		MAX(p1.y, p2.y) < MIN(p3.y, p4.y) ||
		MIN(p1.y, p2.y) > MAX(p3.y, p4.y)
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

void line_line_intersection(
	vertex_t p1, vertex_t p2,
	vertex_t p3, vertex_t p4,
	vertex_t *p_out
) {
	double numer_a = (p4.x-p3.x)*(p1.y-p3.y) - (p4.y-p3.y)*(p1.x-p3.x);
	//double numer_b = (p2.x-p1.x)*(p1.y-p3.y) - (p2.y-p1.y)*(p1.x-p3.x);
	double denom   = (p4.y-p3.y)*(p2.x-p1.x) - (p4.x-p3.x)*(p2.y-p1.y);
	if(!denom) fatal_error("lines are parallel");
	double ua = numer_a / denom;
	//double ub = numer_b / denom;
	(*p_out).x = p1.x + ua*(p2.x-p1.x);
	(*p_out).y = p1.y + ua*(p2.y-p1.y);
}

double ring_oriented_area(ring_t *c) {
	double accum = 0;
	int i;
	for(i=0; i<c->npts; i++) {
		double x0 = c->pts[i].x;
		double y0 = c->pts[i].y;
		double x1 = c->pts[(i+1)%c->npts].x;
		double y1 = c->pts[(i+1)%c->npts].y;
		accum += x0*y1 - x1*y0;
	}
	return accum / 2.0;
}

int ring_is_ccw(ring_t *c) {
	return ring_oriented_area(c) > 0;
}

double ring_area(ring_t *c) {
	return fabs(ring_oriented_area(c));
}

static int ring_contains_point(ring_t *ring, double px, double py) {
	int npts = ring->npts;
	vertex_t *pts = ring->pts;
	int num_crossings = 0;
	int i;
	for(i=0; i<npts; i++) {
		double x0 = pts[i].x;
		double y0 = pts[i].y;
		double x1 = pts[(i+1) % npts].x;
		double y1 = pts[(i+1) % npts].y;
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

int polygon_contains_point(mpoly_t *mp, double px, double py) {
	int num_crossings = 0;
	for(int r_idx=0; r_idx<mp->num_rings; r_idx++) {
		if(ring_contains_point(mp->rings+r_idx, px, py)) num_crossings++;
	}
	// if it is within an odd number of rings it is not in a hole
	return num_crossings & 1;
}

int ring_ring_relation(ring_t *r1, ring_t *r2) {
	bbox_t bb1 = get_ring_bbox(r1);
	bbox_t bb2 = get_ring_bbox(r2);
	if(bboxes_disjoint(&bb1, &bb2)) return RING_DISJOINT;

	int n1 = r1->npts;
	int n2 = r2->npts;
	// test for crossings
	for(int i1=0; i1<n1; i1++) {
		vertex_t p1a = r1->pts[i1];
		vertex_t p1b = r1->pts[(i1+1) % n1];
		if(
			(p1a.x < bb2.min_x && p1b.x < bb2.min_x) ||
			(p1a.y < bb2.min_y && p1b.y < bb2.min_y) ||
			(p1a.x > bb2.max_x && p1b.x > bb2.max_x) ||
			(p1a.y > bb2.max_y && p1b.y > bb2.max_y)
		) continue;
		for(int i2=0; i2<n2; i2++) {
			vertex_t p2a = r2->pts[i2];
			vertex_t p2b = r2->pts[(i2+1) % n2];
			if(line_intersects_line(p1a, p1b, p2a, p2b, 0)) return RING_CROSSES;
		}
	}

	if(ring_contains_point(
		r1, r2->pts[0].x, r2->pts[0].y)) return RING_CONTAINS;
	else if(ring_contains_point(
		r2, r1->pts[0].x, r1->pts[0].y)) return RING_CONTAINED_BY;
	else return RING_DISJOINT;
}

/*
void compute_containments(mpoly_t *mp) {
	int i;
	int nrings = mp->num_rings;

	bbox_t *bboxes = make_bboxes(mp);

	int **ancestors = (int **)malloc_or_die(sizeof(int *) * nrings);
	int *num_ancestors = (int *)malloc_or_die(sizeof(int) * nrings);

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
				ancestors[i] = (int *)realloc_or_die(ancestors[i],
					sizeof(int) * (num_ancestors[i]+1));
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

mpoly_t *mpoly_xy2en(georef_t *georef, mpoly_t *xy_poly) {
	mpoly_t *en_poly = (mpoly_t *)malloc_or_die(sizeof(mpoly_t));
	en_poly->num_rings = xy_poly->num_rings;
	en_poly->rings = (ring_t *)malloc_or_die(sizeof(ring_t) * en_poly->num_rings);

	int r_idx;
	for(r_idx=0; r_idx<xy_poly->num_rings; r_idx++) {
		ring_t *xy_ring = xy_poly->rings + r_idx;
		ring_t *en_ring = en_poly->rings + r_idx;

		*en_ring = duplicate_ring(xy_ring);

		int v_idx;
		for(v_idx=0; v_idx<en_ring->npts; v_idx++) {
			double x = xy_ring->pts[v_idx].x;
			double y = xy_ring->pts[v_idx].y;
			double east, north;
			xy2en(georef, x, y, &east, &north);
			en_ring->pts[v_idx].x = east;
			en_ring->pts[v_idx].y = north;
		}
	}

	return en_poly;
}

mpoly_t *mpoly_en2xy(georef_t *georef, mpoly_t *en_poly) {
	mpoly_t *xy_poly = (mpoly_t *)malloc_or_die(sizeof(mpoly_t));
	xy_poly->num_rings = en_poly->num_rings;
	xy_poly->rings = (ring_t *)malloc_or_die(sizeof(ring_t) * xy_poly->num_rings);

	int r_idx;
	for(r_idx=0; r_idx<en_poly->num_rings; r_idx++) {
		ring_t *en_ring = en_poly->rings + r_idx;
		ring_t *xy_ring = xy_poly->rings + r_idx;

		*xy_ring = duplicate_ring(en_ring);

		int v_idx;
		for(v_idx=0; v_idx<xy_ring->npts; v_idx++) {
			double east = en_ring->pts[v_idx].x;
			double north = en_ring->pts[v_idx].y;
			double x, y;
			en2xy(georef, east, north, &x, &y);
			xy_ring->pts[v_idx].x = x;
			xy_ring->pts[v_idx].y = y;
		}
	}

	return xy_poly;
}

// This way is very simple and robust, but doesn't know anything
// about the actual distortion of the projection.
/*
mpoly_t *mpoly_xy2ll_with_interp(georef_t *georef, mpoly_t *xy_poly, double toler_pixels) {
	mpoly_t *ll_poly = (mpoly_t *)malloc_or_die(sizeof(mpoly_t));
	ll_poly->num_rings = xy_poly->num_rings;
	ll_poly->rings = (ring_t *)malloc_or_die(sizeof(ring_t) * ll_poly->num_rings);

	OGRErr err = OGRERR_NONE;
	double earth_radius = OSRGetSemiMajor(georef->spatial_ref, &err);
	if(err != OGRERR_NONE) fatal_error("could not determine globe radius");

	double toler_radians = toler_pixels * 
		MIN(georef->res_meters_x, georef->res_meters_y) / earth_radius;
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
			double dx = (xy1.x - xy2.x) * georef->res_meters_x / earth_radius;
			double dy = (xy1.y - xy2.y) * georef->res_meters_y / earth_radius;
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
static double estimate_canvas_size_sq(georef_t *georef) {
	double e1, n1, e2, n2;
	xy2en(georef, 0, 0, &e1, &n1);
	xy2en(georef, georef->w, 0, &e2, &n2);
	double dx = (e1 - e2);
	double dy = (n1 - n2);
	double size1 = dx*dx + dy*dy;
	xy2en(georef, 0, georef->h, &e2, &n2);
	dx = (e1 - e2);
	dy = (n1 - n2);
	double size2 = dx*dx + dy*dy;

	double size = MAX(size1, size2);
	size *= georef->units_val * georef->units_val;
	if(OSRIsGeographic(georef->spatial_ref)) {
		size *= georef->semi_major * georef->semi_major;
	}
	return size;
}

mpoly_t *mpoly_xy2ll_with_interp(georef_t *georef, mpoly_t *xy_poly, double toler) {
	mpoly_t *ll_poly = (mpoly_t *)malloc_or_die(sizeof(mpoly_t));
	ll_poly->num_rings = xy_poly->num_rings;
	ll_poly->rings = (ring_t *)malloc_or_die(sizeof(ring_t) * ll_poly->num_rings);
	
	double canvas_size_sq = estimate_canvas_size_sq(georef);
	//printf("canvas_size_sq = %lf\n", canvas_size_sq);

	// FIXME - now that we don't use ll2xy this is probably not needed:
	// This is a kludge that shrinks the canvas by a millionth of a pixel
	// to avoid problems with images that span an entire 360 degrees of
	// longitude.  Without this the map xy -> ll -> xy is not single-valued.
	double epsilon = 5e-7;
	double shrink = ((double)georef->w - 2.0*epsilon) / (double)georef->w;

	int r_idx;
	for(r_idx=0; r_idx<xy_poly->num_rings; r_idx++) {
		// make a copy of input - we will modify this
		ring_t xy_ring = duplicate_ring(xy_poly->rings + r_idx);
		// this will be the output
		ring_t ll_ring = duplicate_ring(xy_poly->rings + r_idx);

		int v_idx;
		for(v_idx=0; v_idx<ll_ring.npts; v_idx++) {
			double x = xy_ring.pts[v_idx].x;
			double y = xy_ring.pts[v_idx].y;
			double lon, lat;
			xy2ll_or_die(georef, x*shrink+epsilon, y, &lon, &lat);
			ll_ring.pts[v_idx].x = lon;
			ll_ring.pts[v_idx].y = lat;
		}

		int num_consec = 0;

		for(v_idx=0; v_idx<ll_ring.npts; ) {
			if(xy_ring.npts != ll_ring.npts) fatal_error("xy_ring.npts != ll_ring.npts");

			vertex_t *xy1 = xy_ring.pts + v_idx;
			vertex_t *xy2 = xy_ring.pts + (v_idx + 1) % xy_ring.npts;
			vertex_t xy_m = (vertex_t) { 
				(xy1->x + xy2->x)/2.0,
				(xy1->y + xy2->y)/2.0 };

			vertex_t ll_m_proj;
			xy2ll_or_die(georef, 
				xy_m.x*shrink+epsilon, xy_m.y,
				&ll_m_proj.x, &ll_m_proj.y);

			vertex_t *ll1 = ll_ring.pts + v_idx;
			vertex_t *ll2 = ll_ring.pts + (v_idx + 1) % ll_ring.npts;

			while(ll1->x - ll2->x >  180) ll1->x -= 360;
			while(ll1->x - ll2->x < -180) ll1->x += 360;

			int need_midpt = 0;

			// avoid topological errors by not allowing segments
			// to be longer than 90 degrees
			if(!need_midpt) {
				double dx = fabs(ll1->x - ll2->x);
				double dy = ll1->y - ll2->y;
				double sqr_error = dx*dx + dy*dy;
				double max_error = 90; // degrees

				need_midpt = sqr_error > max_error*max_error;

				if(VERBOSE && need_midpt) {
					printf("  inserting midpoint at %d,%d (lon/lat difference %g > %g)\n",
						r_idx, v_idx, sqrt(sqr_error), max_error);
				}
			}

			if(!need_midpt) {
				vertex_t ll_m_interp;
				ll_m_interp.x = (ll1->x + ll2->x)/2.0;
				ll_m_interp.y = (ll1->y + ll2->y)/2.0;

				while(ll_m_interp.x - ll_m_proj.x >  180) ll_m_interp.x -= 360;
				while(ll_m_interp.x - ll_m_proj.x < -180) ll_m_interp.x += 360;
				double lonscale = cos(D2R * MAX(fabs(ll_m_interp.y), fabs(ll_m_proj.y)));
				double dx = (ll_m_interp.x - ll_m_proj.x) * lonscale;
				double dy = ll_m_interp.y - ll_m_proj.y;
				dx *= D2R * georef->semi_major;
				dy *= D2R * georef->semi_major;
				double sqr_error = dx*dx + dy*dy;

				// if the midpoint is this far off then something is seriously wrong
				if(sqr_error > canvas_size_sq) {
					fprintf(stderr, "\nInfo on bad point:\n");
					fprintf(stderr, "xy1 = %lf,%lf\n", xy1->x, xy1->y);
					fprintf(stderr, "xy2 = %lf,%lf\n", xy2->x, xy2->y);
					fprintf(stderr, "xy_m = %lf,%lf\n", xy_m.x, xy_m.y);
					fprintf(stderr, "ll1 = %lf,%lf\n", ll1->x, ll1->y);
					fprintf(stderr, "ll2 = %lf,%lf\n", ll2->x, ll2->y);
					fprintf(stderr, "ll_m_proj = %lf,%lf\n", ll_m_proj.x, ll_m_proj.y);
					fprintf(stderr, "ll_m_interp = %lf,%lf\n", ll_m_interp.x, ll_m_interp.y);
					fprintf(stderr, "error = %lf > %lf\n", sqr_error, canvas_size_sq);
					fatal_error("projection error in mpoly_xy2ll_with_interp");
				}

				need_midpt = toler && sqr_error > toler*toler;

				if(VERBOSE && need_midpt) {
					printf("  inserting midpoint at %d,%d (delta=%lf,%lf > %lf)\n",
						r_idx, v_idx, dx, dy, toler);
					printf("    testll=[%lf,%lf] projll=[%lf,%lf]\n", 
						ll_m_interp.x, ll_m_interp.y, ll_m_proj.x, ll_m_proj.y);
				}
			}

			if(need_midpt) {
				if(num_consec++ > 10) fatal_error("convergence error in mpoly_xy2ll_with_interp");

				if(VERBOSE) {
					printf("    xy=[%lf,%lf]:[%lf,%lf]\n", xy1->x, xy1->y, xy2->x, xy2->y);
					printf("    ll=[%lf,%lf]:[%lf,%lf]\n", ll1->x, ll1->y, ll2->x, ll2->y);
					printf("    midxy=[%lf,%lf] midll=[%lf,%lf]\n", 
						xy_m.x, xy_m.y, ll_m_proj.x, ll_m_proj.y);
				}
				insert_point_into_ring(&xy_ring, v_idx+1);
				insert_point_into_ring(&ll_ring, v_idx+1);
				xy_ring.pts[v_idx+1] = xy_m;
				ll_ring.pts[v_idx+1] = ll_m_proj;
			} else {
				v_idx++;
				num_consec = 0;
			}
		}

		// Now reproject everything again, without using the shrink
		// kludge.  We no longer care whether the projection is
		// single-valued and we don't want the loss of accuracy
		// that comes from multiplying by shrink.
		for(v_idx=0; v_idx<ll_ring.npts; v_idx++) {
			double x = xy_ring.pts[v_idx].x;
			double y = xy_ring.pts[v_idx].y;

			double lon, lat;
			xy2ll_or_die(georef, x, y, &lon, &lat);

			ll_ring.pts[v_idx].x = lon;
			ll_ring.pts[v_idx].y = lat;
		}

		free_ring(&xy_ring);
		ll_poly->rings[r_idx] = ll_ring;
	}

	return ll_poly;
}

static char *read_whole_file(FILE *fin) {
	int chunk_size = 1024;
	int data_len = 0;
	int buf_len = chunk_size+1;
	char *buffer = (char *)malloc_or_die(buf_len);
	int num_read;
	while(0 < (num_read = fread(buffer+data_len, 1, chunk_size, fin))) {
		data_len += num_read;
		if(data_len+chunk_size+1 > buf_len) {
			buf_len += chunk_size;
			buffer = (char *)realloc_or_die(buffer, buf_len);
		}
	}
	buffer[data_len++] = 0;
	return buffer;
}

mpoly_t mpoly_from_wktfile(const char *fn) {
	FILE *fh = fopen(fn, "r");
	if(!fh) fatal_error("cannot read file [%s]", fn);
	char *wkt_in = read_whole_file(fh);
	int i;
	for(i=0; wkt_in[i]; i++) {
		if(wkt_in[i] == '\r') wkt_in[i] = ' ';
		if(wkt_in[i] == '\n') wkt_in[i] = ' ';
		if(wkt_in[i] == '\t') wkt_in[i] = ' ';
	}

	OGRGeometryH geom;
	OGRErr err = OGR_G_CreateFromWkt(&wkt_in, NULL, &geom);
	if(OGRERR_NONE != err) {
		fatal_error("OGR_G_CreateFromWkt failed: %d", err);
	}
	
	mpoly_t mp = ogr_to_mpoly(geom);
	return mp;
}
