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



#ifndef POLYGON_H
#define POLYGON_H

#include <ogr_api.h>
#include "common.h"
#include "georef.h"

typedef struct {
	double x, y;
} vertex_t;

typedef struct {
	int npts;
	vertex_t *pts;
	int is_hole;
	int parent_id;
} ring_t;

typedef struct {
	int num_rings;
	ring_t *rings;
} mpoly_t;


mpoly_t empty_polygon();
ring_t duplicate_ring(ring_t *in_ring);
void free_ring(ring_t *ring);
void free_mpoly(mpoly_t *mpoly);
void insert_point_into_ring(ring_t *ring, int idx);
OGRGeometryH mpoly_to_ogr(mpoly_t *mpoly_in);
void split_mpoly_to_polys(mpoly_t *mpoly, int *num_polys, mpoly_t **polys);
mpoly_t compute_reduced_pointset(mpoly_t *in_mpoly, double tolerance);
int polygon_contains(ring_t *c1, ring_t *c2);
double polygon_area(ring_t *c);
void compute_containments(mpoly_t *mp);
void mask_from_mpoly(mpoly_t *mpoly, int w, int h, char *fn);
int line_intersects_line(
	vertex_t p1, vertex_t p2,
	vertex_t p3, vertex_t p4,
	int fail_on_coincident
);
void line_line_intersection(
	vertex_t p1, vertex_t p2,
	vertex_t p3, vertex_t p4,
	vertex_t *p_out
);
void bevel_self_intersections(mpoly_t *mp, double amount);
mpoly_t *mpoly_xy2en(georef_t *georef, mpoly_t *xy_poly);
mpoly_t *mpoly_xy2ll_with_interp(georef_t *georef, mpoly_t *xy_poly, double toler);

#endif // ifndef POLYGON_H
