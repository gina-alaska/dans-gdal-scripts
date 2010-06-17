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



#ifndef POLYGON_H
#define POLYGON_H

#include <vector>
#include <algorithm>
#include <ogr_api.h>
#include "common.h"
#include "georef.h"

namespace dangdal {

// returned by ring_ring_relation
enum RingRelation {
	CONTAINS,
	CONTAINED_BY,
	CROSSES,
	DISJOINT
};

struct Vertex {
	Vertex() : x(0), y(0) { }
	Vertex(double _x, double _y) : x(_x), y(_y) { }

	double x, y;
};

class Bbox {
public:
	Bbox() :
		min_x(0), max_x(0),
		min_y(0), max_y(0),
		empty(true)
	{ }

	void expand(const Vertex v) {
		if(empty) {
			empty = false;
			min_x = max_x = v.x;
			min_y = max_y = v.y;
		} else {
			if(v.x < min_x) min_x = v.x;
			if(v.y < min_y) min_y = v.y;
			if(v.x > max_x) max_x = v.x;
			if(v.y > max_y) max_y = v.y;
		}
	}

	void expand(const Bbox bb);

	double min_x, max_x, min_y, max_y;
	bool empty;
};

Bbox box_union(const Bbox bb1, const Bbox bb2);
bool is_disjoint(const Bbox bb1, const Bbox bb2);

class Ring {
public:
	Ring() : is_hole(false), parent_id(-1) { }

	Bbox getBbox() const;
	double orientedArea() const;
	double area() const;
	bool isCCW() const;
	bool contains(Vertex p) const;
	void reverse() { std::reverse(pts.begin(), pts.end()); }

	std::vector<Vertex> pts;
	bool is_hole;
	int parent_id;
};

class Mpoly {
public:
	Bbox getBbox() const;
	std::vector<Bbox> getRingBboxes() const;
	bool contains(Vertex p) const;

	void xy2en(georef_t *georef);
	void en2xy(georef_t *georef);
	void xy2ll_with_interp(georef_t *georef, double toler);

	std::vector<Ring> rings;
};

OGRGeometryH ring_to_ogr(const Ring &ring);
Ring ogr_to_ring(OGRGeometryH ogr);
OGRGeometryH mpoly_to_ogr(const Mpoly mpoly_in);
Mpoly ogr_to_mpoly(OGRGeometryH geom_in);
std::vector<Mpoly> split_mpoly_to_polys(const Mpoly &mpoly);
bool line_intersects_line(
	Vertex p1, Vertex p2,
	Vertex p3, Vertex p4,
	bool fail_on_coincident
);
Vertex line_line_intersection(
	Vertex p1, Vertex p2,
	Vertex p3, Vertex p4
);
RingRelation ring_ring_relation(const Ring &r1, const Ring &r2);
Mpoly mpoly_from_wktfile(const char *fn);

} // namespace dangdal

#endif // ifndef POLYGON_H
