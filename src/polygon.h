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



#ifndef DANGDAL_POLYGON_H
#define DANGDAL_POLYGON_H

#include <vector>
#include <string>
#include <algorithm>

#include <ogr_api.h>

#include "common.h"
#include "georef.h"

namespace dangdal {

// returned by ring_ring_relation
enum RingRelation {
	RINGREL_CONTAINS,
	RINGREL_CONTAINED_BY,
	RINGREL_CROSSES,
	RINGREL_DISJOINT
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

	Bbox(
		double _min_x, double _max_x,
		double _min_y, double _max_y
	) :
		min_x(_min_x), max_x(_max_x),
		min_y(_min_y), max_y(_max_y),
		empty(false)
	{ }

	void expand(const Vertex &v) {
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

	void expand(const Bbox &bb);

	double min_x, max_x, min_y, max_y;
	bool empty;
};

Bbox box_union(const Bbox &bb1, const Bbox &bb2);

static inline bool is_disjoint(const Bbox &bb1, const Bbox &bb2) {
	return
		bb1.empty || bb2.empty ||
		bb1.min_x >  bb2.max_x ||
		bb1.min_y >  bb2.max_y ||
		bb2.min_x >  bb1.max_x ||
		bb2.min_y >  bb1.max_y;
}

template <typename T>
class BboxBinarySpacePartition {
public:
	BboxBinarySpacePartition(
		std::vector<std::pair<Bbox, T> > items,
		size_t max_leaf_size = 10,
		bool axis = 0
	) :
		left(NULL),
		right(NULL)
	{
		leaf_items = items;

		// compute bbox of all items
		for(size_t i=0; i<items.size(); i++) {
			bbox.expand(items[i].first);
		}

		if(items.size() > max_leaf_size) {
			subdivide(max_leaf_size, axis);
		}
	}

	~BboxBinarySpacePartition() {
		delete(left);
		delete(right);
	}

	std::vector<T> get_intersecting_items(Bbox needle) const {
		std::vector<T> ret;
		append_intersecting_items(ret, needle);
		return ret;
	}

private:
	void append_intersecting_items(
		std::vector<T> &out,
		Bbox needle
	) const {
		if(left) {
			assert(right);
			if(!is_disjoint(needle, left->bbox)) {
				left ->append_intersecting_items(out, needle);
			}
			if(!is_disjoint(needle, right->bbox)) {
				right->append_intersecting_items(out, needle);
			}
		} else {
			for(size_t i=0; i<leaf_items.size(); i++) {
				if(!is_disjoint(needle, leaf_items[i].first)) {
					out.push_back(leaf_items[i].second);
				}
			}
		}
	}

	void subdivide(size_t max_leaf_size, bool axis) {
		// Subdivide space into two halves.  Note: this is used to decide which branch of the
		// tree each item goes into.  However, each branch will then compute the exact Bbox
		// that contains each of its items (it won't use these ones we compute here).
		Bbox left_bbox, right_bbox;
		if(axis) {
			double mid_y = (bbox.min_y + bbox.max_y) / 2.0;
			left_bbox = Bbox(
				bbox.min_x, bbox.max_x,
				bbox.min_y, mid_y);
			right_bbox = Bbox(
				bbox.min_x, bbox.max_x,
				mid_y, bbox.max_y);
		} else {
			double mid_x = (bbox.min_x + bbox.max_x) / 2.0;
			left_bbox = Bbox(
				bbox.min_x, mid_x,
				bbox.min_y, bbox.max_y);
			right_bbox = Bbox(
				mid_x, bbox.max_x,
				bbox.min_y, bbox.max_y);
		}

		// find which items go in each box
		std::vector<std::pair<Bbox, T> > left_items, right_items;
		for(size_t i=0; i<leaf_items.size(); i++) {
			bool hit_l = !is_disjoint(left_bbox,  leaf_items[i].first);
			bool hit_r = !is_disjoint(right_bbox, leaf_items[i].first);
			bool goes_on_right;
			if(hit_l && hit_r) {
				// If it overlaps both, then choose randomly (well, sort of).  The random
				// choice leads to a more balanced tree.
				goes_on_right = i % 2;
			} else if(hit_l) {
				goes_on_right = 0;
			} else if(hit_r) {
				goes_on_right = 1;
			} else {
				// It's gotta be in one of the halfs!
				assert(0);
			}
			if(goes_on_right) {
				right_items.push_back(leaf_items[i]);
			} else {
				left_items.push_back(leaf_items[i]);
			}
		}

		// free memory
		leaf_items.clear();

		left  = new BboxBinarySpacePartition<T>(left_items,  max_leaf_size, !axis);
		right = new BboxBinarySpacePartition<T>(right_items, max_leaf_size, !axis);
	}

private:
	Bbox bbox;
	BboxBinarySpacePartition<T> *left, *right;
	// only used by leafs
	std::vector<std::pair<Bbox, T> > leaf_items;
};

class Ring {
public:
	Ring() : is_hole(false), parent_id(-1) { }

	Bbox getBbox() const;
	double orientedArea() const;
	double area() const;
	bool isCCW() const;
	bool contains(Vertex p) const;
	void reverse() { std::reverse(pts.begin(), pts.end()); }

	// copy everything except for pts
	Ring copyMetadata() const {
		Ring ret;
		ret.is_hole = is_hole;
		ret.parent_id = parent_id;
		return ret;
	}

	std::vector<Vertex> pts;
	bool is_hole;
	int parent_id;
};

class Mpoly {
public:
	Bbox getBbox() const;
	std::vector<Bbox> getRingBboxes() const;
	bool contains(Vertex p) const;
	// Returns true if the given point is contained in the specified ring, but
	// not contained in any of that ring's holes.
	bool component_contains(Vertex p, int outer_ring_id) const;
	void deleteRing(size_t idx);

	void xy2en(const GeoRef &georef);
	void en2xy(const GeoRef &georef);
	void xy2ll_with_interp(const GeoRef &georef, double toler);

	std::vector<Ring> rings;
};

OGRGeometryH ring_to_ogr(const Ring &ring);
Ring ogr_to_ring(OGRGeometryH ogr);
OGRGeometryH mpoly_to_ogr(const Mpoly &mpoly_in);
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
Mpoly mpoly_from_wktfile(const std::string &fn);

} // namespace dangdal

#endif // ifndef DANGDAL_POLYGON_H
