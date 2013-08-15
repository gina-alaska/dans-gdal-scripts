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
#include <utility>
#include <cassert>

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

	bool contains(const Vertex &v) const {
		return !empty &&
			v.x >= min_x && v.x <= max_x &&
			v.y >= min_y && v.y <= max_y;
	}

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

	double width()  { return max_x - min_x; }
	double height() { return max_y - min_y; }

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

// BboxBinarySpacePartition helps you quickly find which of a list of items (that have bounding
// boxes) intersect a given bounding box.  Think of it as a std::map whose keys are bounding
// boxes.  Since several items may intersect a given query box, query returns a list of
// matches.
//
// When many object have overlapping bbox, this tree may not be so efficient and maybe
// something better should be cooked up if it turns out to be slow.
template <typename T>
class BboxBinarySpacePartition {
public:
	BboxBinarySpacePartition(
		std::vector<std::pair<Bbox, T> > items,
		size_t max_leaf_size = 20,
		bool axis = 0
	) :
		left(NULL),
		right(NULL)
	{
		leaf_items = items;

		// Compute bbox containing all items.
		for(size_t i=0; i<items.size(); i++) {
			node_bbox.expand(items[i].first);
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
			if(!is_disjoint(needle, left->node_bbox)) {
				left ->append_intersecting_items(out, needle);
			}
			if(!is_disjoint(needle, right->node_bbox)) {
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
		std::vector<double> midpts(leaf_items.size());
		double cmin = 0, cmax = 0, cavg = 0;
		size_t num_nonempty = 0;
		for(size_t i=0; i<leaf_items.size(); i++) {
			const Bbox &ibb = leaf_items[i].first;
			if(ibb.empty) continue;
			num_nonempty++;
			node_bbox.expand(ibb);
			double item_mid = axis ?
				(ibb.min_y + ibb.max_y) / 2.0 :
				(ibb.min_x + ibb.max_x) / 2.0;
			if(i == 0) cmin = cmax = item_mid;
			cmin = std::min(cmin, item_mid);
			cmax = std::max(cmax, item_mid);
			cavg += item_mid;
		}
		cavg /= num_nonempty;

		// Subdivide space into two halves.  Note: this is used to decide which branch of the
		// tree each item goes into.  However, each branch will then compute the exact Bbox
		// that contains each of its items (it won't use these ones we compute here).
		Bbox left_bbox, right_bbox;
		if(axis) {
			left_bbox = Bbox(
				node_bbox.min_x, node_bbox.max_x,
				cmin, cavg);
			right_bbox = Bbox(
				node_bbox.min_x, node_bbox.max_x,
				cavg, cmax);
		} else {
			left_bbox = Bbox(
				cmin, cavg,
				node_bbox.min_y, node_bbox.max_y);
			right_bbox = Bbox(
				cavg, cmax,
				node_bbox.min_y, node_bbox.max_y);
		}

		// Find which items go in each box.
		//
		// Note: it would be possible to just move items to the left and right side of the
		// original array like with quicksort, and never have to allocate all of these little
		// vectors for every node.
		std::vector<std::pair<Bbox, T> > left_items, right_items;
		for(size_t i=0; i<leaf_items.size(); i++) {
			const Bbox &ibb = leaf_items[i].first;
			Vertex item_mid = Vertex(
				(ibb.min_x + ibb.max_x) / 2.0,
				(ibb.min_y + ibb.max_y) / 2.0);
			if(left_bbox.contains(item_mid)) {
				left_items.push_back(leaf_items[i]);
			} else {
				right_items.push_back(leaf_items[i]);
			}
		}

		// free memory
		leaf_items.clear();

		left  = new BboxBinarySpacePartition<T>(left_items,  max_leaf_size, !axis);
		right = new BboxBinarySpacePartition<T>(right_items, max_leaf_size, !axis);
	}

private:
	Bbox node_bbox;
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

	void debug_dump_binary(FILE *fh) const;
	static Ring debug_load_binary(FILE *fh);

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

	void debug_dump_binary(FILE *fh) const;
	static Mpoly debug_load_binary(FILE *fh);

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
