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



#include <cassert>

#include "common.h"
#include "polygon.h"

namespace dangdal {

struct VertRef {
	VertRef() : ring_idx(0), vert_idx(0) { }
	VertRef(size_t _r, size_t _v) : ring_idx(_r), vert_idx(_v) { }

	const Vertex &getVert(const Mpoly &mp) const {
		return mp.rings[ring_idx].pts[vert_idx];
	}

	size_t ring_idx;
	size_t vert_idx;
};

class CoordsComparator {
public:
	explicit CoordsComparator(Mpoly *_mp) : mp(_mp) { }

	bool operator()(const VertRef &a, const VertRef &b) const {
		const Vertex &va = a.getVert(*mp);
		const Vertex &vb = b.getVert(*mp);

		return 
			va.x < vb.x ? true :
			va.x > vb.x ? false :
			va.y < vb.y ? true :
			va.y > vb.y ? false :
			// It is not strictly necessary to have this tiebreaker here, but
			// it ensures predictability with regard to which of two
			// intersecting corners get shaved.  Consistency is important for
			// unit tests.
			a.ring_idx < b.ring_idx ? true :
			a.ring_idx > b.ring_idx ? false :
			a.vert_idx < b.vert_idx ? true :
			false;
	}

private:
	const Mpoly *mp;
};

class RingsComparator {
public:
	explicit RingsComparator(Mpoly *_mp) : mp(_mp) { }

	bool operator()(const VertRef &a, const VertRef &b) const {
		return
			a.ring_idx < b.ring_idx ? true :
			a.ring_idx > b.ring_idx ? false :
			a.vert_idx < b.vert_idx ? true :
			false;
	}

private:
	const Mpoly *mp;
};

static inline double sgn(double v) {
	return v<0 ? -1 : v>0 ? 1 : 0;
}

// This function is only meant to be called on polygons
// that have orthogonal sides on an integer lattice.
void bevel_self_intersections(Mpoly &mp, double amount) {
	if(VERBOSE) {
		printf("Beveling\n");
	} else {
		printf("Beveling: ");
		GDALTermProgress(0, NULL, NULL);
	}

	size_t total_pts = 0;
	for(size_t i=0; i<mp.rings.size(); i++) {
		total_pts += mp.rings[i].pts.size();
	}

	if(VERBOSE) printf("allocating %zd megs for beveler\n",
		(total_pts*sizeof(VertRef)) >> 20);
	std::vector<VertRef> entries;
	entries.reserve(total_pts);
	for(size_t r_idx=0; r_idx<mp.rings.size(); r_idx++) {
		const Ring &ring = mp.rings[r_idx];
		for(size_t v_idx=0; v_idx<ring.pts.size(); v_idx++) {
			if(VERBOSE >= 2) {
				printf("mp[%zd][%zd] = %g, %g\n", r_idx, v_idx,
					ring.pts[v_idx].x, ring.pts[v_idx].y);
			}
			entries.push_back(VertRef(r_idx, v_idx));
		}
	}
	assert(total_pts == entries.size());

	if(VERBOSE) printf("finding self-intersections\n");
	GDALTermProgress(0.1, NULL, NULL);
	// sort by x,y
	std::sort(entries.begin(), entries.end(), CoordsComparator(&mp));
	GDALTermProgress(0.8, NULL, NULL);

	if(VERBOSE >= 2) {
		printf("\nbefore grep:\n");
		for(size_t i=0; i<total_pts; i++) {
			printf("entry[%zd] = %03zd, %03zd (%g, %g)\n", 
				i, entries[i].ring_idx, entries[i].vert_idx,
				entries[i].getVert(mp).x,
				entries[i].getVert(mp).y
				);
		}
	}

	size_t total_num_touch = 0;
	bool prev_was_same = 0;
	for(size_t i=0; i<total_pts-1; i++) {
		//printf("%lf %lf %zd %zd\n",
		//	entries[i].x, entries[i].y,
		//	entries[i].ring_idx, entries[i].vert_idx);
		const Vertex &va = entries[i  ].getVert(mp);
		const Vertex &vb = entries[i+1].getVert(mp);
		if(
			va.x == vb.x &&
			va.y == vb.y
		) {
			if(prev_was_same) {
				fatal_error("should not have triple intersections in beveler");
			}
			entries[total_num_touch++] = entries[i];
			prev_was_same = 1;
		} else {
			prev_was_same = 0;
		}
	}

	if(VERBOSE) printf("found %zd self-intersections\n", total_num_touch);
	if(!total_num_touch) {
		GDALTermProgress(1, NULL, NULL);
		if(VERBOSE) printf("beveler finish\n");
		return;
	}

	entries.resize(total_num_touch);
	// sort by ring_idx,vert_idx
	std::sort(entries.begin(), entries.end(), RingsComparator(&mp));

	if(VERBOSE >= 2) {
		printf("\nafter sort:\n");
		for(size_t i=0; i<total_num_touch; i++) {
			printf("entry[%zd] = %zd, %zd\n", 
				i, entries[i].ring_idx, entries[i].vert_idx);
		}
	}

	if(VERBOSE) printf("shaving corners\n");

	for(size_t entry_idx=0; entry_idx<total_num_touch; ) {
		const Ring &ring = mp.rings[entries[entry_idx].ring_idx];
		const size_t ring_idx = entries[entry_idx].ring_idx;
		size_t ring_num_touch = 0;
		while(
			entry_idx+ring_num_touch < total_num_touch &&
			entries[entry_idx+ring_num_touch].ring_idx == ring_idx
		) {
			ring_num_touch++;
		}

		if(VERBOSE >= 2) printf("ring %zd: num_touch=%zd\n", ring_idx, ring_num_touch);

		size_t new_numpts = ring.pts.size() + ring_num_touch;
		Ring new_ring;
		new_ring.pts.resize(new_numpts);

		size_t vin_idx = 0;
		size_t vout_idx = 0;
		for(size_t subent=0; subent<ring_num_touch; subent++) {
			size_t touch_vidx = entries[entry_idx+subent].vert_idx;
			if(VERBOSE >= 2) printf("touch at %zd\n", touch_vidx);
			if(touch_vidx < vin_idx) {
				fatal_error("verts out of sequence");
			}
			size_t numcp = touch_vidx - vin_idx;
			if(numcp > 0) {
				if(vin_idx + numcp > ring.pts.size()) {
					fatal_error("index out of bounds");
				}
				std::copy(
					ring.pts.begin() + vin_idx,
					ring.pts.begin() + vin_idx + numcp,
					new_ring.pts.begin() + vout_idx
				);
				vout_idx += numcp;
				vin_idx += numcp;
			}

			if(vin_idx != touch_vidx) {
				fatal_error("didn't copy the right amount of points");
			}

			Vertex this_v = ring.pts[vin_idx];
			Vertex prev_v = ring.pts[(vin_idx+ring.pts.size()-1) % ring.pts.size()];
			Vertex next_v = ring.pts[(vin_idx+1) % ring.pts.size()];
			double sx_prev = sgn(prev_v.x - this_v.x);
			double sy_prev = sgn(prev_v.y - this_v.y);
			double sx_next = sgn(next_v.x - this_v.x);
			double sy_next = sgn(next_v.y - this_v.y);
			new_ring.pts[vout_idx++] = Vertex(this_v.x + sx_prev*amount, this_v.y + sy_prev*amount);
			new_ring.pts[vout_idx++] = Vertex(this_v.x + sx_next*amount, this_v.y + sy_next*amount);
			vin_idx++;
		}

		if(vin_idx < ring.pts.size()) {
			size_t numcp = ring.pts.size() - vin_idx;
			std::copy(
				ring.pts.begin() + vin_idx,
				ring.pts.begin() + vin_idx + numcp,
				new_ring.pts.begin() + vout_idx
			);
			vout_idx += numcp;
			vin_idx += numcp;
		}

		if(vout_idx != new_numpts) {
			fatal_error("wrong number of points in beveled ring (%zd vs. %zd)", vout_idx, new_numpts);
		}

		// "ring" variable was const
		Ring &mutable_ring = mp.rings[entries[entry_idx].ring_idx];
		mutable_ring.pts = new_ring.pts;

		entry_idx += ring_num_touch;
	}

	GDALTermProgress(1, NULL, NULL);
	if(VERBOSE) printf("beveler finish\n");
}

} // namespace dangdal
