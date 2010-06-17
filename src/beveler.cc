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

#ifdef SGN
#undef SGN
#endif
#define SGN(a) ((a)<0?-1:(a)>0?1:0)

namespace dangdal {

typedef struct {
	Ring *ring; // FIXME! use ring_idx here
	size_t vert_idx;
} vertsort_entry;

static int compare_coords(const void *ap, const void *bp) {
	vertsort_entry *a = (vertsort_entry *)ap;
	vertsort_entry *b = (vertsort_entry *)bp;
	const Vertex &va = a->ring->pts[a->vert_idx];
	const Vertex &vb = b->ring->pts[b->vert_idx];
	return 
		va.x < vb.x ? -1 :
		va.x > vb.x ?  1 :
		va.y < vb.y ? -1 :
		va.y > vb.y ?  1 :
		0;
}

static int compare_rings(const void *ap, const void *bp) {
	vertsort_entry *a = (vertsort_entry *)ap;
	vertsort_entry *b = (vertsort_entry *)bp;
	return
		a->ring < b->ring ? -1 :
		a->ring > b->ring ?  1 :
		a->vert_idx < b->vert_idx ? -1 :
		a->vert_idx > b->vert_idx ?  1 :
		0;
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
		(total_pts*sizeof(vertsort_entry)) >> 20);
	vertsort_entry *entries = MYALLOC(vertsort_entry, total_pts);
	size_t entry_idx = 0;
	for(size_t r_idx=0; r_idx<mp.rings.size(); r_idx++) {
		Ring &ring = mp.rings[r_idx]; // FIXME! const
		for(size_t v_idx=0; v_idx<ring.pts.size(); v_idx++) {
			vertsort_entry *entry = entries + (entry_idx++);
			entry->ring = &ring;
			entry->vert_idx = v_idx;
		}
	}

	if(VERBOSE) printf("finding self-intersections\n");
	GDALTermProgress(0.1, NULL, NULL);
	// sort by x,y
	qsort(entries, total_pts, sizeof(vertsort_entry), compare_coords);
	GDALTermProgress(0.8, NULL, NULL);

	size_t total_num_touch = 0;
	bool prev_was_same = 0;
	for(size_t i=0; i<total_pts-1; i++) {
		//printf("%lf %lf %zd %zd\n",
		//	entries[i].x, entries[i].y,
		//	entries[i].ring_idx, entries[i].vert_idx);
		const Vertex &va = entries[i  ].ring->pts[entries[i  ].vert_idx];
		const Vertex &vb = entries[i+1].ring->pts[entries[i+1].vert_idx];
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
		free(entries);
		GDALTermProgress(1, NULL, NULL);
		if(VERBOSE) printf("beveler finish\n");
		return;
	}

	// decrease memory usage
	entries = REMYALLOC(vertsort_entry, entries, total_num_touch);
	// sort by ring_idx,vert_idx
	qsort(entries, total_num_touch, sizeof(vertsort_entry), compare_rings);

	if(VERBOSE) printf("shaving corners\n");

	// FIXME! this goes very slow but could probably be made better
	for(entry_idx=0; entry_idx<total_num_touch; ) {
		Ring *ring = entries[entry_idx].ring;
		size_t ring_num_touch = 0;
		while(
			entry_idx+ring_num_touch < total_num_touch &&
			entries[entry_idx+ring_num_touch].ring == ring
		) {
			ring_num_touch++;
		}

		if(VERBOSE >= 2) printf("ring %zd: num_touch=%zd\n", ring - &mp.rings[0], ring_num_touch);

		size_t new_numpts = ring->pts.size() + ring_num_touch;
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
				if(vin_idx + numcp > ring->pts.size()) {
					fatal_error("index out of bounds");
				}
				std::copy(
					ring->pts.begin() + vin_idx,
					ring->pts.begin() + vin_idx + numcp,
					new_ring.pts.begin() + vout_idx
				);
				vout_idx += numcp;
				vin_idx += numcp;
			}

			if(vin_idx != touch_vidx) {
				fatal_error("didn't copy the right amount of points");
			}

			Vertex this_v = ring->pts[vin_idx];
			Vertex prev_v = ring->pts[(vin_idx+ring->pts.size()-1) % ring->pts.size()];
			Vertex next_v = ring->pts[(vin_idx+1) % ring->pts.size()];
			double sx_prev = SGN(prev_v.x - this_v.x);
			double sy_prev = SGN(prev_v.y - this_v.y);
			double sx_next = SGN(next_v.x - this_v.x);
			double sy_next = SGN(next_v.y - this_v.y);
			new_ring.pts[vout_idx++] = Vertex(this_v.x + sx_prev*amount, this_v.y + sy_prev*amount);
			new_ring.pts[vout_idx++] = Vertex(this_v.x + sx_next*amount, this_v.y + sy_next*amount);
			vin_idx++;
		}

		if(vin_idx < ring->pts.size()) {
			size_t numcp = ring->pts.size() - vin_idx;
			std::copy(
				ring->pts.begin() + vin_idx,
				ring->pts.begin() + vin_idx + numcp,
				new_ring.pts.begin() + vout_idx
			);
			vout_idx += numcp;
			vin_idx += numcp;
		}

		if(vout_idx != new_numpts) {
			fatal_error("wrong number of points in beveled ring (%zd vs. %zd)", vout_idx, new_numpts);
		}
		ring->pts = new_ring.pts;

		entry_idx += ring_num_touch;
	}

	free(entries);

	GDALTermProgress(1, NULL, NULL);
	if(VERBOSE) printf("beveler finish\n");
}

} // namespace dangdal
