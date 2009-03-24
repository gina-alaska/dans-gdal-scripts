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
#include "dp.h"

static inline double sqlen(vertex_t v0, vertex_t v1) {
	double dx = v1.x - v0.x;
	double dy = v1.y - v0.y;
	return dx*dx + dy*dy;
}

static void keep_long_segs(ring_t *ring, reduced_ring_t *reduced) {
	int npts = ring->npts;
	vertex_t *pts = ring->pts;
	bool *keep = MYALLOC(bool, npts);
	for(int i=0; i<npts; i++) keep[i] = false;
	double max_sqlen = 400; // FIXME
	for(int i=0; i<reduced->num_segs; i++) {
		int begin = reduced->segs[i].begin;
		int end = reduced->segs[i].end;
		double seg_sqlen = sqlen(pts[begin], pts[end]);
		if(seg_sqlen >= max_sqlen) {
			for(int j=begin; j!=end; j=(j+1)%npts) {
				keep[j] = true;
			}
			keep[end] = true;
		}
	}
	int nkeep = 0;
	for(int i=0; i<npts; i++) if(keep[i]) nkeep++;
	reduced->num_segs = nkeep;
	reduced->segs = MYALLOC(segment_t, nkeep);
	int last_keep = -1, first_keep = -1;
	int s_idx = 0;
	for(int i=0; i<npts; i++) {
		if(!keep[i]) continue;
		if(first_keep < 0) first_keep = i;
		if(last_keep >= 0) {
			reduced->segs[s_idx].begin = last_keep;
			reduced->segs[s_idx].end = i;
			reduced->segs[s_idx].is_problem = 0;
			s_idx++;
		}
		last_keep = i;
	}
	if(last_keep) {
		reduced->segs[s_idx].begin = last_keep;
		reduced->segs[s_idx].end = first_keep;
		reduced->segs[s_idx].is_problem = 0;
		s_idx++;
	}
	if(s_idx != nkeep) {
		fatal_error("s_idx != nkeep");
	}
}

mpoly_t pinch_excursions(mpoly_t *in_mpoly, report_image_t *dbuf) {
	if(!in_mpoly->num_rings) {
		return empty_polygon();
	}

	reduced_ring_t *reduced_rings = MYALLOC(reduced_ring_t, in_mpoly->num_rings);

	double tolerance = 1;
	int c_idx;
	for(c_idx=0; c_idx<in_mpoly->num_rings; c_idx++) {
		reduced_rings[c_idx] = compute_reduced_ring(&in_mpoly->rings[c_idx], tolerance);
		keep_long_segs(&in_mpoly->rings[c_idx], &reduced_rings[c_idx]);
	}

	//fix_topology(in_mpoly, reduced_rings);

	mpoly_t reduced_mpoly = reduction_to_mpoly(in_mpoly, reduced_rings);

	for(int i=0; i<reduced_mpoly.num_rings; i++) {
		debug_plot_ring(dbuf, reduced_mpoly.rings+i, 255, 0, 0);
	}

	return reduced_mpoly;
}
