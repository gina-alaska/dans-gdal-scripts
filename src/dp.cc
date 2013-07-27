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



#include <vector>
#include <cassert>

#include "common.h"
#include "polygon.h"
#include "dp.h"

namespace dangdal {

// Implementation of Douglas-Peucker polyline reduction algorithm
// rewrite of code from http://www.3dsoftware.com/Cartography/Programming/PolyLineReduction
// (and adapted from src/linework/dp.c in SwathViewer)
// Modifications were made to ensure valid topology of output

static const double EPSILON = 10E-10;

// FIXME - it is possible to factor out the sqrt call
static inline double veclen(double x, double y) {
	return sqrt(x*x + y*y);
}

Mpoly compute_reduced_pointset(const Mpoly &in_mpoly, double tolerance) {
	if(VERBOSE) printf("reducing...\n");

	if(!in_mpoly.rings.size()) {
		return Mpoly();
	}

	std::vector<ReducedRing> reduced_rings(in_mpoly.rings.size());

	for(size_t c_idx=0; c_idx<in_mpoly.rings.size(); c_idx++) {
		reduced_rings[c_idx] = compute_reduced_ring(
			in_mpoly.rings[c_idx], tolerance);
	}

	fix_topology(in_mpoly, reduced_rings);

	return reduction_to_mpoly(in_mpoly, reduced_rings);
}

static inline double get_dist_to_seg(
	double seg_vec_x, double seg_vec_y, 
	Vertex seg_vert1, Vertex seg_vert2, Vertex test_vert
) {
	double vert_vec_x = test_vert.x - seg_vert1.x;
	double vert_vec_y = test_vert.y - seg_vert1.y;
	double scalar_prod = vert_vec_x*seg_vec_x + vert_vec_y*seg_vec_y;
	if(scalar_prod < 0.0) {
		// past beginning of segment - distance from segment
		// is equal to distance from first endpoint of segment
		return veclen(vert_vec_x, vert_vec_y);
	} else {
		vert_vec_x = seg_vert2.x - test_vert.x;
		vert_vec_y = seg_vert2.y - test_vert.y;
		scalar_prod = vert_vec_x*seg_vec_x + vert_vec_y*seg_vec_y;
		if(scalar_prod < 0.0) {
			// past end of segment - distance from segment
			// is equal to distance from second endpoint of segment
			return veclen(vert_vec_x, vert_vec_y);
		} else {
			double c_squared = vert_vec_x*vert_vec_x + vert_vec_y*vert_vec_y;
			double a_squared = scalar_prod * scalar_prod;
			double b_squared = c_squared - a_squared;
			if(b_squared < 0) {
				if(fabs(b_squared) < a_squared * EPSILON) {
					b_squared = 0; // correct for round-off error
				} else {
					fatal_error("a_squared > c_squared");
				}
			}
			// use pythagorean theorem to find distance from
			// vertex to segment
			return sqrt(b_squared);
		}
	}
}

ReducedRing compute_reduced_ring(const Ring &orig_string, double res) {
//printf("enter dp\n");

	const std::vector<Vertex> &pts_in = orig_string.pts;
	const size_t num_in = pts_in.size();

	double tolerance = res;
	int i;

	std::vector<segment_t> stack(num_in);
	size_t stack_ptr = 0;

	ReducedRing keep;
	keep.segs.reserve(num_in);

	// must keep closure segment
	keep.segs.push_back(segment_t(num_in-1, 0));

	stack[stack_ptr].begin = 0;
	stack[stack_ptr].end = num_in-1;
	stack_ptr++;
	while(stack_ptr) {
		stack_ptr--;
		int seg_begin = stack[stack_ptr].begin;
		int seg_end = stack[stack_ptr].end;
//printf("stack_ptr=%d / range=[%d,%d]\n", stack_ptr, seg_begin, seg_end);

		double max_dist = -1.0;
		int idx_of_max = -1; // to prevent compiler warning

		double seg_vec_x = pts_in[seg_end].x - pts_in[seg_begin].x;
		double seg_vec_y = pts_in[seg_end].y - pts_in[seg_begin].y;
		double seg_vec_len = veclen(seg_vec_x, seg_vec_y);
		if(seg_vec_len > 0.0) {
			// normalize vector
			seg_vec_x /= seg_vec_len;
			seg_vec_y /= seg_vec_len;
			for(i=seg_begin+1; i<seg_end; i++) {
				double dist_to_seg = get_dist_to_seg(seg_vec_x, seg_vec_y,
					pts_in[seg_begin], pts_in[seg_end], pts_in[i]);
				if(dist_to_seg < 0.0) fatal_error("dist_to_seg < 0.0");
				if(std::isnan(dist_to_seg)) fatal_error("dist_to_seg == NaN");

				if(dist_to_seg > max_dist) {
					max_dist = dist_to_seg;
					idx_of_max = i;
				}
			}
		} else {
			// Segment is length zero, so we can't use get_dist_to_seg.
			// Instead, just use cartesian distance
			for(i=seg_begin+1; i<seg_end; i++) {
				double dx = pts_in[i].x - pts_in[seg_begin].x;
				double dy = pts_in[i].y - pts_in[seg_begin].y;
				double dist_to_seg = veclen(dx, dy);

				if(dist_to_seg > max_dist) {
					max_dist = dist_to_seg;
					idx_of_max = i;
				}
			}
		}

//printf("max=%.15f, toler=%.15f, idx=%i\n", max_dist, tolerance, idx_of_max);
		if(max_dist >= tolerance) {
			if(idx_of_max <= 0) fatal_error(
				"idx_of_max out of range (perhaps it wasn't set?)");

			// add point and recursively divide subsegments
			stack[stack_ptr].begin = seg_begin;
			stack[stack_ptr].end = idx_of_max;
			stack_ptr++;
			if(stack_ptr >= num_in) fatal_error("stack overflow in dp.c");

			stack[stack_ptr].begin = idx_of_max;
			stack[stack_ptr].end = seg_end;
			stack_ptr++;
			if(stack_ptr >= num_in) fatal_error("stack overflow in dp.c");
		} else {
			// segment doesn't need subdivision - tag
			// endpoint for inclusion
			if(keep.segs.size() == num_in) fatal_error("output stack overflow in dp.c");
			keep.segs.push_back(segment_t(seg_begin, seg_end));
		}
	}

//printf("keeping %d of %d\n", keep.segs.size(), num_in);
//printf("exit dp\n");

	return keep;
}

static Ring make_ring_from_segs(const Ring &c_in, const ReducedRing &r_in) {
	size_t in_npts = c_in.pts.size();
	std::vector<uint8_t> keep_pts(in_npts, 0);

	for(size_t i=0; i<r_in.segs.size(); i++) {
		keep_pts[r_in.segs[i].begin]++;
		keep_pts[r_in.segs[i].end]++;
	}

	size_t num_to_keep = 0;
	for(size_t i=0; i<in_npts; i++) {
		if(keep_pts[i] && keep_pts[i] != 2) {
			fatal_error("point must be in 0 or 2 segs");
		}
		if(keep_pts[i]) num_to_keep++;
	}
	// number of points == number of line segments
	if(num_to_keep != r_in.segs.size()) fatal_error("num_to_keep != r_in.segs.size()");

//printf("keeping %d of %d\n", num_to_keep, num_in);

	Ring new_string = c_in.copyMetadata();
	new_string.pts.reserve(num_to_keep);
	for(size_t i=0; i<in_npts; i++) {
		if(keep_pts[i]) {
			new_string.pts.push_back(c_in.pts[i]);
		}
	}
	if(new_string.pts.size() != num_to_keep) {
		fatal_error("count mismatch after point copy in dp.c");
	}

	return new_string;
}

Mpoly reduction_to_mpoly(const Mpoly &in_mpoly, const std::vector<ReducedRing> &reduced_rings) {
	// This function takes only rings with more than two points.
	// (a polygon with two or less points has no area)
	// Care is taken to make sure that the containment info
	// (specifically the parent_id pointer) is still proper.
	// It is difficult to come up with a test case so I just
	// have to hope that is works.

	// First, find which rings have at least three points.
	std::vector<bool> keep_rings(in_mpoly.rings.size());
	size_t total_npts_in=0, total_npts_out=0;
	for(size_t c_idx=0; c_idx<in_mpoly.rings.size(); c_idx++) {
		const Ring &c_in = in_mpoly.rings[c_idx];
		const ReducedRing &r_in = reduced_rings[c_idx];

		if(VERBOSE >= 2) {
			printf("ring %zd: %zd => %zd pts\n", c_idx, c_in.pts.size(), r_in.segs.size());
		}
		total_npts_in += c_in.pts.size();

		if(r_in.segs.size() > 2) {
			total_npts_out += r_in.segs.size();
			keep_rings[c_idx] = 1;
		} else {
			keep_rings[c_idx] = 0;
		}
	}

	// If an outer ring has been removed, remove its holes
	// also.  Extremely unlikely but just in case...
	for(size_t c_idx=0; c_idx<in_mpoly.rings.size(); c_idx++) {
		const Ring &c_in = in_mpoly.rings[c_idx];
		if(c_in.parent_id >= 0) {
			if(!keep_rings[c_in.parent_id]) {
				keep_rings[c_idx] = 0;
			}
		}
	}

	// Map index of old rings to index of new rings.
	std::vector<int> new_idx_map(in_mpoly.rings.size());
	int out_idx = 0;
	for(size_t c_idx=0; c_idx<in_mpoly.rings.size(); c_idx++) {
		if(keep_rings[c_idx]) {
			//printf("map ring %d => %d\n", c_idx, out_idx);
			new_idx_map[c_idx] = out_idx++;
		} else {
			new_idx_map[c_idx] = -1;
		}
	}
	size_t num_out_rings = out_idx;

	if(!num_out_rings) return Mpoly();

	// Now create the new rings and update parent_id using the
	// remapped index values.
	Mpoly out_mp;
	out_mp.rings.resize(num_out_rings);
	for(size_t c_idx=0; c_idx<in_mpoly.rings.size(); c_idx++) {
		const Ring &c_in = in_mpoly.rings[c_idx];
		const ReducedRing &r_in = reduced_rings[c_idx];

		out_idx = new_idx_map[c_idx];
		if(out_idx < 0) continue;

		Ring new_string = make_ring_from_segs(c_in, r_in);

		int old_parent_id = in_mpoly.rings[c_idx].parent_id;
		if(old_parent_id >= 0) {
			int new_parent_id = new_idx_map[old_parent_id];
			//printf("map parent_id %d => %d\n", old_parent_id, new_parent_id);
			new_string.parent_id = new_parent_id;
		}

		out_mp.rings[out_idx] = new_string;
	}

	if(VERBOSE) printf("reduced %zd => %zd rings, %zd => %zd pts\n",
		in_mpoly.rings.size(), num_out_rings, total_npts_in, total_npts_out);

	return out_mp;
}

static inline int segs_cross(
	bool same_ring, 
	const Ring &c1, const segment_t &s1,
	const Ring &c2, const segment_t &s2
) {
	if(same_ring) {
		// don't test crossing if segments are identical or neighbors
		if(
			(s1.begin == s2.begin) ||
			(s1.begin == s2.end) ||
			(s1.end == s2.begin) ||
			(s1.end == s2.end)
		) return 0;
	}
	
	return line_intersects_line(
		c1.pts[s1.begin], c1.pts[s1.end],
		c2.pts[s2.begin], c2.pts[s2.end],
		1);
}

// The std::pair consists of the ring and segment indices.
BboxBinarySpacePartition<std::pair<size_t, size_t> > get_bsp_for_reduced_rings(
	const Mpoly &mpoly, std::vector<ReducedRing> &reduced_rings
) {
	typedef std::pair<size_t, size_t> segptr_t;
	std::vector<std::pair<Bbox, segptr_t> > items;

	for(size_t ring_idx=0; ring_idx < mpoly.rings.size(); ring_idx++) {
		const Ring &ring = mpoly.rings[ring_idx];
		const ReducedRing &reduced = reduced_rings[ring_idx];
		for(size_t seg_idx=0; seg_idx < reduced.segs.size(); seg_idx++) {
			Bbox bbox = reduced.segs[seg_idx].get_bbox(ring);
			items.push_back(std::make_pair(bbox, segptr_t(ring_idx, seg_idx)));
		}
	}

	return BboxBinarySpacePartition<segptr_t>(items);
}

void fix_topology(const Mpoly &mpoly, std::vector<ReducedRing> &reduced_rings) {
	const double firsthalf_progress = 0.5;
	printf("Fixing topology: ");
	fflush(stdout);

	assert(mpoly.rings.size() == reduced_rings.size());

	// initialize problem arrays
	std::vector<std::vector<bool> > mp_problems(mpoly.rings.size());
	for(size_t r1_idx=0; r1_idx < mpoly.rings.size(); r1_idx++) {
		const ReducedRing &rring = reduced_rings[r1_idx];
		mp_problems[r1_idx].resize(rring.segs.size(), 0);
	}

	// FIXME - no longer needed
	std::vector<Bbox> bboxes = mpoly.getRingBboxes();

	printf("** make tree\n"); fflush(stdout);
	const BboxBinarySpacePartition<std::pair<size_t, size_t> > bsp =
		get_bsp_for_reduced_rings(mpoly, reduced_rings);
	printf("** done making tree\n"); fflush(stdout);

	// flag segments that cross
	int have_problems = 0;
	for(size_t r1_idx=0; r1_idx < mpoly.rings.size(); r1_idx++) {
		GDALTermProgress(firsthalf_progress*
			pow((double)r1_idx / (double)mpoly.rings.size(), 2), NULL, NULL);
		const Ring &c1 = mpoly.rings[r1_idx];
		const ReducedRing &r1 = reduced_rings[r1_idx];
		std::vector<bool> &p1 = mp_problems[r1_idx];

		for(size_t seg1_idx=0; seg1_idx < r1.segs.size(); seg1_idx++) {
			Bbox seg1_bbox = r1.segs[seg1_idx].get_bbox(c1);
			std::vector<std::pair<size_t, size_t> > intersecting_segments =
				bsp.get_intersecting_items(seg1_bbox);
			for(size_t i_s_idx=0; i_s_idx < intersecting_segments.size(); i_s_idx++) {
				size_t r2_idx = intersecting_segments[i_s_idx].first;
				size_t seg2_idx = intersecting_segments[i_s_idx].second;
				if(r2_idx > r1_idx) continue; // symmetry optimization
				if(r2_idx == r1_idx && seg2_idx > seg1_idx) continue; // symmetry optimization

				const Ring &c2 = mpoly.rings[r2_idx];
				const ReducedRing &r2 = reduced_rings[r2_idx];
				std::vector<bool> &p2 = mp_problems[r2_idx];

				int crosses = segs_cross(r1_idx==r2_idx, 
					c1, r1.segs[seg1_idx], c2, r2.segs[seg2_idx]);
				if(crosses) {
					//printf("found a crossing: %d,%d,%d,%d\n",
					//	r1_idx, seg1_idx, r2_idx, seg2_idx);
					p1[seg1_idx] = 1;
					p2[seg2_idx] = 1;
					have_problems += 2;
				}
			} // ring2/seg2 loop
		} // seg1 loop
	} // ring1 loop

	double progress = firsthalf_progress;
	GDALTermProgress(progress, NULL, NULL);

	if(have_problems) {
		if(VERBOSE) printf("fixing %d crossed segments from reduction\n", have_problems/2);
	}

	int did_something = 1;
	while(have_problems && did_something) {
		//printf("%d crossings to fix\n", have_problems/2);
		did_something = 0;
		// subdivide problem segments
		for(size_t r1_idx=0; r1_idx < mpoly.rings.size(); r1_idx++) {
			ReducedRing &r1 = reduced_rings[r1_idx];
			std::vector<bool> &p1 = mp_problems[r1_idx];
			size_t orig_num_segs = r1.segs.size(); // this number will change as we go, so copy it
			for(size_t seg1_idx=0; seg1_idx < orig_num_segs; seg1_idx++) {
				if(!p1[seg1_idx]) continue;
				int begin = r1.segs[seg1_idx].begin;
				int end = r1.segs[seg1_idx].end;
				if(end < begin) continue; // this is a closure segment (and so doesn't have a midpoint)
				if(end == begin+1) continue; // no midpoint

				// subdivide this segment
				int mid = (begin + end) / 2;
				r1.segs[seg1_idx].end = mid;
				r1.segs.push_back(segment_t(mid, end));
				p1.push_back(1);
				did_something = 1;
			} // seg loop
		} // ring loop

		have_problems = 0;
		// now test for resolved problems and new problems
		for(size_t r1_idx=0; r1_idx < mpoly.rings.size(); r1_idx++) {
			{
				double alpha = double(r1_idx) / mpoly.rings.size();
				double p = progress + (1.0-progress)/2*alpha;
				GDALTermProgress(p, NULL, NULL);
			}

			const Ring &c1 = mpoly.rings[r1_idx];
			const ReducedRing &r1 = reduced_rings[r1_idx];
			const Bbox bbox1 = bboxes[r1_idx];
			std::vector<bool> &p1 = mp_problems[r1_idx];
			for(size_t seg1_idx=0; seg1_idx < r1.segs.size(); seg1_idx++) {
				if(!p1[seg1_idx]) continue;
				p1[seg1_idx] = 0;
				// FIXME - use BSP here as well
				for(size_t r2_idx=0; r2_idx < mpoly.rings.size(); r2_idx++) {
					const Ring &c2 = mpoly.rings[r2_idx];
					const ReducedRing &r2 = reduced_rings[r2_idx];
					const Bbox bbox2 = bboxes[r2_idx];
					if(is_disjoint(bbox1, bbox2)) continue;
					for(size_t seg2_idx=0; seg2_idx < r2.segs.size(); seg2_idx++) {
						int crosses = segs_cross(r1_idx==r2_idx,
							c1, r1.segs[seg1_idx], c2, r2.segs[seg2_idx]);
						if(crosses) {
							if(VERBOSE) {
								printf("found a crossing (still): %zd,%zd,%zd,%zd (%f,%f)-(%f,%f) (%f,%f)-(%f,%f)\n",
									r1_idx, seg1_idx, r2_idx, seg2_idx,
									c1.pts[r1.segs[seg1_idx].begin].x,
									c1.pts[r1.segs[seg1_idx].begin].y,
									c1.pts[r1.segs[seg1_idx].end].x,
									c1.pts[r1.segs[seg1_idx].end].y,
									c2.pts[r2.segs[seg2_idx].begin].x,
									c2.pts[r2.segs[seg2_idx].begin].y,
									c2.pts[r2.segs[seg2_idx].end].x,
									c2.pts[r2.segs[seg2_idx].end].y);
							}
							p1[seg1_idx] = 1;
							std::vector<bool> &p2 = mp_problems[r2_idx];
							p2[seg2_idx] = 1;
							have_problems++;
						}
					} // seg loop
				} // ring loop
			} // seg loop
		} // ring loop

		progress += (1.0-progress)/2;
	} // while problems

	GDALTermProgress(1, NULL, NULL);

	if(have_problems) {
		printf("WARNING: Could not fix all topology problems.\n  Please inspect output shapefile manually.\n");
	}
}

} // namespace dangdal
