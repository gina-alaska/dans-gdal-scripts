typedef struct {
	int top_y;
	int bottom_y;
	int *pts;
	int top_linkage;
	int bottom_linkage;
	int *leftray_ids;
	int leftray_count;
	int *lf_refer_ids;
	int lf_refer_count;
	int ring_id;
} descender_t;

typedef struct {
	int num_transitions;
	int *openings;
	int *closings;
	int *descender_ids;
} rowstat_t;

int create_descender_pair(int *num_descenders, descender_t **descenders, int y, int max_pts,
int *leftray_ids, int leftray_count) {
	int d1_idx = *num_descenders;
	int d2_idx = d1_idx+1;
	*descenders = (descender_t *)realloc_or_die(*descenders, sizeof(descender_t)*(d2_idx+1));
	descender_t *d1 = (*descenders) + d1_idx;
	descender_t *d2 = (*descenders) + d2_idx;
	d1->top_y = d2->top_y = y;
	d1->top_linkage = d2_idx;
	d2->top_linkage = d1_idx;
	d1->bottom_y = d2->bottom_y = -1;
	d1->bottom_linkage = -1;
	d2->bottom_linkage = -1;
	d1->pts = (int *)malloc_or_die(sizeof(int) * max_pts);
	d2->pts = (int *)malloc_or_die(sizeof(int) * max_pts);

	// cast a ray to the left - this will be used to sort out containments
	d1->leftray_count = leftray_count;
	if(leftray_count) {
		if(VERBOSE >= 3) {
			int i;
			printf("leftray is");
			for(i=0; i<leftray_count; i++) {
				printf(" %d", leftray_ids[i]);
			}
			printf("\n");
		}
		d1->leftray_ids = (int *)malloc_or_die(sizeof(int) * leftray_count);
		memcpy(d1->leftray_ids, leftray_ids, sizeof(int) * leftray_count);
		int i;
		for(i=0; i<leftray_count; i++) {
			descender_t *ref = (*descenders) + leftray_ids[i];
			ref->lf_refer_ids = (int *)realloc_or_die(
				ref->lf_refer_ids, sizeof(int) * (ref->lf_refer_count+1));
			ref->lf_refer_ids[ref->lf_refer_count++] = d1_idx;
		}
	} else {
		d1->leftray_ids = NULL;
	}
	d2->leftray_count = -1;
	d2->leftray_ids = NULL;

	d1->lf_refer_ids = NULL;
	d1->lf_refer_count = 0;
	d2->lf_refer_ids = NULL;
	d2->lf_refer_count = 0;

	d1->ring_id = -1;
	d2->ring_id = -1;

	(*num_descenders) += 2;
	if(VERBOSE >= 2) printf("num_descenders = %d (y=%d)\n", *num_descenders, y);
	return d1_idx;
}

void close_descender_pair(descender_t *descenders, int d1, int d2, int y, int *ring_id_p) {
	descenders[d1].bottom_linkage = d2;
	descenders[d1].bottom_y = y-1;
	descenders[d1].pts = (int *)realloc_or_die(descenders[d1].pts,
		sizeof(int) * (descenders[d1].bottom_y - descenders[d1].top_y + 1));
	descenders[d2].bottom_linkage = d1;
	descenders[d2].bottom_y = y-1;
	descenders[d2].pts = (int *)realloc_or_die(descenders[d2].pts,
		sizeof(int) * (descenders[d2].bottom_y - descenders[d2].top_y + 1));

	int d_idx;
	d_idx = d1;
	do {
		d_idx = descenders[d_idx].top_linkage;
		d_idx = descenders[d_idx].bottom_linkage;
	} while(d_idx>=0 && d_idx != d1);

	// ring not closed yet - return
	if(d_idx != d1) return;

	int ring_id = (*ring_id_p)++;
	if(VERBOSE >= 3) printf("ring closed: ring_id=%d\n", ring_id);
	int got_rep = 0;
	int flipflop = 1;
	d_idx = d1;
	do {
		descender_t *d = descenders + d_idx;
		if(d->leftray_count >= 0) {
			if(got_rep) {
				if(d->leftray_ids) free(d->leftray_ids);
				d->leftray_ids = NULL;
				d->leftray_count = -1;
			} else {
				got_rep++;
			}
		}

		d->ring_id = ring_id;

		if(flipflop) {
			d_idx = d->top_linkage;
		} else {
			d_idx = d->bottom_linkage;
		}
		flipflop = !flipflop;
	} while(d_idx>=0 && d_idx != d1);
	if(!got_rep) fatal_error("no descenders in ring had leftray_ids");
}

mpoly_t calc_ring_from_mask(unsigned char *mask, int w, int h,
report_image_t *dbuf, int major_ring_only, int no_donuts, 
double min_ring_area, double bevel_size) {
	int x, y;

	if(VERBOSE) printf("finding rings: begin\n");

	int num_rings = 0;

	// up_row is previous row, down_row is current row
	rowstat_t up_row, down_row;
	down_row.num_transitions = 0; // prevent compiler warning;

	// Openings and closings come in pairs having the same index.
	// An opening/closing pair represents entering and then leaving
	// a region of marked pixels in a given row.
	up_row.openings        = (int *)malloc_or_die(sizeof(int) * (w+1)/2);
	up_row.closings        = (int *)malloc_or_die(sizeof(int) * (w+1)/2);
	// Descenders trace the border between marked/unmarked pixels.
	up_row.descender_ids   = (int *)malloc_or_die(sizeof(int) * (w+1));
	down_row.openings      = (int *)malloc_or_die(sizeof(int) * (w+1)/2);
	down_row.closings      = (int *)malloc_or_die(sizeof(int) * (w+1)/2);
	down_row.descender_ids = (int *)malloc_or_die(sizeof(int) * (w+1));

	int num_descenders = 0;
	descender_t *descenders = NULL;

	int show_progress = 1;
	if(show_progress) printf("Tracing: ");

	for(y=0; y<=h; y++) {
		if(show_progress) GDALTermProgress((double)y/(double)(h+1), NULL, NULL);
		if(y) {
			// the previous down_row becomes the new up_row
			up_row.num_transitions = down_row.num_transitions;
			memcpy(up_row.openings, down_row.openings, sizeof(int)*down_row.num_transitions);
			memcpy(up_row.closings, down_row.closings, sizeof(int)*down_row.num_transitions);
			memcpy(up_row.descender_ids, down_row.descender_ids, sizeof(int)*down_row.num_transitions*2);
		} else {
			up_row.num_transitions = 0;
		}

		if(y==h) {
			// This is below the last row of the image.
			// This area is considered to be all unmarked pixels.
			down_row.num_transitions = 0;
		} else {
			// Tabulate the opening/closing pairs for this row.
			down_row.num_transitions = 0;

			int mask_rowlen = (w+7)/8;
			unsigned char mask_bitp = 1;
			unsigned char *mask_bytep = mask + mask_rowlen*y;

			int mask_left = 0;
			for(x=0; x<=w; x++) {
				int mask_right;
				if(x == w) {
					mask_right = 0;
				} else {
					mask_right = (*mask_bytep & mask_bitp) ? 1 : 0;
					mask_bitp <<= 1;
					if(!mask_bitp) {
						mask_bitp = 1;
						mask_bytep++;
					}
				}
				if(mask_right && !mask_left) {
					down_row.openings[down_row.num_transitions] = x;
					//if(VERBOSE) printf("[%d,", x);
				}
				if(mask_left && !mask_right) {
					down_row.closings[down_row.num_transitions++] = x;
					//if(VERBOSE) printf("%d] ", x);
				}
				mask_left = mask_right;
			}
			//if(VERBOSE) printf("\n");
		}

		if(VERBOSE >= 3) {
			int i;
			printf("\nup_row:");
			for(i=0; i<up_row.num_transitions; i++) {
				printf(" [%d-%d,%d-%d]", up_row.openings[i], up_row.closings[i], 
					up_row.descender_ids[i*2], up_row.descender_ids[i*2+1]);
			}
			printf("\ndown_row:");
			for(i=0; i<up_row.num_transitions; i++) {
				printf(" [%d-%d]", down_row.openings[i], down_row.closings[i]);
			}
			printf("\n");
		}

		int up_tid=0, down_tid=0;
		while(up_tid < up_row.num_transitions || down_tid < down_row.num_transitions) {
			//if(VERBOSE) printf("%d/%d:[%d,%d]  %d/%d:[%d,%d]\n",
			//	up_tid, up_row.num_transitions, up_row.openings[up_tid], up_row.closings[up_tid],
			//	down_tid, down_row.num_transitions, down_row.openings[down_tid], down_row.closings[down_tid]);
			//
			if((up_tid < up_row.num_transitions) && // have more transitions in up row and...
				// ran out of transitions in down row or...
				(down_tid == down_row.num_transitions ||
				// the down row opening does not touch the current up row transition
				up_row.closings[up_tid] <= down_row.openings[down_tid])
			) {
				// Close off a pair of descenders:     \____/ 
				if(dbuf && dbuf->mode == PLOT_DESCENDERS) {
					for(x=up_row.openings[up_tid]; x<=up_row.closings[up_tid]; x++)
						plot_point(dbuf, x, y, 255, 0, 0);
				}
				int d1 = up_row.descender_ids[up_tid*2];
				int d2 = up_row.descender_ids[up_tid*2+1];
				if(VERBOSE >= 3) printf("closing pair %d,%d\n", d1, d2);
				close_descender_pair(descenders, d1, d2, y, &num_rings);
				up_tid++;
			} else if((down_tid < down_row.num_transitions) && // have more transitions in down row and...
				// ran out of transitions in up row or...
				(up_tid == up_row.num_transitions ||
				// the up row opening does not touch the current down row transition
				down_row.closings[down_tid] <= up_row.openings[up_tid])
			) {
				// Create a new pair of descenders:     /^^^^\  .
				if(dbuf && dbuf->mode == PLOT_DESCENDERS) {
					for(x=down_row.openings[down_tid]; x<=down_row.closings[down_tid]; x++)
						plot_point(dbuf, x, y, 0, 255, 0);
				}
				int d = create_descender_pair(&num_descenders, &descenders, y, h-y+1,
					down_row.descender_ids, down_tid*2);
				if(VERBOSE >= 3) {
					printf("opening pair %d,%d at %d-%d\n", d, d+1,
						down_row.openings[down_tid], down_row.closings[down_tid]);
				}
				descenders[d  ].pts[0] = down_row.openings[down_tid];
				descenders[d+1].pts[0] = down_row.closings[down_tid];
				down_row.descender_ids[down_tid*2  ] = d;
				down_row.descender_ids[down_tid*2+1] = d+1;

				down_tid++;
			} else if(up_tid < up_row.num_transitions && down_tid < down_row.num_transitions) {
				if(dbuf && dbuf->mode == PLOT_DESCENDERS) {
					plot_point(dbuf, up_row.openings[up_tid], y, 255, 255, 0);
					plot_point(dbuf, down_row.openings[down_tid], y, 255, 255, 0);
				}
				// Link the opening descender in the down row to the one in the up row
				// and add the new boundary X value.
				int dl = up_row.descender_ids[up_tid*2];
//if(VERBOSE) printf("dl=%d\n", dl);
				down_row.descender_ids[down_tid*2] = dl;
				descenders[dl].pts[y-descenders[dl].top_y] = down_row.openings[down_tid];
				for(;;) {
					if(
						(up_tid < up_row.num_transitions-1) &&
						(up_row.openings[up_tid+1] < down_row.closings[down_tid])
					) {
						// Close off a pair of descenders:   | \____/ |
						if(dbuf && dbuf->mode == PLOT_DESCENDERS) {
							for(x=up_row.closings[up_tid]; x<=up_row.openings[up_tid+1]; x++)
								plot_point(dbuf, x, y, 255, 0, 255);
						}
						int d1 = up_row.descender_ids[up_tid*2+1];
						int d2 = up_row.descender_ids[up_tid*2+2];
						if(VERBOSE >= 3) printf("closing inner pair %d,%d\n", d1, d2);
						close_descender_pair(descenders, d1, d2, y, &num_rings);
						up_tid++;
					} else if(
						(down_tid < down_row.num_transitions-1) &&
						(down_row.openings[down_tid+1] < up_row.closings[up_tid])
					) {
						// Create a new pair of descenders:   | /^^^^\ |
						if(dbuf && dbuf->mode == PLOT_DESCENDERS) {
							for(x=down_row.closings[down_tid]; x<=down_row.openings[down_tid+1]; x++)
								plot_point(dbuf, x, y, 0, 255, 255);
						}
						int d = create_descender_pair(&num_descenders, &descenders, y, h-y+1,
							down_row.descender_ids, down_tid*2+1);
						if(VERBOSE >= 3) {
							printf("opening inner pair %d,%d at %d-%d\n", d, d+1,
								down_row.closings[down_tid], down_row.openings[down_tid+1]);
						}
						descenders[d  ].pts[0] = down_row.closings[down_tid];
						descenders[d+1].pts[0] = down_row.openings[down_tid+1];
						down_row.descender_ids[down_tid*2+1] = d;
						down_row.descender_ids[down_tid*2+2] = d+1;
						down_tid++;
					} else break;
				}
				if(dbuf && dbuf->mode == PLOT_DESCENDERS) {
					plot_point(dbuf, up_row.closings[up_tid], y, 255, 255, 0);
					plot_point(dbuf, down_row.closings[down_tid], y, 255, 255, 0);
				}
				// Link the closing descender in the down row to the one in the up row
				// and add the new boundary X value.
				int dr = up_row.descender_ids[up_tid*2+1];
				down_row.descender_ids[down_tid*2+1] = dr;
				descenders[dr].pts[y-descenders[dr].top_y] = down_row.closings[down_tid];
				up_tid++;
				down_tid++;
			}
		}
	}
	if(show_progress) GDALTermProgress(1, NULL, NULL);

	int i;

	if(VERBOSE >= 2) {
		for(i=0; i<num_descenders; i++) printf("%d: top_link=%d bottom_link=%d\n",
			i, descenders[i].top_linkage, descenders[i].bottom_linkage);
	}

	if(!num_descenders) {
		printf("Mask was completely blank - therefore there is no bounding polygon\n");
		return empty_polygon();
	}

	unsigned char *used_desc = (unsigned char *)malloc_or_die(num_descenders);
	for(i=0; i<num_descenders; i++) used_desc[i] = 0;
	int num_desc_used = 0;

	mpoly_t mp = empty_polygon();
	descender_t **representative_descenders = NULL;

	if(show_progress) printf("Joining segments: ");
	int start_d = 0;
	int ring_id = 0;
	for(;;) {
		if(show_progress) GDALTermProgress((double)num_desc_used/(double)num_descenders, NULL, NULL);
		for(; start_d<num_descenders; start_d++) {
			if(!used_desc[start_d]) {
				break;
			}
		}
		if(start_d == num_descenders) break;
		ring_t ring;
		ring.npts = 0;
		ring.pts = NULL;
		// prevent compiler warning - these fields are filled in by compute_containments
		ring.parent_id = ring.is_hole = 0;
		descender_t *ring_representative_descender;

		int cur_d = start_d;
		do {
			if(VERBOSE >= 2) printf("d:%d ring=%d ", cur_d, ring_id);
			descender_t *d = descenders + cur_d;
			if(d->leftray_count >= 0) ring_representative_descender = d;

			if(used_desc[cur_d]) fatal_error("descender used twice");
			used_desc[cur_d]++; num_desc_used++;
			d->ring_id = ring_id;

			int n = d->bottom_y - d->top_y + 1;
			if(n <= 0) fatal_error("npts <= 0 in ring segment");
			ring.pts = (vertex_t *)realloc_or_die(ring.pts, 
				sizeof(vertex_t)*(ring.npts+n*2));
			ring.pts[ring.npts++] = (vertex_t){ d->pts[0], d->top_y };
			for(i=1; i<n; i++) {
				if(d->pts[i] != d->pts[i-1]) {
					ring.pts[ring.npts++] = 
						(vertex_t){ d->pts[i-1], d->top_y + i };
					ring.pts[ring.npts++] =
						(vertex_t){ d->pts[i  ], d->top_y + i };
				}
			}
			ring.pts[ring.npts++] = (vertex_t){ d->pts[n-1], d->bottom_y+1 };

			cur_d = d->bottom_linkage;
			if(cur_d < 0) fatal_error("uninitialized val in descender linkage");

			if(VERBOSE >= 2) printf("u:%d ", cur_d);
			d = descenders + cur_d;
			if(d->leftray_count >= 0) ring_representative_descender = d;

			if(used_desc[cur_d]) fatal_error("descender used twice");
			used_desc[cur_d]++; num_desc_used++;
			d->ring_id = ring_id;

			n = d->bottom_y - d->top_y + 1;
			if(n <= 0) fatal_error("npts <= 0 in ring segment");
			ring.pts = (vertex_t *)realloc_or_die(ring.pts, 
				sizeof(vertex_t)*(ring.npts+n*2));
			ring.pts[ring.npts++] = (vertex_t){ d->pts[n-1], d->bottom_y + 1 };
			for(i=n-2; i>=0; i--) {
				if(d->pts[i] != d->pts[i+1]) {
					ring.pts[ring.npts++] = 
						(vertex_t){ d->pts[i+1], d->top_y + i+1 };
					ring.pts[ring.npts++] = 
						(vertex_t){ d->pts[i  ], d->top_y + i+1 };
				}
			}
			ring.pts[ring.npts++] = (vertex_t){ d->pts[0], d->top_y };
			cur_d = d->top_linkage;
			if(cur_d < 0) fatal_error("uninitialized val in descender linkage");
		} while(cur_d != start_d);
		if(VERBOSE >= 2) printf("\n");

		if(VERBOSE >= 2) printf("ring %d: %d pts\n", mp.num_rings, ring.npts);

		mp.rings = (ring_t *)realloc_or_die(mp.rings, sizeof(ring_t)*(ring_id+1));
		mp.rings[ring_id] = ring;
		representative_descenders = (descender_t **)realloc_or_die(
			representative_descenders, sizeof(descender_t *)*(ring_id+1));
		representative_descenders[ring_id] = ring_representative_descender;
		mp.num_rings = ++ring_id;
	}
	for(i=0; i<num_descenders; i++) {
		free(descenders[i].pts);
	}
	if(show_progress) GDALTermProgress(1, NULL, NULL);

	if(VERBOSE) printf("finding rings: end\n");

	if(VERBOSE) {
		int total_pts = 0;
		int r_idx;
		for(r_idx=0; r_idx<mp.num_rings; r_idx++) {
			total_pts += mp.rings[r_idx].npts;
		}
		printf("tracer produced %d rings with a total of %d points\n",
			mp.num_rings, total_pts);
	}

	if(min_ring_area > 0) {
		if(VERBOSE) printf("removing small rings...\n");

		ring_t *filtered_rings = (ring_t *)malloc_or_die(sizeof(ring_t)*mp.num_rings);
		int num_filtered_rings = 0;
		for(i=0; i<mp.num_rings; i++) {
			double area = polygon_area(mp.rings+i);
			if(VERBOSE) if(area > 10) printf("ring %d has area %.15f\n", i, area);
			if(area >= min_ring_area) {
				filtered_rings[num_filtered_rings++] = mp.rings[i];
			} else {
				free_ring(mp.rings + i);
			}
		}
		printf("filtered by area %d => %d rings\n",
			mp.num_rings, num_filtered_rings);

		free(mp.rings);
		mp.rings = filtered_rings;
		mp.num_rings = num_filtered_rings;
	}

	if(major_ring_only) {
		double biggest_area = 0;
		int best_idx = 0;
		for(i=0; i<mp.num_rings; i++) {
			double area = polygon_area(mp.rings+i);
			if(area > biggest_area) {
				biggest_area = area;
				best_idx = i;
			}
		}
		if(VERBOSE) printf("major ring was %d with %d pts, %.1f area\n",
			best_idx, mp.rings[best_idx].npts, biggest_area);
		mp.rings = mp.rings+best_idx;
		mp.num_rings = 1;
	}

	//printf("computing containments: begin\n");
#if 0
	printf("old containment method\n");
	compute_containments(&mp);
#else
	if(show_progress) printf("Computing containments: ");
	int *crossing_counts = (int *)malloc_or_die(sizeof(int) * mp.num_rings);
	int r_idx;
	for(r_idx=0; r_idx<mp.num_rings; r_idx++) crossing_counts[r_idx] = 0;
	for(r_idx=0; r_idx<mp.num_rings; r_idx++) {
		ring_t *ring = mp.rings + r_idx;

		if(show_progress) GDALTermProgress((double)r_idx/(double)mp.num_rings, NULL, NULL);

		int total_crossings = 0;

		if(VERBOSE >= 3) printf("ring %d: ", r_idx);
		descender_t *rep_desc = representative_descenders[r_idx];
		if(!rep_desc) fatal_error("no representative descender for ring %d\n", r_idx);
		int lr_idx;
		for(lr_idx=0; lr_idx<rep_desc->leftray_count; lr_idx++) {
			int d_idx = rep_desc->leftray_ids[lr_idx];
			int d_ring_id = descenders[d_idx].ring_id;
			if(VERBOSE >= 3) printf(" [d=%d,r=%d]", d_idx, d_ring_id);
			if(d_ring_id != r_idx) {
				crossing_counts[d_ring_id]++;
				total_crossings++;
			}
		}
		if(VERBOSE >= 3) printf("\n");

		ring->is_hole = total_crossings % 2;
		ring->parent_id = -1;

		for(lr_idx=rep_desc->leftray_count-1; lr_idx>=0; lr_idx--) {
			int d_idx = rep_desc->leftray_ids[lr_idx];
			int d_ring_id = descenders[d_idx].ring_id;
			int ncros = crossing_counts[d_ring_id];
			if(ncros % 2) {
				ring->parent_id = d_ring_id;
				break;
			}
		}

		// fast way to clear array
		for(lr_idx=0; lr_idx<rep_desc->leftray_count; lr_idx++) {
			int d_idx = rep_desc->leftray_ids[lr_idx];
			int d_ring_id = descenders[d_idx].ring_id;
			crossing_counts[d_ring_id] = 0;
		}

		if(VERBOSE >= 2) {
			printf("ring %d is_hole=%d parent_id=%d\n", r_idx, ring->is_hole, ring->parent_id);
		}
	}
	for(r_idx=0; r_idx<mp.num_rings; r_idx++) {
		ring_t *ring = mp.rings + r_idx;
		if(ring->parent_id >= 0) {
			ring_t *parent = mp.rings + ring->parent_id;
			if(ring->is_hole == parent->is_hole) {
				fatal_error("topology error in containment calculation");
			}
		} else {
			if(ring->is_hole) {
				fatal_error("topology error in containment calculation");
			}
		}
	}
	free(crossing_counts);
	if(show_progress) GDALTermProgress(1, NULL, NULL);
#endif
	//printf("computing containments: end\n");

	free(used_desc);
	for(i=0; i<num_descenders; i++) {
		if(descenders[i].leftray_ids) {
			free(descenders[i].leftray_ids);
		}
	}
	free(descenders);

	if(no_donuts) {
		// we don't have to worry about remapping parent_id in this
		// case because we only take rings with no parent
		int out_idx = 0;
		for(i=0; i<mp.num_rings; i++) {
			if(mp.rings[i].parent_id < 0) {
				mp.rings[out_idx++] = mp.rings[i];
			}
		}
		mp.num_rings = out_idx;
		if(VERBOSE) printf("number of non-donut rings is %d", mp.num_rings);
	}

	if(bevel_size > 0) {
		// the topology cannot be resolved by us or by geos/jump/postgis if
		// there are self-intersections
		bevel_self_intersections(&mp, bevel_size);
	}

	return mp;
}
