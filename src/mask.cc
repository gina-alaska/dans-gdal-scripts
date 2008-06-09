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



#include "common.h"
#include "mask.h"
#include "polygon.h"
#include "debugplot.h"
#include "ndv.h"

uint8_t *get_mask_for_dataset(GDALDatasetH ds, int bandlist_size, int *bandlist, 
ndv_def_t *ndv_def, report_image_t *dbuf) {
	int i, j;

	int w = GDALGetRasterXSize(ds);
	int h = GDALGetRasterYSize(ds);
	int band_count = GDALGetRasterCount(ds);
	if(VERBOSE) printf("input is %d x %d x %d\n", w, h, band_count);

	if(VERBOSE) printf("mask array is %.1f megabytes\n", (double)(w+2)*(h+2)/1024.0/1024.0);
	uint8_t *mask = (uint8_t *)malloc_or_die((w+2)*(h+2));
	memset(mask, 0, (w+2)*(h+2));

	printf("Reading %d bands of size %d x %d\n", bandlist_size, w, h);

	int bandlist_idx;
	for(bandlist_idx=0; bandlist_idx<bandlist_size; bandlist_idx++) {
		int band_idx = bandlist[bandlist_idx];
		if(band_idx < 1 || band_idx > band_count) fatal_error("bandid out of range");

		GDALRasterBandH band = GDALGetRasterBand(ds, band_idx);

		int blocksize_x, blocksize_y;
		GDALGetBlockSize(band, &blocksize_x, &blocksize_y);

		// gain speed in common 8-bit case
		GDALDataType gdt = GDALGetRasterDataType(band);
		int use_8bit = (gdt == GDT_Byte);

		if(VERBOSE) printf("band %d: block size = %d,%d, use_8bit=%d\n",
			band_idx, blocksize_x, blocksize_y, use_8bit);

		void *block_buf;
		if(use_8bit) {
			block_buf = malloc_or_die(blocksize_x*blocksize_y);
		} else {
			block_buf = malloc_or_die(blocksize_x*blocksize_y*sizeof(double));
		}

		uint8_t *row_ndv = (uint8_t *)malloc(blocksize_x);

		int boff_x, boff_y;
		for(boff_y=0; boff_y<h; boff_y+=blocksize_y) {
			int bsize_y = blocksize_y;
			if(bsize_y + boff_y > h) bsize_y = h - boff_y;
			for(boff_x=0; boff_x<w; boff_x+=blocksize_x) {
				int bsize_x = blocksize_x;
				if(bsize_x + boff_x > w) bsize_x = w - boff_x;

				double progress = 
					((double)bandlist_idx * (double)w * (double)h +
					(double)boff_y * (double)w +
					(double)boff_x * (double)bsize_y) /
					((double)bandlist_size * (double)w * (double)h);
				GDALTermProgress(progress, NULL, NULL);

				GDALRasterIO(band, GF_Read, boff_x, boff_y, bsize_x, bsize_y, 
					block_buf, bsize_x, bsize_y, 
					use_8bit ? GDT_Byte : GDT_Float64,
					0, 0);

				double *p_dbl = NULL;
				uint8_t *p_8bit = NULL;
				if(use_8bit) {
					p_8bit = (uint8_t *)block_buf;
				} else {
					p_dbl = (double *)block_buf;
				}
				for(j=0; j<bsize_y; j++) {
					int y = j + boff_y;
					int is_dbuf_stride_y = dbuf && (bandlist_idx==0) && ((y % dbuf->stride_y) == 0);
					uint8_t *mp = mask + (w+2)*(y+1) + boff_x+1;

					array_check_ndv(ndv_def, bandlist_idx, p_dbl, p_8bit, row_ndv, bsize_x);

					for(i=0; i<bsize_x; i++) {
						int is_dbuf_stride = is_dbuf_stride_y && ((i % dbuf->stride_x) == 0);
						if(is_dbuf_stride) {
							int debug_color = use_8bit ? (int)(*p_8bit) : (int)(*p_dbl);
							int x = i + boff_x;
							int db_v = 100;
							db_v += debug_color / 2;
							if(db_v < 100) db_v = 100;
							if(db_v > 254) db_v = 254;
							uint8_t r = (uint8_t)(db_v*.75);
							plot_point(dbuf, x, y, r, db_v, db_v);
						}

						if(use_8bit) p_8bit++;
						else         p_dbl++;
					}

					if(!bandlist_idx) {
						for(i=0; i<bsize_x; i++) {
							*(mp++) = row_ndv[i] ? 0 : 1;
						}
					} else if(ndv_def->invert) {
						for(i=0; i<bsize_x; i++) {
							if(row_ndv[i]) *mp = 0;
							mp++;
						}
					} else {
						for(i=0; i<bsize_x; i++) {
							if(!row_ndv[i]) *mp = 1;
							mp++;
						}
					}
				}
			}
		}

		free(block_buf);
		free(row_ndv);
	}

	if(dbuf) {
		for(int y=0; y<h; y+=dbuf->stride_y) {
		for(int x=0; x<w; x+=dbuf->stride_x) {
			uint8_t is_valid = mask[(w+2)*(y+1) + x+1];
			if(!is_valid) {
				plot_point(dbuf, x, y, 0, 0, 0);
			}
		}
		}
	}

	GDALTermProgress(1, NULL, NULL);

	return mask;
}

uint8_t *read_dataset_8bit(GDALDatasetH ds, int band_idx, uint8_t *usage_array, report_image_t *dbuf) {
	for(int i=0; i<256; i++) usage_array[i] = 0;

	int w = GDALGetRasterXSize(ds);
	int h = GDALGetRasterYSize(ds);
	int band_count = GDALGetRasterCount(ds);
	if(VERBOSE) printf("input is %d x %d x %d\n", w, h, band_count);

	if(band_idx < 1 || band_idx > band_count) fatal_error("bandid out of range");

	GDALRasterBandH band = GDALGetRasterBand(ds, band_idx);

	int blocksize_x, blocksize_y;
	GDALGetBlockSize(band, &blocksize_x, &blocksize_y);

	GDALDataType gdt = GDALGetRasterDataType(band);
	if(gdt != GDT_Byte) {
		printf("Warning: input is not of type Byte, there may be loss while downsampling!\n");
	}

	if(VERBOSE) printf("band %d: block size = %d,%d\n",
		band_idx, blocksize_x, blocksize_y);

	printf("Reading one band of size %d x %d\n", w, h);

	uint8_t *outbuf = (uint8_t *)malloc_or_die(w*h);
	uint8_t *inbuf = (uint8_t *)malloc_or_die(blocksize_x*blocksize_y);
	for(int boff_y=0; boff_y<h; boff_y+=blocksize_y) {
		int bsize_y = blocksize_y;
		if(bsize_y + boff_y > h) bsize_y = h - boff_y;
		for(int boff_x=0; boff_x<w; boff_x+=blocksize_x) {
			int bsize_x = blocksize_x;
			if(bsize_x + boff_x > w) bsize_x = w - boff_x;

			double progress = 
				((double)boff_y * (double)w +
				(double)boff_x * (double)bsize_y) /
				((double)w * (double)h);
			GDALTermProgress(progress, NULL, NULL);

			GDALRasterIO(band, GF_Read, boff_x, boff_y, bsize_x, bsize_y, 
				inbuf, bsize_x, bsize_y, GDT_Byte, 0, 0);

			uint8_t *p_in = inbuf;
			for(int j=0; j<bsize_y; j++) {
				int y = j + boff_y;
				int is_dbuf_stride_y = dbuf && ((y % dbuf->stride_y) == 0);
				uint8_t *p_out = outbuf + w*y + boff_x;
				for(int i=0; i<bsize_x; i++) {
					uint8_t val = *(p_in++);
					*(p_out++) = val;
					usage_array[val] = 1;

					int is_dbuf_stride = is_dbuf_stride_y && ((i % dbuf->stride_x) == 0);
					if(is_dbuf_stride) {
						int x = i + boff_x;
						uint8_t db_v = 100 + (val/2);
						if(db_v < 100) db_v = 100;
						if(db_v > 254) db_v = 254;
						uint8_t r = (uint8_t)(db_v*.75);
						plot_point(dbuf, x, y, r, db_v, db_v);
					}
				}
			}
		}
	}

	free(inbuf);

	GDALTermProgress(1, NULL, NULL);

	return outbuf;
}

uint8_t *get_mask_for_8bit_raster(int w, int h, const uint8_t *raster, uint8_t wanted) {
	if(VERBOSE) printf("mask array is %.1f megabytes\n", (double)(w+2)*(h+2)/1024.0/1024.0);
	uint8_t *mask = (uint8_t *)malloc_or_die((w+2)*(h+2));
	memset(mask, 0, (w+2)*(h+2));

	for(int y=0; y<h; y++) {
		uint8_t *mp = mask + (w+2)*(y+1) + 1;
		const uint8_t *p_in = raster + y*w;
		for(int x=0; x<w; x++) {
			uint8_t val = *(p_in++);
			if(val == wanted) {
				*mp = 1;
			}
			mp++;
		}
	}

	return mask;
}

void erode_mask(uint8_t *mask, int w, int h) {
	w += 2;
	h += 2;
	uint8_t *rowu = (uint8_t *)malloc_or_die(w);
	uint8_t *rowm = (uint8_t *)malloc_or_die(w);
	uint8_t *rowl = (uint8_t *)malloc_or_die(w);
	memset(rowm, 0, w);
	memcpy(rowl, mask, w);

	uint8_t *mp = mask;

	for(int y=0; y<h; y++) {
		uint8_t *tmp = rowu;
		rowu = rowm; rowm = rowl; rowl = tmp;
		if(y+1 < h) {
			memcpy(rowl, mask+w*(y+1), w);
		} else {
			memset(rowl, 0, w);
		}

		uint8_t ul = 0, um = rowu[0];
		uint8_t ml = 0, mm = rowm[0];
		uint8_t ll = 0, lm = rowl[0];

		for(int x=0; x<w; x++) {
			uint8_t ur = (x+1<w) ? rowu[x+1] : 0;
			uint8_t mr = (x+1<w) ? rowm[x+1] : 0;
			uint8_t lr = (x+1<w) ? rowl[x+1] : 0;

			// remove pixels that don't have two consecutive filled neighbors
			if(!(
				(ul&&um) || (um&&ur) || (ur&&mr) || (mr&&lr) ||
				(lr&&lm) || (lm&&ll) || (ll&&ml) || (ml&&ul)
			)) *mp = 0;

			mp++;
			
			ul=um; ml=mm; ll=lm;
			um=ur; mm=mr; lm=lr;
		}
	}

	free(rowu);
	free(rowm);
	free(rowl);
}

void invert_mask(uint8_t *mask, int w, int h) {
	int len = (w+2)*(h+2);
	for(int i=0; i<len; i++) {
		mask[i] = mask[i] ? 0 : 1;
	}
}

vertex_t calc_centroid_from_mask(const uint8_t *mask, int w, int h) {
	long weight_x=0, weight_y=0, num_datavals=0;
	for(int j=0; j<h; j++) {
		const uint8_t *mp = mask + (w+2)*(j+1) + 1;
		for(int i=0; i<w; i++) {
			if(*(mp++)) {
				weight_x += i;
				weight_y += j;
				num_datavals++;
			}
		}
	}

	return (vertex_t){
		(double)weight_x / (double)num_datavals,
		(double)weight_y / (double)num_datavals
	};
}
