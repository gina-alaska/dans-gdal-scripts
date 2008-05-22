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

unsigned char *get_mask_for_dataset(GDALDatasetH ds, int bandlist_size, int *bandlist, 
int num_ndv, double *ndv_list, double ndv_tolerance, report_image_t *dbuf) {
	int i, j;

	int w = GDALGetRasterXSize(ds);
	int h = GDALGetRasterYSize(ds);
	int band_count = GDALGetRasterCount(ds);
	if(VERBOSE) printf("input is %d x %d x %d\n", w, h, band_count);

	if(VERBOSE) printf("mask array is %.1f megabytes\n", (double)(w+2)*(h+2)/1024.0/1024.0);
	unsigned char *mask = (unsigned char *)malloc_or_die((w+2)*(h+2));
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

		double ndv_lo_dbl = ndv_list[bandlist_idx] - ndv_tolerance;
		double ndv_hi_dbl = ndv_list[bandlist_idx] + ndv_tolerance;
		unsigned char ndv_lo_8bit = (unsigned char)MAX(ceil (ndv_lo_dbl), 0);
		unsigned char ndv_hi_8bit = (unsigned char)MIN(floor(ndv_hi_dbl), 255);

		if(VERBOSE) {
			if(use_8bit) {
				printf("ndv_range = [%d,%d]\n", (int)ndv_lo_8bit, (int)ndv_hi_8bit);
			} else {
				printf("ndv_range = [%.15f,%.15f]\n", ndv_lo_dbl, ndv_hi_dbl);
			}
		}

		void *block_buf;
		if(use_8bit) {
			block_buf = malloc_or_die(blocksize_x*blocksize_y);
		} else {
			block_buf = malloc_or_die(blocksize_x*blocksize_y*sizeof(double));
		}
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
				unsigned char *p_8bit = NULL;
				if(use_8bit) {
					p_8bit = (unsigned char *)block_buf;
				} else {
					p_dbl = (double *)block_buf;
				}
				for(j=0; j<bsize_y; j++) {
					int y = j + boff_y;
					int is_dbuf_stride_y = dbuf && ((y % dbuf->stride_y) == 0);
					unsigned char *mp = mask + (w+2)*(y+1) + boff_x+1;
					for(i=0; i<bsize_x; i++) {
						int debug_color;
						int is_ndv;
						if(use_8bit) {
							unsigned char val = *(p_8bit++);
							is_ndv = (val >= ndv_lo_8bit && val <= ndv_hi_8bit);
							debug_color = (int)val;
						} else {
							double val = *(p_dbl++);
							is_ndv = (val >= ndv_lo_dbl && val <= ndv_hi_dbl);
							debug_color = (int)val;
						}
						if(!is_ndv) {
							*mp = 1;

							int is_dbuf_stride = is_dbuf_stride_y && ((i % dbuf->stride_x) == 0);
							if(is_dbuf_stride) {
								int x = i + boff_x;
								int db_v = 100;
								db_v += debug_color / 2;
								if(db_v < 100) db_v = 100;
								if(db_v > 254) db_v = 254;
								unsigned char r = (unsigned char)(db_v*.75);
								plot_point(dbuf, x, y, r, db_v, db_v);
							}
						}
						mp++;
					}
				}
			}
		}

		free(block_buf);
	}

	GDALTermProgress(1, NULL, NULL);

	return mask;
}

unsigned char *read_dataset_8bit(GDALDatasetH ds, int band_idx, unsigned char *usage_array, report_image_t *dbuf) {
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

	unsigned char *outbuf = (unsigned char *)malloc_or_die(w*h);
	unsigned char *inbuf = (unsigned char *)malloc_or_die(blocksize_x*blocksize_y);
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

			unsigned char *p_in = inbuf;
			for(int j=0; j<bsize_y; j++) {
				int y = j + boff_y;
				int is_dbuf_stride_y = dbuf && ((y % dbuf->stride_y) == 0);
				unsigned char *p_out = outbuf + w*y + boff_x;
				for(int i=0; i<bsize_x; i++) {
					unsigned char val = *(p_in++);
					*(p_out++) = val;
					usage_array[val] = 1;

					int is_dbuf_stride = is_dbuf_stride_y && ((i % dbuf->stride_x) == 0);
					if(is_dbuf_stride) {
						int x = i + boff_x;
						unsigned char db_v = 100 + (val/2);
						if(db_v < 100) db_v = 100;
						if(db_v > 254) db_v = 254;
						unsigned char r = (unsigned char)(db_v*.75);
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

unsigned char *get_mask_for_8bit_raster(int w, int h, const unsigned char *raster, unsigned char wanted) {
	if(VERBOSE) printf("mask array is %.1f megabytes\n", (double)(w+2)*(h+2)/1024.0/1024.0);
	unsigned char *mask = (unsigned char *)malloc_or_die((w+2)*(h+2));
	memset(mask, 0, (w+2)*(h+2));

	for(int y=0; y<h; y++) {
		unsigned char *mp = mask + (w+2)*(y+1) + 1;
		const unsigned char *p_in = raster + y*w;
		for(int x=0; x<w; x++) {
			unsigned char val = *(p_in++);
			if(val == wanted) {
				*mp = 1;
			}
			mp++;
		}
	}

	return mask;
}

void erode_mask(unsigned char *mask, int w, int h) {
	w += 2;
	h += 2;
	unsigned char *rowu = (unsigned char *)malloc_or_die(w);
	unsigned char *rowm = (unsigned char *)malloc_or_die(w);
	unsigned char *rowl = (unsigned char *)malloc_or_die(w);
	memset(rowm, 0, w);
	memcpy(rowl, mask, w);

	unsigned char *mp = mask;

	for(int y=0; y<h; y++) {
		unsigned char *tmp = rowu;
		rowu = rowm; rowm = rowl; rowl = tmp;
		if(y+1 < h) {
			memcpy(rowl, mask+w*(y+1), w);
		} else {
			memset(rowl, 0, w);
		}

		unsigned char ul = 0, um = rowu[0];
		unsigned char ml = 0, mm = rowm[0];
		unsigned char ll = 0, lm = rowl[0];

		for(int x=0; x<w; x++) {
			unsigned char ur = (x+1<w) ? rowu[x+1] : 0;
			unsigned char mr = (x+1<w) ? rowm[x+1] : 0;
			unsigned char lr = (x+1<w) ? rowl[x+1] : 0;

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

void invert_mask(unsigned char *mask, int w, int h) {
	int len = (w+2)*(h+2);
	for(int i=0; i<len; i++) {
		mask[i] = mask[i] ? 0 : 1;
	}
}

vertex_t calc_centroid_from_mask(const unsigned char *mask, int w, int h) {
	long weight_x=0, weight_y=0, num_datavals=0;
	for(int j=0; j<h; j++) {
		const unsigned char *mp = mask + (w+2)*(j+1) + 1;
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
