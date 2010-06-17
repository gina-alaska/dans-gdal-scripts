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
#include "mask.h"
#include "polygon.h"
#include "debugplot.h"
#include "ndv.h"

namespace dangdal {

uint8_t *read_dataset_8bit(GDALDatasetH ds, int band_idx, uint8_t *usage_array, DebugPlot *dbuf) {
	for(int i=0; i<256; i++) usage_array[i] = 0;

	size_t w = GDALGetRasterXSize(ds);
	size_t h = GDALGetRasterYSize(ds);
	int band_count = GDALGetRasterCount(ds);
	if(VERBOSE) printf("input is %zd x %zd x %d\n", w, h, band_count);

	if(band_idx < 1 || band_idx > band_count) fatal_error("bandid out of range");

	GDALRasterBandH band = GDALGetRasterBand(ds, band_idx);

	int blocksize_x_int, blocksize_y_int;
	GDALGetBlockSize(band, &blocksize_x_int, &blocksize_y_int);
	size_t blocksize_x = blocksize_x_int;
	size_t blocksize_y = blocksize_y_int;

	GDALDataType gdt = GDALGetRasterDataType(band);
	if(gdt != GDT_Byte) {
		printf("Warning: input is not of type Byte, there may be loss while downsampling!\n");
	}

	if(VERBOSE) printf("band %d: block size = %zd,%zd\n",
		band_idx, blocksize_x, blocksize_y);

	printf("Reading one band of size %zd x %zd\n", w, h);

	uint8_t *outbuf = MYALLOC(uint8_t, w*h);
	uint8_t *inbuf = MYALLOC(uint8_t, blocksize_x*blocksize_y);
	for(size_t boff_y=0; boff_y<h; boff_y+=blocksize_y) {
		size_t bsize_y = blocksize_y;
		if(bsize_y + boff_y > h) bsize_y = h - boff_y;
		for(size_t boff_x=0; boff_x<w; boff_x+=blocksize_x) {
			size_t bsize_x = blocksize_x;
			if(bsize_x + boff_x > w) bsize_x = w - boff_x;

			double progress = 
				double(
					boff_y * w +
					boff_x * bsize_y
				) / (w * h);
			GDALTermProgress(progress, NULL, NULL);

			GDALRasterIO(band, GF_Read, boff_x, boff_y, bsize_x, bsize_y, 
				inbuf, bsize_x, bsize_y, GDT_Byte, 0, 0);

			uint8_t *p_in = inbuf;
			for(size_t j=0; j<bsize_y; j++) {
				size_t y = j + boff_y;
				bool is_dbuf_stride_y = dbuf && ((y % dbuf->stride_y) == 0);
				uint8_t *p_out = outbuf + w*y + boff_x;
				for(size_t i=0; i<bsize_x; i++) {
					uint8_t val = *(p_in++);
					*(p_out++) = val;
					usage_array[val] = 1;

					bool is_dbuf_stride = is_dbuf_stride_y && ((i % dbuf->stride_x) == 0);
					if(is_dbuf_stride) {
						size_t x = i + boff_x;
						uint8_t db_v = 50 + val/3;
						if(db_v < 50) db_v = 50;
						if(db_v > 254) db_v = 254;
						uint8_t r = (uint8_t)(db_v*.75);
						dbuf->plotPoint(x, y, r, db_v, db_v);
					}
				}
			}
		}
	}

	free(inbuf);

	GDALTermProgress(1, NULL, NULL);

	return outbuf;
}

BitGrid get_bitgrid_for_dataset(GDALDatasetH ds, int bandlist_size, int *bandlist, 
ndv_def_t *ndv_def, DebugPlot *dbuf) {
	size_t w = GDALGetRasterXSize(ds);
	size_t h = GDALGetRasterYSize(ds);
	int band_count = GDALGetRasterCount(ds);
	if(VERBOSE) printf("input is %zd x %zd x %d\n", w, h, band_count);

	BitGrid mask(w, h);
	mask.zero();

	printf("Reading %d bands of size %zd x %zd\n", bandlist_size, w, h);

	for(int bandlist_idx=0; bandlist_idx<bandlist_size; bandlist_idx++) {
		int band_idx = bandlist[bandlist_idx];
		if(band_idx < 1 || band_idx > band_count) fatal_error("bandid out of range");

		GDALRasterBandH band = GDALGetRasterBand(ds, band_idx);

		int blocksize_x_int, blocksize_y_int;
		GDALGetBlockSize(band, &blocksize_x_int, &blocksize_y_int);
		size_t blocksize_x = blocksize_x_int;
		size_t blocksize_y = blocksize_y_int;

		// gain speed in common 8-bit case
		GDALDataType gdt = GDALGetRasterDataType(band);
		bool use_8bit = (gdt == GDT_Byte);

		if(VERBOSE) printf("band %d: block size = %zd,%zd, use_8bit=%d\n",
			band_idx, blocksize_x, blocksize_y, use_8bit?1:0);

		void *block_buf;
		if(use_8bit) {
			block_buf = MYALLOC(uint8_t, blocksize_x*blocksize_y);
		} else {
			block_buf = MYALLOC(double, blocksize_x*blocksize_y);
		}

		uint8_t *row_ndv = MYALLOC(uint8_t, blocksize_x);

		for(size_t boff_y=0; boff_y<h; boff_y+=blocksize_y) {
			size_t bsize_y = blocksize_y;
			if(bsize_y + boff_y > h) bsize_y = h - boff_y;
			for(size_t boff_x=0; boff_x<w; boff_x+=blocksize_x) {
				size_t bsize_x = blocksize_x;
				if(bsize_x + boff_x > w) bsize_x = w - boff_x;

				double progress = 
					double(
						bandlist_idx * w * h +
						boff_y * w +
						boff_x * bsize_y
					) / (bandlist_size * w * h);
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
				for(size_t j=0; j<bsize_y; j++) {
					size_t y = j + boff_y;
					bool is_dbuf_stride_y = dbuf && (bandlist_idx==0) && ((y % dbuf->stride_y) == 0);

					array_check_ndv(ndv_def, bandlist_idx, p_dbl, p_8bit, row_ndv, bsize_x);

					for(size_t i=0; i<bsize_x; i++) {
						bool is_dbuf_stride = is_dbuf_stride_y && ((i % dbuf->stride_x) == 0);
						if(is_dbuf_stride) {
							int val = use_8bit ? (int)(*p_8bit) : (int)(*p_dbl);
							size_t x = i + boff_x;
							int db_v = 50 + val/3;
							if(db_v < 50) db_v = 50;
							if(db_v > 254) db_v = 254;
							uint8_t r = (uint8_t)(db_v*.75);
							dbuf->plotPoint(x, y, r, (uint8_t)db_v, (uint8_t)db_v);
						}

						if(use_8bit) p_8bit++;
						else         p_dbl++;
					}

					if(!bandlist_idx) {
						for(size_t i=0; i<bsize_x; i++) {
							mask.set(boff_x+i, y, !row_ndv[i]);
						}
					} else if(ndv_def->invert) {
						for(size_t i=0; i<bsize_x; i++) {
							if(row_ndv[i]) mask.set(boff_x+i, y, false);
						}
					} else {
						for(size_t i=0; i<bsize_x; i++) {
							if(!row_ndv[i]) mask.set(boff_x+i, y, true);
						}
					}
				}
			}
		}

		free(block_buf);
		free(row_ndv);
	}

	if(dbuf) {
		for(size_t y=0; y<h; y+=dbuf->stride_y) {
		for(size_t x=0; x<w; x+=dbuf->stride_x) {
			if(!mask(x, y)) {
				dbuf->plotPoint(x, y, 0, 0, 0);
			}
		}
		}
	}

	GDALTermProgress(1, NULL, NULL);

	return mask;
}

BitGrid get_bitgrid_for_8bit_raster(size_t w, size_t h, const uint8_t *raster, uint8_t wanted) {
	BitGrid mask(w, h);

	const uint8_t *p = raster;
	for(size_t y=0; y<h; y++) {
		for(size_t x=0; x<w; x++) {
			mask.set(x, y, *(p++) == wanted);
		}
	}

	return mask;
}

void BitGrid::erode() { // FIXME - untested
	bool *rowu = new bool[w];
	bool *rowm = new bool[w];
	bool *rowl = new bool[w];
	for(int i=0; i<w; i++) {
		rowm[i] = 0;
		rowl[i] = get(i, 0);
	}

	for(int y=0; y<h; y++) {
		bool *tmp = rowu;
		rowu = rowm; rowm = rowl; rowl = tmp;
		for(int i=0; i<w; i++) {
			rowl[i] = get(i, y+1);
		}

		bool ul = 0, um = rowu[0];
		bool ml = 0, mm = rowm[0];
		bool ll = 0, lm = rowl[0];

		for(int x=0; x<w; x++) {
			bool ur = (x+1<w) ? rowu[x+1] : 0;
			bool mr = (x+1<w) ? rowm[x+1] : 0;
			bool lr = (x+1<w) ? rowl[x+1] : 0;

			// remove pixels that don't have two consecutive filled neighbors
			if(!(
				(ul&&um) || (um&&ur) || (ur&&mr) || (mr&&lr) ||
				(lr&&lm) || (lm&&ll) || (ll&&ml) || (ml&&ul)
			)) set(x, y, false);

			ul=um; ml=mm; ll=lm;
			um=ur; mm=mr; lm=lr;
		}
	}

	delete[] rowu;
	delete[] rowm;
	delete[] rowl;
}

Vertex BitGrid::centroid() {
	long accum_x=0, accum_y=0, cnt=0;
	for(int y=0; y<h; y++) {
		for(int x=0; x<w; x++) {
			if(get(x, y)) {
				accum_x += x;
				accum_y += y;
				cnt++;
			}
		}
	}
	return Vertex(
		double(accum_x) / cnt,
		double(accum_y) / cnt
	);
}

} // namespace dangdal
