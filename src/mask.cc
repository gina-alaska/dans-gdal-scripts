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

#include "common.h"
#include "mask.h"
#include "polygon.h"
#include "debugplot.h"
#include "ndv.h"
#include "datatype_conversion.h"

namespace dangdal {

BitGrid get_bitgrid_for_dataset(
	GDALDatasetH ds, const std::vector<size_t> &band_ids,
	const NdvDef &ndv_def, DebugPlot *dbuf
) {
	assert(!band_ids.empty());

	size_t w = GDALGetRasterXSize(ds);
	size_t h = GDALGetRasterYSize(ds);
	size_t band_count = GDALGetRasterCount(ds);
	if(VERBOSE) printf("input is %zd x %zd x %zd\n", w, h, band_count);

	std::vector<GDALRasterBandH> bands;
	std::vector<GDALDataType> datatypes;
	for(const size_t band_id : band_ids) {
		if(VERBOSE) printf("opening band %zd\n", band_id);
		GDALRasterBandH band = GDALGetRasterBand(ds, band_id);
		if(!band) fatal_error("Could not open band %zd.", band_id);
		bands.push_back(band);
		GDALDataType dt = GDALGetRasterDataType(band);
		datatypes.push_back(dt);
	}

	int blocksize_x_int, blocksize_y_int;
	GDALGetBlockSize(bands[0], &blocksize_x_int, &blocksize_y_int);
	// Out of laziness, I am hoping that images always have the same block size for each band.
	for(const GDALRasterBandH band : bands) {
		int bx, by;
		GDALGetBlockSize(band, &bx, &by);
		if(bx != blocksize_x_int || by != blocksize_y_int) {
			fatal_error(
				"Bands have different block sizes.  Not currently implemented.  Please contact developer");
		}
	}
	size_t blocksize_x = blocksize_x_int;
	size_t blocksize_y = blocksize_y_int;
	size_t blocksize_xy = blocksize_x * blocksize_y;

	std::vector<std::vector<uint8_t>> band_buf(bands.size());
	for(size_t i=0; i<bands.size(); i++) {
		size_t dt_size = GDALGetDataTypeSize(datatypes[i]) / 8;
		size_t num_bytes = blocksize_xy * dt_size;
		band_buf[i].resize(num_bytes);
	}

	std::vector<uint8_t> block_mask(blocksize_xy);

	BitGrid mask(w, h);
	mask.zero();

	size_t num_valid = 0;
	size_t num_ndv = 0;

	printf("Reading input...\n");

	size_t num_blocks_x = (w + blocksize_x - 1) / blocksize_x;
	size_t num_blocks_y = (h + blocksize_y - 1) / blocksize_y;
	for(size_t block_y=0; block_y<num_blocks_y; block_y++) {
		size_t boff_y = blocksize_y * block_y;
		size_t bsize_y = blocksize_y;
		if(bsize_y + boff_y > h) bsize_y = h - boff_y;
		for(size_t block_x=0; block_x<num_blocks_x; block_x++) {
			size_t boff_x = blocksize_x * block_x;
			size_t bsize_x = blocksize_x;
			if(bsize_x + boff_x > w) bsize_x = w - boff_x;

			double progress =
				double(
					boff_y * w +
					boff_x * bsize_y
				) / (w * h);
			GDALTermProgress(progress, NULL, NULL);

			for(size_t band_idx=0; band_idx<bands.size(); band_idx++) {
				GDALReadBlock(bands[band_idx], block_x, block_y, &band_buf[band_idx][0]);
			}
			ndv_def.getNdvMask(band_buf, datatypes, &block_mask[0], blocksize_xy);

			for(size_t sub_y=0; sub_y<bsize_y; sub_y++) {
				size_t y = sub_y + boff_y;
				bool is_dbuf_stride_y = dbuf && ((y % dbuf->stride_y) == 0);
				for(size_t sub_x=0; sub_x<bsize_x; sub_x++) {
					size_t x = sub_x + boff_x;
					bool is_dbuf_stride = is_dbuf_stride_y && ((sub_x % dbuf->stride_x) == 0);

					bool is_ndv = block_mask[sub_y*blocksize_x + sub_x];
					mask.set(x, y, !is_ndv);
					if(is_ndv) {
						num_ndv++;
					} else {
						num_valid++;
					}

					if(is_dbuf_stride) {
						uint8_t val[3] = { 0, 0, 0 };
						if(!is_ndv) {
							for(size_t rgb_idx=0; rgb_idx<3; rgb_idx++) {
								size_t band_idx = std::min(rgb_idx, bands.size()-1);
								double dbl_val = gdal_scalar_to_double(
									&band_buf[band_idx][sub_y*blocksize_x + sub_x], datatypes[band_idx]);
								// valid pixels have texture of the image, but with a cyanish hue
								if(rgb_idx==0) {
									val[rgb_idx] = std::max(0.0, std::min(127.0, dbl_val*0.5));
								} else {
									val[rgb_idx] = std::max(64.0, std::min(191.0, dbl_val*0.5+64));
								}
							}
						}
						dbuf->plotPoint(x, y, val[0], val[1], val[2]);

						// Old color scheme:
						//int val = gdal_scalar_to_int32(
						//	&band_buf[0][sub_y*blocksize_x + sub_x], datatypes[0]);
						//int db_v = 50 + val/3;
						//if(db_v < 50) db_v = 50;
						//if(db_v > 254) db_v = 254;
						//uint8_t r = (uint8_t)(db_v*.75);
						//dbuf->plotPoint(x, y, r, (uint8_t)db_v, (uint8_t)db_v);
					}
				}
			}
		}
	}

	GDALTermProgress(1, NULL, NULL);

	printf("Found %zd valid and %zd NDV pixels.\n", num_valid, num_ndv);

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
	int64_t accum_x=0, accum_y=0, cnt=0;
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
