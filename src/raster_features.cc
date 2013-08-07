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



#include "raster_features.h"

namespace dangdal {

FeatureBitmap *read_raster_features(
	GDALDatasetH ds, std::vector<size_t> band_ids, const NdvDef &ndv_def, DebugPlot *dbuf
) {
	assert(!band_ids.empty());

	size_t w = GDALGetRasterXSize(ds);
	size_t h = GDALGetRasterYSize(ds);
	size_t band_count = GDALGetRasterCount(ds);
	if(VERBOSE) printf("input is %zd x %zd x %zd\n", w, h, band_count);

	std::vector<GDALRasterBandH> bands;
	for(const size_t band_id : band_ids) {
		if(VERBOSE) printf("opening band %zd\n", band_id);
		GDALRasterBandH band = GDALGetRasterBand(ds, band_id);
		if(!band) fatal_error("Could not open band %zd.", band_id);
		bands.push_back(band);
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

	std::vector<GDALDataType> datatypes;
	std::vector<size_t> dt_sizes;
	size_t dt_total_size = 0;
	if(VERBOSE >= 2) printf("datatype sizes:");
	for(const GDALRasterBandH band : bands) {
		GDALDataType dt = GDALGetRasterDataType(band);
		datatypes.push_back(dt);
		size_t s = GDALGetDataTypeSize(dt) / 8;
		dt_sizes.push_back(s);
		dt_total_size += s;
		if(VERBOSE >= 2) printf(" %zd", s);
	}
	if(VERBOSE >= 2) printf("\n");

	std::vector<std::vector<uint8_t>> band_buf;
	for(size_t i=0; i<bands.size(); i++) {
		size_t num_bytes = blocksize_xy * dt_sizes[i];
		band_buf.push_back(std::vector<uint8_t>(num_bytes));
	}

	std::vector<uint8_t> ndv_mask(blocksize_xy);
	std::vector<uint8_t> band_mask(blocksize_xy);
	// we need to convert to a known datatype in order to interpret NDV values
	std::vector<double> band_buf_dbl(blocksize_xy);

	FeatureBitmap *fbm = new FeatureBitmap(w, h, dt_total_size);

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

				if(!ndv_def.empty()) {
					GDALCopyWords(
						&band_buf[band_idx][0], datatypes[band_idx], dt_sizes[band_idx],
						&band_buf_dbl[0], GDT_Float64, sizeof(double),
						blocksize_xy);

					if(band_idx == 0) {
						ndv_def.arrayCheckNdv(band_idx, &band_buf_dbl[0], &ndv_mask[0], blocksize_xy);
					} else {
						ndv_def.arrayCheckNdv(band_idx, &band_buf_dbl[0], &band_mask[0], blocksize_xy);
						ndv_def.aggregateMask(&ndv_mask[0], &band_mask[0], blocksize_xy);
					}
				}
			}

			FeatureRawVal pixel;
			pixel.resize(dt_total_size);

			std::vector<uint8_t *> buf_p(bands.size());
			for(size_t band_id=0; band_id<bands.size(); band_id++) {
				buf_p[band_id] = &band_buf[band_id][0];
			}
			uint8_t *ndv_mask_p = &ndv_mask[0];

			for(size_t sub_y=0; sub_y<bsize_y; sub_y++) {
				size_t y = sub_y + boff_y;
				bool is_dbuf_stride_y = dbuf && ((y % dbuf->stride_y) == 0);
				for(size_t sub_x=0; sub_x<bsize_x; sub_x++) {
					size_t x = sub_x + boff_x;
					bool is_dbuf_stride = is_dbuf_stride_y && ((sub_x % dbuf->stride_x) == 0);

					// Read a pixel from the buffers.  This is called even for NDV pixels,
					// because it is needed to increment the pointers.
					{
						size_t p = 0;
						for(size_t band_id=0; band_id<bands.size(); band_id++) {
							for(size_t i=0; i<dt_sizes[band_id]; i++) {
								pixel[p++] = *(buf_p[band_id]++);
							}
						}
						assert(p == dt_total_size);
					}

					if(*(ndv_mask_p++)) {
						num_ndv++;

						if(is_dbuf_stride) {
							dbuf->plotPoint(x, y, 0, 0, 0);
						}
					} else {
						num_valid++;

						FeatureBitmap::Index index_val = fbm->get_index(pixel);
						fbm->raster[y*w + x] = index_val;

						if(is_dbuf_stride) {
							size_t x = sub_x + boff_x;
							// assign a random palette for the debug report
							uint8_t r = index_val * 100 + 100;
							uint8_t g = index_val * 173 + 36;
							uint8_t b = index_val * 47  + 202;
							dbuf->plotPoint(x, y, r, g, b);
						}
					}
				}
			}
		}
	}

	GDALTermProgress(1, NULL, NULL);

	printf("Found %zd valid and %zd NDV pixels.\n", num_valid, num_ndv);

	return fbm;
}

} // namespace dangdal

//using namespace dangdal;
//int main(int argc, char **argv) {
//	if(argc != 2) fatal_error("give filename");
//
//	VERBOSE = 2;
//
//	GDALAllRegister();
//
//	std::string input_raster_fn(argv[1]);
//	GDALDatasetH ds = GDALOpen(input_raster_fn.c_str(), GA_ReadOnly);
//	if(!ds) fatal_error("open failed");
//	size_t w = GDALGetRasterXSize(ds);
//	size_t h = GDALGetRasterYSize(ds);
//
//	std::vector<size_t> band_ids;
//	size_t nbands = GDALGetRasterCount(ds);
//	for(size_t i=0; i<nbands; i++) band_ids.push_back(i+1);
//
//	DebugPlot *dbuf = new DebugPlot(w, h, PLOT_CONTOURS);
//
//	FeatureBitmap *fbm = read_raster_features(ds, band_ids, dbuf);
//	FeatureInterpreter interp(ds, band_ids);
//
//	dbuf->writePlot("zz.ppm");
//
//	for(const auto &f : fbm->get_sorted_feature_table()) {
//		printf("feature %d: %s\n", f.first, interp.pixel_to_string(f.second).c_str());
//	}
//
//	GDALClose(ds);
//
//	return 0;
//}
