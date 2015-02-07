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



#include <boost/foreach.hpp>

#include "raster_features.h"

namespace dangdal {

FeatureInterpreter::BandInfo::BandInfo() :
	raw_val_offset(0),
	raw_val_size(0),
	is_int(0),
	is_double(0),
	is_complex(0),
	color_table(NULL),
	ogr_fld(NULL)
{ }

// GDALCopyWords doesn't allow const input
int32_t FeatureInterpreter::BandInfo::val_to_int(FeatureRawVal rawval) const {
	return gdal_scalar_to_int32(&rawval[raw_val_offset], dt);
}

// GDALCopyWords doesn't allow const input
double FeatureInterpreter::BandInfo::val_to_double(FeatureRawVal rawval) const {
	return gdal_scalar_to_double(&rawval[raw_val_offset], dt);
}

std::string FeatureInterpreter::BandInfo::val_to_str(const FeatureRawVal &rawval) const {
	if(color_table) {
		const GDALColorEntry *color = GDALGetColorEntry(color_table,
				val_to_int(rawval));
		return str(boost::format("%d, %d, %d, %d") %
				color->c1 % color->c2 % color->c3 % color->c4);
	} else if(is_int) {
		return str(boost::format("%d") % val_to_int(rawval));
	} else if(is_double) {
		return str(boost::format("%g") % val_to_double(rawval));
	} else {
		std::string ret = "0x";
		for(size_t i=0; i<raw_val_size; i++) {
			uint8_t v = rawval[raw_val_offset + i];
			ret += str(boost::format("%02x") % int(v));
		}
		return ret;
	}
}

FeatureInterpreter::FeatureInterpreter(GDALDatasetH ds, std::vector<size_t> band_ids) {
	size_t offset = 0;
	BOOST_FOREACH(const size_t band_id, band_ids) {
		BandInfo bi;

		GDALRasterBandH band = GDALGetRasterBand(ds, band_id);
		if(!band) fatal_error("Could not open band %zd.", band_id);

		if(band_ids.size() == 1) {
			bi.label = "value";
		} else {
			bi.label = str(boost::format("band%d") % band_id);
		}

		bi.dt = GDALGetRasterDataType(band);
		bi.raw_val_size = GDALGetDataTypeSize(bi.dt) / 8;
		bi.raw_val_offset = offset;
		offset += bi.raw_val_size;

		bi.is_int = bi.is_double = bi.is_complex = false;
		switch(bi.dt) {
			case GDT_Byte:
			case GDT_UInt16:
			case GDT_Int16:
			case GDT_UInt32:
			case GDT_Int32:
				bi.is_int = true; break;
			case GDT_Float32:
			case GDT_Float64:
				bi.is_double = true; break;
			case GDT_CInt16:
			case GDT_CInt32:
			case GDT_CFloat32:
			case GDT_CFloat64:
				bi.is_complex = true; break;
			case GDT_Unknown:
			default:
				break; // no-op
		}

		switch(bi.dt) {
			case GDT_Byte:
			case GDT_UInt16:
			case GDT_Int16:
			case GDT_UInt32:
			case GDT_Int32:
				bi.ogr_ft = OFTInteger; break;
			case GDT_Float32:
			case GDT_Float64:
				bi.ogr_ft = OFTReal; break;
			case GDT_CInt16:
			case GDT_CInt32:
			case GDT_CFloat32:
			case GDT_CFloat64:
				// FIXME - currently complex goes to string, but a pair of numeric fields
				// would be better.
				bi.ogr_ft = OFTString; break;
			case GDT_Unknown:
			default:
				bi.ogr_ft = OFTString; break;
		}
		bi.ogr_fld = OGR_Fld_Create(bi.label.c_str(), bi.ogr_ft);

		if(GDALGetRasterColorInterpretation(band) == GCI_PaletteIndex) {
			bi.color_table = GDALGetRasterColorTable(band);
			for(int i=0; i<gdal_color_table_size; i++) {
				std::string label;
				if(band_ids.size() > 1) {
					label = str(boost::format("band%d.") % band_id);
				}
				label += str(boost::format("c%d") % (i+1));

				bi.ogr_pal_fld[i].first = label;
				bi.ogr_pal_fld[i].second = OGR_Fld_Create(label.c_str(), OFTInteger);
			}
		} else {
			bi.color_table = NULL;
		}

		band_info_list.push_back(bi);
	}
}

std::string FeatureInterpreter::pixel_to_string(const FeatureRawVal &rawval) const {
	std::string ret;
	for(size_t i=0; i<band_info_list.size(); i++) {
		if(i) ret += ", ";
		ret += band_info_list[i].val_to_str(rawval);
	}
	return ret;
}

void FeatureInterpreter::create_ogr_fields(OGRLayerH ogr_layer) const {
	BOOST_FOREACH(const BandInfo &bi, band_info_list) {
		OGR_L_CreateField(ogr_layer, bi.ogr_fld, TRUE);
		if(bi.color_table) {
			typedef std::map<std::string, OGRFieldDefnH>::value_type fld_pair_t;
			BOOST_FOREACH(const fld_pair_t &f, bi.ogr_pal_fld) {
				OGR_L_CreateField(ogr_layer, f.second, TRUE);
			}
		}
	}
}

void FeatureInterpreter::set_ogr_fields(OGRLayerH ogr_layer, OGRFeatureH ogr_feat, const FeatureRawVal &rawval) const {
	BOOST_FOREACH(const BandInfo &bi, band_info_list) {
		int fld_idx = OGR_FD_GetFieldIndex(OGR_L_GetLayerDefn(ogr_layer), bi.label.c_str());
		// Did you remember to call create_ogr_fields?
		if(fld_idx < 0) fatal_error("cannot get OGR field %s", bi.label.c_str());
		switch(bi.ogr_ft) {
			case OFTInteger:
				OGR_F_SetFieldInteger(ogr_feat, fld_idx, bi.val_to_int(rawval));
				break;
			case OFTReal:
				OGR_F_SetFieldDouble(ogr_feat, fld_idx, bi.val_to_double(rawval));
				break;
			case OFTString:
				OGR_F_SetFieldString(ogr_feat, fld_idx, bi.val_to_str(rawval).c_str());
				break;
			default:
				// This should not happen, since ogr_ft is only set by this class.
				fatal_error("not implemented");
		}

		if(bi.color_table) {
			const GDALColorEntry *color = GDALGetColorEntry(bi.color_table, bi.val_to_int(rawval));
			for(int i=0; i<gdal_color_table_size; i++) {
				const std::pair<std::string, OGRFieldDefnH> &f = bi.ogr_pal_fld[i];
				int fld_idx = OGR_FD_GetFieldIndex(OGR_L_GetLayerDefn(ogr_layer), f.first.c_str());
				if(fld_idx < 0) fatal_error("cannot get OGR field %s", f.first.c_str());
				int val =
					i==0 ? color->c1 :
					i==1 ? color->c2 :
					i==2 ? color->c3 : color->c4;
				OGR_F_SetFieldInteger(ogr_feat, fld_idx, val);
			}
		}
	}
}

FeatureBitmap::FeatureBitmap(const size_t _w, const size_t _h, const size_t _raw_vals_size) :
	w(_w), h(_h),
	raw_vals_size(_raw_vals_size),
	raster(w, h)
{ }

FeatureBitmap::Index FeatureBitmap::get_index(const FeatureRawVal &pixel) {
	std::map<FeatureRawVal, Index>::iterator it = table.find(pixel);
	if(it == table.end()) {
		if(table.size() >= std::numeric_limits<FeatureBitmap::Index>::max()) {
			fatal_error("Input had too many feature values (max is %zd)",
				size_t(std::numeric_limits<FeatureBitmap::Index>::max()));
		}
		FeatureBitmap::Index v = table.size();
		table[pixel] = v;
		return v;
	} else {
		return it->second;
	}
}

void FeatureBitmap::dump_feature_table() const {
	typedef std::map<FeatureRawVal, Index>::value_type table_pair_t;
	BOOST_FOREACH(const table_pair_t &f, table) {
		printf("feature %d:", f.second);
		for(size_t i=0; i<f.first.size(); i++) {
			if(i) printf(",");
			printf(" %d", f.first[i]);
		}
		printf("\n");
	}
}

BitGrid FeatureBitmap::get_mask_for_feature(FeatureBitmap::Index wanted) const {
	BitGrid mask(w, h);

	for(size_t y=0; y<h; y++) {
		for(size_t x=0; x<w; x++) {
			mask.set(x, y, raster(x, y) == wanted);
		}
	}

	return mask;
}

FeatureBitmap *FeatureBitmap::from_raster(
	GDALDatasetH ds, std::vector<size_t> band_ids, const NdvDef &ndv_def, DebugPlot *dbuf
) {
	assert(!band_ids.empty());

	size_t w = GDALGetRasterXSize(ds);
	size_t h = GDALGetRasterYSize(ds);
	size_t band_count = GDALGetRasterCount(ds);
	if(VERBOSE) printf("input is %zd x %zd x %zd\n", w, h, band_count);

	std::vector<GDALRasterBandH> bands;
	BOOST_FOREACH(const size_t band_id, band_ids) {
		if(VERBOSE) printf("opening band %zd\n", band_id);
		GDALRasterBandH band = GDALGetRasterBand(ds, band_id);
		if(!band) fatal_error("Could not open band %zd.", band_id);
		bands.push_back(band);
	}

	int blocksize_x_int, blocksize_y_int;
	GDALGetBlockSize(bands[0], &blocksize_x_int, &blocksize_y_int);
	// Out of laziness, I am hoping that images always have the same block size for each band.
	BOOST_FOREACH(const GDALRasterBandH band, bands) {
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
	BOOST_FOREACH(const GDALRasterBandH band, bands) {
		GDALDataType dt = GDALGetRasterDataType(band);
		datatypes.push_back(dt);
		size_t s = GDALGetDataTypeSize(dt) / 8;
		dt_sizes.push_back(s);
		dt_total_size += s;
		if(VERBOSE >= 2) printf(" %zd", s);
	}
	if(VERBOSE >= 2) printf("\n");

	std::vector<std::vector<uint8_t> > band_buf(bands.size());
	for(size_t i=0; i<bands.size(); i++) {
		size_t num_bytes = blocksize_xy * dt_sizes[i];
		band_buf[i].resize(num_bytes);
	}

	std::vector<uint8_t> ndv_mask(blocksize_xy);

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
			}
			if(!ndv_def.empty()) {
				ndv_def.getNdvMask(band_buf, datatypes, &ndv_mask[0], blocksize_xy);
			}

			FeatureRawVal pixel;
			pixel.resize(dt_total_size);

			for(size_t sub_y=0; sub_y<bsize_y; sub_y++) {
				size_t y = sub_y + boff_y;
				bool is_dbuf_stride_y = dbuf && ((y % dbuf->stride_y) == 0);

				for(size_t sub_x=0; sub_x<bsize_x; sub_x++) {
					size_t x = sub_x + boff_x;
					bool is_dbuf_stride = is_dbuf_stride_y && ((sub_x % dbuf->stride_x) == 0);

					size_t in_idx = blocksize_x*sub_y + sub_x;

					if(ndv_mask[in_idx]) {
						num_ndv++;

						if(is_dbuf_stride) {
							dbuf->plotPoint(x, y, 0, 0, 0);
						}
					} else {
						num_valid++;

						size_t j = 0;
						for(size_t band_id=0; band_id<bands.size(); band_id++) {
							uint8_t *p = &band_buf[band_id][in_idx*dt_sizes[band_id]];
							for(size_t i=0; i<dt_sizes[band_id]; i++) {
								pixel[j++] = *(p++);
							}
						}
						assert(j == dt_total_size);

						FeatureBitmap::Index index_val = fbm->get_index(pixel);
						fbm->raster(x, y) = index_val;

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
//void usage(const std::string &cmdname) {}
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
//	NdvDef ndv_def = NdvDef(ds, band_ids);
//
//	FeatureBitmap *fbm = FeatureBitmap::from_raster(ds, band_ids, ndv_def, dbuf);
//	FeatureInterpreter interp(ds, band_ids);
//
//	dbuf->writePlot("zz.ppm");
//
//	typedef std::map<FeatureRawVal, FeatureBitmap::Index>::value_type feature_pair_t;
//	BOOST_FOREACH(const feature_pair_t &f, fbm->feature_table()) {
//		printf("feature %d: %s\n", f.second, interp.pixel_to_string(f.first).c_str());
//	}
//
//	GDALClose(ds);
//
//	return 0;
//}
