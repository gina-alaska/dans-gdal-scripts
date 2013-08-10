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



// This module supports the gdal_trace_outline -classify option.  A raster is read and pixel
// values are converted to unique indices (kind of like a palette, except for arbitrary data
// types).  With this indexed array, it is possible to pull out a bitmask for each individual
// pixel value.  Also, the pixel values can be formatted as strings or as fields in an OGR
// file.

#include <vector>
#include <map>
#include <utility>
#include <limits>
#include <string>

#include <boost/format.hpp>

#include <ogr_api.h>

#include "common.h"
#include "mask.h"
#include "debugplot.h"
#include "datatype_conversion.h"

namespace dangdal {

// Represents a raw pixel value.  For convenience, it inherits directly from std::vector so as
// to inherit the comparison operators.
struct FeatureRawVal : public std::vector<uint8_t> { };

// Converts pixel values to strings or OGR fields.
struct FeatureInterpreter {
	static const int gdal_color_table_size = 4;

	struct BandInfo {
		BandInfo() :
			raw_val_offset(0),
			raw_val_size(0),
			is_int(0),
			is_double(0),
			is_complex(0),
			color_table(NULL),
			ogr_fld(NULL)
		{ }

		// GDALCopyWords doesn't allow const input
		int32_t val_to_int(FeatureRawVal rawval) const {
			return gdal_scalar_to_int32(&rawval[raw_val_offset], dt);
		}

		// GDALCopyWords doesn't allow const input
		double val_to_double(FeatureRawVal rawval) const {
			return gdal_scalar_to_double(&rawval[raw_val_offset], dt);
		}

		std::string val_to_str(const FeatureRawVal &rawval) const {
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

		std::string label;
		size_t raw_val_offset;
		size_t raw_val_size;
		GDALDataType dt;
		bool is_int, is_double, is_complex;
		GDALColorTableH color_table;
		OGRFieldType ogr_ft;
		OGRFieldDefnH ogr_fld;
		// For paletted bands, this gives the label and field for each color component.
		std::pair<std::string, OGRFieldDefnH> ogr_pal_fld[gdal_color_table_size];
	};

	FeatureInterpreter(GDALDatasetH ds, std::vector<size_t> band_ids) {
		size_t offset = 0;
		for(const size_t band_id : band_ids) {
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

	std::string pixel_to_string(const FeatureRawVal &rawval) const {
		std::string ret;
		for(size_t i=0; i<band_info_list.size(); i++) {
			if(i) ret += ", ";
			ret += band_info_list[i].val_to_str(rawval);
		}
		return ret;
	}

	void create_ogr_fields(OGRLayerH ogr_layer) const {
		for(const BandInfo &bi : band_info_list) {
			OGR_L_CreateField(ogr_layer, bi.ogr_fld, TRUE);
			if(bi.color_table) {
				for(const auto &f : bi.ogr_pal_fld) {
					OGR_L_CreateField(ogr_layer, f.second, TRUE);
				}
			}
		}
	}

	void set_ogr_fields(OGRLayerH ogr_layer, OGRFeatureH ogr_feat, const FeatureRawVal &rawval) const {
		for(const BandInfo &bi : band_info_list) {
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
					const auto &f = bi.ogr_pal_fld[i];
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

	std::vector<BandInfo> band_info_list;
};

// A bitmap of features.  Features are stored as type Index in a bitmap, and can be mapped to
// FeatureRawVal.
struct FeatureBitmap {
	typedef uint16_t Index;

	FeatureBitmap(const size_t _w, const size_t _h, const size_t _raw_vals_size) :
		w(_w), h(_h),
		raw_vals_size(_raw_vals_size),
		raster(w * h)
	{ }

	Index get_index(const FeatureRawVal &pixel) {
		auto it = table.find(pixel);
		if(it == table.end()) {
			if(table.size() >= std::numeric_limits<Index>::max()) {
				fatal_error("Input had too many feature values (max is %zd)",
					size_t(std::numeric_limits<Index>::max()));
			}
			Index v = table.size();
			table[pixel] = v;
			return v;
		} else {
			return it->second;
		}
	}

	const std::map<FeatureRawVal, Index> &feature_table() const {
		return table;
	}

	void dump_feature_table() const {
		for(const auto &f : table) {
			printf("feature %d:", f.second);
			for(size_t i=0; i<f.first.size(); i++) {
				if(i) printf(",");
				printf(" %d", f.first[i]);
			}
			printf("\n");
		}
	}

	BitGrid get_mask_for_feature(Index wanted) const {
		BitGrid mask(w, h);

		const Index *p = &raster[0];
		for(size_t y=0; y<h; y++) {
			for(size_t x=0; x<w; x++) {
				mask.set(x, y, *(p++) == wanted);
			}
		}

		return mask;
	}

	const size_t w, h;
	const size_t raw_vals_size;
	std::vector<Index> raster;
	std::map<FeatureRawVal, Index> table;
};

FeatureBitmap *read_raster_features(
	GDALDatasetH ds, std::vector<size_t> band_ids, const NdvDef &ndv_def, DebugPlot *dbuf);

} // namespace dangdal
