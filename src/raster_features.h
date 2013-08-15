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
private:
	static const int gdal_color_table_size = 4;

	struct BandInfo {
		BandInfo();

		// GDALCopyWords doesn't allow const input
		int32_t val_to_int(FeatureRawVal rawval) const;
		// GDALCopyWords doesn't allow const input
		double val_to_double(FeatureRawVal rawval) const;

		std::string val_to_str(const FeatureRawVal &rawval) const;

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

public:
	FeatureInterpreter(GDALDatasetH ds, std::vector<size_t> band_ids);
	std::string pixel_to_string(const FeatureRawVal &rawval) const;
	void create_ogr_fields(OGRLayerH ogr_layer) const;
	void set_ogr_fields(OGRLayerH ogr_layer, OGRFeatureH ogr_feat, const FeatureRawVal &rawval) const;

private:
	std::vector<BandInfo> band_info_list;
};

// A bitmap of features.  Features are stored as type Index in a bitmap, and can be mapped to
// FeatureRawVal.
struct FeatureBitmap {
	typedef uint16_t Index;

	FeatureBitmap(const size_t _w, const size_t _h, const size_t _raw_vals_size);

	static FeatureBitmap *from_raster(
		GDALDatasetH ds, std::vector<size_t> band_ids, const NdvDef &ndv_def, DebugPlot *dbuf);

	const std::map<FeatureRawVal, Index> &feature_table() const {
		return table;
	}

	Index get_index(const FeatureRawVal &pixel);
	void dump_feature_table() const;
	BitGrid get_mask_for_feature(Index wanted) const;

private:
	const size_t w, h;
	const size_t raw_vals_size;
	GridArray<Index> raster;
	std::map<FeatureRawVal, Index> table;
};

} // namespace dangdal
