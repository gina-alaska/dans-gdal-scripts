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



#include "common.h"

#include <boost/lexical_cast.hpp>

#include <vector>

using namespace dangdal;

struct ScaledBand {
	ScaledBand() :
		oversample(0), lo_w(0), lo_h(0), hi_w(0), hi_h(0),
		delta_x(0), delta_y(0), band(NULL), line_buf_idx(0)
	{ }

	int oversample;
	size_t lo_w, lo_h;
	size_t hi_w, hi_h;
	int delta_x, delta_y;
	std::vector<std::vector<double> > kernel_x;
	std::vector<std::vector<double> > kernel_y;
	GDALRasterBandH band;
	std::vector<std::vector<double> > lines_buf;
	int line_buf_idx;
};

void copyGeoCode(GDALDatasetH dst_ds, GDALDatasetH src_ds);
ScaledBand getScaledBand(GDALDatasetH lores_ds, int band_id, GDALDatasetH hires_ds);
void readLineScaled(ScaledBand &sb, int row, double *hires_buf);
double avoidNDV(double in, double ndv, GDALDataType out_dt);

void usage(const std::string &cmdname) {
	printf("Usage:\n %s\n", cmdname.c_str());
	printf(
"      -rgb <src_rgb.tif> [ -rgb <src.tif> ... ]\n"
"      [ -lum <lum.tif> <weight> ... ] -pan <pan.tif>\n"
"      [ -ndv <nodataval> ] -o <out-rgb.tif>\n"
"\nWhere:\n"
"    rgb.tif    Source bands that are to be enhanced\n"
"    lum.tif    Bands used to simulate lo-res pan band\n"
"    pan.tif    Hi-res panchromatic band\n"
"\nExamples, basic usage:\n"
"    gdal_landsat_pansharp -rgb quickbird_rgb.tif -pan quickbird_pan.tif -o out.tif\n"
"\nExamples, using simulated pan band (gives better results):\n"
"(coefficients from http://www.dpi.inpe.br/~leila/publications/boggione2003spb.pdf)\n"
"    gdal_landsat_pansharp -rgb lansat321.tif -lum landsat234.tif 0.25 0.23 0.52 \\\n"
"      -pan landsat8.tif -ndv 0 -o out.tif\n\n"
"    gdal_landsat_pansharp -rgb landsat3.tif -rgb landsat2.tif -rgb landsat1.tif \\\n"
"      -lum landsat2.tif 0.25 -lum landsat3.tif 0.23 -lum landsat4.tif 0.52 \\\n"
"      -pan landsat8.tif -ndv 0 -o out.tif\n\n"
);
	exit(1);
}

int main(int argc, char *argv[]) {
	const std::string cmdname = argv[0];
	if(argc == 1) usage(cmdname);
	std::vector<std::string> arg_list = argv_to_list(argc, argv);

	std::vector<GDALDatasetH> rgb_ds;
	std::vector<GDALDatasetH> lum_ds;
	std::vector<double> lum_weights;

	std::string pan_fn;
	std::string dst_fn;
	std::string output_format;
	double ndv = 0;
	bool use_ndv = 0;

	GDALAllRegister();

	size_t argp = 1;
	while(argp < arg_list.size()) {
		const std::string &arg = arg_list[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			try {
				if(arg == "-ndv") {
					if(argp == arg_list.size()) usage(cmdname);
					ndv = boost::lexical_cast<double>(arg_list[argp++].c_str());
					use_ndv = true;
				} 
				else if(arg == "-of" ) { if(argp == arg_list.size()) usage(cmdname); output_format = arg_list[argp++]; }
				else if(arg == "-o"  ) { if(argp == arg_list.size()) usage(cmdname); dst_fn = arg_list[argp++]; }
				else if(arg == "-pan") { if(argp == arg_list.size()) usage(cmdname); pan_fn = arg_list[argp++]; }
				else if(arg == "-rgb") {
					if(argp == arg_list.size()) usage(cmdname);
					std::string fn = arg_list[argp++];
					GDALDatasetH ds = GDALOpen(fn.c_str(), GA_ReadOnly);
					if(!ds) fatal_error("open failed");
					rgb_ds.push_back(ds); 
				}
				else if(arg == "-lum") {
					if(argp == arg_list.size()) usage(cmdname);
					std::string fn = arg_list[argp++];
					GDALDatasetH ds = GDALOpen(fn.c_str(), GA_ReadOnly);
					if(!ds) fatal_error("open failed");
					lum_ds.push_back(ds); 
					int nb = GDALGetRasterCount(ds);
					while(nb) {
						if(argp == arg_list.size()) usage(cmdname);
						double w = boost::lexical_cast<double>(arg_list[argp++]);
						lum_weights.push_back(w);
						nb--;
					}
				} else {
					usage(cmdname);
				}
			} catch(boost::bad_lexical_cast &e) {
				fatal_error("cannot parse number given on command line");
			}
		} else {
			usage(cmdname);
		}
	}

	if(pan_fn.empty() || dst_fn.empty()) usage(cmdname);

	if(output_format.empty()) output_format = "GTiff";

	//////// open source ////////

	GDALDatasetH pan_ds = GDALOpen(pan_fn.c_str(), GA_ReadOnly);
	if(!pan_ds) fatal_error("open failed");

	size_t w = GDALGetRasterXSize(pan_ds);
	size_t h = GDALGetRasterYSize(pan_ds);
	if(!w || !h) fatal_error("missing width/height");

	if(GDALGetRasterCount(pan_ds) != 1) fatal_error("Pan input must be only one band");

	GDALRasterBandH pan_band = GDALGetRasterBand(pan_ds, 1);

	// this will be updated via GDALDataTypeUnion as RGB bands are opened
	GDALDataType out_dt = GDALGetRasterDataType(pan_band);

	/////

	std::vector<ScaledBand> rgb_bands;
	for(size_t ds_idx=0; ds_idx<rgb_ds.size(); ds_idx++) {
		int nb = GDALGetRasterCount(rgb_ds[ds_idx]);
		for(int i=0; i<nb; i++) {
			rgb_bands.push_back(getScaledBand(rgb_ds[ds_idx], i+1, pan_ds));
			out_dt = GDALDataTypeUnion(out_dt, GDALGetRasterDataType(
				GDALGetRasterBand(rgb_ds[ds_idx], i+1)));
		}
	}
	size_t rgb_band_count = rgb_bands.size();
	if(!rgb_band_count) usage(cmdname);

	std::vector<ScaledBand> lum_bands;
	if(!lum_ds.empty()) {
		for(size_t ds_idx=0; ds_idx<lum_ds.size(); ds_idx++) {
			int nb = GDALGetRasterCount(lum_ds[ds_idx]);
			for(int i=0; i<nb; i++) {
				lum_bands.push_back(getScaledBand(lum_ds[ds_idx], i+1, pan_ds));
			}
		}
	} else {
		lum_bands = rgb_bands;
		lum_weights.assign(rgb_band_count, 1);
	}
	size_t lum_band_count = lum_bands.size();

	double lum_weight_total = 0;
	for(size_t i=0; i<lum_band_count; i++) {
		lum_weight_total += lum_weights[i];
	}
	//printf("lum weights:");
	for(size_t i=0; i<lum_band_count; i++) {
		lum_weights[i] /= lum_weight_total;
		//printf(" %lf", lum_weights[i]);
	}
	//printf("\n");

	//////// open output ////////

	printf("Output size is %zd x %zd x %zd\n", w, h, rgb_band_count);
	printf("Output datatype is %s\n", GDALGetDataTypeName(out_dt));

	GDALDriverH dst_driver = GDALGetDriverByName(output_format.c_str());
	if(!dst_driver) fatal_error("unrecognized output format (%s)", output_format.c_str());
	GDALDatasetH dst_ds = GDALCreate(dst_driver, dst_fn.c_str(), w, h, rgb_band_count, out_dt, NULL);
	if(!dst_ds) fatal_error("could not create output");
	copyGeoCode(dst_ds, pan_ds);

	std::vector<GDALRasterBandH> dst_bands;
	for(size_t i=0; i<rgb_band_count; i++) {
		dst_bands.push_back(GDALGetRasterBand(dst_ds, i+1));
		if(use_ndv) {
			GDALSetRasterNoDataValue(dst_bands[i], ndv);
		}
	}

	//////// process data ////////

	std::vector<std::vector<double> > lum_buf(lum_band_count);
	for(size_t band_idx=0; band_idx<lum_band_count; band_idx++) {
		lum_buf[band_idx].resize(w);
	}
	std::vector<double> pan_buf(w);
	std::vector<double> rgb_buf(w);
	std::vector<double> out_buf(w);
	std::vector<double> scale_buf(w);

	for(size_t row=0; row<h; row++) {
		GDALTermProgress((double)row/h, NULL, NULL);

		GDALRasterIO(pan_band, GF_Read, 0, row, w, 1, &pan_buf[0], w, 1, GDT_Float64, 0, 0);
		for(size_t band_idx=0; band_idx<lum_band_count; band_idx++) {
			readLineScaled(lum_bands[band_idx], row, &lum_buf[band_idx][0]);
		}

		for(size_t col=0; col<w; col++) {
			bool skip = 0;

			if(use_ndv) {
				if(pan_buf[col] == ndv) skip = 1;
				for(size_t band_idx=0; band_idx<lum_band_count; band_idx++) {
					if(lum_buf[band_idx][col] == ndv) {
						skip = 1;
					}
				}
			}

			if(skip) {
				scale_buf[col] = 1;
			} else {
				double lum_out = (double)pan_buf[col];
				double lum_in = 0;
				for(size_t i=0; i<lum_band_count; i++) {
					lum_in += lum_buf[i][col] * lum_weights[i];
				}

				scale_buf[col] = lum_in>0 ? lum_out/lum_in : 0;
			} // skip
		} // col

		for(size_t band_idx=0; band_idx<rgb_band_count; band_idx++) {
			readLineScaled(rgb_bands[band_idx], row, &rgb_buf[0]);

			for(size_t col=0; col<w; col++) {
				if(use_ndv && rgb_buf[col] == ndv) {
					out_buf[col] = ndv;
				} else {
					double dbl_val = rgb_buf[col] * scale_buf[col];

					if(use_ndv) {
						dbl_val = avoidNDV(dbl_val, ndv, out_dt);
					}

					out_buf[col] = dbl_val;
				}
			}

			GDALRasterIO(dst_bands[band_idx], GF_Write, 0, row, w, 1,
				&out_buf[0], w, 1, GDT_Float64, 0, 0);
		}
	} // row

	for(size_t i=0; i<rgb_ds.size(); i++) {
		GDALClose(rgb_ds[i]);
	}
	for(size_t i=0; i<lum_ds.size(); i++) {
		GDALClose(lum_ds[i]);
	}
	GDALClose(pan_ds);
	GDALClose(dst_ds);

	GDALTermProgress(1, NULL, NULL);

	return 0;
}

void copyGeoCode(GDALDatasetH dst_ds, GDALDatasetH src_ds) {
	double affine[6];
	if(GDALGetGeoTransform(src_ds, affine) == CE_None) {
		GDALSetGeoTransform(dst_ds, affine);
	}
	GDALSetProjection(dst_ds, GDALGetProjectionRef(src_ds));
}

ScaledBand getScaledBand(GDALDatasetH lores_ds, int band_id, GDALDatasetH hires_ds) {
	double lores_affine[6];
	double hires_affine[6];
	if(GDALGetGeoTransform(lores_ds, lores_affine) != CE_None) {
		fatal_error("cannot determine affine");
	}
	if(GDALGetGeoTransform(hires_ds, hires_affine) != CE_None) {
		fatal_error("cannot determine affine");
	}
	// FIXME - we could probably work it out, as long as all the inputs
	// are rotated the same way
	if(
		lores_affine[2] || lores_affine[4] ||
		hires_affine[2] || hires_affine[4]
	) fatal_error("input must not be rotated");
	// again, it could probably be made to work...
	if(
		lores_affine[1] != -lores_affine[5] ||
		hires_affine[1] != -hires_affine[5]
	) fatal_error("input must have square pixels");

	ScaledBand sb;

	double band_res = lores_affine[1];
	double base_res = hires_affine[1];
	double offset_x =  (hires_affine[0] - lores_affine[0]) / band_res;
	double offset_y = -(hires_affine[3] - lores_affine[3]) / band_res;
	double scale = band_res / base_res;
	sb.oversample = (int)round(scale);
	if(fabs(scale - (double)sb.oversample) > 1e+6) 
		fatal_error("scales of input images must differ by an integer factor");

	sb.lo_w = GDALGetRasterXSize(lores_ds);
	sb.lo_h = GDALGetRasterYSize(lores_ds);
	if(!sb.lo_w || !sb.lo_h) fatal_error("missing width/height");

	sb.hi_w = GDALGetRasterXSize(hires_ds);
	sb.hi_h = GDALGetRasterYSize(hires_ds);
	if(!sb.hi_w || !sb.hi_h) fatal_error("missing width/height");

	//printf("%d: %f %f %d %dx%d\n", band_id, 
	//	offset_x, offset_y, sb.oversample, sb.lo_w, sb.lo_h);

	sb.band = GDALGetRasterBand(lores_ds, band_id);

	sb.delta_x = (int)floor(offset_x);
	offset_x -= (double)sb.delta_x;
	sb.delta_y = (int)floor(offset_y);
	offset_y -= (double)sb.delta_y;

	// Cubic convolution resampling.  For more info see:
	// http://www.imgfsr.com/ResamplingCVPR.pdf
	sb.kernel_x.resize(sb.oversample);
	sb.kernel_y.resize(sb.oversample);
	for(int mod=0; mod<sb.oversample; mod++) {
		sb.kernel_x[mod].resize(4);
		sb.kernel_y[mod].resize(4);
		double t = offset_x + (double)mod / (double)sb.oversample;
		sb.kernel_x[mod][0] = -0.5*t*t*t + 1.0*t*t - 0.5*t;
		sb.kernel_x[mod][1] =  1.5*t*t*t - 2.5*t*t + 1;
		sb.kernel_x[mod][2] = -1.5*t*t*t + 2.0*t*t + 0.5*t;
		sb.kernel_x[mod][3] =  0.5*t*t*t - 0.5*t*t;
		t = offset_y + (double)mod / (double)sb.oversample;
		sb.kernel_y[mod][0] = -0.5*t*t*t + 1.0*t*t - 0.5*t;
		sb.kernel_y[mod][1] =  1.5*t*t*t - 2.5*t*t + 1;
		sb.kernel_y[mod][2] = -1.5*t*t*t + 2.0*t*t + 0.5*t;
		sb.kernel_y[mod][3] =  0.5*t*t*t - 0.5*t*t;
	}

	sb.lines_buf.resize(4);
	for(int j=0; j<4; j++) {
		sb.lines_buf[j].resize(sb.hi_w);
	}
	sb.line_buf_idx = -1000000;

	return sb;
}

void readLineScaled1D(ScaledBand &sb, int row, double *hires_buf) {
	std::vector<double> lores_buf(sb.lo_w);

	if(row < 0 || size_t(row) >= sb.lo_h) {
		for(size_t col=0; col<sb.hi_w; col++) {
			hires_buf[col] = 0;
		}
	} else {
		GDALRasterIO(sb.band, GF_Read, 0, row, sb.lo_w, 1, &lores_buf[0], sb.lo_w, 1, GDT_Float64, 0, 0);
		for(int mx=0; mx<sb.oversample; mx++) {
			double *kernel = &sb.kernel_x[mx][0];
			for(int x0=0; ; x0++) {
				size_t col = x0 * sb.oversample + mx;
				if(col >= sb.hi_w) break;
				double accum = 0;
				for(int i=0; i<4; i++) {
					int x = x0 - 1 + i + sb.delta_x;
					double v = (x<0 || size_t(x) >= sb.lo_w) ? 0 : lores_buf[x];
					accum += v * kernel[i];
				}
				hires_buf[col] = accum;
			}
		}
	}
}

void readLineScaled(ScaledBand &sb, int row, double *hires_buf) {
	int y0 = row / sb.oversample;
	int my = row % sb.oversample;
	double *kernel = &sb.kernel_y[my][0];

	int top_y = y0 - 1 + sb.delta_y;
	if(top_y == sb.line_buf_idx) {
		// no action
	} else if(top_y == sb.line_buf_idx+1) {
		std::vector<double> tmp;
		std::swap(tmp, sb.lines_buf[0]);
		for(int j=0; j<3; j++) {
			std::swap(sb.lines_buf[j], sb.lines_buf[j+1]);
		}
		std::swap(sb.lines_buf[3], tmp);
		readLineScaled1D(sb, top_y+3, &sb.lines_buf[3][0]);
		sb.line_buf_idx = top_y;
	} else {
		for(int j=0; j<4; j++) {
			readLineScaled1D(sb, top_y+j, &sb.lines_buf[j][0]);
		}
		sb.line_buf_idx = top_y;
	}

	for(size_t col=0; col<sb.hi_w; col++) {
		double accum = 0;
		for(int j=0; j<4; j++) {
			accum += sb.lines_buf[j][col] * kernel[j];
		}
		hires_buf[col] = accum;
	}
}

// Helper for avoidNDV, used when datatype is an integer type.
template <typename T>
double avoidNDV_int(double in, double ndv) {
	int64_t in_int = int64_t(round(in));
	int64_t ndv_int = int64_t(round(ndv));
	if(in_int != ndv_int) return in_int;

	int64_t valid_low  = std::numeric_limits<T>::min();
	int64_t valid_high = std::numeric_limits<T>::max();

	// Compute center of valid range.
	// NOTE: this addition doesn't overflow.
	int64_t mid = (valid_low+valid_high)/2;
	// Avoid NDV, by moving toward center of valid range.
	// This way, we don't fall off the edge.
	return (ndv < mid) ? in_int+1 : in_int-1;
}

// If in==ndv, then perturb the value to avoid NDV.
double avoidNDV(double in, double ndv, GDALDataType out_dt) {
	// First, clip to valid range.  This is the only way to know whether the pixel will end up
	// being NDV after being written.  This is probably not the fastest way to do it though...
	{
		// size of largest possible datatype (GDT_CFloat64), doubled (just in case).
		char tmp[32];
		GDALCopyWords(&in, GDT_Float64, 0, tmp, out_dt, 0, 1);
		GDALCopyWords(tmp, out_dt, 0, &in, GDT_Float64, 0, 1);
	}

	switch(out_dt) {
		case GDT_Byte:   return avoidNDV_int< uint8_t>(in, ndv);
		case GDT_UInt16: return avoidNDV_int<uint16_t>(in, ndv);
		case GDT_Int16:  return avoidNDV_int< int16_t>(in, ndv);
		case GDT_UInt32: return avoidNDV_int<uint32_t>(in, ndv);
		case GDT_Int32:  return avoidNDV_int< int32_t>(in, ndv);
		default:         return (in==ndv) ? in+1 : in;
	}
}
