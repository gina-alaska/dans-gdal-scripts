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



#include <cassert>

#include "common.h"
#include "ndv.h"

using namespace dangdal;

std::vector<std::vector<size_t> > compute_histogram(
	std::vector<GDALRasterBandH> src_bands, const NdvDef &ndv_def, 
	size_t w, size_t h, int input_range
);
void print_band_stats(size_t band_idx, const std::vector<size_t> &histogram, size_t npix);
void get_scale_from_stddev(const std::vector<size_t> &histogram,
	double dst_avg, double dst_stddev, double *scale_out, double *offset_out);
void get_scale_from_percentile(const std::vector<size_t> &histogram, int output_range,
	double from_percentile, double to_percentile, double *scale_out, double *offset_out);
std::vector<size_t> invert_histogram_to_gaussian(const std::vector<size_t> &histogram_in, double variance, 
	int output_range, double max_rel_freq);
void copyGeoCode(GDALDatasetH dst_ds, GDALDatasetH src_ds);

void usage(const std::string &cmdname) {
	// FIXME - layout/grammar
	printf("Usage: %s <options> src.tif dst.tif\n\n", cmdname.c_str());
	NdvDef::printUsage();
	printf("\
  -outndv <output_nodata_val>             Output no-data value\n\
\n\
Stretch mode:\n\
  -linear-stretch <target_avg> <target_stddev>      Linear stretch to a target range\n\
  -percentile-range <from: 0.0-1.0> <to: 0.0-1.0>   Linear stretch using a percentile range of input\n\
  -histeq <target_stddev>                           Histogram normalize to a target bell curve\n\
\n\
Input must be either 8-bit or 16-bit (unsigned).  Output is 8-bit.\n\
");
	exit(1);
}

int main(int argc, char *argv[]) {
	const std::string cmdname = argv[0];
	if(argc == 1) usage(cmdname);
	std::vector<std::string> arg_list = argv_to_list(argc, argv);

	std::string src_fn;
	std::string dst_fn;
	std::string output_format;

	char mode_histeq = 0;
	char mode_stddev = 0;
	double dst_avg = -1;
	double dst_stddev = -1;
	char mode_percentile = 0;
	double from_percentile = -1;
	double to_percentile = -1;
	int out_ndv = 0, set_out_ndv = 0;

	NdvDef ndv_def = NdvDef(arg_list);

	size_t argp = 1;
	while(argp < arg_list.size()) {
		const std::string &arg = arg_list[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(arg == "-of") {
				if(argp == arg_list.size()) usage(cmdname);
				output_format = arg_list[argp++];
			} else if(arg == "-linear-stretch") {
				if(argp == arg_list.size()) usage(cmdname);
				char *endptr;
				dst_avg = strtod(arg_list[argp++].c_str(), &endptr);
				if(*endptr) usage(cmdname);

				if(argp == arg_list.size()) usage(cmdname);
				dst_stddev = strtod(arg_list[argp++].c_str(), &endptr);
				if(*endptr) usage(cmdname);

				mode_stddev = 1;
			} else if(arg == "-percentile-range") {
				if(argp == arg_list.size()) usage(cmdname);
				char *endptr;
				from_percentile = strtod(arg_list[argp++].c_str(), &endptr);
				if(*endptr) usage(cmdname);

				if(argp == arg_list.size()) usage(cmdname);
				to_percentile = strtod(arg_list[argp++].c_str(), &endptr);
				if(*endptr) usage(cmdname);

				mode_percentile = 1;
			} else if(arg == "-histeq") {
				if(argp == arg_list.size()) usage(cmdname);
				char *endptr;
				dst_stddev = strtod(arg_list[argp++].c_str(), &endptr);
				if(*endptr) usage(cmdname);

				mode_histeq = 1;
			} else if(arg == "-outndv") {
 				if(argp == arg_list.size()) usage(cmdname);
				char *endptr;
				long ndv_long = strtol(arg_list[argp++].c_str(), &endptr, 10);
				out_ndv = (uint8_t)ndv_long;
				if(ndv_long != (long)out_ndv) fatal_error("ndv must be in the range 0..255");
				set_out_ndv++;
				if(*endptr) usage(cmdname);
			} else usage(cmdname);
		} else {
			if(src_fn.empty()) {
				src_fn = arg;
			} else if(dst_fn.empty()) {
				dst_fn = arg;
			} else {
				usage(cmdname);
			}
		}
	}

	if(src_fn.empty() || dst_fn.empty()) usage(cmdname);
	if(mode_percentile + mode_stddev + mode_histeq > 1) usage(cmdname);
	if(mode_stddev && (dst_avg < 0 || dst_stddev < 0)) usage(cmdname);
	if(mode_percentile && !(
		0 <= from_percentile && 
		from_percentile < to_percentile &&
		to_percentile <= 1)) usage(cmdname);

	if(output_format.empty()) output_format = "GTiff";

	GDALAllRegister();

	//////// open source ////////

	GDALDatasetH src_ds = GDALOpen(src_fn.c_str(), GA_ReadOnly);
	if(!src_ds) fatal_error("open failed");

	size_t w = GDALGetRasterXSize(src_ds);
	size_t h = GDALGetRasterYSize(src_ds);
	if(!w || !h) fatal_error("missing width/height");
	size_t src_band_count = GDALGetRasterCount(src_ds);
	printf("Input size is %zd, %zd, %zd\n", w, h, src_band_count);

	std::vector<size_t> bandlist;
	for(size_t i=0; i<src_band_count; i++) {
		bandlist.push_back(i+1);
	}
	size_t dst_band_count = bandlist.size();

	if(ndv_def.empty()) {
		ndv_def = NdvDef(src_ds, bandlist);
	}

	bool use_ndv = !ndv_def.empty();
	if(use_ndv && !set_out_ndv) {
		if(ndv_def.slabs.size() == 1) {
			const NdvSlab &slab = ndv_def.slabs[0];
			assert(slab.range_by_band.size());
			const NdvInterval &range = slab.range_by_band[0];
			if(range.first == range.second) {
				double v = range.first;
				out_ndv = (uint8_t)v;
				if((double)out_ndv == v) {
					set_out_ndv++;
				} else {
					printf("Cannot use %g as an NDV value because it cannot be cast to an 8-bit number.\n", v);
				}
			}
		}
		if(!set_out_ndv) {
			out_ndv = 0;
			set_out_ndv++;
		}
		printf("Output NDV is defaulting to %d.\n", out_ndv);
	}

	//////// open output ////////

	printf("Output size is %zd x %zd x %zd\n", w, h, dst_band_count);

	GDALDriverH dst_driver = GDALGetDriverByName(output_format.c_str());
	if(!dst_driver) fatal_error("unrecognized output format (%s)", output_format.c_str());
	GDALDatasetH dst_ds = GDALCreate(dst_driver, dst_fn.c_str(), w, h, dst_band_count, GDT_Byte, NULL);
	if(!dst_ds) fatal_error("couldn't create output");
	copyGeoCode(dst_ds, src_ds);

	//////// open bands ////////

	std::vector<GDALRasterBandH> src_bands;
	std::vector<GDALRasterBandH> dst_bands;

	for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
		src_bands.push_back(GDALGetRasterBand(src_ds, bandlist[band_idx]));
		dst_bands.push_back(GDALGetRasterBand(dst_ds, band_idx+1));
	}

	//////// compute lookup table ////////

	int input_range = 65536;
	int output_range = 256;

	int mode_noeq = !(mode_percentile || mode_stddev || mode_histeq);

	std::vector<std::vector<size_t> > histograms;
	if(!mode_noeq) {
		histograms = compute_histogram(src_bands, ndv_def, w, h, input_range);
		printf("\n");
	}

	std::vector<std::vector<size_t> > xform_table(dst_band_count);

	for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
		if(mode_noeq) {
			xform_table[band_idx].resize(input_range);
			for(int i=0; i<input_range; i++) xform_table[band_idx][i] = i;
		} else {
			print_band_stats(band_idx, histograms[band_idx], w*h);
			if(mode_stddev || mode_percentile) {
				double scale, offset;
				if(mode_stddev) {
					get_scale_from_stddev(histograms[band_idx],
						dst_avg, dst_stddev, &scale, &offset);
				} else if(mode_percentile) {
					get_scale_from_percentile(histograms[band_idx], output_range,
						from_percentile, to_percentile, &scale, &offset);
				}
				printf("  scale=%f, offset=%f, src_range=[%f, %f]\n",
					scale, offset, -offset/scale, ((double)(output_range-1)-offset)/scale);
				xform_table[band_idx].resize(input_range);
				for(int i=0; i<input_range; i++) {
					xform_table[band_idx][i] = (int)( (double)i * scale + offset );
				}
			} else if(mode_histeq) {
				xform_table[band_idx] = invert_histogram_to_gaussian(
					histograms[band_idx], dst_stddev, output_range, 5.0);
			} else {
				fatal_error("unrecognized mode");
			}
		}
	}
	histograms.clear(); // free up memory

	// clipping
	for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
		for(int i=0; i<input_range; i++) {
			int v = xform_table[band_idx][i];
			if(v < 0) v = 0;
			if(v > output_range-1) v = output_range-1;
			xform_table[band_idx][i] = v;
		}
	}

	// avoid ndv in output for good pixels
	if(use_ndv) {
		for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
			for(int i=0; i<input_range; i++) {
				int v = xform_table[band_idx][i];
				if(v == out_ndv) {
					if(out_ndv < output_range/2) v++;
					else v--;
				}
				xform_table[band_idx][i] = v;
			}
		}
	}

	//////// do transformation ////////

	printf("\nComputing output...\n");

	int blocksize_x_int, blocksize_y_int;
	GDALGetBlockSize(src_bands[0], &blocksize_x_int, &blocksize_y_int);
	size_t blocksize_x = blocksize_x_int;
	size_t blocksize_y = blocksize_y_int;
	size_t block_len = blocksize_x*blocksize_y;

	std::vector<std::vector<GUInt16> > buf_in(dst_band_count);
	std::vector<std::vector<uint8_t> > buf_out(dst_band_count);
	for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
		buf_in[band_idx].resize(block_len);
		buf_out[band_idx].resize(block_len);
	}
	std::vector<double> buf_in_dbl(block_len);
	std::vector<uint8_t> ndv_mask(block_len);
	std::vector<uint8_t> band_mask(block_len);

	for(size_t boff_y=0; boff_y<h; boff_y+=blocksize_y) {
		size_t bsize_y = blocksize_y;
		if(bsize_y + boff_y > h) bsize_y = h - boff_y;
		for(size_t boff_x=0; boff_x<w; boff_x+=blocksize_x) {
			size_t bsize_x = blocksize_x;
			if(bsize_x + boff_x > w) bsize_x = w - boff_x;

			block_len = bsize_x*bsize_y;

			double progress = 
				((double)boff_y * (double)w +
				(double)boff_x * (double)bsize_y) /
				((double)w * (double)h);
			GDALTermProgress(progress, NULL, NULL);

			for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
				GDALRasterIO(src_bands[band_idx], GF_Read, boff_x, boff_y, bsize_x, bsize_y, 
					&buf_in[band_idx][0], bsize_x, bsize_y, GDT_UInt16, 0, 0);

				for(size_t i=0; i<block_len; i++) {
					buf_in_dbl[i] = buf_in[band_idx][i];
				}
				if(band_idx == 0) {
					ndv_def.arrayCheckNdv(band_idx, &buf_in_dbl[0], &ndv_mask[0], block_len);
				} else {
					ndv_def.arrayCheckNdv(band_idx, &buf_in_dbl[0], &band_mask[0], block_len);
					ndv_def.aggregateMask(&ndv_mask[0], &band_mask[0], block_len);
				}
			}

			for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
				size_t *xform = &xform_table[band_idx][0];
				GUInt16 *p_in = &buf_in[band_idx][0];
				uint8_t *p_out = &buf_out[band_idx][0];
				uint8_t *p_ndv = &ndv_mask[0];

				for(size_t i=0; i<block_len; i++) {
					if(*p_ndv) {
						*p_out = out_ndv;
					} else {
						*p_out = (uint8_t)xform[*p_in];
					}
					p_in++; p_out++; p_ndv++;
				}

				GDALRasterIO(dst_bands[band_idx], GF_Write, boff_x, boff_y, bsize_x, bsize_y, 
					&buf_out[band_idx][0], bsize_x, bsize_y, GDT_Byte, 0, 0);
			} // band
		} // block x
	} // block y

	GDALClose(src_ds);
	GDALClose(dst_ds);

	GDALTermProgress(1, NULL, NULL);

	return 0;
}

std::vector<std::vector<size_t> > compute_histogram(
	std::vector<GDALRasterBandH> src_bands, const NdvDef &ndv_def, 
	size_t w, size_t h, int input_range
) {
	printf("\nComputing histogram...\n");

	size_t band_count = src_bands.size();
	std::vector<std::vector<size_t> > histograms(band_count);

	for(size_t band_idx=0; band_idx<band_count; band_idx++) {
		histograms[band_idx].assign(input_range, 0);
	}

	int blocksize_x_int, blocksize_y_int;
	GDALGetBlockSize(src_bands[0], &blocksize_x_int, &blocksize_y_int);
	size_t blocksize_x = blocksize_x_int;
	size_t blocksize_y = blocksize_y_int;
	size_t block_len = blocksize_x*blocksize_y;

	std::vector<std::vector<GUInt16> > buf_in(band_count);
	for(size_t band_idx=0; band_idx<band_count; band_idx++) {
		buf_in[band_idx].resize(block_len);
	}
	std::vector<double> buf_in_dbl(block_len);
	std::vector<uint8_t> ndv_mask(block_len);
	std::vector<uint8_t> band_mask(block_len);

	for(size_t boff_y=0; boff_y<h; boff_y+=blocksize_y) {
		size_t bsize_y = blocksize_y;
		if(bsize_y + boff_y > h) bsize_y = h - boff_y;
		for(size_t boff_x=0; boff_x<w; boff_x+=blocksize_x) {
			size_t bsize_x = blocksize_x;
			if(bsize_x + boff_x > w) bsize_x = w - boff_x;

			block_len = bsize_x*bsize_y;

			double progress = 
				((double)boff_y * (double)w +
				(double)boff_x * (double)bsize_y) /
				((double)w * (double)h);
			GDALTermProgress(progress, NULL, NULL);

			for(size_t band_idx=0; band_idx<band_count; band_idx++) {
				GDALRasterIO(src_bands[band_idx], GF_Read, boff_x, boff_y, bsize_x, bsize_y, 
					&buf_in[band_idx][0], bsize_x, bsize_y, GDT_UInt16, 0, 0);

				for(size_t i=0; i<block_len; i++) {
					buf_in_dbl[i] = buf_in[band_idx][i];
				}
				if(band_idx == 0) {
					ndv_def.arrayCheckNdv(band_idx, &buf_in_dbl[0], &ndv_mask[0], block_len);
				} else {
					ndv_def.arrayCheckNdv(band_idx, &buf_in_dbl[0], &band_mask[0], block_len);
					ndv_def.aggregateMask(&ndv_mask[0], &band_mask[0], block_len);
				}
			}
			for(size_t band_idx=0; band_idx<band_count; band_idx++) {
				GUInt16 *p = &buf_in[band_idx][0];
				size_t *hg = &histograms[band_idx][0];
				
				for(size_t i=0; i<block_len; i++) {
					if(!ndv_mask[i]) hg[p[i]]++;
				}
			}
		}
	}
	GDALTermProgress(1, NULL, NULL);

	return histograms;
}

void print_band_stats(size_t band_idx, const std::vector<size_t> &histogram, size_t npix) {
	size_t ngood=0;
	size_t min=0, max=0;
	for(size_t i=0; i<histogram.size(); i++) {
		if(histogram[i]) {
			if(!ngood || i < min) min = i;
			if(!ngood || i > max) max = i;
			ngood += histogram[i];
		}
	}
	printf("Band %zd:\n  valid_pixels=%zd, ndv_pixels=%zd, min=%zd, max=%zd\n", 
		band_idx, ngood, npix-ngood, min, max);
}

void get_scale_from_stddev(
	const std::vector<size_t> &histogram,
	double dst_avg, double dst_stddev,
	double *scale_out, double *offset_out
) {
	size_t num_pixels = 0;
	double pixel_total = 0;
	for(size_t i=0; i<histogram.size(); i++) {
		num_pixels += histogram[i];
		pixel_total += i*histogram[i];
	}
	double src_avg = pixel_total / (double)num_pixels;
	double error_total = 0;
	for(size_t i=0; i<histogram.size(); i++) {
		double diff = (double)i - src_avg;
		error_total += diff*diff*(double)histogram[i];
	}
	double src_stddev = sqrt(error_total / (double)num_pixels);

	printf("  avg=%f, std_dev=%f\n", src_avg, src_stddev);

	*scale_out = src_stddev ? dst_stddev / src_stddev : 0;
	*offset_out = dst_avg - src_avg * (*scale_out);
}

void get_scale_from_percentile(
	const std::vector<size_t> &histogram, int output_range,
	double from_percentile, double to_percentile,
	double *scale_out, double *offset_out
) {
	size_t num_pixels = 0;
	for(size_t i=0; i<histogram.size(); i++) num_pixels += histogram[i];

	size_t start_count = (size_t)(num_pixels * from_percentile);
	size_t end_count = (size_t)(num_pixels * to_percentile);

	size_t cnt = 0;
	int from_val = -1;
	int to_val = -1;
	for(size_t i=0; i<histogram.size(); i++) {
		if(cnt <= start_count) from_val = i;
		cnt += histogram[i];
		if(cnt <= end_count) to_val = i; 
		else break;
	}
	if(from_val<0 || to_val<0) fatal_error("impossible: could not find window");
	if(from_val == to_val) { from_val=0; to_val=histogram.size()-1; } // FIXME

	*scale_out = (double)(output_range-1) / (double)(to_val-from_val);
	*offset_out = -(double)from_val * (*scale_out);
}

std::vector<double> gen_gaussian(double variance, int bin_count) {
	std::vector<double> arr(bin_count);
	double total = 0;
	for(int i=0; i<bin_count; i++) {
		if(variance == 0.0) {
			// if variance==0.0, just use a uniform distribution
			arr[i] = 1;
		} else {
			// Gaussian distribution.  Large variance gives
			// a more level curve, small variance gives
			// a curved peaked at center.
			double x = (double)(i-bin_count/2) / variance;
			arr[i] = exp(-x*x);
		}
		total += arr[i];
	}

	for(int i=0; i<bin_count; i++) arr[i] /= total;

	return arr;
}

std::vector<size_t> invert_histogram(
	const std::vector<size_t> &src_h_in,
	const std::vector<double> dst_h,
	size_t output_range, double max_rel_freq
) {
	std::vector<size_t> src_h(src_h_in);
	size_t pixel_count = 0;
	for(size_t i=0; i<src_h.size(); i++) pixel_count += src_h[i];

	if(max_rel_freq) {
		size_t bin_max = size_t(max_rel_freq * pixel_count / src_h.size());
		for(size_t i=0; i<src_h.size(); i++) {
			if(src_h[i] > bin_max) {
				pixel_count -= src_h[i] - bin_max;
				src_h[i] = bin_max;
			}
		}
	}

	std::vector<size_t> out_h(src_h.size());
	double src_total = 0;
	double dst_total = 0;
	size_t j = 0;
	for(size_t i=0; i<src_h.size(); i++) {
		out_h[i] = j;
		src_total += src_h[i];
		while(j<output_range && dst_total < src_total) {
			dst_total += dst_h[j++] * (double)pixel_count;
		}
	}

	return out_h;
}

std::vector<size_t> invert_histogram_to_gaussian(
	const std::vector<size_t> &histogram_in, double variance, 
	int output_range, double max_rel_freq
) {
	std::vector<double> gaussian = gen_gaussian(variance, output_range);
	return invert_histogram(histogram_in, gaussian, output_range, max_rel_freq);
}

void copyGeoCode(GDALDatasetH dst_ds, GDALDatasetH src_ds) {
	double affine[6];
	if(GDALGetGeoTransform(src_ds, affine) == CE_None) {
		GDALSetGeoTransform(dst_ds, affine);
	}
	GDALSetProjection(dst_ds, GDALGetProjectionRef(src_ds));
}
