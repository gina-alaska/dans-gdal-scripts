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

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "common.h"
#include "ndv.h"

using namespace dangdal;

struct Binning {
	int nbins;
	double offset;
	double scale;

	int to_bin(double v) const {
		if(isinf(v) == -1) return 0;
		if(isinf(v) ==  1) return nbins-1;
		double bin_dbl = round((v-offset)/scale);
		if(isnan(bin_dbl)) fatal_error("nan in to_bin");
		if(bin_dbl < 0) return 0;
		if(bin_dbl > nbins-1) return nbins-1;
		int bin_int = int(bin_dbl);
		assert(bin_int >= 0 && bin_int < nbins);
		return bin_int;
	}

	double from_bin(int i) const {
		return double(i) * scale + offset;
	}
};

struct Histogram {
	Binning binning;
	double min, max, mean, stddev;
	size_t data_count;
	size_t ndv_count;
	std::vector<size_t> counts;
};

std::vector<std::pair<double, double> > compute_minmax(
	const std::vector<GDALRasterBandH> src_bands, const NdvDef &ndv_def, 
	size_t w, size_t h
);
std::vector<Histogram> compute_histogram(
	const std::vector<GDALRasterBandH> src_bands, const NdvDef &ndv_def, 
	size_t w, size_t h, const std::vector<Binning> binnings
);
void get_scale_from_percentile(
	const Histogram &histogram, int output_range,
	double from_percentile, double to_percentile,
	double *scale_out, double *offset_out
);
std::vector<uint8_t> invert_histogram_to_gaussian(const Histogram &histogram_in, double variance, 
	int output_range);
void copyGeoCode(GDALDatasetH dst_ds, GDALDatasetH src_ds);

void usage(const std::string &cmdname) {
	printf("Usage: %s <options> src.tif dst.tif\n\n", cmdname.c_str());
	NdvDef::printUsage();
	printf(
"  -outndv <output_nodata_val>        Output no-data value\n"
"\n"
"Operation:\n"
"  -linear-stretch <target_avg> <target_stddev>      Linear stretch to a target range\n"
"  -percentile-range <from: 0.0-1.0> <to: 0.0-1.0>   Linear stretch using a percentile range of input\n"
"  -histeq <target_stddev>                           Histogram normalize to a target bell curve\n"
"  -dump-histogram                                   Just print the histogram to console\n"
"\n"
"Input can be any integer or floating type (but not complex).  Output is 8-bit.\n"
);
	exit(1);
}

int main(int argc, char *argv[]) {
	const std::string cmdname = argv[0];
	if(argc == 1) usage(cmdname);
	std::vector<std::string> arg_list = argv_to_list(argc, argv);

	std::string src_fn;
	std::string dst_fn;
	std::string output_format;

	int mode_histeq = 0;
	int mode_stddev = 0;
	double dst_avg = -1;
	double dst_stddev = -1;
	int mode_percentile = 0;
	int mode_dump_histogram = 0;
	double from_percentile = -1;
	double to_percentile = -1;
	int out_ndv = 0, set_out_ndv = 0;

	NdvDef ndv_def = NdvDef(arg_list);

	size_t argp = 1;
	while(argp < arg_list.size()) {
		const std::string &arg = arg_list[argp++];
		// FIXME - check for duplicate values
		if(arg[0] == '-') {
			try {
				if(arg == "-of") {
					if(argp == arg_list.size()) usage(cmdname);
					output_format = arg_list[argp++];
				} else if(arg == "-linear-stretch") {
					if(argp == arg_list.size()) usage(cmdname);
					dst_avg = boost::lexical_cast<double>(arg_list[argp++]);

					if(argp == arg_list.size()) usage(cmdname);
					dst_stddev = boost::lexical_cast<double>(arg_list[argp++]);

					mode_stddev = 1;
				} else if(arg == "-percentile-range") {
					if(argp == arg_list.size()) usage(cmdname);
					from_percentile = boost::lexical_cast<double>(arg_list[argp++]);

					if(argp == arg_list.size()) usage(cmdname);
					to_percentile = boost::lexical_cast<double>(arg_list[argp++]);

					mode_percentile = 1;
				} else if(arg == "-histeq") {
					if(argp == arg_list.size()) usage(cmdname);
					dst_stddev = boost::lexical_cast<double>(arg_list[argp++]);

					mode_histeq = 1;
				} else if(arg == "-dump-histogram") {
					mode_dump_histogram = 1;
				} else if(arg == "-outndv") {
					if(argp == arg_list.size()) usage(cmdname);
					int64_t ndv_long = boost::lexical_cast<int64_t>(arg_list[argp++]);
					if(ndv_long < 0 || ndv_long > 255) fatal_error("ndv must be in the range 0..255");
					out_ndv = boost::numeric_cast<uint8_t>(ndv_long);
					set_out_ndv++;
				} else {
					usage(cmdname);
				}
			} catch(boost::bad_lexical_cast &e) {
				fatal_error("cannot parse number given on command line");
			} catch(boost::bad_numeric_cast &e) {
				fatal_error("number given on command line out of range: %s", e.what());
			}
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

	if(src_fn.empty()) usage(cmdname);
	if(dst_fn.empty() != (mode_dump_histogram > 0)) usage(cmdname);
	if(mode_percentile + mode_stddev + mode_histeq + mode_dump_histogram > 1) usage(cmdname);
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

	//////// open bands ////////

	std::vector<GDALRasterBandH> src_bands;
	std::vector<GDALRasterBandH> dst_bands;

	for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
		src_bands.push_back(GDALGetRasterBand(src_ds, bandlist[band_idx]));
	}

	//////// find optimal binning ////////

	std::vector<Binning> binnings(dst_band_count);
	{
		// computed on demand
		std::vector<std::pair<double, double> > minmax;

		for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
			Binning &binning = binnings[band_idx];
			GDALDataType dt = GDALGetRasterDataType(src_bands[band_idx]);
			switch(dt) {
				case GDT_Byte:
					binning.nbins = 256;
					binning.offset = 0;
					binning.scale = 1;
					break;
				case GDT_UInt16:
					binning.nbins = 65536;
					binning.offset = 0;
					binning.scale = 1;
					break;
				case GDT_Int16:
					binning.nbins = 65536;
					binning.offset = -32768;
					binning.scale = 1;
					break;
				default:
					if(!minmax.size()) {
						printf("Computing min/max values...\n");
						minmax = compute_minmax(src_bands, ndv_def, w, h);
					}
					double min = minmax[band_idx].first;
					double max = minmax[band_idx].second;
					// a compromise between memory usage and datavalue resolution
					binning.nbins = 10000000;
					binning.offset = min;
					binning.scale = (max - min) / double(binning.nbins-1);
			}
		}
	}

	//////// compute lookup table ////////

	printf("\nComputing histogram...\n");
	std::vector<Histogram> histograms =
		compute_histogram(src_bands, ndv_def, w, h, binnings);

	for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
		Histogram &hg = histograms[band_idx];
		printf("band %zd: min=%g, max=%g, mean=%g, stddev=%g, valid_count=%zd, ndv_count=%zd\n",
			band_idx, hg.min, hg.max, hg.mean, hg.stddev, hg.data_count, hg.ndv_count);
		if(mode_dump_histogram) {
			for(int i=0; i<hg.binning.nbins; i++) {
				printf("bin %d: val=%g cnt=%zd\n",
					i, hg.binning.from_bin(i), hg.counts[i]);
			}
		}
	}
	if(mode_dump_histogram) {
		return 0;
	}

	//////// open output ////////

	GDALDatasetH dst_ds = GDALCreate(dst_driver, dst_fn.c_str(), w, h, dst_band_count, GDT_Byte, NULL);
	if(!dst_ds) fatal_error("couldn't create output");
	copyGeoCode(dst_ds, src_ds);

	for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
		dst_bands.push_back(GDALGetRasterBand(dst_ds, band_idx+1));
	}

	//////// compute tranformation parameters ////////

	const int output_range = 256;
	bool use_table; // otherwise, use linear
	std::vector<std::vector<uint8_t> > xform_table(dst_band_count);
	std::vector<double> lin_scales(dst_band_count);
	std::vector<double> lin_offsets(dst_band_count);

	if(mode_histeq) {
		use_table = true;
		for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
			xform_table[band_idx] = invert_histogram_to_gaussian(
				histograms[band_idx], dst_stddev, output_range);
		}
	} else{
		use_table = false;
		if(mode_percentile) {
			for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
				get_scale_from_percentile(
					histograms[band_idx], output_range, from_percentile, to_percentile,
					&lin_scales[band_idx], &lin_offsets[band_idx]);
			}
		} else if(mode_stddev) {
			for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
				Histogram &hg = histograms[band_idx];
				lin_scales[band_idx] = hg.stddev ? dst_stddev / hg.stddev : 0;
				lin_offsets[band_idx] = hg.mean - dst_avg / lin_scales[band_idx];
			}
		} else { // no transformation
			for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
				lin_scales[band_idx] = 1;
				lin_offsets[band_idx] = 0;
			}
		}
		printf("Linear stretch:\n");
		for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
			double scale = lin_scales[band_idx];
			double offset = lin_offsets[band_idx];
			printf("band %zd: scale=%f, offset=%f, src_range=[%f, %f]\n",
				band_idx, scale, offset, -offset/scale, ((double)(output_range-1)-offset)/scale);
		}
	}

	histograms.clear(); // free up memory

	if(use_table) {
		// avoid ndv in output for good pixels
		if(use_ndv) {
			for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
				for(size_t i=0; i<xform_table[band_idx].size(); i++) {
					uint8_t v = xform_table[band_idx][i];
					if(v == out_ndv) {
						if(out_ndv < output_range/2) v++;
						else v--;
					}
					xform_table[band_idx][i] = v;
				}
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

	std::vector<std::vector<double> > buf_in(dst_band_count);
	std::vector<std::vector<uint8_t> > buf_out(dst_band_count);
	for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
		buf_in[band_idx].resize(block_len);
		buf_out[band_idx].resize(block_len);
	}
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
					&buf_in[band_idx][0], bsize_x, bsize_y, GDT_Float64, 0, 0);

				if(band_idx == 0) {
					ndv_def.arrayCheckNdv(band_idx, &buf_in[band_idx][0], &ndv_mask[0], block_len);
				} else {
					ndv_def.arrayCheckNdv(band_idx, &buf_in[band_idx][0], &band_mask[0], block_len);
					ndv_def.aggregateMask(&ndv_mask[0], &band_mask[0], block_len);
				}
			}

			for(size_t band_idx=0; band_idx<dst_band_count; band_idx++) {
				double *p_in = &buf_in[band_idx][0];
				uint8_t *p_out = &buf_out[band_idx][0];
				uint8_t *p_ndv = &ndv_mask[0];
				if(use_table) {
					uint8_t *xform = &xform_table[band_idx][0];
					Binning binning = binnings[band_idx];

					for(size_t i=0; i<block_len; i++) {
						if(*p_ndv) {
							*p_out = out_ndv;
						} else {
							*p_out = xform[binning.to_bin(*p_in)];
							//printf("%g %d %d\n", *p_in, binning.to_bin(*p_in), *p_out);
						}
						p_in++; p_out++; p_ndv++;
					}
				} else {
					double scale = lin_scales[band_idx];
					double offset = lin_offsets[band_idx];
					for(size_t i=0; i<block_len; i++) {
						if(*p_ndv) {
							*p_out = out_ndv;
						} else {
							double out_dbl = (*p_in - offset) * scale;
							uint8_t v =
								(out_dbl < 0) ? 0 :
								(out_dbl > output_range-1) ? output_range-1 :
								uint8_t(out_dbl);
							if(v == out_ndv) {
								if(out_ndv < output_range/2) v++;
								else v--;
							}
							*p_out = v;
						}
						p_in++; p_out++; p_ndv++;
					}
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

std::vector<std::pair<double, double> > compute_minmax(
	const std::vector<GDALRasterBandH> src_bands, const NdvDef &ndv_def, 
	size_t w, size_t h
) {
	size_t band_count = src_bands.size();
	std::vector<std::pair<double, double> > minmax(band_count);

	int blocksize_x_int, blocksize_y_int;
	GDALGetBlockSize(src_bands[0], &blocksize_x_int, &blocksize_y_int);
	size_t blocksize_x = blocksize_x_int;
	size_t blocksize_y = blocksize_y_int;
	size_t block_len = blocksize_x*blocksize_y;

	std::vector<std::vector<double> > buf_in(band_count);
	for(size_t band_idx=0; band_idx<band_count; band_idx++) {
		buf_in[band_idx].resize(block_len);
	}
	std::vector<uint8_t> ndv_mask(block_len);
	std::vector<uint8_t> band_mask(block_len);

	std::vector<bool> got_data(band_count);

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
					&buf_in[band_idx][0], bsize_x, bsize_y, GDT_Float64, 0, 0);

				if(band_idx == 0) {
					ndv_def.arrayCheckNdv(band_idx, &buf_in[band_idx][0], &ndv_mask[0], block_len);
				} else {
					ndv_def.arrayCheckNdv(band_idx, &buf_in[band_idx][0], &band_mask[0], block_len);
					ndv_def.aggregateMask(&ndv_mask[0], &band_mask[0], block_len);
				}
			}
			for(size_t band_idx=0; band_idx<band_count; band_idx++) {
				for(size_t i=0; i<block_len; i++) {
					if(ndv_mask[i]) continue;
					double v = buf_in[band_idx][i];
					if(isnan(v) || isinf(v)) continue;

					double &min = minmax[band_idx].first;
					double &max = minmax[band_idx].second;
					if(!got_data[band_idx]) {
						min = v;
						max = v;
						got_data[band_idx] = true;
					}
					if(v < min) min = v;
					if(v > max) max = v;
				}
			}
		}
	}
	GDALTermProgress(1, NULL, NULL);

	return minmax;
}

std::vector<Histogram> compute_histogram(
	const std::vector<GDALRasterBandH> src_bands, const NdvDef &ndv_def, 
	size_t w, size_t h, const std::vector<Binning> binnings
) {
	size_t band_count = src_bands.size();
	std::vector<Histogram> histograms(band_count);

	for(size_t band_idx=0; band_idx<band_count; band_idx++) {
		histograms[band_idx].binning = binnings[band_idx];
		histograms[band_idx].counts.assign(binnings[band_idx].nbins, 0);
	}

	int blocksize_x_int, blocksize_y_int;
	GDALGetBlockSize(src_bands[0], &blocksize_x_int, &blocksize_y_int);
	size_t blocksize_x = blocksize_x_int;
	size_t blocksize_y = blocksize_y_int;
	size_t block_len = blocksize_x*blocksize_y;

	std::vector<std::vector<double> > buf_in(band_count);
	for(size_t band_idx=0; band_idx<band_count; band_idx++) {
		buf_in[band_idx].resize(block_len);
	}
	std::vector<uint8_t> ndv_mask(block_len);
	std::vector<uint8_t> band_mask(block_len);

	bool first_valid_pixel = true;

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
					&buf_in[band_idx][0], bsize_x, bsize_y, GDT_Float64, 0, 0);

				if(band_idx == 0) {
					ndv_def.arrayCheckNdv(band_idx, &buf_in[band_idx][0], &ndv_mask[0], block_len);
				} else {
					ndv_def.arrayCheckNdv(band_idx, &buf_in[band_idx][0], &band_mask[0], block_len);
					ndv_def.aggregateMask(&ndv_mask[0], &band_mask[0], block_len);
				}
			}
			for(size_t band_idx=0; band_idx<band_count; band_idx++) {
				Histogram &hg = histograms[band_idx];
				double *p = &buf_in[band_idx][0];
				
				for(size_t i=0; i<block_len; i++) {
					if(ndv_mask[i]) {
						hg.ndv_count++;
					} else {
						double v = p[i];
						hg.counts[hg.binning.to_bin(v)]++;
						if(first_valid_pixel) {
							hg.min = hg.max = v;
							first_valid_pixel = false;
						}
						if(v < hg.min) hg.min = v;
						if(v > hg.max) hg.max = v;
					}
				}
			}
		}
	}
	GDALTermProgress(1, NULL, NULL);

	for(size_t band_idx=0; band_idx<band_count; band_idx++) {
		Histogram &hg = histograms[band_idx];
		double accum = 0;
		for(int i=0; i<hg.binning.nbins; i++) {
			size_t cnt = hg.counts[i];
			double v = hg.binning.from_bin(i);
			hg.data_count += cnt;
			accum += v * cnt;
		}
		hg.mean = accum / hg.data_count;

		double var_accum = 0;
		for(int i=0; i<hg.binning.nbins; i++) {
			size_t cnt = hg.counts[i];
			double v = hg.binning.from_bin(i);
			var_accum += (v-hg.mean) * (v-hg.mean) * cnt;
		}
		hg.stddev = sqrt(var_accum / hg.data_count);
	}

	return histograms;
}

void get_scale_from_percentile(
	const Histogram &histogram, int output_range,
	double from_percentile, double to_percentile,
	double *scale_out, double *offset_out
) {
	size_t start_count = (size_t)(histogram.data_count * from_percentile);
	size_t end_count = (size_t)(histogram.data_count * to_percentile);

	size_t cnt = 0;
	int from_idx = -1;
	int to_idx = -1;
	for(int i=0; i<histogram.binning.nbins; i++) {
		if(cnt <= start_count) from_idx = i;
		cnt += histogram.counts[i];
		if(cnt <= end_count) to_idx = i; 
		else break;
	}
	if(from_idx<0 || to_idx<0) fatal_error("impossible: could not find window");
	if(from_idx == to_idx) { from_idx=0; to_idx=histogram.binning.nbins-1; }

	double from_val = histogram.binning.from_bin(from_idx);
	double to_val = histogram.binning.from_bin(to_idx);

	*scale_out = (double)(output_range-1) / (double)(to_val-from_val);
	*offset_out = from_val;
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

std::vector<uint8_t> invert_histogram(
	const Histogram &src_h_in,
	const std::vector<double> dst_h,
	size_t output_range
) {
	std::vector<size_t> src_h(src_h_in.counts);
	size_t pixel_count = 0;
	for(size_t i=0; i<src_h.size(); i++) pixel_count += src_h[i];

	std::vector<uint8_t> out_h(src_h.size());
	double src_total = 0;
	double dst_total = 0;
	uint8_t j = 0;
	for(size_t i=0; i<src_h.size(); i++) {
		out_h[i] = j;
		src_total += src_h[i];
		while(j<output_range-1 && dst_total < src_total) {
			dst_total += dst_h[j++] * (double)pixel_count;
		}
	}

	return out_h;
}

std::vector<uint8_t> invert_histogram_to_gaussian(
	const Histogram &histogram_in, double variance, 
	int output_range
) {
	std::vector<double> gaussian = gen_gaussian(variance, output_range);
	return invert_histogram(histogram_in, gaussian, output_range);
}

void copyGeoCode(GDALDatasetH dst_ds, GDALDatasetH src_ds) {
	double affine[6];
	if(GDALGetGeoTransform(src_ds, affine) == CE_None) {
		GDALSetGeoTransform(dst_ds, affine);
	}
	GDALSetProjection(dst_ds, GDALGetProjectionRef(src_ds));
}
