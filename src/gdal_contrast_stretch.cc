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
#include "ndv.h"

void compute_histogram(GDALRasterBandH *src_bands, int band_count, ndv_def_t *ndv_def, 
	int w, int h, int **histograms, int input_range);
void print_band_stats(int band_idx, int *histogram, int input_range, int npix);
void get_scale_from_stddev(int *histogram, int input_range,
	double dst_avg, double dst_stddev, double *scale_out, double *offset_out);
void get_scale_from_percentile(int *histogram, int input_range, int output_range,
	double from_percentile, double to_percentile, double *scale_out, double *offset_out);
void invert_histogram_to_gaussian(int *histogram_in, double variance, 
	int *table_out, int input_range, int output_range, double max_rel_freq);
void copyGeoCode(GDALDatasetH dst_ds, GDALDatasetH src_ds);

void usage(const char *cmdname) {
	// FIXME - layout/grammar
	printf("Usage: %s <options> src.tif dst.tif\n\n", cmdname);
	print_ndv_usage();
	printf("\
  -outndv <output_nodata_val>             Output no-data value\n\
\n\
Stretch mode:\n\
  -linear-stretch <target_avg> <target_stddev>      Linear stretch to a target range\n\
  -percentile-range <from: 0.0-1.0> <to: 0.0-1.0>   Linear stretch using a percentile range of input\n\
  -histeq <target_stddev>                           Histogram normalize to a target bell curve\n\
\n\
Input must be either 8-bit or 16-bit.  Output is 8-bit.\n\
");
	exit(1);
}

int main(int argc, char *argv[]) {
	const char *src_fn = NULL;
	const char *dst_fn = NULL;
	const char *output_format = NULL;

	char mode_histeq = 0;
	char mode_stddev = 0;
	double dst_avg = -1;
	double dst_stddev = -1;
	char mode_percentile = 0;
	double from_percentile = -1;
	double to_percentile = -1;
	int out_ndv = 0, set_out_ndv = 0;

	ndv_def_t ndv_def = init_ndv_options(&argc, &argv);

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(!strcmp(arg, "-of")) {
				if(argp == argc) usage(argv[0]);
				output_format = argv[argp++];
			} else if(!strcmp(arg, "-linear-stretch")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				dst_avg = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);

				if(argp == argc) usage(argv[0]);
				dst_stddev = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);

				mode_stddev = 1;
			} else if(!strcmp(arg, "-percentile-range")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				from_percentile = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);

				if(argp == argc) usage(argv[0]);
				to_percentile = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);

				mode_percentile = 1;
			} else if(!strcmp(arg, "-histeq")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				dst_stddev = strtod(argv[argp++], &endptr);
				if(*endptr) usage(argv[0]);

				mode_histeq = 1;
			} else if(!strcmp(arg, "-outndv")) {
 				if(argp == argc) usage(argv[0]);
				char *endptr;
				long ndv_long = strtol(argv[argp++], &endptr, 10);
				out_ndv = (uint8_t)ndv_long;
				if(ndv_long != (long)out_ndv) fatal_error("ndv must be in the range 0..255");
				set_out_ndv++;
				if(*endptr) usage(argv[0]);
			} else usage(argv[0]);
		} else {
			if(!src_fn) {
				src_fn = arg;
			} else if(!dst_fn) {
				dst_fn = arg;
			} else {
				usage(argv[0]);
			}
		}
	}

	if(!src_fn || !dst_fn) usage(argv[0]);
	if(mode_percentile + mode_stddev + mode_histeq > 1) usage(argv[0]);
	if(mode_stddev && (dst_avg < 0 || dst_stddev < 0)) usage(argv[0]);
	if(mode_percentile && !(
		0 <= from_percentile && 
		from_percentile < to_percentile &&
		to_percentile <= 1)) usage(argv[0]);

	if(!output_format) output_format = "GTiff";

	GDALAllRegister();

	//////// open source ////////

	GDALDatasetH src_ds = GDALOpen(src_fn, GA_ReadOnly);
	if(!src_ds) fatal_error("open failed");

	int w = GDALGetRasterXSize(src_ds);
	int h = GDALGetRasterYSize(src_ds);
	if(!w || !h) fatal_error("missing width/height");
	int band_count = GDALGetRasterCount(src_ds);
	printf("Input size is %d, %d, %d\n", w, h, band_count);

	if(!ndv_def.nranges) {
		int *bandlist = MYALLOC(int, band_count);
		for(int i=0; i<band_count; i++) bandlist[i] = i+1;
		add_ndv_from_raster(&ndv_def, src_ds, band_count, bandlist);
	}

	int use_ndv = (ndv_def.nranges > 0);
	if(use_ndv && !set_out_ndv) {
		if(ndv_def.nranges == 1 && ndv_def.ranges[0].min[0] == ndv_def.ranges[0].max[0]) {
			double v = ndv_def.ranges[0].min[0];
			out_ndv = (uint8_t)v;
			if((double)out_ndv == v) {
				set_out_ndv++;
			} else {
				printf("Cannot use %g as an NDV value because it cannot be cast to an 8-bit number.\n", v);
			}
		}
		if(!set_out_ndv) {
			out_ndv = 0;
			set_out_ndv++;
		}
		printf("Output NDV is defaulting to %d.\n", out_ndv);
	}

	//////// open output ////////

	printf("Output size is %d x %d x %d\n", w, h, band_count);

	GDALDriverH dst_driver = GDALGetDriverByName(output_format);
	if(!dst_driver) fatal_error("unrecognized output format (%s)", output_format);
	GDALDatasetH dst_ds = GDALCreate(dst_driver, dst_fn, w, h, band_count, GDT_Byte, NULL);
	if(!dst_ds) fatal_error("could create output");
	copyGeoCode(dst_ds, src_ds);

	//////// open bands ////////

	GDALRasterBandH *src_bands = MYALLOC(GDALRasterBandH, band_count);
	GDALRasterBandH *dst_bands = MYALLOC(GDALRasterBandH, band_count);

	for(int band_idx=0; band_idx<band_count; band_idx++) {
		src_bands[band_idx] = GDALGetRasterBand(src_ds, band_idx+1);
		dst_bands[band_idx] = GDALGetRasterBand(dst_ds, band_idx+1);
	}

	//////// compute lookup table ////////

	int input_range = 65536;
	int output_range = 256;

	int mode_noeq = !(mode_percentile || mode_stddev || mode_histeq);

	int **histograms = NULL;
	if(!mode_noeq) {
		histograms = MYALLOC(int *, band_count);
		for(int band_idx=0; band_idx<band_count; band_idx++) {
			histograms[band_idx] = MYALLOC(int, input_range);
		}
		compute_histogram(src_bands, band_count, &ndv_def, w, h, histograms, input_range);
		printf("\n");
	}

	int **xform_table = MYALLOC(int *, band_count);

	for(int band_idx=0; band_idx<band_count; band_idx++) {
		xform_table[band_idx] = MYALLOC(int, input_range);

		if(mode_noeq) {
			for(int i=0; i<input_range; i++) xform_table[band_idx][i] = i;
		} else {
			print_band_stats(band_idx, histograms[band_idx], input_range, w*h);
			if(mode_stddev || mode_percentile) {
				double scale, offset;
				if(mode_stddev) {
					get_scale_from_stddev(histograms[band_idx], input_range,
						dst_avg, dst_stddev, &scale, &offset);
				} else if(mode_percentile) {
					get_scale_from_percentile(histograms[band_idx], input_range, output_range,
						from_percentile, to_percentile, &scale, &offset);
				}
				printf("  scale=%f, offset=%f, src_range=[%f, %f]\n",
					scale, offset, -offset/scale, ((double)(output_range-1)-offset)/scale);
				for(int i=0; i<input_range; i++) {
					xform_table[band_idx][i] = (int)( (double)i * scale + offset );
				}
			} else if(mode_histeq) {
				invert_histogram_to_gaussian(histograms[band_idx], dst_stddev,
					xform_table[band_idx], input_range, output_range, 5.0);
			} else {
				fatal_error("unrecognized mode");
			}
			free(histograms[band_idx]);
		}
	}
	free(histograms);

	// clipping
	for(int band_idx=0; band_idx<band_count; band_idx++) {
		for(int i=0; i<input_range; i++) {
			int v = xform_table[band_idx][i];
			if(v < 0) v = 0;
			if(v > output_range-1) v = output_range-1;
			xform_table[band_idx][i] = v;
		}
	}

	// avoid ndv in output for good pixels
	if(use_ndv) {
		for(int band_idx=0; band_idx<band_count; band_idx++) {
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

	int blocksize_x, blocksize_y;
	GDALGetBlockSize(src_bands[0], &blocksize_x, &blocksize_y);
	int block_len = blocksize_x*blocksize_y;

	GUInt16 **buf_in = MYALLOC(GUInt16 *, band_count);
	uint8_t **buf_out = MYALLOC(uint8_t *, band_count);
	for(int band_idx=0; band_idx<band_count; band_idx++) {
		buf_in[band_idx] = MYALLOC(GUInt16, block_len);
		buf_out[band_idx] = MYALLOC(uint8_t, block_len);
	}
	double *buf_in_dbl = MYALLOC(double, block_len);
	uint8_t *ndv_mask = MYALLOC(uint8_t, block_len);
	uint8_t *band_mask = MYALLOC(uint8_t, block_len);

	for(int boff_y=0; boff_y<h; boff_y+=blocksize_y) {
		int bsize_y = blocksize_y;
		if(bsize_y + boff_y > h) bsize_y = h - boff_y;
		for(int boff_x=0; boff_x<w; boff_x+=blocksize_x) {
			int bsize_x = blocksize_x;
			if(bsize_x + boff_x > w) bsize_x = w - boff_x;

			block_len = bsize_x*bsize_y;

			double progress = 
				((double)boff_y * (double)w +
				(double)boff_x * (double)bsize_y) /
				((double)w * (double)h);
			GDALTermProgress(progress, NULL, NULL);

			for(int band_idx=0; band_idx<band_count; band_idx++) {
				GDALRasterIO(src_bands[band_idx], GF_Read, boff_x, boff_y, bsize_x, bsize_y, 
					buf_in[band_idx], bsize_x, bsize_y, GDT_UInt16, 0, 0);

				for(int i=0; i<block_len; i++) {
					buf_in_dbl[i] = buf_in[band_idx][i];
				}
				if(band_idx == 0) {
					array_check_ndv(&ndv_def, band_idx, buf_in_dbl, NULL, ndv_mask, block_len);
				} else {
					array_check_ndv(&ndv_def, band_idx, buf_in_dbl, NULL, band_mask, block_len);
					aggregate_ndv_mask(&ndv_def, ndv_mask, band_mask, block_len);
				}
			}

			for(int band_idx=0; band_idx<band_count; band_idx++) {
				int *xform = xform_table[band_idx];
				GUInt16 *p_in = buf_in[band_idx];
				uint8_t *p_out = buf_out[band_idx];
				uint8_t *p_ndv = ndv_mask;

				for(int i=0; i<block_len; i++) {
					if(*p_ndv) {
						*p_out = out_ndv;
					} else {
						*p_out = (uint8_t)xform[*p_in];
					}
					p_in++; p_out++; p_ndv++;
				}

				GDALRasterIO(dst_bands[band_idx], GF_Write, boff_x, boff_y, bsize_x, bsize_y, 
					buf_out[band_idx], bsize_x, bsize_y, GDT_Byte, 0, 0);
			} // band
		} // block x
	} // block y

	GDALClose(src_ds);
	GDALClose(dst_ds);

	GDALTermProgress(1, NULL, NULL);

	return 0;
}

void compute_histogram(
	GDALRasterBandH *src_bands, int band_count, ndv_def_t *ndv_def, 
	int w, int h, int **histograms, int input_range
) {
	printf("\nComputing histogram...\n");

	for(int band_idx=0; band_idx<band_count; band_idx++) {
		for(int i=0; i<input_range; i++) histograms[band_idx][i] = 0;
	}

	int blocksize_x, blocksize_y;
	GDALGetBlockSize(src_bands[0], &blocksize_x, &blocksize_y);
	int block_len = blocksize_x*blocksize_y;

	GUInt16 **buf_in = MYALLOC(GUInt16 *, band_count);
	for(int band_idx=0; band_idx<band_count; band_idx++) {
		buf_in[band_idx] = MYALLOC(GUInt16, block_len);
	}
	double *buf_in_dbl = MYALLOC(double, block_len);
	uint8_t *ndv_mask = MYALLOC(uint8_t, block_len);
	uint8_t *band_mask = MYALLOC(uint8_t, block_len);

	for(int boff_y=0; boff_y<h; boff_y+=blocksize_y) {
		int bsize_y = blocksize_y;
		if(bsize_y + boff_y > h) bsize_y = h - boff_y;
		for(int boff_x=0; boff_x<w; boff_x+=blocksize_x) {
			int bsize_x = blocksize_x;
			if(bsize_x + boff_x > w) bsize_x = w - boff_x;

			block_len = bsize_x*bsize_y;

			double progress = 
				((double)boff_y * (double)w +
				(double)boff_x * (double)bsize_y) /
				((double)w * (double)h);
			GDALTermProgress(progress, NULL, NULL);

			for(int band_idx=0; band_idx<band_count; band_idx++) {
				GDALRasterIO(src_bands[band_idx], GF_Read, boff_x, boff_y, bsize_x, bsize_y, 
					buf_in[band_idx], bsize_x, bsize_y, GDT_UInt16, 0, 0);

				for(int i=0; i<block_len; i++) {
					buf_in_dbl[i] = buf_in[band_idx][i];
				}
				if(band_idx == 0) {
					array_check_ndv(ndv_def, band_idx, buf_in_dbl, NULL, ndv_mask, block_len);
				} else {
					array_check_ndv(ndv_def, band_idx, buf_in_dbl, NULL, band_mask, block_len);
					aggregate_ndv_mask(ndv_def, ndv_mask, band_mask, block_len);
				}
			}
			for(int band_idx=0; band_idx<band_count; band_idx++) {
				GUInt16 *p = buf_in[band_idx];
				int *hg = histograms[band_idx];
				
				for(int i=0; i<block_len; i++) {
					if(!ndv_mask[i]) hg[p[i]]++;
				}
			}
		}
	}
	GDALTermProgress(1, NULL, NULL);

	for(int band_idx=0; band_idx<band_count; band_idx++) {
		free(buf_in[band_idx]);
	}
	free(buf_in);
	free(buf_in_dbl);
	free(ndv_mask);
	free(band_mask);
}

void print_band_stats(int band_idx, int *histogram, int input_range, int npix) {
	int ngood=0;
	int min=0, max=0;
	for(int i=0; i<input_range; i++) {
		if(histogram[i]) {
			if(!ngood || i < min) min = i;
			if(!ngood || i > max) max = i;
			ngood += histogram[i];
		}
	}
	printf("Band %d:\n  valid_pixels=%d, ndv_pixels=%d, min=%d, max=%d\n", 
		band_idx, ngood, npix-ngood, min, max);
}

void get_scale_from_stddev(
	int *histogram, int input_range,
	double dst_avg, double dst_stddev,
	double *scale_out, double *offset_out
) {
	int num_pixels = 0;
	double pixel_total = 0;
	for(int i=0; i<input_range; i++) {
		num_pixels += histogram[i];
		pixel_total += i*histogram[i];
	}
	double src_avg = pixel_total / (double)num_pixels;
	double error_total = 0;
	for(int i=0; i<input_range; i++) {
		double diff = (double)i - src_avg;
		error_total += diff*diff*(double)histogram[i];
	}
	double src_stddev = sqrt(error_total / (double)num_pixels);

	printf("  avg=%f, std_dev=%f\n", src_avg, src_stddev);

	*scale_out = src_stddev ? dst_stddev / src_stddev : 0;
	*offset_out = dst_avg - src_avg * (*scale_out);
}

void get_scale_from_percentile(
	int *histogram, int input_range, int output_range,
	double from_percentile, double to_percentile,
	double *scale_out, double *offset_out
) {
	int num_pixels = 0;
	for(int i=0; i<input_range; i++) num_pixels += histogram[i];

	int start_count = (int)(num_pixels * from_percentile);
	int end_count = (int)(num_pixels * to_percentile);

	int cnt = 0;
	int from_val = -1;
	int to_val = -1;
	for(int i=0; i<input_range; i++) {
		if(cnt <= start_count) from_val = i;
		cnt += histogram[i];
		if(cnt <= end_count) to_val = i; 
		else break;
	}
	if(from_val<0 || to_val<0) fatal_error("impossible: could not find window");
	if(from_val == to_val) { from_val=0; to_val=input_range-1; } // FIXME

	*scale_out = (double)(output_range-1) / (double)(to_val-from_val);
	*offset_out = -(double)from_val * (*scale_out);
}

double *gen_gaussian(double variance, int bin_count) {
	double *arr = MYALLOC(double, bin_count);
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

void invert_histogram(int *src_h, double *dst_h, int *out_h, 
int input_range, int output_range, double max_rel_freq) {
	int pixel_count = 0;
	for(int i=0; i<input_range; i++) pixel_count += src_h[i];

	if(max_rel_freq) {
		int bin_max = (int)(max_rel_freq * (double)pixel_count / (double)input_range);
		for(int i=0; i<input_range; i++) {
			if(src_h[i] > bin_max) {
				pixel_count -= src_h[i] - bin_max;
				src_h[i] = bin_max;
			}
		}
	}

	double src_total = 0;
	double dst_total = 0;
	int j = 0;
	for(int i=0; i<input_range; i++) {
		out_h[i] = j;
		src_total += src_h[i];
		while(j<output_range && dst_total < src_total) {
			dst_total += dst_h[j++] * (double)pixel_count;
		}
	}
}

void invert_histogram_to_gaussian(int *histogram_in, double variance, 
int *table_out, int input_range, int output_range, double max_rel_freq) {
	double *gaussian = gen_gaussian(variance, output_range);
	invert_histogram(histogram_in, gaussian, table_out, input_range, output_range, max_rel_freq);
	free(gaussian);
}

void copyGeoCode(GDALDatasetH dst_ds, GDALDatasetH src_ds) {
	double affine[6];
	if(GDALGetGeoTransform(src_ds, affine) == CE_None) {
		GDALSetGeoTransform(dst_ds, affine);
	}
	GDALSetProjection(dst_ds, GDALGetProjectionRef(src_ds));
}
