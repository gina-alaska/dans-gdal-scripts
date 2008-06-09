/*
Copyright (c) 2007, Regents of the University of Alaska

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
*/



#include <common.h>

typedef struct {
	int oversample;
	int lo_w, lo_h;
	int hi_w, hi_h;
	int delta_x, delta_y;
	double **kernel_x;
	double **kernel_y;
	GDALRasterBandH band;
	double **lines_buf;
	int line_buf_idx;
} scaled_band_t;

void copyGeoCode(GDALDatasetH dst_ds, GDALDatasetH src_ds);
scaled_band_t getScaledBand(GDALDatasetH lores_ds, int band_id, GDALDatasetH hires_ds);
void readLineScaled(scaled_band_t *sb, int row, double *hires_buf);

void usage(const char *cmdname) {
	printf("Usage:\n %s\n", cmdname);
	printf("\
      -rgb <src_rgb.tif> [ -rgb <src.tif> ... ]\n\
      [ -lum <lum.tif> <weight> ... ] -pan <pan.tif>\n\
      [ -ndv <nodataval> ] -o <out-rgb.tif>\n\
Where:\n\
    rgb.tif    Source bands that are to be enhanced\n\
    lum.tif    Bands used to simulate lo-res pan band\n\
    pan.tif    Hi-res panchromatic band\n\
Examples:\n\
    gdal_landsat_pansharp -rgb lansat321.tif -lum landsat234.tif 0.25 0.23 0.52 \\\n\
      -pan landsat8.tif -ndv 0 -o out.tif\n\n\
    gdal_landsat_pansharp -rgb landsat3.tif -rgb landsat2.tif -rgb landsat1.tif \\\n\
      -lum landsat2.tif 0.25 -lum landsat3.tif 0.23 -lum landsat4.tif 0.52 \\\n\
      -pan landsat8.tif -ndv 0 -out.tif\n\n\
    gdal_landsat_pansharp -rgb quickbird_rgb.tif -pan quickbird_pan.tif -o out.tif\n\
");
	exit(1);
}

int main(int argc, char *argv[]) {
	int i;

	int lum_ds_count = 0;
	int lum_band_count = 0;
	GDALDatasetH *lum_ds = NULL;
	double *lum_weights = NULL;

	int rgb_ds_count = 0;
	int rgb_band_count = 0;
	GDALDatasetH *rgb_ds = NULL;

	const char *pan_fn = NULL;
	const char *dst_fn = NULL;
	const char *output_format = NULL;
	double ndv = 0;
	char use_ndv = 0;

	GDALAllRegister();

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(!strcmp(arg, "-ndv")) {
				if(argp == argc) usage(argv[0]);
				char *endptr;
				ndv = strtol(argv[argp++], &endptr, 10);
				use_ndv++;
				if(*endptr) usage(argv[0]);
				if(ndv < 0 || ndv > 255) fatal_error("no_data_val must be in the range 0-255");
			} 
			else if(!strcmp(arg, "-of" )) { if(argp == argc) usage(argv[0]); output_format = argv[argp++]; }
			else if(!strcmp(arg, "-o"  )) { if(argp == argc) usage(argv[0]); dst_fn = argv[argp++]; }
			else if(!strcmp(arg, "-pan")) { if(argp == argc) usage(argv[0]); pan_fn = argv[argp++]; }
			else if(!strcmp(arg, "-rgb")) {
				if(argp == argc) usage(argv[0]);
				char *fn = argv[argp++];
				rgb_ds = (GDALDatasetH *)realloc_or_die(rgb_ds, sizeof(GDALDatasetH) * (rgb_ds_count+1));
				GDALDatasetH ds = GDALOpen(fn, GA_ReadOnly);
				if(!ds) fatal_error("open failed");
				rgb_ds[rgb_ds_count++] = ds; 
				rgb_band_count += GDALGetRasterCount(ds);
			}
			else if(!strcmp(arg, "-lum")) {
				if(argp == argc) usage(argv[0]);
				char *fn = argv[argp++];
				lum_ds = (GDALDatasetH *)realloc_or_die(lum_ds, sizeof(GDALDatasetH) * (lum_ds_count+1));
				GDALDatasetH ds = GDALOpen(fn, GA_ReadOnly);
				if(!ds) fatal_error("open failed");
				lum_ds[lum_ds_count++] = ds; 
				int nb = GDALGetRasterCount(ds);
				lum_weights = (double *)realloc_or_die(lum_weights, sizeof(double) * (lum_band_count+nb));
				while(nb) {
					if(argp == argc) usage(argv[0]);
					char *endptr;
					lum_weights[lum_band_count] = strtod(argv[argp++], &endptr);
					if(*endptr) usage(argv[0]);
					lum_band_count++;
					nb--;
				}
			}
			else usage(argv[0]);
		} else {
			usage(argv[0]);
		}
	}

	if(!rgb_band_count) usage(argv[0]);
	if(!pan_fn || !dst_fn) usage(argv[0]);

	if(!output_format) output_format = "GTiff";

	//////// open source ////////

	GDALDatasetH pan_ds = GDALOpen(pan_fn, GA_ReadOnly);
	if(!pan_ds) fatal_error("open failed");

	int w = GDALGetRasterXSize(pan_ds);
	int h = GDALGetRasterYSize(pan_ds);
	if(!w || !h) fatal_error("missing width/height");

	if(GDALGetRasterCount(pan_ds) != 1) fatal_error("Pan input must be only one band");

	GDALRasterBandH pan_band = GDALGetRasterBand(pan_ds, 1);

	/////

	int ds_idx;
	int band_idx;

	scaled_band_t *rgb_bands = (scaled_band_t *)malloc_or_die(sizeof(scaled_band_t) * rgb_band_count);
	for(ds_idx=0, band_idx=0; ds_idx<rgb_ds_count; ds_idx++) {
		int nb = GDALGetRasterCount(rgb_ds[ds_idx]);
		int i;
		for(i=0; i<nb; i++) {
			rgb_bands[band_idx++] = getScaledBand(rgb_ds[ds_idx], i+1, pan_ds);
		}
	}

	scaled_band_t *lum_bands;
	if(lum_ds_count) {
		lum_bands = (scaled_band_t *)malloc_or_die(sizeof(scaled_band_t) * lum_band_count);
		for(ds_idx=0, band_idx=0; ds_idx<lum_ds_count; ds_idx++) {
			int nb = GDALGetRasterCount(lum_ds[ds_idx]);
			for(i=0; i<nb; i++) {
				lum_bands[band_idx++] = getScaledBand(lum_ds[ds_idx], i+1, pan_ds);
			}
		}
	} else {
		lum_band_count = rgb_band_count;
		lum_bands = rgb_bands;
		lum_weights = (double *)malloc_or_die(sizeof(double) * lum_band_count);
		for(i=0; i<lum_band_count; i++) lum_weights[i] = 1;
	}

	double lum_weight_total = 0;
	for(i=0; i<lum_band_count; i++) {
		lum_weight_total += lum_weights[i];
	}
	//printf("lum weights:");
	for(i=0; i<lum_band_count; i++) {
		lum_weights[i] /= lum_weight_total;
		//printf(" %lf", lum_weights[i]);
	}
	//printf("\n");

	//////// open output ////////

	printf("Output size is %d x %d x %d\n", w, h, rgb_band_count);

	GDALDriverH dst_driver = GDALGetDriverByName(output_format);
	if(!dst_driver) fatal_error("unrecognized output format (%s)", output_format);
	GDALDatasetH dst_ds = GDALCreate(dst_driver, dst_fn, w, h, rgb_band_count, GDT_Byte, NULL);
	if(!dst_ds) fatal_error("could not create output");
	copyGeoCode(dst_ds, pan_ds);

	GDALRasterBandH *dst_bands = (GDALRasterBandH *)malloc_or_die(sizeof(GDALRasterBandH) * rgb_band_count);
	for(i=0; i<rgb_band_count; i++) {
		dst_bands[i] = GDALGetRasterBand(dst_ds, i+1);
	}

	//////// process data ////////

	double **lum_buf = (double **)malloc_or_die(sizeof(double *) * lum_band_count);
	for(band_idx=0; band_idx<lum_band_count; band_idx++) {
		lum_buf[band_idx] = (double *)malloc_or_die(sizeof(double) * w);
	}
	double *pan_buf = (double *)malloc_or_die(sizeof(double) * w);
	double *rgb_buf = (double *)malloc_or_die(sizeof(double) * w);
	uint8_t *out_buf = (uint8_t *)malloc_or_die(sizeof(uint8_t) * w);
	double *scale_buf = (double *)malloc_or_die(sizeof(double) * w);

	int row;
	for(row=0; row<h; row++) {
		int col;

		GDALTermProgress((double)row/(double)h, NULL, NULL);

		GDALRasterIO(pan_band, GF_Read, 0, row, w, 1, pan_buf, w, 1, GDT_Float64, 0, 0);
		for(band_idx=0; band_idx<lum_band_count; band_idx++) {
			readLineScaled(lum_bands+band_idx, row, lum_buf[band_idx]);
		}

		for(col=0; col<w; col++) {
			int skip = 0;

			if(use_ndv) {
				if(pan_buf[col] == ndv) skip = 1;
				for(band_idx=0; band_idx<lum_band_count; band_idx++) {
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
				for(i=0; i<lum_band_count; i++) {
					lum_in += lum_buf[i][col] * lum_weights[i];
				}

				scale_buf[col] = lum_in>0 ? lum_out/lum_in : 0;
			} // skip
		} // col

		for(band_idx=0; band_idx<rgb_band_count; band_idx++) {
			readLineScaled(rgb_bands+band_idx, row, rgb_buf);

			for(col=0; col<w; col++) {
				if(use_ndv && rgb_buf[col] == ndv) {
					out_buf[col] = (uint8_t)ndv;
				} else {
					double dbl_val = rgb_buf[col] * scale_buf[col];

					uint8_t byte_val = 
						dbl_val < 0 ? 0 :
						dbl_val > 255.0 ? 255 :
						(uint8_t)dbl_val;

					// avoid ndv in output
					if(use_ndv && byte_val == ndv) {
						if(ndv < 128) byte_val++;
						else byte_val--;
					}

					out_buf[col] = byte_val;
				}
			}

			GDALRasterIO(dst_bands[band_idx], GF_Write, 0, row, w, 1, out_buf, w, 1, GDT_Byte, 0, 0);
		}
	} // row

	for(ds_idx=0; ds_idx<rgb_ds_count; ds_idx++) {
		GDALClose(rgb_ds[ds_idx]);
	}
	for(ds_idx=0; ds_idx<lum_ds_count; ds_idx++) {
		GDALClose(lum_ds[ds_idx]);
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

scaled_band_t getScaledBand(GDALDatasetH lores_ds, int band_id, GDALDatasetH hires_ds) {
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

	scaled_band_t sb;

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
	sb.kernel_x = (double **)malloc_or_die(sizeof(double *) * sb.oversample);
	sb.kernel_y = (double **)malloc_or_die(sizeof(double *) * sb.oversample);
	int mod;
	for(mod=0; mod<sb.oversample; mod++) {
		sb.kernel_x[mod] = (double *)malloc_or_die(sizeof(double) * 4);
		sb.kernel_y[mod] = (double *)malloc_or_die(sizeof(double) * 4);
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

	sb.lines_buf = (double **)malloc_or_die(sizeof(double *) * 4);
	int j;
	for(j=0; j<4; j++) {
		sb.lines_buf[j] = (double *)malloc_or_die(sizeof(double) * sb.hi_w);
	}
	sb.line_buf_idx = -1000000;

	return sb;
}

void readLineScaled1D(scaled_band_t *sb, int row, double *hires_buf) {
	double *lores_buf = (double *)malloc_or_die(sizeof(double) * sb->lo_w);

	if(row<0 || row>=sb->lo_h) {
		int col;
		for(col=0; col<sb->hi_w; col++) {
			hires_buf[col] = 0;
		}
	} else {
		GDALRasterIO(sb->band, GF_Read, 0, row, sb->lo_w, 1, lores_buf, sb->lo_w, 1, GDT_Float64, 0, 0);
		int mx;
		for(mx=0; mx<sb->oversample; mx++) {
			double *kernel = sb->kernel_x[mx];
			int x0;
			for(x0=0; ; x0++) {
				int col = x0 * sb->oversample + mx;
				if(col >= sb->hi_w) break;
				int i;
				double accum = 0;
				for(i=0; i<4; i++) {
					int x = x0 - 1 + i + sb->delta_x;
					double v = (x<0 || x>=sb->lo_w) ? 0 : lores_buf[x];
					accum += v * kernel[i];
				}
				hires_buf[col] = accum;
			}
		}
	}

	free(lores_buf);
}

void readLineScaled(scaled_band_t *sb, int row, double *hires_buf) {
	int j;
	int y0 = row / sb->oversample;
	int my = row % sb->oversample;
	double *kernel = sb->kernel_y[my];

	int top_y = y0 - 1 + sb->delta_y;
	if(top_y == sb->line_buf_idx) {
		// no action
	} else if(top_y == sb->line_buf_idx+1) {
		double *tmp = sb->lines_buf[0];
		for(j=0; j<3; j++) {
			sb->lines_buf[j] = sb->lines_buf[j+1];
		}
		sb->lines_buf[3] = tmp;
		readLineScaled1D(sb, top_y+3, sb->lines_buf[3]);
		sb->line_buf_idx = top_y;
	} else {
		for(j=0; j<4; j++) {
			readLineScaled1D(sb, top_y+j, sb->lines_buf[j]);
		}
		sb->line_buf_idx = top_y;
	}

	int col;
	for(col=0; col<sb->hi_w; col++) {
		double accum = 0;
		for(j=0; j<4; j++) {
			accum += sb->lines_buf[j][col] * kernel[j];
		}
		hires_buf[col] = accum;
	}
}
