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




#ifndef GEOCODE_H
#define GEOCODE_H

typedef struct {
	char *s_srs;
	char *geo_srs;
	int w, h;
	int got_ll_en;
	int got_ul_en;
	double given_left_e, given_lower_n, given_upper_n;
	double res_x, res_y;
} geo_opts_t;

typedef struct {
	char *s_srs;
	char *geo_srs;
	double res_x, res_y; // zero if there is rotation
	double res_meters_x, res_meters_y; // zero if rotated/unknown
	char *units_name;
	double units_val;
	double semi_major, semi_minor;
	int w, h;
	OGRSpatialReferenceH spatial_ref;
	OGRSpatialReferenceH geo_sref;
	OGRCoordinateTransformationH fwd_xform;
	OGRCoordinateTransformationH inv_xform;
	double *fwd_affine;
	double *inv_affine;
	double lon_range1, lon_range2, lon_loopsize;
} georef_t;

void print_georef_usage();
geo_opts_t init_geo_options(int *argc_ptr, char ***argv_ptr);
georef_t init_georef(geo_opts_t *opt, GDALDatasetH ds);

void xy2en(georef_t *georef, double xpos, double ypos, double *e_out, double *n_out);
void en2xy(georef_t *georef, double east, double north, double *x_out, double *y_out);
int en2ll(georef_t *georef, double east, double north, double *lon_out, double *lat_out);
int ll2en(georef_t *georef, double lon, double lat, double *e_out, double *n_out);
int xy2ll(georef_t *georef, double x, double y, double *lon_out, double *lat_out);
int ll2xy(georef_t *georef, double lon, double lat, double *x_out, double *y_out);
void en2ll_or_die(georef_t *georef, double east, double north, double *lon_out, double *lat_out);
void ll2en_or_die(georef_t *georef, double lon, double lat, double *e_out, double *n_out);
void xy2ll_or_die(georef_t *georef, double x, double y, double *lon_out, double *lat_out);
void ll2xy_or_die(georef_t *georef, double lon, double lat, double *x_out, double *y_out);

#endif // ifndef GEOCODE_H
