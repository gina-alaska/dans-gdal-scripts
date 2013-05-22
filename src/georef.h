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




#ifndef DANGDAL_GEOCODE_H
#define DANGDAL_GEOCODE_H

#include <string>
#include <vector>

namespace dangdal {

struct GeoOpts {
	static void printUsage();
	explicit GeoOpts(std::vector<std::string> &arg_list);

	std::string s_srs;
	std::string geo_srs;
	size_t w, h;
	bool got_ll_en;
	bool got_ul_en;
	double given_left_e, given_lower_n, given_upper_n;
	double res_x, res_y;
};

class GeoRef {
public:
	GeoRef(GeoOpts opt, const GDALDatasetH ds);

	bool hasAffine() const { return !fwd_affine.empty(); }

	void xy2en(double xpos, double ypos, double *e_out, double *n_out) const;
	void en2xy(double east, double north, double *x_out, double *y_out) const;
	bool en2ll(double east, double north, double *lon_out, double *lat_out) const;
	bool ll2en(double lon, double lat, double *e_out, double *n_out) const;
	bool xy2ll(double x, double y, double *lon_out, double *lat_out) const;
	bool ll2xy(double lon, double lat, double *x_out, double *y_out) const;
	void en2ll_or_die(double east, double north, double *lon_out, double *lat_out) const;
	void ll2en_or_die(double lon, double lat, double *e_out, double *n_out) const;
	void xy2ll_or_die(double x, double y, double *lon_out, double *lat_out) const;
	void ll2xy_or_die(double lon, double lat, double *x_out, double *y_out) const;

	std::string s_srs;
	std::string geo_srs;
	double res_x, res_y; // zero if there is rotation
	double res_meters_x, res_meters_y; // zero if rotated/unknown
	std::string units_name;
	double units_val;
	bool have_semi_major;
	double semi_major;
	size_t w, h;
	OGRSpatialReferenceH spatial_ref;
	OGRSpatialReferenceH geo_sref;
	OGRCoordinateTransformationH fwd_xform;
	OGRCoordinateTransformationH inv_xform;
	std::vector<double> fwd_affine;
	std::vector<double> inv_affine;
	double lon_range1, lon_range2, lon_loopsize;
};

} // namespace dangdal

#endif // ifndef DANGDAL_GEOCODE_H
