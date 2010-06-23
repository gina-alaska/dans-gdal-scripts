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



#include <string>

#include <boost/lexical_cast.hpp>

#include "common.h"
#include "georef.h"

static const double EPSILON = 1e-9;

void usage(const std::string &cmdname); // externally defined

namespace dangdal {

void GeoOpts::printUsage() {
	printf("\
Geocoding:\n\
  -s_srs 'proj4 def'              Set or override source SRS\n\
  -geo_srs 'proj4 def'            Set or override geographic SRS (datum, etc.) used for output\n\
  -ll_en left_east lower_north    Set or override lower-left coordinate\n\
  -ul_en left_east lower_north    Set or override upper-left coordinate\n\
                                  (don't use both ll_en and ul_en)\n\
  -wh width height                Set or override image size\n\
  -res res_x res_y                Set or override resolution\n\
");
}

GeoOpts::GeoOpts(std::vector<std::string> &arg_list) :
	w(0),
	h(0),
	got_ll_en(0),
	got_ul_en(0),
	given_left_e(0),
	given_lower_n(0),
	given_upper_n(0),
	res_x(0),
	res_y(0)
{
	std::vector<std::string> args_out;
	const std::string cmdname = arg_list[0];
	args_out.push_back(cmdname);

	size_t argp = 1;
	while(argp < arg_list.size()) {
		const std::string arg = arg_list[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			try {
				if(arg == "-s_srs") {
					if(argp == arg_list.size()) usage(cmdname);
					s_srs = arg_list[argp++];
				} else if(arg == "-geo_srs") {
					if(argp == arg_list.size()) usage(cmdname);
					geo_srs = arg_list[argp++];
				} else if(arg == "-ll_en") {
					if(argp+2 > arg_list.size()) usage(cmdname);
					given_left_e = boost::lexical_cast<double>(arg_list[argp++]);
					given_lower_n = boost::lexical_cast<double>(arg_list[argp++]);
					got_ll_en++;
				} else if(arg == "-ul_en") {
					if(argp+2 > arg_list.size()) usage(cmdname);
					given_left_e = boost::lexical_cast<double>(arg_list[argp++]);
					given_upper_n = boost::lexical_cast<double>(arg_list[argp++]);
					got_ul_en++;
				} else if(arg == "-wh") {
					if(argp+2 > arg_list.size()) usage(cmdname);
					w = boost::lexical_cast<size_t>(arg_list[argp++]);
					h = boost::lexical_cast<size_t>(arg_list[argp++]);
				} else if(arg == "-res") {
					if(argp == arg_list.size()) usage(cmdname);
					res_x = boost::lexical_cast<double>(arg_list[argp++]);

					if(argp == arg_list.size()) usage(cmdname);
					res_y = boost::lexical_cast<double>(arg_list[argp++]);
				} else {
					args_out.push_back(arg);
				}
			} catch(boost::bad_lexical_cast &e) {
				fatal_error("cannot parse number given on command line");
			}
		} else {
			args_out.push_back(arg);
		}
	}

	if(got_ll_en && got_ul_en) fatal_error("don't specify both -ll_en and -ul_en");

	arg_list = args_out;
}

GeoRef::GeoRef(GeoOpts opt, const GDALDatasetH ds) {
	if(!ds && !(opt.s_srs.size() && (opt.got_ll_en || opt.got_ul_en) && 
		opt.w && opt.h && opt.res_x && opt.res_y)) fatal_error("not enough information to determine geolocation");

	spatial_ref = NULL;
	if(opt.s_srs.size()) {
		spatial_ref = OSRNewSpatialReference(NULL);
		if(OSRImportFromProj4(spatial_ref, opt.s_srs.c_str())
			!= OGRERR_NONE) fatal_error("cannot parse proj4 definition");
	} else if(ds) {
		const char *wkt = GDALGetProjectionRef(ds);
		if(wkt && strlen(wkt)) {
			//if(VERBOSE) printf("%s\n", wkt);
			spatial_ref = OSRNewSpatialReference(wkt);
		}
	}

	if(spatial_ref) {
		char *s_srs_str = NULL;
		OSRExportToProj4(spatial_ref, &s_srs_str);
		s_srs = s_srs_str;

		if(opt.geo_srs.size()) {
			geo_sref = OSRNewSpatialReference(NULL);
			if(OSRImportFromProj4(geo_sref, opt.geo_srs.c_str())
				!= OGRERR_NONE) fatal_error("cannot parse proj4 definition");
			// take only the geographic part of the definition
			geo_sref = OSRCloneGeogCS(geo_sref);
		} else {
			geo_sref = OSRCloneGeogCS(spatial_ref);
		}

		char *geo_srs_str = NULL;
		OSRExportToProj4(geo_sref, &geo_srs_str);
		geo_srs = geo_srs_str;

		fwd_xform = OCTNewCoordinateTransformation(spatial_ref, geo_sref);
		inv_xform = OCTNewCoordinateTransformation(geo_sref, spatial_ref);
	} else {
		fwd_xform = NULL;
		inv_xform = NULL;
		s_srs.erase();
		geo_sref = NULL;
	}

	if(ds) {
		if(!opt.w) opt.w = GDALGetRasterXSize(ds);
		if(!opt.h) opt.h = GDALGetRasterYSize(ds);
	}
	if(!opt.w || !opt.h) fatal_error("missing width/height");

	opt.res_x = fabs(opt.res_x);
	opt.res_y = fabs(opt.res_y);

	if((opt.got_ul_en || opt.got_ll_en) && opt.res_x && opt.res_y) {
		if(opt.got_ll_en) {
			opt.given_upper_n = opt.given_lower_n + (double)opt.h*opt.res_y;
			opt.got_ul_en = 1;
		}

		if(!opt.got_ul_en) fatal_error("impossibility");
		fwd_affine.resize(6);
		fwd_affine[0] = opt.given_left_e;  
		fwd_affine[3] = opt.given_upper_n; 
		fwd_affine[1] = opt.res_x; fwd_affine[2] =          0;
		fwd_affine[4] =         0; fwd_affine[5] = -opt.res_y;
	} else if(ds) {
		bool has_rotation = 0;
		fwd_affine.resize(6);
		if(GDALGetGeoTransform(ds, &fwd_affine[0]) == CE_None) {
			has_rotation = fwd_affine[2] || fwd_affine[4];
		} else {
			fwd_affine.clear();
		}

		if(opt.res_x && opt.res_y) {
			// if corner coordinate were specified, the first branch of the outer if
			// statement would have been taken
			if(!hasAffine()) fatal_error("missing ll_en/ul_en parameter");

			if(has_rotation) fatal_error("cannot override resolution if source image has rotation");

			fwd_affine[1] =  opt.res_x;
			fwd_affine[5] = -opt.res_y;
		} else {
			if(has_rotation || !hasAffine()) {
				opt.res_x = opt.res_y = 0;
			} else {
				opt.res_x = fabs(fwd_affine[1]);
				opt.res_y = fabs(fwd_affine[5]);
			}
		}

		if(opt.got_ul_en || opt.got_ll_en) {
			if(has_rotation) fatal_error("cannot override ll_en/ul_en if source image has rotation");

			if(!hasAffine() || !opt.res_x || !opt.res_y) fatal_error("missing -res parameter");

			if(opt.got_ll_en) {
				opt.given_upper_n = opt.given_lower_n + (double)opt.h*opt.res_y;
				opt.got_ul_en = 1;
			}

			if(!opt.got_ul_en) fatal_error("impossibility");
			fwd_affine[0] = opt.given_left_e;
			fwd_affine[3] = opt.given_upper_n;
		}
	}

	if(hasAffine()) {
		inv_affine.resize(6);
		if(!GDALInvGeoTransform(&fwd_affine[0], &inv_affine[0])) {
			fatal_error("affine is not invertible");
		}
	}

	res_x = opt.res_x;
	res_y = opt.res_y;
	w = opt.w;
	h = opt.h;

	res_meters_x = 0;
	res_meters_y = 0;
	units_val = 0;
	have_semi_major = false;
	if(spatial_ref) {
		OGRErr err = OGRERR_NONE;
		semi_major = OSRGetSemiMajor(spatial_ref, &err);
		if(err == OGRERR_NONE) {
			have_semi_major = true;
		}

		if(OSRIsProjected(spatial_ref)) {
			char *units_name_str = NULL;
			units_val = OSRGetLinearUnits(spatial_ref, &units_name_str);
			//printf("units: %s, %lf\n", units_name_str?units_name_str:"null", units_val);
			if(units_name_str) units_name = units_name_str;
			res_meters_x = units_val * res_x;
			res_meters_y = units_val * res_y;
		} else if(OSRIsGeographic(spatial_ref)) {
			char *units_name_str = NULL;
			units_val = OSRGetAngularUnits(spatial_ref, &units_name_str);
			//printf("ang units: %s, %lf\n", units_name_str?units_name_str:"null", units_val);
			if(units_name_str) units_name = units_name_str;

			// FIXME - what is the best way to convert degrees to meters on ellipsoid?
			// The X-resolution will be fictional anyway since the size of a degree
			// of longitude varies depending on latitude.

			if(have_semi_major) {
				res_meters_x = semi_major * units_val * res_x;
				res_meters_y = semi_major * units_val * res_y;
			}
		}
	}

	// For geographic projections, compute the maximum/minimum longitude.  Coordinates
	// will be mapped into this range by adding a multiple of lon_loopsize, which
	// corresponds to 360 degrees.
	if(spatial_ref && OSRIsGeographic(spatial_ref) && units_val) {
		double east, north;
		double min_lon, max_lon;
		xy2en(0, 0, &east, &north);
		min_lon = max_lon = east;
		xy2en(w, 0, &east, &north);
		if(east < min_lon) min_lon = east;
		if(east > max_lon) max_lon = east;
		xy2en(0, h, &east, &north);
		if(east < min_lon) min_lon = east;
		if(east > max_lon) max_lon = east;
		xy2en(w, h, &east, &north);
		if(east < min_lon) min_lon = east;
		if(east > max_lon) max_lon = east;
		lon_range1 = min_lon;
		lon_range2 = max_lon;
		lon_loopsize = 2.0 * M_PI / units_val;
		double span = lon_range2 - lon_range1;
		if(span - lon_loopsize > 1e-12) {
			fprintf(stderr, "WARNING: input spans more than 360 degrees - projection is not single-valued\n");
		}
	} else {
		lon_range1 = 0;
		lon_range2 = 0;
		lon_loopsize = 0;
	}

	if(VERBOSE) {
		printf("lon_range1 = %lf\n", lon_range1);
		printf("lon_range2 = %lf\n", lon_range2);
		printf("lon_loopsize = %lf\n", lon_loopsize);
	}
}

void GeoRef::xy2en(
	double xpos, double ypos,
	double *e_out, double *n_out
) const {
	if(!hasAffine()) fatal_error("missing affine");
	*e_out = fwd_affine[0] + fwd_affine[1] * xpos + fwd_affine[2] * ypos;
	*n_out = fwd_affine[3] + fwd_affine[4] * xpos + fwd_affine[5] * ypos;
}

void GeoRef::en2xy(
	double east, double north,
	double *x_out, double *y_out
) const {
	if(!hasAffine()) fatal_error("missing affine");
	*x_out = inv_affine[0] + inv_affine[1] * east + inv_affine[2] * north;
	*y_out = inv_affine[3] + inv_affine[4] * east + inv_affine[5] * north;
}

bool GeoRef::en2ll(
	double east, double north,
	double *lon_out, double *lat_out
) const {
	if(!fwd_xform) fatal_error("missing xform");

	double u = east;
	double v = north;
	if(!OCTTransform(fwd_xform, 1, &u, &v, NULL)) {
		return 1;
	}
	double lon = u;
	double lat = v;

	if(lat < -90.0-EPSILON || lat > 90.0+EPSILON) return 1; //fatal_error("latitude out of range (%lf)", lat);
	// images in latlong projection that cross the dateline can
	// have numbers outside of this range...
	//if(lon < -180.0 || lon > 180.0) return 1; //fatal_error("longitude out of range");
	// but it shouldn't be outside of *this* range no matter what!
	if(lon < -360.0-EPSILON || lon > 540.0+EPSILON) return 1; //fatal_error("longitude out of range (%lf)", lon);

	*lon_out = lon;
	*lat_out = lat;
	return 0;
}

void GeoRef::en2ll_or_die(
	double east, double north,
	double *lon_out, double *lat_out
) const {
	if(en2ll(east, north, lon_out, lat_out)) {
		fatal_error("en2ll transform failed [%g,%g]", east, north);
	}
}

bool GeoRef::ll2en(
	double lon, double lat,
	double *e_out, double *n_out
) const {
	if(!inv_xform) fatal_error("missing xform");

	if(lat < -90.0-EPSILON || lat > 90.0+EPSILON) return 1; //fatal_error("latitude out of range (%lf)", lat);
	// images in latlong projection that cross the dateline can
	// have numbers outside of this range...
	//if(lon < -180.0 || lon > 180.0) return 1; //fatal_error("longitude out of range");
	// but it shouldn't be outside of *this* range no matter what!
	if(lon < -360.0-EPSILON || lon > 540.0+EPSILON) return 1; //fatal_error("longitude out of range (%lf)", lon);

	double u = lon;
	double v = lat;
	if(!OCTTransform(inv_xform, 1, &u, &v, NULL)) {
		return 1;
	}
	double east = u;
	double north = v;

	// This will add a multiple of 360 degrees in order to bring the
	// coordinate into the proper range.  This is needed because
	// OCTTransform will usually return a number in the -180..180
	// range, but the raster may be defined on, for example, a range
	// of 0..360.
	if(lon_loopsize) {
		double east_orig = east;
		while(east < lon_range1) east += lon_loopsize;
		while(east > lon_range2) east -= lon_loopsize;
		if(east < lon_range1) east = east_orig;
		//printf("%g => %g\n", east_orig, east);
	}

	*e_out = east;
	*n_out = north;
	return 0;
}

void GeoRef::ll2en_or_die(
	double lon, double lat,
	double *e_out, double *n_out
) const {
	if(ll2en(lon, lat, e_out, n_out)) {
		fatal_error("ll2en transform failed [%g,%g]", lon, lat);
	}
}

bool GeoRef::xy2ll(
	double x, double y,
	double *lon_out, double *lat_out
) const {
	double east, north;
	xy2en(x, y, &east, &north);
	return en2ll(east, north, lon_out, lat_out);
}

void GeoRef::xy2ll_or_die(
	double x, double y,
	double *lon_out, double *lat_out
) const {
	double east, north;
	xy2en(x, y, &east, &north);
	en2ll_or_die(east, north, lon_out, lat_out);
}

bool GeoRef::ll2xy(
	double lon, double lat,
	double *x_out, double *y_out
) const {
	double east, north;
	bool ret = ll2en(lon, lat, &east, &north);
	if(ret) return ret;
	en2xy(east, north, x_out, y_out);
	return 0;
}

void GeoRef::ll2xy_or_die(
	double lon, double lat,
	double *x_out, double *y_out
) const {
	double east, north;
	ll2en_or_die(lon, lat, &east, &north);
	en2xy(east, north, x_out, y_out);
}

} // namespace dangdal
