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



#include "common.h"
#include "geocode.h"

void xy2en(
	double *affine,
	double xpos, double ypos,
	double *e_out, double *n_out
) {
	*e_out = affine[0] + affine[1] * xpos + affine[2] * ypos;
	*n_out = affine[3] + affine[4] * xpos + affine[5] * ypos;
}

void en2ll(
	OGRCoordinateTransformationH xform,
	double east, double north,
	double *lon_out, double *lat_out
) {
	if(xform) {
		if(!OCTTransform(xform, 1, &east, &north, NULL)) {
			fatal_error("OCTTransform failed");
		}
	}

	if(north < -90.0 || north > 90.0) fatal_error("latitude out of range");
	// images in latlong projection that cross the dateline can
	// have numbers outside of this range...
	//if(east < -180.0 || east > 180.0) fatal_error("longitude out of range");
	// but it shouldn't be outside of *this* range no matter what!
	if(east < -360.0 || east > 540.0) fatal_error("longitude out of range");

	*lon_out = east;
	*lat_out = north;
}

void xy2ll(
	georef_t *georef,
	double x, double y,
	double *lon_out, double *lat_out
) {
	double east, north;
	xy2en(georef->fwd_affine, x, y, &east, &north);
	en2ll(georef->fwd_xform, east, north, lon_out, lat_out);
}
