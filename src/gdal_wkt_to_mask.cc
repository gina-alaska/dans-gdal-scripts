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
#include "polygon.h"
#include "polygon-rasterizer.h"
#include "debugplot.h"
#include "georef.h"

#include <ogrsf_frmts.h>
#include <cpl_string.h>
#include <cpl_conv.h>
#include <cpl_port.h>

using namespace dangdal;

void usage(const char *cmdname) {
	printf("Usage:\n  %s [options] \n", cmdname);
	printf("\n");
	
	GeoOpts::printUsage();
	printf("  -geo-from <fn>                  Get georeference from this raster file\n");
	printf("\nOptions:\n");
	printf("  -wkt <fn>                       File containing WKT def in easting/northing units\n");
	printf("  -mask-out <fn.pbm>              Filename for mask output in PBM format\n");

	exit(1);
}

int main(int argc, char **argv) {
	char *wkt_fn = NULL;
	char *mask_fn = NULL;
	char *geo_fn = NULL;

	if(argc == 1) usage(argv[0]);

	GeoOpts geo_opts = GeoOpts(&argc, &argv);

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(!strcmp(arg, "-v")) {
				VERBOSE++;
			} else if(!strcmp(arg, "-wkt")) {
				if(argp == argc) usage(argv[0]);
				wkt_fn = argv[argp++];
			} else if(!strcmp(arg, "-mask-out")) {
				if(argp == argc) usage(argv[0]);
				mask_fn = argv[argp++];
			} else if(!strcmp(arg, "-geo-from")) {
				if(argp == argc) usage(argv[0]);
				geo_fn = argv[argp++];
			} else usage(argv[0]);
		} else {
			usage(argv[0]);
		}
	}

	if(!wkt_fn || !mask_fn) usage(argv[0]);

	GDALAllRegister();

	GDALDatasetH ds = NULL;
	if(geo_fn) {
		ds = GDALOpen(geo_fn, GA_ReadOnly);
		if(!ds) fatal_error("open failed");
	}

	CPLPushErrorHandler(CPLQuietErrorHandler);

	GeoRef georef = GeoRef(geo_opts, ds);
	if(!georef.hasAffine()) fatal_error("missing affine transform");

	Mpoly mp = mpoly_from_wktfile(wkt_fn);

	mp.en2xy(georef);

	mask_from_mpoly(mp, georef.w, georef.h, mask_fn);

	return 0;
}
