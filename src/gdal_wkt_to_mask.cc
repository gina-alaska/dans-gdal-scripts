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
#include "polygon.h"
#include "debugplot.h"
#include "georef.h"
#include "mask.h"
#include "mask-tracer.h"
#include "dp.h"

#include <ogrsf_frmts.h>
#include <cpl_string.h>
#include <cpl_conv.h>
#include <cpl_port.h>

char *read_whole_file(FILE *fin);

void usage(char *cmdname) {
	// FIXME
	printf("Usage:\n  %s [options] \n", cmdname);
	printf("\n");
	
	print_georef_usage();

	exit(1);
}

int main(int argc, char **argv) {
	char *wkt_fn = NULL;
	char *mask_fn = NULL;
	char *geo_fn = NULL;

	if(argc == 1) usage(argv[0]);

	geo_opts_t geo_opts = init_geo_options(&argc, &argv);

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

	georef_t georef = init_georef(&geo_opts, ds);
	if(!georef.inv_affine) fatal_error("missing affine transform");

	///////////////////////
	
	FILE *fh = fopen(wkt_fn, "r");
	char *wkt_in = read_whole_file(fh);
	int i;
	for(i=0; wkt_in[i]; i++) {
		if(wkt_in[i] == '\r') wkt_in[i] = ' ';
		if(wkt_in[i] == '\n') wkt_in[i] = ' ';
		if(wkt_in[i] == '\t') wkt_in[i] = ' ';
	}

	OGRGeometryH geom;
	OGRErr err = OGR_G_CreateFromWkt(&wkt_in, NULL, &geom);
	if(OGRERR_NONE != err) {
		fatal_error("OGR_G_CreateFromWkt failed: %d", err);
	}
	
	///////////////////////
	
	mpoly_t mp_en = ogr_to_mpoly(geom);
	mpoly_t *mp_xy = mpoly_en2xy(&georef, &mp_en);

	mask_from_mpoly(mp_xy, georef.w, georef.h, mask_fn);

	return 0;
}

char *read_whole_file(FILE *fin) {
	int chunk_size = 1024;
	int data_len = 0;
	int buf_len = chunk_size+1;
	char *buffer = (char *)malloc_or_die(buf_len);
	int num_read;
	while(0 < (num_read = fread(buffer+data_len, 1, chunk_size, fin))) {
		data_len += num_read;
		if(data_len+chunk_size+1 > buf_len) {
			buf_len += chunk_size;
			buffer = (char *)realloc_or_die(buffer, buf_len);
		}
	}
	buffer[data_len++] = 0;
	return buffer;
}
