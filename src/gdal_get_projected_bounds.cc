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
#include "mask.h"
#include "mask-tracer.h"
#include "dp.h"

#include <ogrsf_frmts.h>
#include <cpl_string.h>
#include <cpl_conv.h>
#include <cpl_port.h>

void usage(const char *cmdname) {
	// FIXME
	printf("Usage:\n  %s [options] \n", cmdname);
	printf("\n");
	
	exit(1);
}

int main(int argc, char **argv) {
	const char *src_wkt_fn = NULL;
	const char *dst_wkt_fn = NULL;

	if(argc == 1) usage(argv[0]);

	int argp = 1;
	while(argp < argc) {
		char *arg = argv[argp++];
		// FIXME - check duplicate values
		if(arg[0] == '-') {
			if(!strcmp(arg, "-v")) {
				VERBOSE++;
			} else if(!strcmp(arg, "-src-wkt")) {
				if(argp == argc) usage(argv[0]);
				src_wkt_fn = argv[argp++];
			} else if(!strcmp(arg, "-dst-wkt")) {
				if(argp == argc) usage(argv[0]);
				dst_wkt_fn = argv[argp++];
			} else usage(argv[0]);
		} else {
			usage(argv[0]);
		}
	}

	GDALAllRegister();

	CPLPushErrorHandler(CPLQuietErrorHandler);

	///////////////////////
	
	mpoly_t mp = mpoly_from_wktfile(src_wkt_fn);

	return 0;
}
