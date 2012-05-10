/*
Copyright (c) 2012, Regents of the University of Alaska

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



#include <vector>
#include <string>

#include "common.h"

int VERBOSE = 0;

void fatal_error(const std::string &s) {
	fprintf(stderr, "\n\nerror:\n%s\n\n", s.c_str());
	exit(1);
}

void fatal_error(const char *fmt, ...) {
	va_list argp;
	
	fprintf(stderr, "\n\nerror:\n");
	
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);

	fprintf(stderr, "\n\n");

	exit(1);
}

std::vector<std::string> argv_to_list(int argc, char **argv) {
	std::vector<std::string> ret;
	for(int i=0; i<argc; i++) {
		ret.push_back(argv[i]);
	}
	return ret;
}
