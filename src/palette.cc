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



#include "common.h"
#include "palette.h"
#include "default_palette.h"

namespace dangdal {

Palette Palette::fromLines(const char * const *lines) {
	std::vector<std::string> lv;
	while(*lines) {
		lv.push_back(*lines);
		lines++;
	}
	return fromLines(lv);
}

Palette Palette::fromLines(const std::vector<std::string> &lines) {
	Palette p;

	p.nan_color = RGB(0, 0, 0);

	for(size_t line_num=0; line_num<lines.size(); line_num++) {
		const std::string &line = lines[line_num];
		if(line.empty() || line[0] == '#') continue;
		RGB rgb;
		double val;
		if(
			4 != sscanf(line.c_str(), "%100lf %100hhd %100hhd %100hhd\n",
				&val, &rgb.r, &rgb.g, &rgb.b)
		) {
			fatal_error("cannot parse line in palette file: [%s]", line.c_str());
		}
		if(std::isnan(val)) {
			p.nan_color = rgb;
		} else {
			p.vals.push_back(val);
			p.colors.push_back(rgb);
		}
	}
	
	if(p.vals.size() < 2) fatal_error("not enough entries in palette");

	return p;
}

Palette Palette::fromFile(const std::string &fn) {
	FILE *fh = fopen(fn.c_str(), "r");
	if(!fh) fatal_error("cannot open palette file");
	std::vector<std::string> lines;
	char *line;
	char buf[1000]; // FIXME - read it the c++ way
	while((line = fgets(buf, 1000, fh))) {
		lines.push_back(line);
	}
	fclose(fh);

	if(lines.empty()) fatal_error("palette file was empty");

	return fromLines(lines);
}

Palette Palette::createDefault() {
	return fromLines(DEFAULT_PALETTE);
}

RGB Palette::get(double val) const {
	if(val < vals[0]) val = vals[0];
	if(val > vals[vals.size()-1]) val = vals[vals.size()-1];

	// FIXME - it would be better to use std::lower_bound
	for(size_t i=0; i<vals.size()-1; i++) {
		double v1 = vals[i];
		double v2 = vals[i+1];
		if(val >= v1 && val <= v2) {
			double alpha = (val - v1) / (v2 - v1);
			RGB c1 = colors[i];
			RGB c2 = colors[i+1];
			return RGB(
				uint8_t(c1.r*(1.0-alpha) + c2.r*alpha + .5),
				uint8_t(c1.g*(1.0-alpha) + c2.g*alpha + .5),
				uint8_t(c1.b*(1.0-alpha) + c2.b*alpha + .5)
			);
		}
	}

	fatal_error("palette file out of sequence, can't find val %g\n", val);
}

} // namespace dangdal
