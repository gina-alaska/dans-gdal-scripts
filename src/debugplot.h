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



#ifndef DEBUGPLOT_H
#define DEBUGPLOT_H

#include "common.h"
#include "polygon.h"

#define PLOT_RECT4 1
#define PLOT_DESCENDERS 2
#define PLOT_CONTOURS 3
typedef struct {
	unsigned char *img;
	double canvas_w, canvas_h;
	int img_w, img_h;
	int stride_x, stride_y;
	int mode;
} report_image_t;

report_image_t *create_plot(double w, double h);
void write_plot(report_image_t *dbuf, char *fn);
void plot_point_big(report_image_t *dbuf, double x, double y, unsigned char r, unsigned char g, unsigned char b);
void plot_point(report_image_t *dbuf, double x, double y, unsigned char r, unsigned char g, unsigned char b);
void plot_line(report_image_t *dbuf, vertex_t p0, vertex_t p1, 
	unsigned char r, unsigned char g, unsigned char b);
void debug_plot_contours(mpoly_t *mpoly, report_image_t *dbuf);

#endif // ifndef DEBUGPLOT_H
