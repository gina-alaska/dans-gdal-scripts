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




#ifndef MASK_H
#define MASK_H

#include "common.h"
#include "polygon.h"
#include "debugplot.h"
#include "ndv.h"

unsigned char *get_mask_for_dataset(GDALDatasetH ds, int bandlist_size, int *bandlist, 
	ndv_def_t *ndv_def, report_image_t *dbuf);
unsigned char *read_dataset_8bit(GDALDatasetH ds, int band_idx, unsigned char *usage_array, report_image_t *dbuf);
unsigned char *get_mask_for_8bit_raster(int w, int h, const unsigned char *raster, unsigned char wanted);
void erode_mask(unsigned char *in_mask, int w, int h);
void invert_mask(unsigned char *in_mask, int w, int h);
vertex_t calc_centroid_from_mask(const unsigned char *mask, int w, int h);

#endif // ifndef MASK_H
