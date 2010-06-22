# Copyright (c) 2009, Regents of the University of Alaska
# 
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
#     * Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the Geographic Information Network of Alaska nor
#       the names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



AM_CPPFLAGS = @GDALCFLAGS@ -Wall -Wextra -Werror -O2
LIBS = @GDALLIBS@

bin_PROGRAMS = gdal_raw2geotiff gdal_dem2rgb gdal_list_corners gdal_trace_outline gdal_contrast_stretch gdal_landsat_pansharp gdal_wkt_to_mask gdal_merge_simple gdal_merge_vrt gdal_get_projected_bounds gdal_make_ndv_mask

default_palette.h: default_palette.pal
	(echo "const char *DEFAULT_PALETTE[] = {"; sed 's/^/\t\"/;s/$$/",/' < default_palette.pal; echo "	NULL };") > default_palette.h

gdal_raw2geotiff_SOURCES = gdal_raw2geotiff.cc common.cc

gdal_dem2rgb.o: default_palette.h
gdal_dem2rgb_SOURCES = gdal_dem2rgb.cc common.cc georef.cc ndv.cc palette.cc

gdal_list_corners_SOURCES = gdal_list_corners.cc common.cc polygon.cc polygon-rasterizer.cc debugplot.cc georef.cc mask.cc rectangle_finder.cc ndv.cc

gdal_trace_outline_SOURCES = gdal_trace_outline.cc common.cc polygon.cc polygon-rasterizer.cc debugplot.cc georef.cc mask.cc mask-tracer.cc beveler.cc dp.cc ndv.cc excursion_pincher2.cc

gdal_contrast_stretch_SOURCES = gdal_contrast_stretch.cc common.cc ndv.cc

gdal_landsat_pansharp_SOURCES = gdal_landsat_pansharp.cc common.cc

gdal_wkt_to_mask_SOURCES = gdal_wkt_to_mask.cc common.cc polygon.cc polygon-rasterizer.cc georef.cc

gdal_get_projected_bounds_SOURCES = gdal_get_projected_bounds.cc common.cc polygon.cc georef.cc debugplot.cc

gdal_merge_simple_SOURCES = gdal_merge_simple.cc common.cc

gdal_merge_vrt_SOURCES = gdal_merge_vrt.cc common.cc

gdal_make_ndv_mask_SOURCES = gdal_make_ndv_mask.cc common.cc ndv.cc mask.cc debugplot.cc

lint:
	cpplint.py --filter=-whitespace,-readability/streams,-build/header_guard,-build/include_order,-readability/multiline_string \
	*.cc *.h 2>&1 \
	| grep -v 'Include the directory when naming .h files' \
	| grep -v 'Is this a non-const reference' \
	| grep -v 'Do not use namespace using-directives.' \
	| grep -v 'Done processing '

noinst_HEADERS = common.h debugplot.h default_palette.h dp.h excursion_pincher.h georef.h mask-tracer.h mask.h ndv.h polygon-rasterizer.h polygon.h rectangle_finder.h
EXTRA_DIST = default_palette.pal tests/good-tc1.wkt tests/good-tc2.wkt tests/good-tc3.wkt tests/good-tc4.wkt tests/good-tc5.wkt tests/pal.tif tests/test1.sh tests/test2.sh tests/testcase.tif tests/testcase2.tif tests/testcase3.tif tests/testcase4.png tests/testcase5.png attic/excursion_pincher.cc attic/excursion_pincher3.cc attic/hough.cc attic/hough.h
