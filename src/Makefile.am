AM_CFLAGS   = @GDAL_CFLAGS@ -Wall
AM_CXXFLAGS = @GDAL_CFLAGS@ -Wall
LIBS = @GDAL_LIBS@

bin_PROGRAMS = gdal_raw2geotiff gdal_dem2rgb gdal_list_corners gdal_trace_outline gdal_contrast_stretch gdal_landsat_pansharp gdal_wkt_to_mask gdal_merge_simple gdal_merge_vrt gdal_get_projected_bounds

default_palette.h: default_palette.pal
	(echo "const char *DEFAULT_PALETTE[] = {"; sed 's/^/\t\"/;s/$$/",/' < default_palette.pal; echo "	NULL };") > default_palette.h

gdal_raw2geotiff_SOURCES = gdal_raw2geotiff.cc common.cc
gdal_raw2geotiff_CXXFLAGS = $(AM_CXXFLAGS)

gdal_dem2rgb_SOURCES = gdal_dem2rgb.cc common.cc georef.cc ndv.cc
gdal_dem2rgb_CXXFLAGS = $(AM_CXXFLAGS)

gdal_list_corners_SOURCES = gdal_list_corners.cc common.cc polygon.cc polygon-rasterizer.cc debugplot.cc georef.cc mask.cc rectangle_finder.cc ndv.cc
gdal_list_corners_CXXFLAGS = $(AM_CXXFLAGS)

gdal_trace_outline_SOURCES = gdal_trace_outline.cc common.cc polygon.cc polygon-rasterizer.cc debugplot.cc georef.cc mask.cc mask-tracer.cc beveler.cc dp.cc ndv.cc excursion_pincher2.cc
gdal_trace_outline_CXXFLAGS = $(AM_CXXFLAGS)

gdal_contrast_stretch_SOURCES = gdal_contrast_stretch.cc common.cc
gdal_contrast_stretch_CXXFLAGS = $(AM_CXXFLAGS)

gdal_landsat_pansharp_SOURCES = gdal_landsat_pansharp.cc common.cc
gdal_landsat_pansharp_CXXFLAGS = $(AM_CXXFLAGS)

gdal_wkt_to_mask_SOURCES = gdal_wkt_to_mask.cc common.cc polygon.cc polygon-rasterizer.cc georef.cc
gdal_wkt_to_mask_CXXFLAGS = $(AM_CXXFLAGS)

gdal_get_projected_bounds_SOURCES = gdal_get_projected_bounds.cc common.cc polygon.cc georef.cc debugplot.cc
gdal_get_projected_bounds_CXXFLAGS = $(AM_CXXFLAGS)

gdal_merge_simple_SOURCES = gdal_merge_simple.cc common.cc
gdal_merge_simple_CXXFLAGS = $(AM_CXXFLAGS)

gdal_merge_vrt_SOURCES = gdal_merge_vrt.cc common.cc
gdal_merge_vrt_CXXFLAGS = $(AM_CXXFLAGS)

noinst_HEADERS = common.h debugplot.h default_palette.h dp.h excursion_pincher.h georef.h mask-tracer.h mask.h ndv.h polygon-rasterizer.h polygon.h rectangle_finder.h
EXTRA_DIST = default_palette.pal tests/good-tc1.wkt tests/good-tc2.wkt tests/good-tc3.wkt tests/good-tc4.wkt tests/good-tc5.wkt tests/pal.tif tests/test1.sh tests/test2.sh tests/testcase.tif tests/testcase2.tif tests/testcase3.tif tests/testcase4.png tests/testcase5.png
