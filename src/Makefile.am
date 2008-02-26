AM_CFLAGS = -Wall -g -O3 @GDAL_CFLAGS@
AM_CXXFLAGS = -Wall -g -O3 @GDAL_CFLAGS@
LIBS = @GDAL_LIBS@

bin_PROGRAMS = gdal_raw2geotiff gdal_dem2rgb gdal_list_corners gdal_contrast_stretch gdal_landsat_pansharp gdal_merge_simple gdal_merge_vrt

gdal_raw2geotiff_SOURCES = gdal_raw2geotiff.cc common.cc
gdal_raw2geotiff_CXXFLAGS = $(AM_CXXFLAGS)

gdal_dem2rgb_SOURCES = gdal_dem2rgb.cc common.cc
gdal_dem2rgb_CXXFLAGS = $(AM_CXXFLAGS)

gdal_list_corners_SOURCES = gdal_list_corners.cc common.cc polygon.cc debugplot.cc georef.cc mask.cc
gdal_list_corners_CXXFLAGS = $(AM_CXXFLAGS)

gdal_contrast_stretch_SOURCES = gdal_contrast_stretch.cc common.cc
gdal_contrast_stretch_CXXFLAGS = $(AM_CXXFLAGS)

gdal_landsat_pansharp_SOURCES = gdal_landsat_pansharp.cc common.cc
gdal_landsat_pansharp_CXXFLAGS = $(AM_CXXFLAGS)

gdal_merge_simple_SOURCES = gdal_merge_simple.cc common.cc
gdal_merge_simple_CXXFLAGS = $(AM_CXXFLAGS)

gdal_merge_vrt_SOURCES = gdal_merge_vrt.cc common.cc
gdal_merge_vrt_CXXFLAGS = $(AM_CXXFLAGS)

noinst_HEADERS = common.h polygon.h debugplot.h georef.h mask.h
EXTRA_DIST = dem.pal test2.sh testcase.tif testcase2.tif
