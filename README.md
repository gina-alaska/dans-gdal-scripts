dans-gdal-scripts
=================

This package consists of a number of utilities for use in conjunction with
GDAL.  Currently these programs are only supported on Unix systems, but it
should be possible to compile them on Windows if you know how to make the
corresponding Makefiles.

## Programs

### gdal_contrast_stretch      
Contrast stretch and conversion from 8,16-bit or floating point to 8-bit

### gdal_dem2rgb               
Generate hillshaded images from DEMs

### gdal_get_projected_bounds  
Project a polygon and return its bounding rectangle

### gdal_landsat_pansharp      
Pansharpening - works best for Landsat 7 images

### gdal_list_corners          
Prints raster geocode information in YAML format (similar to gdalinfo but gives YAML)

### gdal_merge_simple          
Merge individual bands into a single GeoTIFF image (8-bit only)

### gdal_merge_vrt             
Merge individual bands into a single VRT image

### gdal_raw2geotiff           
Convert raw binary files into GeoTIFFs

### gdal_trace_outline         
Trace the outline of an image and generate WKT or Shapefile (options exist for cleaning up raggedy edges, can also do feature classification)

### gdal_wkt_to_mask           
Generate a bitmap of the area covered by a polygon

## Install

Debian/ubuntu: apt-get install dans-gdal-scripts

from source:

    ./autogen.sh
    ./configure  # --prefix=/your/favorite/prefix
    make
    make install

### OSX (homebrew)

To get dans-gdal-scripts working on OSX you'll need to get a few deps:

* [Homebrew](http://mxcl.github.io/homebrew/)
 * brew install gdal
 * brew install boost

## Credits

These programs were written by Dan Stahlke of the Geographic Information Network of Alaska.
Send questions or comments to [dan@stahlke.org](mailto:dan@stahlke.org).

