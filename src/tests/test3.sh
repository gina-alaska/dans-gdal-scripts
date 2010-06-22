#!/bin/sh
rm -f out_*

../gdal_dem2rgb nedcut.tif out_dem.tif -default-palette && tifftopnm out_dem.tif >out_dem.pnm
../gdal_contrast_stretch -ndv '0..255 0..255 0..255 0' -histeq 50 testcase_4.png out_histeq.tif && tifftopnm out_histeq.tif >out_histeq.pnm

echo '####################'

for i in dem.pnm histeq.pnm ; do
	if diff --brief good_$i out_$i ; then
		echo "GOOD $i"
	else
		echo "BAD $i"
	fi
done
