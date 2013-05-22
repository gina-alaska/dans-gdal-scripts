#!/bin/sh

rm -f out_test3_test3_*

#BINDIR="valgrind -q .."
BINDIR=..

$BINDIR/gdal_dem2rgb nedcut.tif out_test3_dem.tif -default-palette && tifftopnm out_test3_dem.tif >out_test3_dem.pnm
$BINDIR/gdal_contrast_stretch -ndv '0..255 0..255 0..255 0' -histeq 50 testcase_4.png out_test3_histeq.tif && tifftopnm out_test3_histeq.tif >out_test3_histeq.pnm

echo '####################'

for i in good_test3_* ; do
	if diff --brief $i ${i/good/out} ; then
		echo "GOOD ${i/good_/}"
	else
		echo "BAD ${i/good_/}"
	fi
done
