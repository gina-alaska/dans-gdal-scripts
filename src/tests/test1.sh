#!/bin/sh
rm -f out_*

BINDIR=..

$BINDIR/gdal_trace_outline testcase_1.tif -ndv 255 -out-cs xy -wkt-out out_1.wkt    -report out_1.ppm -split-polys -dp-toler 0
$BINDIR/gdal_trace_outline testcase_1.tif -ndv 255 -out-cs en -wkt-out out_1_en.wkt -dp-toler 0
$BINDIR/gdal_trace_outline testcase_2.tif -ndv 255 -out-cs xy -wkt-out out_2.wkt    -report out_2.ppm -split-polys -dp-toler 0
$BINDIR/gdal_trace_outline testcase_3.tif -ndv 255 -out-cs xy -wkt-out out_3.wkt    -report out_3.ppm -split-polys -dp-toler 0
$BINDIR/gdal_trace_outline testcase_4.png -ndv '0..255 0..255 0..255 0' -out-cs xy -wkt-out out_4.wkt    -report out_4.ppm -split-polys -dp-toler 0
$BINDIR/gdal_trace_outline testcase_5.png -ndv 255 -out-cs xy -wkt-out out_5.wkt    -report out_5.ppm -split-polys -dp-toler 0
$BINDIR/gdal_trace_outline testcase_maze.png  -ndv 255 -out-cs xy -wkt-out out_maze.wkt  -report out_maze.ppm  -split-polys -dp-toler 0
$BINDIR/gdal_trace_outline testcase_noise.png -b 1 -ndv   0 -out-cs xy -wkt-out out_noise.wkt -report out_noise.ppm -split-polys -dp-toler 0
$BINDIR/gdal_trace_outline testcase_noise.png -b 1 -ndv   0 -out-cs xy -wkt-out out_noise_dp3.wkt -report out_noise_dp3.ppm -split-polys -dp-toler 3

$BINDIR/gdal_trace_outline testcase_3.tif -out-cs xy -wkt-out out_3_classify.wkt -dp-toler 0 -classify
$BINDIR/gdal_trace_outline pal.tif -out-cs xy -wkt-out out_3_classify_pal.wkt -dp-toler 0 -classify

$BINDIR/gdal_list_corners -inspect-rect4 -erosion -ndv 0 testcase_4.png -report out_4-rect.ppm > out_4-rect.wkt

$BINDIR/gdal_wkt_to_mask -wkt good_1_en.wkt -geo-from testcase_1.tif -mask-out out_1_mask.ppm

$BINDIR/gdal_get_projected_bounds -s_wkt good_1_en.wkt -s_srs '+proj=utm +zone=6 +ellps=WGS84 +units=m +no_defs ' -t_srs '+proj=stere +lat_ts=80 +lat_0=90 +lon_0=0 +ellps=WGS84' -report out_projbounds_report.ppm > out_projbounds.yml

echo '####################'

#for i in 1 2 3 4 5 ; do md5sum good-tc$i.wkt test-tc$i.wkt ; done
for i in 1 1_en 2 3 3_classify 3_classify_pal 4 5 maze noise noise_dp3 ; do 
	if diff --brief good_$i.wkt out_$i.wkt ; then
		echo "GOOD outline $i"
	else
		echo "BAD outline $i"
	fi
done

for i in 4-rect.wkt 1_mask.ppm projbounds_report.ppm projbounds.yml ; do
	if diff --brief good_$i out_$i ; then
		echo "GOOD $i"
	else
		echo "BAD $i"
	fi
done
