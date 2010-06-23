#!/bin/sh

BINDIR=..

$BINDIR/gdal_trace_outline testcase_1.tif     -ndv 255 -out-cs xy -wkt-out out_1.wkt     -report out_1.ppm     -dp-toler 0
$BINDIR/gdal_trace_outline testcase_2.tif     -ndv 255 -out-cs xy -wkt-out out_2.wkt     -report out_2.ppm     -dp-toler 0
$BINDIR/gdal_trace_outline testcase_3.tif     -ndv 255 -out-cs xy -wkt-out out_3.wkt     -report out_3.ppm     -dp-toler 0
$BINDIR/gdal_trace_outline testcase_4.png     -b 1 -ndv   0 -out-cs xy -wkt-out out_4.wkt     -report out_4.ppm     -dp-toler 0
$BINDIR/gdal_trace_outline testcase_5.png     -b 1 -ndv 255 -out-cs xy -wkt-out out_5.wkt     -report out_5.ppm     -dp-toler 0
$BINDIR/gdal_trace_outline testcase_maze.png  -b 1 -ndv 255 -out-cs xy -wkt-out out_maze.wkt  -report out_maze.ppm  -dp-toler 0
$BINDIR/gdal_trace_outline testcase_noise.png -b 1 -ndv   0 -out-cs xy -wkt-out out_noise.wkt -report out_noise.ppm -dp-toler 0
$BINDIR/gdal_trace_outline testcase_noise.png -b 1 -ndv   0 -out-cs xy -wkt-out out_noise_dp3.wkt -report out_noise_dp3.ppm -dp-toler 3

for i in 1 2 3 4 5 maze noise noise_dp3 ; do 
	fn=out_$i.wkt
	echo -n "$fn: "
	( echo -n "select isvalid(mpolyfromtext('" ; cat $fn ; echo -n "'));" ) |psql $TEST_DB |grep -q '\<t\>' && echo "Good" || echo "Bad"
done
