#!/bin/bash

rm -f out_test2_*

#BINDIR="valgrind -q .."
BINDIR=..

$BINDIR/gdal_trace_outline testcase_1.tif     -ndv 255 -out-cs xy -wkt-out out_test2_1.wkt     -report out_test2_1.ppm     -dp-toler 0
$BINDIR/gdal_trace_outline testcase_2.tif     -ndv 255 -out-cs xy -wkt-out out_test2_2.wkt     -report out_test2_2.ppm     -dp-toler 0
$BINDIR/gdal_trace_outline testcase_3.tif     -ndv 255 -out-cs xy -wkt-out out_test2_3.wkt     -report out_test2_3.ppm     -dp-toler 0
$BINDIR/gdal_trace_outline testcase_4.png     -b 1 -ndv   0 -out-cs xy -wkt-out out_test2_4.wkt     -report out_test2_4.ppm     -dp-toler 0 -min-ring-area 100
$BINDIR/gdal_trace_outline testcase_5.png     -b 1 -ndv 255 -out-cs xy -wkt-out out_test2_5.wkt     -report out_test2_5.ppm     -dp-toler 0
$BINDIR/gdal_trace_outline testcase_maze.png  -b 1 -ndv 255 -out-cs xy -wkt-out out_test2_maze.wkt  -report out_test2_maze.ppm  -dp-toler 0
$BINDIR/gdal_trace_outline testcase_noise.png -b 1 -ndv   0 -out-cs xy -wkt-out out_test2_noise.wkt -report out_test2_noise.ppm -dp-toler 0
$BINDIR/gdal_trace_outline testcase_noise.png -b 1 -ndv   0 -out-cs xy -wkt-out out_test2_noise_dp3.wkt -report out_test2_noise_dp3.ppm -dp-toler 3

echo '####################'

for i in out_test2_* ; do
	if [ ! -e ${i/out/good} ] ; then
		echo "!!! ${i/out/good} doesn't exist"
	fi
done

for i in good_test2_* ; do
	if diff --brief $i ${i/good/out} ; then
		echo "GOOD ${i/good_/}"
	else
		echo "BAD ${i/good_/}"
	fi
done

if [ "z$TEST_DB" != "z" ]; then
	for i in 1 2 3 4 5 maze noise noise_dp3 ; do 
		fn=out_test2_$i.wkt
		echo -n "postgis test of $fn: "
		( echo -n "select isvalid(mpolyfromtext('" ; cat $fn ; echo -n "'));" ) |psql $TEST_DB |grep -q '\<t\>' && echo "Good" || echo "Bad"
	done
fi
