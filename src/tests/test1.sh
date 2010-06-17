#!/bin/sh
rm -f out_*.wkt
../gdal_trace_outline testcase_1.tif     -ndv 255 -out-cs xy -wkt-out out_1.wkt     -report out_1.ppm     -split-polys -dp-toler 0
../gdal_trace_outline testcase_2.tif     -ndv 255 -out-cs xy -wkt-out out_2.wkt     -report out_2.ppm     -split-polys -dp-toler 0
../gdal_trace_outline testcase_3.tif     -ndv 255 -out-cs xy -wkt-out out_3.wkt     -report out_3.ppm     -split-polys -dp-toler 0
../gdal_trace_outline testcase_4.png     -b 1 -ndv   0 -out-cs xy -wkt-out out_4.wkt     -report out_4.ppm     -split-polys -dp-toler 0
../gdal_trace_outline testcase_5.png     -b 1 -ndv 255 -out-cs xy -wkt-out out_5.wkt     -report out_5.ppm     -split-polys -dp-toler 0
../gdal_trace_outline testcase_maze.png  -b 1 -ndv 255 -out-cs xy -wkt-out out_maze.wkt  -report out_maze.ppm  -split-polys -dp-toler 0
../gdal_trace_outline testcase_noise.png -b 1 -ndv   0 -out-cs xy -wkt-out out_noise.wkt -report out_noise.ppm -split-polys -dp-toler 0

../gdal_list_corners -inspect-rect4 -erosion -ndv 0 testcase_4.png -report out_4-rect.ppm > out_4-rect.wkt

#for i in 1 2 3 4 5 ; do md5sum good-tc$i.wkt test-tc$i.wkt ; done
for i in 1 2 3 4 5 maze noise ; do 
	if diff --brief good_$i.wkt out_$i.wkt ; then
		echo "outline $i - Good"
	else
		echo "outline $i - Bad"
	fi
done

if diff --brief good_4-rect.wkt out_4-rect.wkt ; then
	echo "rect - Good"
else
	echo "rect - Bad"
fi
