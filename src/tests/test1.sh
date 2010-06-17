#!/bin/sh
#../gdal_trace_outline testcase.tif  -ndv 255 -out-cs xy -wkt-out good-tc1.wkt -split-polys -dp-toler 0
#../gdal_trace_outline testcase2.tif -ndv 255 -out-cs xy -wkt-out good-tc2.wkt -split-polys -dp-toler 0
#../gdal_trace_outline testcase3.tif -ndv 255 -out-cs xy -wkt-out good-tc3.wkt -split-polys -dp-toler 0
#../gdal_trace_outline testcase4.png -ndv 255 -out-cs xy -wkt-out good-tc4.wkt -split-polys -dp-toler 0
#../gdal_trace_outline testcase5.png -ndv 255 -out-cs xy -wkt-out good-tc5.wkt -split-polys -dp-toler 0
rm -f test-tc*.wkt
../gdal_trace_outline testcase_1.tif    -ndv 255 -out-cs xy -wkt-out out_1.wkt    -report out_1.ppm    -split-polys -dp-toler 0
../gdal_trace_outline testcase_2.tif    -ndv 255 -out-cs xy -wkt-out out_2.wkt    -report out_2.ppm    -split-polys -dp-toler 0
../gdal_trace_outline testcase_3.tif    -ndv 255 -out-cs xy -wkt-out out_3.wkt    -report out_3.ppm    -split-polys -dp-toler 0
../gdal_trace_outline testcase_4.png    -ndv   0 -out-cs xy -wkt-out out_4.wkt    -report out_4.ppm    -split-polys -dp-toler 0
../gdal_trace_outline testcase_5.png    -ndv 255 -out-cs xy -wkt-out out_5.wkt    -report out_5.ppm    -split-polys -dp-toler 0
../gdal_trace_outline testcase_maze.png -ndv 255 -out-cs xy -wkt-out out_maze.wkt -report out_maze.ppm -split-polys -dp-toler 0

../gdal_list_corners -inspect-rect4 -erosion -ndv 0 testcase_4.png -report out_4-rect.ppm > out_4-rect.wkt

#for i in 1 2 3 4 5 ; do md5sum good-tc$i.wkt test-tc$i.wkt ; done
for i in 1 2 3 4 5 maze ; do 
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
