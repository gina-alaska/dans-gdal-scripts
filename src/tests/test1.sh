#!/bin/sh
#../gdal_trace_outline testcase.tif  -ndv 255 -out-cs xy -wkt-out good-tc1.wkt -split-polys -dp-toler 0
#../gdal_trace_outline testcase2.tif -ndv 255 -out-cs xy -wkt-out good-tc2.wkt -split-polys -dp-toler 0
#../gdal_trace_outline testcase3.tif -ndv 255 -out-cs xy -wkt-out good-tc3.wkt -split-polys -dp-toler 0
#../gdal_trace_outline testcase4.png -ndv 255 -out-cs xy -wkt-out good-tc4.wkt -split-polys -dp-toler 0
#../gdal_trace_outline testcase5.png -ndv 255 -out-cs xy -wkt-out good-tc5.wkt -split-polys -dp-toler 0
rm -f test-tc*.wkt
../gdal_trace_outline testcase.tif  -ndv 255 -out-cs xy -wkt-out test-tc1.wkt -split-polys -dp-toler 0
../gdal_trace_outline testcase2.tif -ndv 255 -out-cs xy -wkt-out test-tc2.wkt -split-polys -dp-toler 0
../gdal_trace_outline testcase3.tif -ndv 255 -out-cs xy -wkt-out test-tc3.wkt -split-polys -dp-toler 0
../gdal_trace_outline testcase4.png -ndv 255 -out-cs xy -wkt-out test-tc4.wkt -split-polys -dp-toler 0
../gdal_trace_outline testcase5.png -ndv 255 -out-cs xy -wkt-out test-tc5.wkt -split-polys -dp-toler 0

../gdal_list_corners -inspect-rect4 -erosion -ndv 0 testcase4.png > test-tc4-rect.wkt

#for i in 1 2 3 4 5 ; do md5sum good-tc$i.wkt test-tc$i.wkt ; done
for i in 1 2 3 4 5 ; do 
	if diff --brief good-tc$i.wkt test-tc$i.wkt ; then
		echo "outline $i - Good"
	else
		echo "outline $i - Bad"
	fi
done

if diff --brief good-tc4-rect.wkt test-tc4-rect.wkt ; then
	echo "rect - Good"
else
	echo "rect - Bad"
fi
