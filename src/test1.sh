#./gdal_trace_outline testcase.tif  -nodataval 255 -out-cs xy -wkt-out good-tc1.wkt -split-polys
#./gdal_trace_outline testcase2.tif -nodataval 255 -out-cs xy -wkt-out good-tc2.wkt -split-polys
#./gdal_trace_outline testcase3.tif -nodataval 255 -out-cs xy -wkt-out good-tc3.wkt -split-polys
#./gdal_trace_outline testcase4.tif -nodataval 255 -out-cs xy -wkt-out good-tc4.wkt -split-polys
#./gdal_trace_outline testcase5.png -nodataval 255 -out-cs xy -wkt-out good-tc5.wkt -split-polys
rm -f test-tc*.wkt
./gdal_trace_outline testcase.tif  -nodataval 255 -out-cs xy -wkt-out test-tc1.wkt -split-polys
./gdal_trace_outline testcase2.tif -nodataval 255 -out-cs xy -wkt-out test-tc2.wkt -split-polys
./gdal_trace_outline testcase3.tif -nodataval 255 -out-cs xy -wkt-out test-tc3.wkt -split-polys
./gdal_trace_outline testcase4.tif -nodataval 255 -out-cs xy -wkt-out test-tc4.wkt -split-polys
./gdal_trace_outline testcase5.png -nodataval 255 -out-cs xy -wkt-out test-tc5.wkt -split-polys
for i in 1 2 3 5 ; do md5sum good-tc$i.wkt test-tc$i.wkt ; done
