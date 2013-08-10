#!/bin/sh

rm -f out_test1_*

#BINDIR="valgrind -q .."
BINDIR=..

$BINDIR/gdal_trace_outline testcase_1.tif -ndv 255 -out-cs xy -wkt-out out_test1_1.wkt    -report out_test1_1.ppm -split-polys -dp-toler 0
$BINDIR/gdal_trace_outline testcase_1.tif -ndv 255 -out-cs en -wkt-out out_test1_1_en.wkt -dp-toler 0
$BINDIR/gdal_trace_outline testcase_2.tif -ndv 255 -out-cs xy -wkt-out out_test1_2.wkt    -report out_test1_2.ppm -split-polys -dp-toler 0
$BINDIR/gdal_trace_outline testcase_3.tif -ndv 255 -out-cs xy -wkt-out out_test1_3.wkt    -report out_test1_3.ppm -split-polys -dp-toler 0
$BINDIR/gdal_trace_outline testcase_4.png -ndv '0..255 0..255 0..255 0' -out-cs xy -wkt-out out_test1_4.wkt    -report out_test1_4.ppm -split-polys -dp-toler 0
$BINDIR/gdal_trace_outline testcase_5.png -ndv 255 -out-cs xy -wkt-out out_test1_5.wkt    -report out_test1_5.ppm -split-polys -dp-toler 0
$BINDIR/gdal_trace_outline testcase_maze.png  -ndv 255 -out-cs xy -wkt-out out_test1_maze.wkt  -report out_test1_maze.ppm  -split-polys -dp-toler 0
$BINDIR/gdal_trace_outline testcase_noise.png -b 1 -ndv   0 -out-cs xy -wkt-out out_test1_noise.wkt -report out_test1_noise.ppm -split-polys -dp-toler 0
$BINDIR/gdal_trace_outline testcase_noise.png -b 1 -ndv   0 -out-cs xy -wkt-out out_test1_noise_dp3.wkt -report out_test1_noise_dp3.ppm -split-polys -dp-toler 3

$BINDIR/gdal_trace_outline testcase_3.tif -out-cs xy -wkt-out out_test1_3_classify.wkt -ogr-out out_test1_3_classify.shp -dp-toler 0 -classify
$BINDIR/gdal_trace_outline testcase_3_paletted.tif -out-cs xy -wkt-out out_test1_3_classify_pal.wkt -ogr-out out_test1_3_classify_pal.shp -dp-toler 0 -classify
$BINDIR/gdal_trace_outline testcase_features.png -out-cs xy -wkt-out out_test1_features.wkt -ogr-out out_test1_features.shp -dp-toler 0 -classify
$BINDIR/gdal_trace_outline testcase_double.tif -out-cs xy -wkt-out out_test1_double.wkt -ogr-out out_test1_double.shp -dp-toler 0 -classify
$BINDIR/gdal_trace_outline testcase_double.tif -out-cs xy -wkt-out out_test1_double_clip.wkt -dp-toler 0 -classify -valid-range '3..6'

$BINDIR/gdal_list_corners -inspect-rect4 -erosion -ndv 0 testcase_4.png -report out_test1_4-rect.ppm > out_test1_4-rect.wkt

$BINDIR/gdal_wkt_to_mask -wkt good_test1_1_en.wkt -geo-from testcase_1.tif -mask-out out_test1_1_mask.ppm

$BINDIR/gdal_get_projected_bounds -s_wkt good_test1_1_en.wkt -s_srs '+proj=utm +zone=6 +ellps=WGS84 +units=m +no_defs ' -t_srs '+proj=stere +lat_ts=80 +lat_0=90 +lon_0=0 +ellps=WGS84' -report out_test1_projbounds_report.ppm > out_test1_projbounds.yml

$BINDIR/gdal_list_corners testcase_1.tif > out_test1_1_metdata.yml

$BINDIR/gdal_make_ndv_mask -ndv '155 52 52' -ndv '24 173 79'     testcase_3.tif out_test1_3_ndvmask.pbm
$BINDIR/gdal_make_ndv_mask -ndv '155 52 52' -ndv '24 173 79.9..80.1' testcase_3.tif out_test1_3_ndvmask2.pbm

# Make a gradient image.
python <<END
import numpy as np
import osgeo.gdal as gdal

driver = gdal.GetDriverByName('GTiff')

(w, h) = (256, 256)

x = np.outer(np.ones(h), range(w))
y = np.outer(range(h), np.ones(w))

dst_ds = driver.Create('gradient1.tif', w, h, 2, gdal.GDT_Byte)
dst_ds.GetRasterBand(1).WriteArray(x, 0, 0)
dst_ds.GetRasterBand(2).WriteArray(y, 0, 0)
dst_ds = None

dst_ds = driver.Create('gradient2.tif', w, h, 1, gdal.GDT_Float64)
dst_ds.GetRasterBand(1).WriteArray((y+1)/(x+1), 0, 0)
dst_ds = None
END
gdal_merge_vrt -in gradient1.tif -in gradient2.tif -out gradient3.tif

$BINDIR/gdal_make_ndv_mask \
    -ndv '10..30 30..70 *' \
    -ndv '* * 4..Inf' \
    -ndv '100..140 50..80 0.5..Inf' \
    -ndv '* * 0.8..0.3' \
    -ndv '* * 0.3..0.4' \
    gradient3.tif out_test1_gradient_ndv.pbm

$BINDIR/gdal_make_ndv_mask \
    -ndv '10..30 30..70 *' \
    -ndv '* * 4..Inf' \
    -ndv '100..140 50..80 0.5..Inf' \
    -ndv '* * 0.8..0.3' \
    -ndv '* * 0.3..0.4' \
	-invert \
    gradient3.tif out_test1_gradient_ndv_inv.pbm

$BINDIR/gdal_make_ndv_mask \
    -valid-range '10..30 30..70 *' \
    -valid-range '* * 4..Inf' \
    -valid-range '100..140 50..80 0.5..Inf' \
    -valid-range '* * 0.8..0.3' \
    -valid-range '* * 0.3..0.4' \
    gradient3.tif out_test1_gradient_valid.pbm

python <<END
import numpy as np
import osgeo.gdal as gdal

driver = gdal.GetDriverByName('GTiff')
(w, h) = (256, 256)
dst_ds = driver.Create('has_nan.tif', w, h, 2, gdal.GDT_Float64)

vals = np.zeros((w, h))
vals[30:40,100:110] = np.nan
vals[60:70,100:110] = np.inf
dst_ds.GetRasterBand(1).WriteArray(vals, 0, 0)

vals = np.zeros((w, h))
vals[30:40,150:160] = np.nan
vals[60:70,150:160] = -np.inf
dst_ds.GetRasterBand(2).WriteArray(vals, 0, 0)

dst_ds = None
END

# NDV is specified here just to make the script run.  But the point is to see that NaN values
# (in either band) are considered to be NDV.
$BINDIR/gdal_make_ndv_mask has_nan.tif out_test1_nan1.pbm -ndv 7
# Check detection of infinity
$BINDIR/gdal_make_ndv_mask has_nan.tif out_test1_nan2.pbm -ndv 'Inf *'
$BINDIR/gdal_make_ndv_mask has_nan.tif out_test1_nan3.pbm -ndv '* Inf'
$BINDIR/gdal_make_ndv_mask has_nan.tif out_test1_nan4.pbm -ndv '* -Inf'

echo '####################'

for i in out_test1_* ; do
	if [ ! -e ${i/out/good} ] ; then
		echo "!!! ${i/out/good} doesn't exist"
	fi
done

for i in good_test1_* ; do
	if diff --brief $i ${i/good/out} ; then
		echo "GOOD ${i/good_/}"
	else
		echo "BAD ${i/good_/}"
	fi
done
