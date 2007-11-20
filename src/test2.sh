./gdal_list_corners testcase2.tif -ll_en 0 0 -res 1 -inspect-contour -nodataval 255 -dp-toler .3 -skip-erosion -wkt-en zz.wkt
( echo -n "select isvalid(mpolyfromtext('" ; cat zz.wkt ; echo -n "'));" ) |psql sv_ion
