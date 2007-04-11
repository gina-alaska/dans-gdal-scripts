/*
    Copyright 2006 University of Alaska,
    Geographic Information Network of Alaska
    All Rights Reserved
    UAF Confidential Information
*/


#include <config.h>

#include <math.h>

#if HAVE_INTTYPES_H
#  include <inttypes.h>
#endif

#if HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#if HAVE_STRING_H
#  include <string.h>
#endif

#include <stdio.h>
#include <cpl_string.h>
#include <gdal.h>
#include <ogr_spatialref.h>

void fatal_error(char *s);
void *malloc_or_die(size_t size);
void *realloc_or_die(void *p, size_t size);
