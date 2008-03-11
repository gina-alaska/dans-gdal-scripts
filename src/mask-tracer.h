#ifndef MASK_TRACER_H

#include "polygon.h"

mpoly_t trace_mask(unsigned char *mask_1bit, int w, int h, long min_area, int no_donuts);

#endif // ifndef MASK_TRACER_H
