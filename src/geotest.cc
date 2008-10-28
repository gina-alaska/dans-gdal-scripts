#include "common.h"
#include "polygon.h"
#include "georef.h"

void usage(const char *cmdname) {
	printf("Usage:\n  %s [options] [image_name]\n", cmdname);
	printf("\n");
	print_georef_usage();
}

int main(int argc, char **argv) {
	geo_opts_t geo_opts = init_geo_options(&argc, &argv);

	GDALAllRegister();

	if(argc < 2) usage(argv[0]);
	char *fn = argv[1];
	GDALDatasetH ds = GDALOpen(fn, GA_ReadOnly);
	if(!ds) fatal_error("open failed");

	georef_t georef = init_georef(&geo_opts, ds);
	
	double x=1000, y=10;
	double east, north;
	double lon, lat;
	printf("xy = %g, %g\n", x, y);
	xy2en(&georef, x, y, &east, &north);
	printf("en = %g, %g\n", east, north);
	en2ll_or_die(&georef, east, north, &lon, &lat);
	printf("ll = %g, %g\n", lon, lat);
	ll2en_or_die(&georef, lon, lat, &east, &north);
	printf("en = %g, %g\n", east, north);
	en2xy(&georef, east, north, &x, &y);
	printf("xy = %g, %g\n", x, y);
}
