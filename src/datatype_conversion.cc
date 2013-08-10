#include <gdal.h>

#include "datatype_conversion.h"

namespace dangdal {

int32_t gdal_scalar_to_int32(void *p, GDALDataType dt) {
    int32_t ret;
    GDALCopyWords(p, dt, 0, &ret, GDT_Int32, 0, 1);
    return ret;
}

double gdal_scalar_to_double(void *p, GDALDataType dt) {
    double ret;
    GDALCopyWords(p, dt, 0, &ret, GDT_Float64, 0, 1);
    return ret;
}

template <typename T>
struct NaN_Checker {
	static bool isnan(const void *) { return false; }
};

template <>
struct NaN_Checker<double> {
	static bool isnan(const void *p) {
		return std::isnan(*reinterpret_cast<const double *>(p));
	}
};

template <>
struct NaN_Checker<float> {
	static bool isnan(const void *p) {
		return std::isnan(*reinterpret_cast<const float *>(p));
	}
};

template <>
struct NaN_Checker<std::complex<double> > {
	static bool isnan(const void *p) {
		return
			std::isnan(reinterpret_cast<const std::complex<double> *>(p)->real()) ||
			std::isnan(reinterpret_cast<const std::complex<double> *>(p)->imag());
	}
};

template <>
struct NaN_Checker<std::complex<float> > {
	static bool isnan(const void *p) {
		return
			std::isnan(reinterpret_cast<const std::complex<float> *>(p)->real()) ||
			std::isnan(reinterpret_cast<const std::complex<float> *>(p)->imag());
	}
};

template <typename T>
static bool nan_check(const void *p) {
	return NaN_Checker<T>::isnan(p);
}

bool gdal_scalar_pointer_isnan(const void *p, GDALDataType dt) {
	return DANGDAL_RUNTIME_TEMPLATE(dt, nan_check, p);
}

} // namespace dangdal
