#ifndef DANGDAL_DATATYPE_CONVERSION_H
#define DANGDAL_DATATYPE_CONVERSION_H

#include <stdint.h>

#include <complex>
#include <stdexcept>

namespace dangdal {

// This lets you convert a GDALDataType to a template argument.
#define DANGDAL_RUNTIME_TEMPLATE(dt, fn, ...) ( \
    (dt == GDT_Byte    ) ? fn< uint8_t>(__VA_ARGS__) : \
    (dt == GDT_UInt16  ) ? fn<uint16_t>(__VA_ARGS__) : \
    (dt == GDT_Int16   ) ? fn< int16_t>(__VA_ARGS__) : \
    (dt == GDT_UInt32  ) ? fn<uint32_t>(__VA_ARGS__) : \
    (dt == GDT_Int32   ) ? fn< int32_t>(__VA_ARGS__) : \
    (dt == GDT_Float32 ) ? fn<   float>(__VA_ARGS__) : \
    (dt == GDT_Float64 ) ? fn<  double>(__VA_ARGS__) : \
    (dt == GDT_CInt16  ) ? fn<std::complex<int16_t> >(__VA_ARGS__) : \
    (dt == GDT_CInt32  ) ? fn<std::complex<int32_t> >(__VA_ARGS__) : \
    (dt == GDT_CFloat32) ? fn<std::complex<  float> >(__VA_ARGS__) : \
    (dt == GDT_CFloat64) ? fn<std::complex< double> >(__VA_ARGS__) : \
    throw(std::invalid_argument("unrecognized datatype")) \
)

template <typename T>
struct GetGDALDataTypeFor { static const GDALDataType t = GDT_Unknown; };
template <> struct GetGDALDataTypeFor< uint8_t> { static const GDALDataType t = GDT_Byte   ; };
template <> struct GetGDALDataTypeFor<uint16_t> { static const GDALDataType t = GDT_UInt16 ; };
template <> struct GetGDALDataTypeFor< int16_t> { static const GDALDataType t = GDT_Int16  ; };
template <> struct GetGDALDataTypeFor<uint32_t> { static const GDALDataType t = GDT_UInt32 ; };
template <> struct GetGDALDataTypeFor< int32_t> { static const GDALDataType t = GDT_Int32  ; };
template <> struct GetGDALDataTypeFor<   float> { static const GDALDataType t = GDT_Float32; };
template <> struct GetGDALDataTypeFor<  double> { static const GDALDataType t = GDT_Float64; };
template <> struct GetGDALDataTypeFor<std::complex<int16_t> > { static const GDALDataType t = GDT_CInt16  ; };
template <> struct GetGDALDataTypeFor<std::complex<int32_t> > { static const GDALDataType t = GDT_CInt32  ; };
template <> struct GetGDALDataTypeFor<std::complex<  float> > { static const GDALDataType t = GDT_CFloat32; };
template <> struct GetGDALDataTypeFor<std::complex< double> > { static const GDALDataType t = GDT_CFloat64; };

// unfortunately, GDALCopyWords doesn't allow const input
int32_t gdal_scalar_to_int32(void *p, GDALDataType dt);
double gdal_scalar_to_double(void *p, GDALDataType dt);

bool gdal_scalar_pointer_isnan(const void *p, GDALDataType dt);

} // namespace dangdal

#endif // DANGDAL_DATATYPE_CONVERSION_H
