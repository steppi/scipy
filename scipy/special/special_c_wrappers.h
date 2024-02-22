#pragma once


#ifdef __cplusplus
extern "C" {
#endif

#include <numpy/npy_math.h>

    npy_cdouble hyp2f1_complex_wrap(double a, double b, double c, npy_cdouble zp);
    double cephes_jv_wrap(double v, double x);
    double cephes_airy_wrap(double x, double *ai, double *aip, double *bi, double *bip);

#ifdef __cplusplus
}
#endif
