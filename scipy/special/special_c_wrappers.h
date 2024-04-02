#pragma once


#ifdef __cplusplus
extern "C" {
#endif

#include <numpy/npy_math.h>

    npy_cdouble hyp2f1_complex_wrap(double a, double b, double c, npy_cdouble zp);
    double cephes_airy_wrap(double x, double *ai, double *aip, double *bi, double *bip);
    double cephes_jv_wrap(double v, double x);
    int cephes_ellpj_wrap(double u, double m, double *sn, double *cn, double *dn, double *ph);
    int cephes_fresnl_wrap(double xxa, double *ssa, double *cca);
    int cephes_sici_wrap(double x, double *si, double *ci);

#ifdef __cplusplus
}
#endif
