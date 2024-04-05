#pragma once


#ifdef __cplusplus
extern "C" {
#endif

#include <numpy/npy_math.h>

    double binom_wrap(double n, double k);
    npy_cdouble hyp2f1_complex_wrap(double a, double b, double c, npy_cdouble zp);
    double cephes_hyp2f1_wrap(double a, double b, double c, double x);
    double cephes_airy_wrap(double x, double *ai, double *aip, double *bi, double *bip);
    double cephes_beta_wrap(double a, double b);
    double cephes_lbeta_wrap(double a, double b);
    double cephes_cosm1_wrap(double x);
    double cephes_expm1_wrap(double x);
    double cephes_log1p_wrap(double x);
    double cephes_gamma_wrap(double x);
    double cephes_gammasgn_wrap(double x);
    double cephes_lgam_wrap(double x);
    double cephes_iv_wrap(double v, double x);
    double cephes_jv_wrap(double v, double x);
    double cephes_ellpk_wrap(double x);
    int cephes_ellpj_wrap(double u, double m, double *sn, double *cn, double *dn, double *ph);
    int cephes_fresnl_wrap(double xxa, double *ssa, double *cca);
    double cephes_ndtr_wrap(double x);
    double cephes_ndtri_wrap(double x);
    double cephes_poch_wrap(double x, double m);
    int cephes_sici_wrap(double x, double *si, double *ci);
    int cephes_shichi_wrap(double x, double *si, double *ci);
    double cephes__struve_asymp_large_z(double v, double z, int is_h, double *err);
    double cephes__struve_bessel_series(double v, double z, int is_h, double *err);
    double cephes__struve_power_series(double v, double z, int is_h, double *err);
    

#ifdef __cplusplus
}
#endif
