extern "C" {
#include <numpy/npy_math.h>
#include "special_c_wrappers.h"
}

#include "special/binom.h"
#include "special/hyp2f1.h"
#include "special/cephes/airy.h"
#include "special/cephes/beta.h"
#include "special/cephes/bdtr.h"
#include "special/cephes/ellpj.h"
#include "special/cephes/ellpk.h"
#include "special/cephes/expn.h"
#include "special/cephes/hyp2f1.h"
#include "special/cephes/fresnl.h"
#include "special/cephes/gamma.h"
#include "special/cephes/scipy_iv.h"
#include "special/cephes/jv.h"
#include "special/cephes/kolmogorov.h"
#include "special/cephes/nbdtr.h"
#include "special/cephes/ndtr.h"
#include "special/cephes/ndtri.h"
#include "special/cephes/pdtr.h"
#include "special/cephes/poch.h"
#include "special/cephes/polevl.h"
#include "special/cephes/sici.h"
#include "special/cephes/shichi.h"
#include "special/cephes/struve.h"
#include "special/cephes/unity.h"
#include "special/cephes/yn.h"


extern "C" double binom_wrap(double n, double k) {
    return special::binom(n, k);
}

extern "C" npy_cdouble hyp2f1_complex_wrap(double a, double b, double c, npy_cdouble zp) {
    std::complex<double> z(npy_creal(zp), npy_cimag(zp));
    std::complex<double> w = special::hyp2f1(a, b, c, z);
    return npy_cpack(real(w), imag(w));
}

extern "C" double cephes_hyp2f1_wrap(double a, double b, double c, double x) {
    return special::cephes::hyp2f1(a, b, c, x);
}

extern "C" double cephes_airy_wrap(double x, double *ai, double *aip, double *bi, double *bip) {
    return special::cephes::airy(x, ai, aip, bi, bip);
}

extern "C" double cephes_beta_wrap(double a, double b) {
    return special::cephes::beta(a, b);
}

extern "C" double cephes_lbeta_wrap(double a, double b) {
    return special::cephes::lbeta(a, b);
}

extern "C" double cephes_bdtr_wrap(double k, int n, double p) {
    return special::cephes::bdtr(k, n, p);
}

extern "C" double cephes_bdtri_wrap(double k, int n, double y) {
    return special::cephes::bdtri(k, n, y);
}

extern "C" double cephes_bdtrc_wrap(double k, int n, double p) {
    return special::cephes::bdtrc(k, n, p);
}

extern "C" double cephes_cosm1_wrap(double x) {
    return special::cephes::cosm1(x);
}

extern "C" double cephes_expm1_wrap(double x) {
    return special::cephes::expm1(x);
}

extern "C" double cephes_expn_wrap(int n, double x) {
    return special::cephes::expn(n, x);
}

extern "C" double cephes_log1p_wrap(double x) {
    return special::cephes::log1p(x);
}

extern "C" double cephes_gamma_wrap(double x) {
    return special::cephes::Gamma(x);
}

extern "C" double cephes_gammasgn_wrap(double x) {
    return special::cephes::gammasgn(x);
}

extern "C" double cephes_lgam_wrap(double x) {
    return special::cephes::lgam(x);
}

extern "C" double cephes_iv_wrap(double v, double x) {
    return special::cephes::iv(v, x);
}

extern "C" double cephes_jv_wrap(double v, double x) {
    return special::cephes::jv(v, x);
}

extern "C" int cephes_ellpj_wrap(double u, double m, double *sn, double *cn, double *dn, double *ph) {
    return special::cephes::ellpj(u, m, sn, cn, dn, ph);
}

extern "C" double cephes_ellpk_wrap(double x) {
    return special::cephes::ellpk(x);
}

extern "C" int cephes_fresnl_wrap(double xxa, double *ssa, double *cca) {
    return special::cephes::fresnl(xxa, ssa, cca);
}

extern "C" double cephes_nbdtr_wrap(int k, int n, double p) {
    return special::cephes::nbdtr(k, n, p);
}

extern "C" double cephes_nbdtrc_wrap(int k, int n, double p) {
    return special::cephes::nbdtrc(k, n, p);
}

extern "C" double cephes_nbdtri_wrap(int k, int n, double p) {
    return special::cephes::nbdtri(k, n, p);
}

extern "C" double cephes_ndtr_wrap(double x) {
    return special::cephes::ndtr(x);
}

extern "C" double cephes_ndtri_wrap(double x) {
    return special::cephes::ndtri(x);
}

extern "C" double cephes_pdtri_wrap(int k, double y) {
    return special::cephes::pdtri(k, y);
}

extern "C" double cephes_poch_wrap(double x, double m) {
    return special::cephes::poch(x, m);
}

extern "C" int cephes_sici_wrap(double x, double *si, double *ci){
    return special::cephes::sici(x, si, ci);
}

extern "C" int cephes_shichi_wrap(double x, double *si, double *ci){
    return special::cephes::shichi(x, si, ci);
}

extern "C" double cephes_smirnov_wrap(int n, double x) {
    return special::cephes::smirnov(n, x);
}

extern "C" double cephes_smirnovc_wrap(int n, double x) {
    return special::cephes::smirnovc(n, x);
}

extern "C" double cephes_smirnovi_wrap(int n, double x) {
    return special::cephes::smirnovi(n, x);
}

extern "C" double cephes_smirnovci_wrap(int n, double x) {
    return special::cephes::smirnovci(n, x);
}

extern "C" double cephes_smirnovp_wrap(int n, double x) {
    return special::cephes::smirnovp(n, x);
}

extern "C" double cephes__struve_asymp_large_z(double v, double z, int is_h, double *err) {
    return special::cephes::detail::struve_asymp_large_z(v, z, is_h, err);
}

extern "C" double cephes__struve_bessel_series(double v, double z, int is_h, double *err) {
    return special::cephes::detail::struve_bessel_series(v, z, is_h, err);
}

extern "C" double cephes__struve_power_series(double v, double z, int is_h, double *err) {
    return special::cephes::detail::struve_power_series(v, z, is_h, err);
}

extern "C" double cephes_yn_wrap(int n, double x) {
    return special::cephes::yn(n, x);
}
    
extern "C" double cephes_polevl_wrap(double x, const double coef[], int N) {
    return special::cephes::polevl(x, coef, N);
}

extern "C" double cephes_p1evl_wrap(double x, const double coef[], int N) {
    return special::cephes::p1evl(x, coef, N);
}
