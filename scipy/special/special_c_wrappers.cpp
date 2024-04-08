/* These wrappers exist to allow for calling functions from the under development
 * header only C++ special function library in scipy/special/special within files
 * for special functions which are still implemented in C or Cython. The hope is
 * that these  wrappers will become unnecessary if all of special is translated
 * into C++ */

#include <numpy/npy_math.h>

#include "special_c_wrappers.h"

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


double binom_wrap(double n, double k) {
    return special::binom(n, k);
}

npy_cdouble hyp2f1_complex_wrap(double a, double b, double c, npy_cdouble zp) {
    std::complex<double> z(npy_creal(zp), npy_cimag(zp));
    std::complex<double> w = special::hyp2f1(a, b, c, z);
    return npy_cpack(real(w), imag(w));
}

double cephes_hyp2f1_wrap(double a, double b, double c, double x) {
    return special::cephes::hyp2f1(a, b, c, x);
}

double cephes_airy_wrap(double x, double *ai, double *aip, double *bi, double *bip) {
    return special::cephes::airy(x, ai, aip, bi, bip);
}

double cephes_beta_wrap(double a, double b) {
    return special::cephes::beta(a, b);
}

double cephes_lbeta_wrap(double a, double b) {
    return special::cephes::lbeta(a, b);
}

double cephes_bdtr_wrap(double k, int n, double p) {
    return special::cephes::bdtr(k, n, p);
}

double cephes_bdtri_wrap(double k, int n, double y) {
    return special::cephes::bdtri(k, n, y);
}

double cephes_bdtrc_wrap(double k, int n, double p) {
    return special::cephes::bdtrc(k, n, p);
}

double cephes_cosm1_wrap(double x) {
    return special::cephes::cosm1(x);
}

double cephes_expm1_wrap(double x) {
    return special::cephes::expm1(x);
}

double cephes_expn_wrap(int n, double x) {
    return special::cephes::expn(n, x);
}

double cephes_log1p_wrap(double x) {
    return special::cephes::log1p(x);
}

double cephes_gamma_wrap(double x) {
    return special::cephes::Gamma(x);
}

double cephes_gammasgn_wrap(double x) {
    return special::cephes::gammasgn(x);
}

double cephes_lgam_wrap(double x) {
    return special::cephes::lgam(x);
}

double cephes_iv_wrap(double v, double x) {
    return special::cephes::iv(v, x);
}

double cephes_jv_wrap(double v, double x) {
    return special::cephes::jv(v, x);
}

int cephes_ellpj_wrap(double u, double m, double *sn, double *cn, double *dn, double *ph) {
    return special::cephes::ellpj(u, m, sn, cn, dn, ph);
}

double cephes_ellpk_wrap(double x) {
    return special::cephes::ellpk(x);
}

int cephes_fresnl_wrap(double xxa, double *ssa, double *cca) {
    return special::cephes::fresnl(xxa, ssa, cca);
}

double cephes_nbdtr_wrap(int k, int n, double p) {
    return special::cephes::nbdtr(k, n, p);
}

double cephes_nbdtrc_wrap(int k, int n, double p) {
    return special::cephes::nbdtrc(k, n, p);
}

double cephes_nbdtri_wrap(int k, int n, double p) {
    return special::cephes::nbdtri(k, n, p);
}

double cephes_ndtr_wrap(double x) {
    return special::cephes::ndtr(x);
}

double cephes_ndtri_wrap(double x) {
    return special::cephes::ndtri(x);
}

double cephes_pdtri_wrap(int k, double y) {
    return special::cephes::pdtri(k, y);
}

double cephes_poch_wrap(double x, double m) {
    return special::cephes::poch(x, m);
}

int cephes_sici_wrap(double x, double *si, double *ci){
    return special::cephes::sici(x, si, ci);
}

int cephes_shichi_wrap(double x, double *si, double *ci){
    return special::cephes::shichi(x, si, ci);
}

double cephes_smirnov_wrap(int n, double x) {
    return special::cephes::smirnov(n, x);
}

double cephes_smirnovc_wrap(int n, double x) {
    return special::cephes::smirnovc(n, x);
}

double cephes_smirnovi_wrap(int n, double x) {
    return special::cephes::smirnovi(n, x);
}

double cephes_smirnovci_wrap(int n, double x) {
    return special::cephes::smirnovci(n, x);
}

double cephes_smirnovp_wrap(int n, double x) {
    return special::cephes::smirnovp(n, x);
}

double cephes__struve_asymp_large_z(double v, double z, int is_h, double *err) {
    return special::cephes::detail::struve_asymp_large_z(v, z, is_h, err);
}

double cephes__struve_bessel_series(double v, double z, int is_h, double *err) {
    return special::cephes::detail::struve_bessel_series(v, z, is_h, err);
}

double cephes__struve_power_series(double v, double z, int is_h, double *err) {
    return special::cephes::detail::struve_power_series(v, z, is_h, err);
}

double cephes_yn_wrap(int n, double x) {
    return special::cephes::yn(n, x);
}
    
double cephes_polevl_wrap(double x, const double coef[], int N) {
    return special::cephes::polevl(x, coef, N);
}

double cephes_p1evl_wrap(double x, const double coef[], int N) {
    return special::cephes::p1evl(x, coef, N);
}
