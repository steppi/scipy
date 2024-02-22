extern "C" {
#include <numpy/npy_math.h>
#include "special_c_wrappers.h"
}

#include "special/hyp2f1.h"
#include "special/cephes/airy.h"
#include "special/cephes/jv.h"

extern "C" npy_cdouble hyp2f1_complex_wrap(double a, double b, double c, npy_cdouble zp) {
    std::complex<double> z(npy_creal(zp), npy_cimag(zp));
    std::complex<double> w = special::hyp2f1(a, b, c, z);
    return npy_cpack(real(w), imag(w));
}


extern "C" double cephes_airy_wrap(double x, double *ai, double *aip, double *bi, double *bip) {
    return special::cephes::airy(x, ai, aip, bi, bip);
}

extern "C" double cephes_jv_wrap(double v, double x) {
        return special::cephes::jv(v, x);
}
