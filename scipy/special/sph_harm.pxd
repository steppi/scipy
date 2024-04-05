# cython: cpow=True

from . cimport sf_error

cdef extern from "special_c_wrappers.h" nogil:
    double cephes_poch_wrap(double x, double m)

cdef extern from "special_wrappers.h":
    double pmv_wrap(double, double, double) nogil

from ._complexstuff cimport *
from libc.math cimport cos, sqrt, fabs, NAN, M_PI
from libc.stdlib cimport abs

cdef inline double complex sph_harmonic(int m, int n, double theta, double phi) noexcept nogil:
    cdef double x, prefactor
    cdef double complex val
    cdef int mp
    x = cos(phi)
    if abs(m) > n :
        sf_error.error("sph_harm", sf_error.ARG, "m should not be greater than n")
        return NAN
    if n < 0:
        sf_error.error("sph_harm", sf_error.ARG, "n should not be negative")
        return NAN
    if m < 0:
        mp = -m
        prefactor = (-1)**mp * cephes_poch_wrap(n + mp + 1, -2 * mp)
    else:
        mp = m
    val = pmv_wrap(mp, n, x)
    if  m < 0:
        val *= prefactor
    val *= sqrt((2*n + 1) / 4.0 / M_PI)
    val *= sqrt(cephes_poch_wrap(n + m + 1, -2 * m))
    val *= zexp(1j * m * theta)
    return val
