#pragma once

#include "cephes/beta.h"

namespace special {

template <typename T>
SPECFUN_HOST_DEVICE T beta(T a, T b) {
    return cephes::beta(a, b);
}

template <>
SPECFUN_HOST_DEVICE inline float beta(float af, float bf) {
    double a = af;
    double b = bf;

    return beta(a, b);
}

template <typename T>
SPECFUN_HOST_DEVICE T lbeta(T a, T b) {
    return cephes::lbeta(a, b);
}

template <>
SPECFUN_HOST_DEVICE inline float lbeta(float af, float bf) {
    double a = af;
    double b = bf;

    return lbeta(a, b);
}

} // namespace special