#pragma once

#include "cephes/gamma.h"

namespace special {

template <typename T>
SPECFUN_HOST_DEVICE T gamma(T x) {
    return cephes::Gamma(x);
}

template <>
SPECFUN_HOST_DEVICE inline float gamma(float xf) {
    double x = xf;

    return gamma(x);
}

} // namespace special