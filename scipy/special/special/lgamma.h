#pragma once

#include "cephes/gamma.h"

namespace special {

template <typename T>
SPECFUN_HOST_DEVICE T lgam(T x);

template <>
SPECFUN_HOST_DEVICE double lgam(double x) {
    return cephes::lgam(x);
}

template <>
SPECFUN_HOST_DEVICE float lgam(float x) {
    return static_cast<float>(lgam(static_cast<double>(x)));
}

} // namespace special