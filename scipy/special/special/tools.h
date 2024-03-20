/* Building blocks for implementing special functions */

#pragma once

#include "config.h"
#include "error.h"

namespace special {
namespace detail {

    namespace maybe_complex_numeric_limits {
        // Handle numeric limits when type may be complex.
        template <typename T>
        SPECFUN_HOST_DEVICE inline std::enable_if_t<std::is_floating_point_v<T>, T> quiet_NaN() {
            return std::numeric_limits<T>::quiet_NaN();
        }

        template <typename T>
        SPECFUN_HOST_DEVICE inline std::enable_if_t<!std::is_floating_point_v<T>, T> quiet_NaN() {
            using V = typename T::value_type;
            return {std::numeric_limits<V>::quiet_NaN(), std::numeric_limits<V>::quiet_NaN()};
        }

        template <typename T>
        SPECFUN_HOST_DEVICE inline std::enable_if_t<std::is_floating_point_v<T>, T> min() {
            return std::numeric_limits<T>::min();
        }

        template <typename T>
        SPECFUN_HOST_DEVICE inline std::enable_if_t<!std::is_floating_point_v<T>, T> min() {
            using V = typename T::value_type;
            return std::numeric_limits<V>::min();
        }

    } // namespace maybe_complex_numeric_limits

    template <typename Generator, typename T>
    SPECFUN_HOST_DEVICE T continued_fraction_eval(Generator &g, T init_val, double tol, std::uint64_t max_terms,
                                                  const char *func_name) {
        /* Evaluate continued fraction with modified Lentz's algorithm.
         *
         * Evaluates the continued fraction b0 + a1 / (b1 + a1 / (b2 + a2 / ( ...
         *
         * g : Generator of pairs of values (a1, b1), (a2, b2), (a3, b3), ...
         *
         * init_val : Initial value b0
         *
         * tol : relative tolerance for stopping criterion.
         *
         * max_terms : The maximum number of iterations before giving up and declaring
         *     non-convergence.
         *
         * func_name : The name of the function within SciPy where this call to series_eval
         *     will ultimately be used. This is needed to pass to set_error in case
         *     of non-convergence.
         */
        std::pair<T, T> v;

        T tiny_value = 16.0 * maybe_complex_numeric_limits::min<T>();
        T f = (init_val == 0.0) ? tiny_value : init_val;

        double C = f;
        double D = 0.0;
        double delta;

        for (uint64_t i = 0; i < max_terms; i++) {
            v = g();
            D = v.second + v.first * D;
            if (D == 0.0) {
                D = tiny_value;
            }
            C = v.second + v.first / C;
            if (C == 0.0) {
                C = tiny_value;
            }
            D = 1.0 / D;
            delta = C * D;
            f *= delta;
            if (std::abs(delta - 1.0) <= tol) {
                return f;
            }
        }
        // Exceeded max terms without converging. Return NaN.
        set_error(func_name, SF_ERROR_NO_RESULT, NULL);
        return maybe_complex_numeric_limits::quiet_NaN<T>();
    }

} // namespace detail
} // namespace special
