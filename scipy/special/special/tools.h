/* Building blocks for implementing special functions */

#pragma once

#include "config.h"
#include "error.h"

namespace special {
namespace detail {

    // Series evaluators.
    template <typename Generator, typename Number>
    SPECFUN_HOST_DEVICE Number series_eval(Generator &g, Number init_val, double tol, std::uint64_t max_terms,
                                           const char *func_name) {
        /* Sum an infinite series to a given precision.
         *
         * g : a generator of terms for the series.
         *
         * init_val : A starting value that terms are added to. This argument determines the
         *     type of the result.
         *
         * tol : relative tolerance for stopping criterion.
         *
         * max_terms : The maximum number of terms to add before giving up and declaring
         *     non-convergence.
         *
         * func_name : The name of the function within SciPy where this call to series_eval
         *     will ultimately be used. This is needed to pass to set_error in case
         *     of non-convergence.
         */
        Number result = init_val;
        Number term;
        for (std::uint64_t i = 0; i < max_terms; ++i) {
            term = g();
            result += term;
            if (std::abs(term) < std::abs(result) * tol) {
                return result;
            }
        }
        // Exceeded max terms without converging. Return NaN.
        set_error(func_name, SF_ERROR_NO_RESULT, NULL);
        if constexpr (std::is_floating_point<Number>::value) {
            return std::numeric_limits<Number>::quiet_NaN();
        } else {
            // If result type is not a floating point type, assume it is complex.
            using FloatingType = typename Number::value_type;
            return Number(std::numeric_limits<FloatingType>::quiet_NaN(), std::numeric_limits<FloatingType>::quiet_NaN());
        }
    }


    template <typename Generator, typename Number>
    SPECFUN_HOST_DEVICE Number series_eval_custom_stop(Generator &g, Number init_val, double tol, std::uint64_t max_terms,
                                           const char *func_name) {
        /* Sum an infinite series to a given precision.
         *
         * g : a generator of pairs (xi, ci). The xi are the terms to be summed, and the ci give stopping criterion.
         *     Convergence is declared and the summation is terminated when |ci| < tol * |result|. Where result is
         *     the current partial sum.
         *
         * init_val : A starting value that terms are added to. This argument determines the
         *     type of the result.
         *
         * tol : Terminate when |ci| < tol * |result|.
         *
         * max_terms : The maximum number of terms to add before giving up and declaring
         *     non-convergence.
         *
         * func_name : The name of the function within SciPy where this call to series_eval
         *     will ultimately be used. This is needed to pass to set_error in case
         *     of non-convergence.
         */
        Number result = init_val;
     
        std::pair<Number, Number> v;
        for (std::uint64_t i = 0; i < max_terms; ++i) {
            v = g();
            result += v.first;
            if (std::abs(v.second) < tol * std::abs(result)) {
                return result;
            }
        }
        // Exceeded max terms without converging. Return NaN.
        set_error(func_name, SF_ERROR_NO_RESULT, NULL);
        if constexpr (std::is_floating_point<Number>::value) {
            return std::numeric_limits<Number>::quiet_NaN();
        } else {
            // If result type is not a floating point type, assume it is complex.
            using FloatingType = typename Number::value_type;
            return Number(std::numeric_limits<FloatingType>::quiet_NaN(), std::numeric_limits<FloatingType>::quiet_NaN());
        }
    }

    template <typename Generator, typename Number>
    SPECFUN_HOST_DEVICE Number series_eval_fixed_length(Generator &g, Number init_val, std::uint64_t num_terms) {
        /* Sum a fixed number of terms from a series.
         *
         * g : a generator of terms for the series.
         *
         * init_val : A starting value that terms are added to. This argument determines the
         *     type of the result.
         *
         * max_terms : The number of terms from the series to sum.
         *
         */
        Number result = init_val;
        for (std::uint64_t i = 0; i < num_terms; ++i) {
            result += g();
        }
        return result;
    }

    template <typename Generator, typename Number>
    SPECFUN_HOST_DEVICE Number continued_fraction_eval(Generator &g, Number init_val, double tol, std::uint64_t max_terms,
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
         * max_terms : The maximum number of iterations before giving up and declaring non-convergence.
         *
         * func_name : The name of the function within SciPy where this call to series_eval
         *     will ultimately be used. This is needed to pass to set_error in case
         *     of non-convergence.
         */
        std::pair<Number, Number> v;

        Number tiny_value;

        if constexpr (std::is_floating_point<Number>::value) {
            tiny_value = 16.0 * std::numeric_limits<Number>::min();
        } else {
            // If not a floating point type, assume it is complex.
            using FloatingType = typename Number::value_type;
            tiny_value = 16.0 * std::numeric_limits<FloatingType>::min();
        }

        Number f = (f == 0.0) ? tiny_value : init_val;

        double C = f;
        double D = 0.0;
        double delta;

        for (uint64_t i=0; i < max_terms; i++) {
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
        if constexpr (std::is_floating_point<Number>::value) {
            return std::numeric_limits<Number>::quiet_NaN();
        } else {
            // If result type is not a floating point type, assume it is complex.
            using FloatingType = typename Number::value_type;
            return Number(std::numeric_limits<FloatingType>::quiet_NaN(), std::numeric_limits<FloatingType>::quiet_NaN());
        }
    }

} // namespace detail
} // namespace special
