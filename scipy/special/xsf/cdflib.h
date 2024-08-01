#include "config.h"
#include "error.h"

#include "tools.h"
#include "cephes/gdtr.h"


namespace special {

    SPECFUN_HOST_DEVICE inline double gdtria(double p, double b, double x) {
	if (std::isnan(p) || std::isnan(b) || std::isnan(x)) {
	    return std::numeric_limits<double>::quiet_NaN();
	}
	if (!((0 <= p) && (p <= 1))) {
	    set_error("gdtria", SF_ERROR_ARG, "Input parameter p is out of range");
	    return std::numeric_limits<double>::quiet_NaN();
	}
	if (!(b > 0)) {
	    set_error("gdtria", SF_ERROR_ARG, "Input parameter b is out of range");
	    return std::numeric_limits<double>::quiet_NaN();
	}
	if (!(x >= 0)) {
	    set_error("gdtria", SF_ERROR_ARG, "Input parameter x is out of range");
	    return std::numeric_limits<double>::quiet_NaN();
	}
	auto func = [p, b, x](double a) {
	    return cephes::gdtr(a, b, x) - p;
	};
	double result = detail::find_root_bus_dekker_r(func, 1e-100, 1e100, 10000);
	return result;
    }

    SPECFUN_HOST_DEVICE inline double gdtrib(double a, double p, double x) {
	if (std::isnan(p) || std::isnan(a) || std::isnan(x)) {
	    return std::numeric_limits<double>::quiet_NaN();
	}
	if (!((0 <= p) && (p <= 1))) {
	    set_error("gdtrib", SF_ERROR_ARG, "Input parameter p is out of range");
	    return std::numeric_limits<double>::quiet_NaN();
	}
	if (!(a > 0)) {
	    set_error("gdtrib", SF_ERROR_ARG, "Input parameter a is out of range");
	    return std::numeric_limits<double>::quiet_NaN();
	}
	if (!(x >= 0)) {
	    set_error("gdtrib", SF_ERROR_ARG, "Input parameter x is out of range");
	    return std::numeric_limits<double>::quiet_NaN();
	}
	auto func = [p, a, x](double b) {
	    return cephes::gdtr(a, b, x) - p;
	};
	double result = detail::find_root_bus_dekker_r(func, 1e-100, 1e100, 10000);
	return result;
    }

    SPECFUN_HOST_DEVICE inline double gdtrix(double a, double b, double p) {
	if (std::isnan(p) || std::isnan(a) || std::isnan(b)) {
	    return std::numeric_limits<double>::quiet_NaN();
	}
	if (!((0 <= p) && (p <= 1))) {
	    set_error("gdtria", SF_ERROR_ARG, "Input parameter p is out of range");
	    return std::numeric_limits<double>::quiet_NaN();
	}
	if (!(a > 0)) {
	    set_error("gdtrix", SF_ERROR_ARG, "Input parameter a is out of range");
	    return std::numeric_limits<double>::quiet_NaN();
	}
	if (!(b > 0)) {
	    set_error("gdtrix", SF_ERROR_ARG, "Input parameter b is out of range");
	    return std::numeric_limits<double>::quiet_NaN();
	}
	auto func = [p, a, b](double x) {
	    return cephes::gdtr(a, b, x) - p;
	};
	double result = detail::find_root_bus_dekker_r(func, 0, 1e100, 10000);
	return result;
    }

}
