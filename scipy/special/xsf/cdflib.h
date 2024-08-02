#pragma once

#include "config.h"
#include "error.h"

#include "tools.h"
#include "trig.h"
#include "digamma.h"
#include "cephes/const.h"
#include "cephes/gamma.h"
#include "cephes/gdtr.h"
#include "cephes/unity.h"


namespace xsf {
namespace cdflib {

    XSF_HOST_DEVICE inline double devlpl(const double a[], int n, double x) {
	//            Double precision EVALuate a PoLynomial at X
	//
	//
	//                            Function
	//
	//
	//    returns
	//        A(1) + A(2)*X + ... + A(N)*X**(N-1)
	//
	//
	//                            Arguments
	//
	//
	//    A --> Array of coefficients of the polynomial.
	//                                    A is DOUBLE PRECISION(N)
	//
	//    N --> Length of A, also degree of polynomial - 1.
	//                                    N is INTEGER
	//
	//    X --> Point at which the polynomial is to be evaluated.
	//                                    X is DOUBLE PRECISION
	double temp = a[n-1];
	int i;
	
	for (i = n - 2; i >= 0; i--) {
	    temp = a[i] + temp*x;
	}
	return temp;
    }

    constexpr double algdiv_carr[6] = {
	0.833333333333333e-01, -0.277777777760991e-02,
	0.793650666825390e-03, -0.595202931351870e-03,
	0.837308034031215e-03, -0.165322962780713e-02
    };

    XSF_HOST_DEVICE inline double algdiv(double a, double b) {
	//         Computation of ln(gamma(b)/gamma(a+b)) when b >= 8
	//
	//                             --------
	//
	//         In this algorithm, del(x) is the function defined by
	//         Ln(gamma(x)) = (x - 0.5)*ln(x) - x + 0.5*ln(2*pi) + del(x).
	//

	double c, d, h, s11, s3, s5, s7, s9, t, u, v, w, x, x2;

	if (a > b) {
	    h = b / a;
	    c = 1./(1. + h);
	    x = h/(1. + h);
	    d = a + (b - 0.5);
	} else {
	    h = a / b;
	    c = h/(1. + h);
	    x = 1./(1. + h);
	    d = b + (a - 0.5);
	}
	// Set sn = (1 - x**n)/(1 - x)
	x2 = x*x;
	s3 = 1. + (x + x2);
	s5 = 1. + (x + x2*s3);
	s7 = 1. + (x + x2*s5);
	s9 = 1. + (x + x2*s7);
	s11 = 1. + (x + x2*s9);

	// Set w = del(b) - del(a + b)
	t = std::pow((1. / b), 2);
	w = (((((algdiv_carr[5]*s11
		 )*t + algdiv_carr[4]*s9
		)*t + algdiv_carr[3]*s7
	       )*t + algdiv_carr[2]*s5
	      )*t + algdiv_carr[1]*s3
	     )*t + algdiv_carr[0];
	w *= c / b;
	// Combine the results
	u = d * std::log1p(a / b);
	v = a * (std::log(b) - 1.);
	return (u > v ? (w - v) - u : (w - u) - v);
    }
    
    constexpr double alngam_scoefn[9] = {
	0.62003838007127258804e2, 0.36036772530024836321e2,
	0.20782472531792126786e2, 0.6338067999387272343e1,
	0.215994312846059073e1, 0.3980671310203570498e0,
	0.1093115956710439502e0, 0.92381945590275995e-2,
	0.29737866448101651e-2
    };
    constexpr double alngam_scoefd[4] = {
	0.62003838007126989331e2, 0.9822521104713994894e1,
	-0.8906016659497461257e1, 0.1000000000000000000e1
    };
    constexpr double alngam_coef[5] = {
	0.83333333333333023564e-1, -0.27777777768818808e-2,
	0.79365006754279e-3, -0.594997310889e-3, 0.8065880899e-3
    };

    XSF_HOST_DEVICE inline double alngam(double x) {
	/* Included only for benchmarking. it's not clear why cdflib has
	* two gammaln implementations. */

	//                   Double precision log of the gamma function
	//
	//
	//                                  Function
	//
	//
	//         Returns the natural logarithm of gamma(x).
	//
	//
	//                                  Arguments
	//
	//
	//         X --> value at which scaled log gamma is to be returned
	//                        x is double precision
	//
	//
	//                                  Method
	//
	//
	//         If x <= 6.0, then use recursion to get x below 3
	//         then apply rational approximation number 5236 of
	//         hart et al, computer approximations, john wiley and
	//         sons, ny, 1968.
	//
	//         If x > 6.0, then use recursion to get x to at least 12 and
	//         then use formula 5423 of the same source.
	
	double prod, xx, result, offset;
	int i, n;

	if (x <= 6.0) {
	    prod = 1.0;
	    xx = x;

	    if (x > 3.0) {
		while (xx > 3.0) {
		    xx -= 1.0;
		    prod *= xx;
		}
	    }
	    if (x < 2.0) {
		while (xx < 2.0) {
		    prod /= xx;
		    xx += 1.0;
		}
	    }
	    result = devlpl(alngam_scoefn, 9, xx - 2.) / devlpl(alngam_scoefd, 4, xx - 2.);
	    // Compute rational approximation to gamma(x)
	    return std::log(result * prod);
	}
	offset = 0.5*std::log(2.*M_PI);
	// If necessary make x at least 12 and carry correction in offset
	if (x <= 12.0) {
	    n = static_cast<int>(12. - x);
	    if (n > 0) {
		prod = 1.0;
		for (i = 0; i < n; i++) {
		    prod *= x + i;
		}
		offset -= std::log(prod);
		xx = x + n;
	    } else {
		xx = x;
	    }
	} else {
	    xx = x;
	}
	// Compute power series
	result = devlpl(alngam_coef, 5, std::pow((1./xx), 2)) / xx;
	result += offset + (xx - 0.5)*std::log(xx) - xx;
	return result;
    }

    constexpr double psi_p1[7] = {
	0.895385022981970e-02, 0.477762828042627e+01,
	0.142441585084029e+03, 0.118645200713425e+04,
	0.363351846806499e+04, 0.413810161269013e+04,
	0.130560269827897e+04
    };
    constexpr double psi_q1[6] = {
	0.448452573429826e+02, 0.520752771467162e+03,
	0.221000799247830e+04, 0.364127349079381e+04,
	0.190831076596300e+04, 0.691091682714533e-05
    };
    constexpr double psi_p2[4] = {
	-0.212940445131011e+01, -0.701677227766759e+01,
	-0.448616543918019e+01, -0.648157123766197e+00
    };
    constexpr double psi_q2[4] = {
	0.322703493791143e+02, 0.892920700481861e+02,
	0.546117738103215e+02, 0.777788548522962e+01
    };

    XSF_HOST_DEVICE inline double psi(double xx) {
	/* Included only for benchmarking, prefer xsf::digamma */
 
	//                    Evaluation of the digamma function
	//
	//                          -----------
	//
	//    Psi(xx) is assigned the value 0 when the digamma function cannot
	//    be computed.
	//
	//    The main computation involves evaluation of rational chebyshev
	//    approximations published in math. Comp. 27, 123-127(1973) By
	//    cody, strecok and thacher.
	//
	//    ----------------------------------------------------------------
	//    Psi was written at Argonne National Laboratory for the FUNPACK
	//    package of special function subroutines. Psi was modified by
	//    A.H. Morris (nswc).
	
	double aug, den, dx0, sgn, upper, w, x, xmax1, xmx0, xsmall, z;
	int nq, i;
	dx0 = 1.461632144968362341262659542325721325;
	xmax1 = 4503599627370496.0;
	xsmall = 1e-9;
	x = xx;
	aug = 0.0;

	if (x < 0.5) {
	    if (std::abs(x) <= xsmall) {
		if (x == 0.) {
		    return 0.0;
		}
		aug = -1./x;
	    } else {
		// 10
		w = -x;
		sgn = M_PI / 4;
		if (w <= 0.) {
		    w = -w;
		    sgn = -sgn;
		}
		// 20
		if (w >= xmax1) {
		    return 0.0;
		}
		w -= static_cast<int>(w);
		nq = static_cast<int>(w*4.0);
		w = 4.*(w - 0.25*nq);
		
		if (nq % 2 == 1) {
		    w = 1. - w;
		}
		z = (M_PI / 4.)*w;
		
		if ((nq / 2) % 2 == 1) {
		    sgn = -sgn;
		}
		if ((((nq + 1) / 2) % 2) == 1) {
		    aug = sgn * (std::tan(z)*4.);
		} else {
		    if (z == 0.) {
			return 0.0;
		    }
		    aug = sgn * (4./std::tan(z));
		}
	    }
	    x = 1 - x;
	}

	if (x <= 3.0) {
	    // 50
	    den = x;
	    upper = psi_p1[0]*x;
	    for (i = 0; i < 5; i++)
		{
		    den = (den + psi_q1[i])*x;
		    upper = (upper + psi_p1[i+1])*x;
		}
	    den = (upper + psi_p1[6]) / (den + psi_q1[5]);
	    xmx0 = x - dx0;
	    return (den * xmx0) + aug;
	} else {
	    // 70
	    if (x < xmax1) {
		w = 1. / (x*x);
		den = w;
		upper = psi_p2[0]*w;
		
		for (i = 0; i < 3; i++) {
		    den = (den + psi_q2[i])*w;
		    upper = (upper + psi_p2[i+1])*w;
		}
		aug += upper / (den + psi_q2[3]) - 0.5/x;
	    }
	    return aug + std::log(x);
	}
    }
    
    XSF_HOST_DEVICE inline double apser(double a, double b, double x, double eps) {
	//    Apser yields the incomplete beta ratio I_(1-x))(b,a) for
	//    a <= Min(eps,eps*b), b*x <= 1, And x <= 0.5. Used when
	//    a is very small. Use only if above inequalities are satisfied.
	
	double aj, bx, c, j, s, t, tol;
	double g = 0.577215664901532860606512090082;

	bx = b*x;
	t = x - bx;
	if (b*eps > 0.02) {
	    c = std::log(bx) + g + t;
	} else {
	    c = std::log(x) + xsf::digamma(b) + g + t;
	}

	tol = 5.*eps*std::abs(c);
	j = 1.0;
	s = 0.0;

	while (1) {
	    j += 1;
	    t *= (x - bx/j);
	    aj = t / j;
	    s += aj;
	    if (std::abs(aj) <= tol) { break; }
	}
	return -a * (c + s);
    }

    constexpr double rlog_p[3] = {
	.333333333333333, -.224696413112536, .620886815375787e-02
    };
    constexpr double rlog_q[2] = {-.127408923933623e+01, .354508718369557};

    XSF_HOST_DEVICE inline double rlog(double x) {
	// Computation of  x - 1 - ln(x)

	double r, t, u, w, w1;
	double a = .566749439387324e-01;
	double b = .456512608815524e-01;

	if ((x < 0.61) || (x > 1.57))
	    {
		return ((x - 0.5) - 0.5) - std::log(x);
	    }

	if (x < 0.82)
	    {
		u = (x - 0.7) / 0.7;
		w1 = a - u*0.3;
	    }
	else if (x > 1.18)
	    {
		u = 0.75*x - 1.;
		w1 = b + u/3.;
	    }
	else
	    {
		u = (x - 0.5) - 0.5;
		w1 = 0.0;
	    }

	r = u / (u + 2.);
	t = r*r;
	w = ((rlog_p[2]*t+rlog_p[1])*t+rlog_p[0])/ ((rlog_q[1]*t+rlog_q[0])*t+1.);
	return 2.*t*(1. / (1. - r) - r*w) + w1;
    }

    XSF_HOST_DEVICE inline double rlog1(double x) {
	// Evaluation of the function x - ln(1 + x)

	double a = .566749439387324e-01;
	double b = .456512608815524e-01;
	double p0 = .333333333333333e+00;
	double p1 = -.224696413112536e+00;
	double p2 = .620886815375787e-02;
	double q1 = -.127408923933623e+01;
	double q2 = .354508718369557e+00;
	double h, r, t, w, w1;

	if ((-0.39 <= x) && (x <= 0.57)) {
	    if ((-0.18 <= x) && (x <= 0.18))
		{
		    h = x;
		    w1 = 0.0;
		} else if (x < -0.18)
		{
		    h = (x + 0.3)/0.7;
		    w1 = a - h*0.3;
		} else
		{
		    // 0.57 >= x > 0.18
		    h = 0.75*x - 0.25;
		    w1 = b + h/3.0;
		}

	    r = h / (h + 2);
	    t = r*r;
	    w = ((p2*t + p1)*t + p0) / ((q2*t + q1)*t + 1.);
	    return 2*t*(1./(1.-r) - r*w) + w1;

	} else {

	    return x - std::log((x + 0.5) + 0.5);

	}
    }

    constexpr double erfc1_a[5] = {.771058495001320e-04, -.133733772997339e-02,
	.323076579225834e-01, .479137145607681e-01,
	.128379167095513e+00};
    constexpr double erfc1_b[3] = {.301048631703895e-02, .538971687740286e-01,
	.375795757275549e+00};
    constexpr double erfc1_p[8] = {-1.36864857382717e-07, 5.64195517478974e-01,
	7.21175825088309e+00, 4.31622272220567e+01,
	1.52989285046940e+02, 3.39320816734344e+02,
	4.51918953711873e+02, 3.00459261020162e+02};
    constexpr double erfc1_q[8] = {1.00000000000000e+00, 1.27827273196294e+01,
	7.70001529352295e+01, 2.77585444743988e+02,
	6.38980264465631e+02, 9.31354094850610e+02,
	7.90950925327898e+02, 3.00459260956983e+02};
    constexpr double erfc1_r[5] = {2.10144126479064e+00, 2.62370141675169e+01,
	2.13688200555087e+01, 4.65807828718470e+00,
	2.82094791773523e-01};
    constexpr double erfc1_s[4] = {9.41537750555460e+01, 1.87114811799590e+02,
	9.90191814623914e+01, 1.80124575948747e+01};
    
    XSF_HOST_DEVICE inline double erfc1(int ind, double x) {
	//    Evaluation of the complementary error function
	//
	//    Erfc1(ind,x) = erfc(x)            if ind = 0
	//    Erfc1(ind,x) = exp(x*x)*erfc(x)   otherwise
	
	double ax, bot, t, top, result;
	double c = 0.564189583547756;

	if (x <= -5.6) { return (ind == 0 ? 2.0 : (2*std::exp(x*x))); }

	// sqrt(log(np.finfo(np.float64).max)) ~= 26.64
	if ((ind == 0) && (x > 26.64))  { return 0.0; }
	
	ax = std::abs(x);

	if (ax <= 0.5) {
	    t = x*x;
	    top = (((((erfc1_a[0])*t+erfc1_a[1])*t+erfc1_a[2])*t+erfc1_a[3])*t+erfc1_a[4]) + 1.0;
	    bot = (((erfc1_b[0])*t+erfc1_b[1])*t+erfc1_b[2])*t + 1.0;
	    result = 0.5 + (0.5 - x*(top/bot));
	    return (ind == 0 ? result : result*std::exp(t));
	} else if ((0.5 < ax) && (ax <= 4.0)) {
	    top = (((((((erfc1_p[0]
			 )*ax+erfc1_p[1]
			)*ax+erfc1_p[2]
		       )*ax+erfc1_p[3]
		      )*ax+erfc1_p[4]
		     )*ax+erfc1_p[5]
		    )*ax+erfc1_p[6]
		   )*ax + erfc1_p[7];
	    bot = (((((((erfc1_q[0]
			 )*ax+erfc1_q[1]
			)*ax+erfc1_q[2]
		       )*ax+erfc1_q[3]
		      )*ax+erfc1_q[4]
		     )*ax+erfc1_q[5]
		    )*ax+erfc1_q[6])*ax + erfc1_q[7];
	    result = top / bot;
	} else {
	    t = std::pow(1 / x, 2);
	    top = (((erfc1_r[0]*t+erfc1_r[1])*t+erfc1_r[2])*t+erfc1_r[3])*t + erfc1_r[4];
	    bot = (((erfc1_s[0]*t+erfc1_s[1])*t+erfc1_s[2])*t+erfc1_s[3])*t + 1.0;
	    result = (c - t*(top/bot)) / ax;
	}
	if (ind == 0) {
	    result *= std::exp(-(x*x));
	    return (x < 0 ? (2.0 - result) : result);
	} else {
	    return (x < 0 ? (2.0*std::exp(x*x) - result) : result);
	}
    }
    
    constexpr double bcorr_carr[6] = {
	0.833333333333333e-01, -0.277777777760991e-02,
	0.793650666825390e-03, -0.595202931351870e-03,
	0.837308034031215e-03, -0.165322962780713e-02
    };

    XSF_HOST_DEVICE inline double bcorr(double a0, double b0) {
	//    Evaluation of  del(a0) + del(b0) - del(a0 + b0)  where
	//    ln(gamma(a)) = (a - 0.5)*ln(a) - a + 0.5*ln(2*pi) + del(a).
	//    It is assumed that a0 >= 8 And b0 >= 8.

	double a,b,c,h,s11,s3,s5,s7,s9,t,w,x,x2;

	a = std::fmin(a0, b0);
	b = std::fmax(a0, b0);
	h = a / b;
	c = h/(1. + h);
	x = 1./(1. + h);
	x2 = x*x;
	//  Set sn = (1 - x**n)/(1 - x)
	s3 = 1. + (x + x2);
	s5 = 1. + (x + x2*s3);
	s7 = 1. + (x + x2*s5);
	s9 = 1. + (x + x2*s7);
	s11 = 1. + (x + x2*s9);
	// Set w = del(b) - del(a + b)
	t = std::pow((1. / b), 2);
	w = (((((bcorr_carr[5]*s11
		 )*t + bcorr_carr[4]*s9
		)*t + bcorr_carr[3]*s7
	       )*t + bcorr_carr[2]*s5
	      )*t + bcorr_carr[1]*s3
	     )*t + bcorr_carr[0];
	w *= c / b;
	// Compute  del(a) + w
	t = std::pow((1. / a), 2);
	return ((((((bcorr_carr[5])*t + bcorr_carr[4]
		    )*t + bcorr_carr[3]
		   )*t + bcorr_carr[2]
		  )*t + bcorr_carr[1]
		 )*t + bcorr_carr[0]
		)/a + w;
    }

    XSF_HOST_DEVICE inline double basym(double a, double b, double lmbda, double eps) {
	//    Asymptotic expansion for ix(a,b) for large a and b.
	//    Lambda = (a + b)*y - b  and eps is the tolerance used.
	//    It is assumed that lambda is nonnegative and that
	//    a and b are greater than or equal to 15.
	
	double a0[21] = { 0.0 };
	double b0[21] = { 0.0 };
	double c[21] = { 0.0 };
	double d[21] = { 0.0 };
	double bsum, dsum, f, h, h2, hn, j0, j1, r, r0, r1, s, ssum;
	double t, t0, t1, u, w, w0, z, z0, z2, zn, znm1;
	double e0 = 2. / std::sqrt(M_PI);
	double e1 = std::pow(2.0, (-3./2.));
	int i, imj, j, m, mmj, n, num;

	// ****** Num is the maximum value that n can take in the do looprlog1
	//        ending at statement 50. It is required that num be even.
	//        The arrays a0, b0, c, d have dimension num + 1.
	num = 20;

	if (a < b) {
	    h = a / b;
	    r0 = 1./(1.+h);
	    r1 = (b - a) / b;
	    w0 = 1. / std::sqrt(a * (1. + h));
	} else {
	    h = b / a;
	    r0 = 1./(1.+h);
	    r1 = (b - a) / a;
	    w0 = 1. / std::sqrt(b * (1. + h));
	}
	f = -a*xsf::cephes::log1pmx(-lmbda/a) - b*xsf::cephes::log1pmx(lmbda/b);
	f = -a*xsf::cephes::log1pmx(-lmbda/a) - b*xsf::cephes::log1pmx(lmbda/b);
	t = std::exp(-f);
	if (t == 0.0) { return 0.0; }
	z0 = std::sqrt(f);
	z = 0.5*(z0/e1);
	z2 = f + f;
	
	a0[0] = (2./3.)*r1;
	c[0] = -0.5*a0[0];
	d[0] = -c[0];
	j0 = (0.5/e0)*erfc1(1, z0);
	j1 = e1;
	ssum = j0 + d[0]*w0*j1;

	s = 1.0;
	h2 = h*h;
	hn = 1.0;
	w = w0;
	znm1 = z;
	zn = z2;

	for (n = 2; n <= num; n += 2) {
	    hn *= h2;
	    a0[n-1] = 2.*r0*(1.+ h*hn)/(n + 2.);
	    s += hn;
	    a0[n] = 2.*r1*s/(n + 3.);

	    for (i = n; i <= n + 1; i++) {
		r = -0.5*(i + 1.);
		b0[0] = r*a0[0];
		
		for (m = 2; m <= i; m++) {
		    bsum = 0.0;
		    for (j = 1; j < m; j++) {
			mmj = m - j;
			bsum += (j*r - mmj)*a0[j-1]*b0[mmj-1];
		    }
		    b0[m-1] = r*a0[m-1] + bsum/m;
		}
		c[i-1] = b0[i-1] / (i + 1.);
		dsum = 0.0;
		
		for (j = 1; j < i; j++) {
		    imj = i - j;
		    dsum += d[imj-1]*c[j-1];
		}
		d[i-1] = - (dsum+c[i-1]);
	    }
	    j0 = e1*znm1 + (n-1.)*j0;
	    j1 = e1*zn + n*j1;
	    znm1 *= z2;
	    zn *= z2;
	    w *= w0;
	    t0 = d[n-1]*w*j0;
	    w *= w0;
	    t1 = d[n]*w*j1;
	    ssum += t0 + t1;
	    if ((std::abs(t0) + std::abs(t1)) <= eps*ssum) { break; }
	}
	u = std::exp(-bcorr(a, b));
	return e0*t*u*ssum;
    }

    constexpr double gamln1_p[7] = {
	.577215664901533e+00,  .844203922187225e+00,
	-.168860593646662e+00, -.780427615533591e+00,
	-.402055799310489e+00, -.673562214325671e-01,
	-.271935708322958e-02
    };
    constexpr double gamln1_q[6] = {
	.288743195473681e+01, .312755088914843e+01,
	.156875193295039e+01, .361951990101499e+00,
	.325038868253937e-01, .667465618796164e-03
    };
    constexpr double gamln1_r[6] = {
	.422784335098467e+00, .848044614534529e+00,
	.565221050691933e+00, .156513060486551e+00,
	.170502484022650e-01, .497958207639485e-03
    };
    constexpr double gamln1_s[5] = {
	.124313399877507e+01, .548042109832463e+00,
	.101552187439830e+00, .713309612391000e-02,
	.116165475989616e-03
    };

    XSF_HOST_DEVICE inline double gamln1(double a) {
	//    Evaluation of ln(gamma(1 + a)) for -0.2 <= A <= 1.25

	double bot, top, w, x;
	
	if (a < 0.6) {
	    top = ((((((gamln1_p[6]
			)*a+gamln1_p[5]
		       )*a+gamln1_p[4]
		      )*a+gamln1_p[3]
		     )*a+gamln1_p[2]
		    )*a+gamln1_p[1]
		   )*a+gamln1_p[0];
	    bot = ((((((gamln1_q[5]
			)*a+gamln1_q[4]
		       )*a+gamln1_q[3]
		      )*a+gamln1_q[2]
		     )*a+gamln1_q[1]
		    )*a+gamln1_q[0]
		   )*a+1.;
	    w = top/bot;
	    return -a*w;
	} else {
	    x = (a - 0.5) - 0.5;
	    top = (((((gamln1_r[5]
		       )*x+gamln1_r[4]
		      )*x+gamln1_r[3]
		     )*x+gamln1_r[2]
		    )*x+gamln1_r[1]
		   )*x+gamln1_r[0];
	    bot = (((((gamln1_s[4]
		       )*x+gamln1_s[3]
		      )*x+gamln1_s[2]
		     )*x+gamln1_s[1]
		    )*x+gamln1_s[0]
		   )*x+1.;
	    w = top/bot;
	    return x*w;
	}
    }

    constexpr double gamln_c[6] = {
	.833333333333333e-01, -.277777777760991e-02,
	.793650666825390e-03, -.595202931351870e-03,
	.837308034031215e-03, -.165322962780713e-02
    };

    XSF_HOST_DEVICE inline double gamln(double a) {
	//    Evaluation of ln(gamma(a)) for positive a
	
	double t, w, d = .418938533204673;
	int i,n;


	if (a <= 0.8) { return gamln1(a) - std::log(a); }

	if (a <= 2.25) {
	    t = (a-0.5) - 0.5;
	    return gamln1(t);
	}

	if (a < 10) {
	    n = static_cast<int>(a - 1.25);
	    t = a;
	    w = 1.0;
	    for (i = 0; i < n; i++)
		{
		    t -= 1.0;
		    w *= t;
		}
	    return gamln1(t-1.) + std::log(w);
	}
	t = std::pow(1/a, 2);
 	w = (((((gamln_c[5]*t+gamln_c[4])*t+gamln_c[3])*t+gamln_c[2])*t+gamln_c[1])*t+gamln_c[0])/a;
	return (d + w) + (a-0.5)*(std::log(a) - 1.);
    }

    XSF_HOST_DEVICE double gammaln(double x) {
	/* We determined empirically where cephes and cdflib do best */
	if (
	    (x <= 0) || (x > 1.18 && x < 1.6) || (x > 2 && x < 2.85) ||
	    (x > 10 && x < 13)
	    ) {
	    return xsf::cephes::lgam(x);
	}
	return xsf::cdflib::gamln(x);
    }

    XSF_HOST_DEVICE inline double gsumln(double a, double b) {
	//     Evaluation of the function ln(gamma(a + b))
	//     for 1 <= A <= 2  And  1 <= B <= 2

	double x;

	x = a + b - 2;
	if (x <= 0.25) {
	    return gamln1(1. + x);
	}

	if (x <= 1.25) {
	    return gamln1(x) + std::log1p(x);
	}

	return gamln1(x - 1.) + std::log(x*(1. + x));
    }

    XSF_HOST_DEVICE inline double betaln(double a0, double b0) {
	//    Evaluation of the logarithm of the beta function
	double a, b, c, h, u, v, w, z;
	double e = 0.918938533204673;
	int i, n;

	a = std::fmin(a0, b0);
	b = std::fmax(a0, b0);

	if (a >= 8.0) {
	    w = bcorr(a, b);
	    h = a / b;
	    c = h/(1. + h);
	    u = -(a - 0.5)*std::log(c);
	    v = b*std::log1p(h);
	    if (u > v) {
		return (((-0.5*std::log(b)+e)+w)-v) - u;
	    } else {
		return (((-0.5*std::log(b)+e)+w)-u) - v;
	    }
	}
	if (a < 1) {
	    if (b > 8) {
		return gammaln(a) + algdiv(a,b);
	    } else {
		return gammaln(a) + (gammaln(b) - gammaln(a+b));
	    }
	}
	
	if (a <= 2) {
	    if (b <= 2) {
		return gammaln(a) + gammaln(b) - gsumln(a, b);
	    }
	    if (b >= 8) {
		return gammaln(a) + algdiv(a, b);
	    }
	    w = 0.;
	}

	if (a > 2) {
	    if (b <= 1000) {
		n = static_cast<int>(a - 1.);
		w = 1.;
		for (i = 0; i < n; i++) {
		    a -= 1.0;
		    h = a / b;
		    w *= h/(1.+h);
		}
		w = std::log(w);
		if (b >= 8.0) {
		    return w + gammaln(a) + algdiv(a, b);
		}
	    } else {
		n = static_cast<int>(a - 1.);
		w = 1.0;
		for (i = 0; i < n; i++) {
		    a -= 1.0;
		    w *= a/(1. + (a/b));
		}
		return (std::log(w) - n*std::log(b)) + (gammaln(a) + algdiv(a, b));
	    }
	}
	n = static_cast<int>(b - 1.);
	z = 1.0;
	for (i = 0; i < n; i++) {
	    b -= 1.0;
	    z *= b / (a + b);
	}
	return w + std::log(z) + (gammaln(a) + gammaln(b) - gsumln(a, b));
    }

    constexpr double gam1_p[7] = {
	.577215664901533e+00, -.409078193005776e+00,
	-.230975380857675e+00, .597275330452234e-01,
	.766968181649490e-02, -.514889771323592e-02,
	.589597428611429e-03
    };
    constexpr double gam1_q[5] = {
	.100000000000000e+01, .427569613095214e+00,
	.158451672430138e+00, .261132021441447e-01,
	.423244297896961e-02
    };
    constexpr double gam1_r[9] = {
	-.422784335098468e+00, -.771330383816272e+00,
	-.244757765222226e+00, .118378989872749e+00,
	.930357293360349e-03, -.118290993445146e-01,
	.223047661158249e-02, .266505979058923e-03,
	-.132674909766242e-03
    };
    constexpr double gam1_s[2] = {.273076135303957e+00, .559398236957378e-01};

    XSF_HOST_DEVICE inline double gam1(double a) {
	//    Computation of 1/gamma(a+1) - 1  for -0.5 <= A <= 1.5
	
	double bot, d, t, top, w;

	d = a - 0.5;
	t = (d > 0 ? d - 0.5 : a);
	
	if (t == 0.0) { return 0.0; }

	if (t < 0) {
	    top = ((((((((gam1_r[8]
			  )*t+gam1_r[7]
			 )*t+gam1_r[6]
			)*t+gam1_r[5]
		       )*t+gam1_r[4]
		      )*t+gam1_r[3]
		     )*t+gam1_r[2]
		    )*t+gam1_r[1]
		   )*t + gam1_r[0];
	    bot = (gam1_s[1]*t + gam1_s[0])*t + 1.0;
	    w = top / bot;
	    if (d > 0.0) {
		return t*w/a;
	    } else {
		return a * ((w + 0.5) + 0.5);
	    }
	}
	top = ((((((gam1_p[6]
		    )*t+gam1_p[5]
		   )*t+gam1_p[4]
		  )*t+gam1_p[3]
		 )*t+gam1_p[2]
		)*t+gam1_p[1]
	       )*t + gam1_p[0];
	bot = ((((gam1_q[4]
		  )*t+gam1_q[3]
		 )*t+gam1_q[2]
		)*t+gam1_q[1]
	       )*t + 1.0;
	w = top / bot;
	if (d > 0.0) {
	    return (t/a) * ((w - 0.5) - 0.5);
	} else {
	    return a * w;
	}
    }

    XSF_HOST_DEVICE inline double brcomp(double a, double b, double x, double y) {
	//    Evaluation of x**a*y**b/beta(a,b)

	double a0, apb, b0, c, e, h, lmbda, lnx, lny, t, u, v, x0, y0, z;
	double constexpr r2pi = 0.3989422804014327;  // 1. / std::sqrt(2 * M_PI);
	int i, n;

	if ((x == 0.) || (y == 0.)) {
	    return 0.;
	}
	a0 = std::fmin(a, b);

	if (a0 >= 8.) {
	    if (a > b) {
		h = b / a;
		x0 = 1. / (1. + h);
		y0 = h / (1. + h);
		lmbda = (a + b)*y - b;
	    } else {
		h = a / b;
		x0 = h / (1. + h);
		y0 = 1. / (1. + h);
		lmbda = a - (a + b)*x;
	    }
	    e = -lmbda / a;
	    if (fabs(e) > 0.6) {
		u = e - std::log(x / x0);
	    } else {
		u = -xsf::cephes::log1pmx(e);
	    }
	    e = lmbda / b;
	    if (std::abs(e) > 0.6) {
		v = e - std::log(y / y0);
	    } else {
		v = -xsf::cephes::log1pmx(e);
	    }

	    z = std::exp(-(a*u + b*v));
	    return r2pi*std::sqrt(b*x0)*z*std::exp(-bcorr(a, b));
	}
	if (x <= 0.375) {
	    lnx = std::log(x);
	    lny = std::log1p(-x);
	} else {
	    lnx = (y > 0.375 ? log(x) : std::log1p(-y));
	    lny = std::log(y);
	}
	z = a*lnx + b*lny;
	if (a0 >= 1.0) {
	    z -= betaln(a, b);
	    return std::exp(z);
	}

	b0 = std::fmax(a, b);
	if (b0 >= 8.0) {
	    u = gamln1(a0) + algdiv(a0, b0);
	    return a0*std::exp(z - u);
	}

	if (b0 > 1.0) {
	    u = gamln1(a0);
	    n = static_cast<int>(b0 - 1.);
	    if (n >= 1) {
		c = 1.0;
		for (i = 0; i < n; i++) {
		    b0 -= 1.0;
		    c *= b0 / (a0 + b0);
		}
		u += std::log(c);
	    }
	    z -= u;
	    b0 -= 1.0;
	    apb = a0 + b0;
	    
	    if (apb > 1.0) {
		u = a0 + b0 - 1.0;
		t = (1. + gam1(u)) / apb;
	    } else {
		t = 1. + gam1(apb);
	    }
	    return a0*std::exp(z)*(1. + gam1(b0)) / t;
	}
	if (std::exp(z) == 0.) {
	    return 0.0;
	}

	apb = a + b;
	t = std::exp(z);
	if (apb > 1.0) {
	    u = a + b - 1.0;
	    z = (1. + gam1(u)) / apb;
	} else {
	    z = 1. + gam1(apb);
	}
	c = (1. + gam1(a)) * (1. + gam1(b)) / z;
	return t * (a0*c) / (1. + a0 / b0);
    }

    XSF_HOST_DEVICE inline double bfrac(double a, double b, double x, double y, double lmbda, double eps) {
	//    Continued fraction expansion for ix(a,b) when a,b > 1.
	//    It is assumed that  lambda = (a + b)*y - b.
	
	double alpha, beta, e, r0, t, w, result;
	double c = 1. + lmbda;
	double c0 = b / a;
	double c1 = 1. + (1. / a);
	double yp1 = y + 1.;
	double n = 0.;
	double p = 1.;
	double s = a + 1.;
	double an = 0.;
	double bn = 1.;
	double anp1 = 1.;
	double bnp1 = c / c1;
	double r = c1 / c;
	
	result = brcomp(a, b, x, y);
	
	if (result == 0.0) { return 0; }
	
	// Continued fraction calculation
	while (1) {
	    n += 1.0;
	    t = n / a;
	    w = n * (b - n)*x;
	    e = a / s;
	    alpha = (p*(p + c0)*e*e) * (w*x);
	    e = (1. + t) / (c1 + t + t);
	    beta = n + (w / s) + e*(c + n*yp1);
	    p = 1. + t;
	    s += 2.;
	    // Update an, bn, anp1, and bnp1
	    t = alpha*an + beta*anp1;
	    an = anp1;
	    anp1 = t;
	    t = alpha*bn + beta*bnp1;
	    bn = bnp1;
	    bnp1 = t;
	    r0 = r;
	    r = anp1 / bnp1;
	    if (!(std::abs(r - r0) > eps*r)) { break; }
	    // Rescale an, bn, anp1, and bnp1
	    an /= bnp1;
	    bn /= bnp1;
	    anp1 = r;
	    bnp1 = 1.0;
	}
	return result*r;
    }

    constexpr double rexp_p[2] = {.914041914819518e-09, .238082361044469e-01};
    constexpr double rexp_q[4] = {
	-.499999999085958e+00, .107141568980644e+00,
	-.119041179760821e-01, .595130811860248e-03
    };
    XSF_HOST_DEVICE inline double rexp(double x) {
	// Evaluation of exp(-x)*x**a/gamma(a)
	double w;

	if (std::abs(x) <= 0.15) {
	    return x*(((rexp_p[1]*x+rexp_p[0])*x+1.)/((((rexp_q[3]*x+rexp_q[2])*x+rexp_q[1])*x+rexp_q[0])*x + 1.));
	}
	w = std::exp(x);
	if (x > 0.) {
	    return w*(0.5+ (0.5 - 1./w));
	} 
	return (w - 0.5) - 0.5;
    }

    constexpr double erf_a[5] = {
	.771058495001320e-04, -.133733772997339e-02,
	.323076579225834e-01, .479137145607681e-01,
	.128379167095513e+00
    };
    constexpr double erf_b[3] = {
	.301048631703895e-02, .538971687740286e-01,
	.375795757275549e+00
    };
    constexpr double erf_p[8] = {
	-1.36864857382717e-07, 5.64195517478974e-01,
	7.21175825088309e+00, 4.31622272220567e+01,
	1.52989285046940e+02, 3.39320816734344e+02,
	4.51918953711873e+02, 3.00459261020162e+02
    };
    constexpr double erf_q[8] = {
	1.00000000000000e+00, 1.27827273196294e+01,
	7.70001529352295e+01, 2.77585444743988e+02,
	6.38980264465631e+02, 9.31354094850610e+02,
	7.90950925327898e+02, 3.00459260956983e+02
    };
    constexpr double erf_r[5] = {
	2.10144126479064e+00, 2.62370141675169e+01,
	2.13688200555087e+01, 4.65807828718470e+00,
	2.82094791773523e-01
    };
    constexpr double erf_s[4] = {
	9.41537750555460e+01, 1.87114811799590e+02,
	9.90191814623914e+01, 1.80124575948747e+01
    };

    XSF_HOST_DEVICE inline double erf(double x) {
	//    Evaluation of the real error function

	double ax, bot, t, top;
	double c = .564189583547756;
	ax = std::abs(x);
	if (ax <= 0.5) {
	    t = x*x;
	    top = ((((erf_a[0]*t+erf_a[1])*t+erf_a[2])*t+erf_a[3])*t+erf_a[4]) + 1.0;
	    bot = ((erf_b[0]*t+erf_b[1])*t+erf_b[2])*t + 1.0;
	    return x*(top/bot);
	}
	if (ax <= 4.0) {
	    top = (((((((erf_p[0]
			 )*ax+erf_p[1]
			)*ax+erf_p[2]
		       )*ax+erf_p[3]
		      )*ax+erf_p[4]
		     )*ax+erf_p[5]
		    )*ax+erf_p[6]
		   )*ax + erf_p[7];
	    bot = (((((((erf_q[0]
			 )*ax+erf_q[1]
			)*ax+erf_q[2]
		       )*ax+erf_q[3]
		      )*ax+erf_q[4]
		     )*ax+erf_q[5]
		    )*ax+erf_q[6])*ax + erf_q[7];
	    t = 0.5 + (0.5 - std::exp(-x*x)*(top/bot));
	    return (x < 0 ? -t : t);
	}
	if (ax < 5.8) {
	    t = std::pow(1/x, 2);
	    top = (((erf_r[0]*t+erf_r[1])*t+erf_r[2])*t+erf_r[3])*t + erf_r[4];
	    bot = (((erf_s[0]*t+erf_s[1])*t+erf_s[2])*t+erf_s[3])*t + 1.0;
	    t = 0.5 + (0.5 - std::exp(-x*x) * (c - top/(x*x*bot))/ax);
	    return (x < 0 ? -t : t);
	}
	return (x < 0 ? -1 : 1);
    }
    
    XSF_HOST_DEVICE inline std::pair<double, double> grat1(double a, double x, double r, double eps) {
	//        Evaluation of the incomplete gamma ratio functions
	//                    p(a,x) and q(a,x)
	//
	//    it is assumed that a <= 1.  Eps is the tolerance to be used.
	//    the input argument r has the value e**(-x)*x**a/gamma(a).

	double a2n, a2nm1, am0, an, an0, b2n, b2nm1, c, cma, g, h, j, l;
	double p, q, ssum, t, tol, w, z;

	if (a*x == 0.) {
	    if (x > a) {
		return {1.0, 0.0};
	    } else {
		return {0.0, 1.0};
	    }
	}

	if (a == 0.5) {
	    if (x < 0.25) {
		p = erf(sqrt(x));
		return {p, 0.5 + (0.5 - p)};
	    } else {
		q = erfc1(0, std::sqrt(x));
		return {0.5 + (0.5 - q), q};
	    }
	}

	if (x < 1.1) {
	    //
	    // Taylor series for p(a,x)/x**a
	    //
	    an = 3.0;
	    c = x;
	    ssum = x / (a + 3.);
	    tol = 0.1*eps / (a + 1.);
	    while (1) {
		an += 1;
		c *= -(x / an);
		t = c / (a + an);
		ssum += t;
		if (std::abs(t) <= tol) { break; }
	    }
	    j = a*x*((ssum/6. - 0.5/(a+2.))*x + 1./(a+1.));
	    
	    z = a * std::log(x);
	    h = gam1(a);
	    g = 1. + h;

	    if (((x >= 0.25) && (a >= x /2.59)) || ((x < 0.25) && (z <= -0.13394))) {
		w = std::exp(z);
		p = w*g*(0.5 + (0.5 - j));
		q = 0.5 + (0.5 - p);
		return {p, q};
	    } else {
		l = rexp(z);
		w = 0.5 + (0.5 + l);
		q = (w*j - l)*g - h;
		if (q < 0.0) {
		    return {1.0, 0.0};
		}
		p = 0.5 + (0.5 - q);
		return {p, q};
	    }
	}
	//  Continued fraction expansion
	a2nm1 = 1.0;
	a2n = 1.0;
	b2nm1 = x;
	b2n = x + (1. - a);
	c = 1.0;
	while (1) {
	    a2nm1 = x*a2n + c*a2nm1;
	    b2nm1 = x*b2n + c*b2nm1;
	    am0 = a2nm1/b2nm1;
	    c = c + 1.;
	    cma = c - a;
	    a2n = a2nm1 + cma*a2n;
	    b2n = b2nm1 + cma*b2n;
	    an0 = a2n/b2n;
	    if (!(std::abs(an0-am0) >= eps*an0)) { break; }
	}
	q = r*an0;
	p = 0.5 + (0.5 - q);
	return {p, q};
    }

    constexpr double gamma_p[7] = {
	.539637273585445e-03, .261939260042690e-02,
	.204493667594920e-01, .730981088720487e-01,
	.279648642639792e+00, .553413866010467e+00,
	1.0
    };
    constexpr double gamma_q[7] = {
	-.832979206704073e-03, .470059485860584e-02,
	.225211131035340e-01, -.170458969313360e+00,
	-.567902761974940e-01, .113062953091122e+01,
	1.0
    };
    constexpr double gamma_r[5] = {
	.820756370353826e-03, -.595156336428591e-03,
	.793650663183693e-03, -.277777777770481e-02,
	.833333333333333e-01
    };

    XSF_HOST_DEVICE inline double gamma(double a) {
	//        Evaluation of the gamma function for real arguments
	//
	//    Gamma(a) is assigned the value 0 when the gamma function cannot
	//    be computed.

	double bot, g, lnx, t, top, w, z, result;
	int i, j, m, n;
	double s = 0.0;
	double d = 0.5*(std::log(2.*M_PI) - 1);
	double x = a;


	result = 0.0;
	if (std::abs(a) < 15) {
	    t = 1.0;
	    m = static_cast<int>(a) - 1;
	    if (m > 0) {
		for (j = 0; j < m; j++)
		    {
			x -= 1.0;
			t *= x;
		    }
		x -= 1.0;
	    } else if (m == 0) {
		x -= 1.0;
	    } else {
		t = a;
		if (a <= 0.) {
		    m = -m - 1;
		    if (m != 0.) {
			for (j = 0; j < m; j++)
			    {
				x += 1.0;
				t *= x;
			    }
		    }
		    x += 0.5;
		    x += 0.5;
		    t *= x;
		    
		    if (t == 0.) { return result; }
		}
		if (std::abs(t) < 1e-30) {
		    if (std::abs(t)*std::numeric_limits<double>::max() <= 1.0001) { return result; }
		    return 1./t;
		}
	    }
	    top = gamma_p[0];
	    bot = gamma_q[0];
	    for (i = 1; i < 7; i++)
		{
		    top *= x;
		    top += gamma_p[i];
		    bot *= x;
		    bot += gamma_q[i];
		}
	    result = top / bot;
	    return (a < 1.0 ? result/t : result*t);
	}
	
	if (std::abs(a) >= 1.e3) { return result; }

	if (a <= 0.0) {
	    x = -a;
	    n = static_cast<int>(x);
	    t = x - n;
	    if (t > 0.9) {
		t = 1. - t;
	    }
	    s = xsf::sinpi(t) / M_PI;
	    if (n % 2 == 0) {
		s = -s;
	    }
	    if (s == 0.0) { return result; }
	}
	t = std::pow(1 / x, 2);
	g = ((((gamma_r[0]*t+gamma_r[1])*t+gamma_r[2])*t+gamma_r[3])*t+gamma_r[4]) / x;
	lnx = std::log(x);
	z = x;
	g = (d + g) + (z -0.5)*(lnx - 1.);
	w = g;
	t = g - w;
	if (w > 0.99999*709) { return result; }
	result = std::exp(w)*(1. + t);
	return (a < 0.0 ? (1. / (result * s)) / x : result);
    }

    constexpr double gratio_big[3] = {20., 14., 10.};
    constexpr double gratio_e00[3] = {.00025, .025, .14};
    constexpr double gratio_x00[3] = {31., 17., 9.7};
    constexpr double gratio_d0[13] = {
	.833333333333333e-01, -.148148148148148e-01,
	.115740740740741e-02, .352733686067019e-03,
	-.178755144032922e-03, .391926317852244e-04,
	-.218544851067999e-05, -.185406221071516e-05,
	.829671134095309e-06, -.176659527368261e-06,
	.670785354340150e-08, .102618097842403e-07,
	-.438203601845335e-08
    };
    constexpr double gratio_d1[12] = {
	-.347222222222222e-02, .264550264550265e-02,
	-.990226337448560e-03, .205761316872428e-03,
	-.401877572016461e-06, -.180985503344900e-04,
	.764916091608111e-05, -.161209008945634e-05,
	.464712780280743e-08, .137863344691572e-06,
	-.575254560351770e-07, .119516285997781e-07
    };
    constexpr double gratio_d2[10] = {
	-.268132716049383e-02, .771604938271605e-03,
	.200938786008230e-05, -.107366532263652e-03,
	.529234488291201e-04, -.127606351886187e-04,
	.342357873409614e-07, .137219573090629e-05,
	-.629899213838006e-06, .142806142060642e-06
    };
    constexpr double gratio_d3[8] = {
	.229472093621399e-03, -.469189494395256e-03,
	.267720632062839e-03, -.756180167188398e-04,
	-.239650511386730e-06, .110826541153473e-04,
	-.567495282699160e-05, .142309007324359e-05
    };
    constexpr double gratio_d4[6] = {
	.784039221720067e-03, -.299072480303190e-03,
	-.146384525788434e-05, .664149821546512e-04,
	-.396836504717943e-04, .113757269706784e-04
    };
    constexpr double gratio_d5[4] = {
	-.697281375836586e-04, .277275324495939e-03,
	-.199325705161888e-03, .679778047793721e-04
    };
    constexpr double gratio_d6[2] = {-.592166437353694e-03, .270878209671804e-03};
    constexpr double gratio_acc0[3] = {5.e-15, 5.e-7, 5.e-4};

    XSF_HOST_DEVICE inline std::pair<double, double> gratio(double a, double x, int ind) {
	//    Evaluation of the incomplete gamma ratio functions
	//                    P(a,x) and Q(a,x)
	//
	//                    ----------
	//
	//    It is assumed that a and x are nonnegative, where a and x
	//    Are not both 0.
	//
	//    Ans and qans are variables. Gratio assigns ans the value
	//    P(a,x) and qans the value q(a,x). Ind may be any integer.
	//    If ind = 0 then the user is requesting as much accuracy as
	//    Possible (up to 14 significant digits). Otherwise, if
	//    Ind = 1 then accuracy is requested to within 1 unit of the
	//    6-Th significant digit, and if ind .ne. 0,1 Then accuracy
	//    Is requested to within 1 unit of the 3rd significant digit.
	//
	//    Error return ...
	//    Ans is assigned the value 2 when a or x is negative,
	//    When a*x = 0, or when p(a,x) and q(a,x) are indeterminant.
	//    P(a,x) and q(a,x) are computationally indeterminant when
	//    X is exceedingly close to a and a is extremely large.

	double d10 = -.185185185185185e-02;
	double d20 = .413359788359788e-02;
	double d30 = .649434156378601e-03;
	double d40 = -.861888290916712e-03;
	double d50 = -.336798553366358e-03;
	double d60 = .531307936463992e-03;
	double d70 = .344367606892378e-03;
	double alog10 = std::log(10);
	double rt2pi = std::sqrt(1. / (2.*M_PI));
	double rtpi = xsf::cephes::detail::SQRTPI;
	double eps = std::numeric_limits<double>::epsilon();
	double acc, a2n, a2nm1, am0, amn, an, an0, ans, apn, b2n, b2nm1;
	double c, c0, c1, c2, c3, c4, c5, c6, cma, e0, g, h, j, l, qans, r;
	double rta, rtx, s, ssum, t, t1, tol, twoa, u, w, x0, y, z;
	int i, m, n, last_entry;
	double wk[20] = {0.0};

	if ((!(a >= 0.0)) || (!(x >= 0.0))) {
	    return {2.0, 0.0};
	}
	if ((a == 0.0) && (x == 0.0)) {
	    return {2.0, 0.0};
	}

	if (a*x == 0.0) {
	    if (x > a) {
		return {1.0, 0.0};
	    } else {
		return {0.0, 1.0};
	    }
	}

	if ((!(ind == 0)) && (!(ind == 1))) {
	    ind = 2;
	}

	acc = std::fmax(gratio_acc0[ind], eps);
	e0 = gratio_e00[ind];
	x0 = gratio_x00[ind];

	if (a >= 1.0) {
	    if (a >= gratio_big[ind]) {
		// 30
		l = x / a;
		if (l == 0.0) {
		    // 370
		    return {0.0, 1.0};
		}
		s = 0.5 + (0.5 - l);
		z = rlog(l);
		if (z >= (700. / a)) {
		    // 410
		    if (std::abs(s) <= 2.*eps) {
			return {2.0, 0.0};
		    }
		    if (x > a) {
			return {1.0, 0.0};
		    } else {
			return {0.0, 1.0};
		    }
		}
		y = a*z;
		rta = std::sqrt(a);
		if (std::abs(s) <= (e0 / rta)) {
		    // 330
		    // TEMME EXPANSION FOR L = 1
		    //
		    if (a*eps*eps > 3.28e-3) {
			return {2.0, 0.0};
		    }
		    c = 0.5 + (0.5 - y);
		    w = (0.5 - std::sqrt(y) * (0.5 + (0.5 - y/3.))/rtpi)/c;
		    u = 1. / a;
		    z = std::sqrt(z+z);
		    if (l < 1.0) {
			z = -z;
		    }
		    if (ind == 0) {
			c0 = (((((((gratio_d0[6]
				    )*z+gratio_d0[5]
				   )*z+gratio_d0[4]
				  )*z+gratio_d0[3]
				 )*z+gratio_d0[2]
				)*z+gratio_d0[1]
			       )*z+gratio_d0[0]
			      )*z - (1./3.);
			c1 = ((((((gratio_d1[5]
				   )*z+gratio_d1[4]
				  )*z+gratio_d1[3]
				 )*z+gratio_d1[2]
				)*z+gratio_d1[1]
			       )*z+gratio_d1[0]
			      )*z + d10;
			c2 = (((((gratio_d2[4])*z+gratio_d2[3])*z+gratio_d2[2])*z+gratio_d2[1])*z+gratio_d2[0])*z + d20;
			c3 = ((((gratio_d3[3])*z+gratio_d3[2])*z+gratio_d3[1])*z+gratio_d3[0])*z + d30;
			c4 = (gratio_d4[1]*z+gratio_d4[0])*z + d40;
			c5 = (gratio_d5[1]*z+gratio_d5[0])*z + d50;
			c6 = gratio_d6[0]*z + d60;
			t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0;
		    } else if (ind == 1) {
			c0 = (gratio_d0[1]*z+gratio_d0[0])*z - (1. / 3.);
			c1 = gratio_d1[0]*z + d10;
			t = (d20*u+c1)*u + c0;
		    } else {
			t = gratio_d0[0]*z - (1. / 3.);
		    }
		    if (l < 1.) {
			ans = c * (w - rt2pi*t/rta);
			qans = 0.5 + (0.5 - ans);
		    } else {
			qans = c * (w + rt2pi*t/rta);
			ans = 0.5 + (0.5 - qans);
		    }
		    return {ans, qans};
		}
		if (std::abs(s) <= 0.4) {
		    // 270
		    // GENERAL TEMME EXPANSION
		    //
		    if ((std::abs(s) <= 2.*eps) && (a*eps*eps > 3.28e-3)) {
			return {2.0, 0.0};
		    }
		    c = std::exp(-y);
		    w = 0.5*erfc1(1, std::sqrt(y));
		    u = 1./a;
		    z = std::sqrt(z+z);
		    if (l < 1.) {
			z = -z;
		    }
		    if (ind == 0) {
			if (std::abs(s) <= 1e-3) {
			    c0 = (((((((gratio_d0[6]
					)*z+gratio_d0[5]
				       )*z+gratio_d0[4]
				      )*z+gratio_d0[3]
				     )*z+gratio_d0[2]
				    )*z+gratio_d0[1]
				   )*z+gratio_d0[0]
				  )*z - (1./3.);
                        c1 = ((((((gratio_d1[5]
                                  )*z+gratio_d1[4]
                                 )*z+gratio_d1[3]
                                )*z+gratio_d1[2]
                               )*z+gratio_d1[1]
                              )*z+gratio_d1[0]
                             )*z + d10;
                        c2 = (((((gratio_d2[4])*z+gratio_d2[3])*z+gratio_d2[2])*z+gratio_d2[1])*z+gratio_d2[0])*z + d20;
                        c3 = ((((gratio_d3[3])*z+gratio_d3[2])*z+gratio_d3[1])*z+gratio_d3[0])*z + d30;
                        c4 = (gratio_d4[1]*z+gratio_d4[0])*z + d40;
                        c5 = (gratio_d5[1]*z+gratio_d5[0])*z + d50;
                        c6 = gratio_d6[0]*z + d60;
                        t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0;
			} else {
			    c0 = (((((((((((((gratio_d0[12]
					      )*z+gratio_d0[11]
					     )*z+gratio_d0[10]
					    )*z+gratio_d0[9]
					   )*z+gratio_d0[8]
					  )*z+gratio_d0[7]
					 )*z+gratio_d0[6]
					)*z+gratio_d0[5]
				       )*z+gratio_d0[4]
				      )*z+gratio_d0[3]
				     )*z+gratio_d0[2]
				    )*z+gratio_d0[1]
				   )*z+gratio_d0[0]
				  )*z - (1./3.);
			    c1 = ((((((((((((gratio_d1[11]
					     )*z+gratio_d1[10]
					    )*z+gratio_d1[9]
					   )*z+gratio_d1[8]
					  )*z+gratio_d1[7]
					 )*z+gratio_d1[6]
					)*z+gratio_d1[5]
				       )*z+gratio_d1[4]
				      )*z+gratio_d1[3]
				     )*z+gratio_d1[2]
				    )*z+gratio_d1[1]
				   )*z+gratio_d1[0]
				  )*z + d10;
			    c2 = ((((((((((gratio_d2[9]
					   )*z+gratio_d2[8]
					  )*z+gratio_d2[7]
					 )*z+gratio_d2[6]
					)*z+gratio_d2[5]
				       )*z+gratio_d2[4]
				      )*z+gratio_d2[3]
				     )*z+gratio_d2[2]
				    )*z+gratio_d2[1]
				   )*z+gratio_d2[0]
				  )*z + d20;
			    c3 = ((((((((gratio_d3[7]
					 )*z+gratio_d3[6]
					)*z+gratio_d3[5]
				       )*z+gratio_d3[4]
				      )*z+gratio_d3[3]
				     )*z+gratio_d3[2]
				    )*z+gratio_d3[1]
				   )*z+gratio_d3[0]
				  )*z + d30;
			    c4 = ((((((gratio_d4[5]
				       )*z+gratio_d4[4]
				      )*z+gratio_d4[3]
				     )*z+gratio_d4[2]
				    )*z+gratio_d4[1]
				   )*z+gratio_d4[0]
				  )*z + d40;
			    c5 = (((gratio_d5[3]*z+gratio_d5[2])*z+gratio_d5[1])*z+gratio_d5[0])*z + d50;
			    c6 = (gratio_d6[1]*z+gratio_d6[0])*z + d60;
			    t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0;
			}
		    } else if (ind == 1) {
			c0 = ((((((gratio_d0[5]
				   )*z+gratio_d0[4]
				  )*z+gratio_d0[3]
				 )*z+gratio_d0[2]
				)*z+gratio_d0[1]
			       )*z+gratio_d0[0]
			      )*z - (1./3.);
			c1 = (((gratio_d1[3]*z+gratio_d1[2])*z+gratio_d1[1])*z+gratio_d1[0])*z + d10;
			c2 = gratio_d2[0]*z + d20;
			t = (c2*u+c1)*u + c0;
		    } else {
			t = ((gratio_d0[2]*z+gratio_d0[1])*z+gratio_d0[0])*z - (1./3.);
		    }
		    if (l < 1.) {
			ans = c * (w - rt2pi*t/rta);
			qans = 0.5 + (0.5 - ans);
		    } else {
			qans = c * (w + rt2pi*t/rta);
			ans = 0.5 + (0.5 - qans);
		    }
		    return {ans, qans};
		}
		
		t = std::pow(1 / a, 2);
		t1 = (((0.75*t-1.)*t+3.5)*t-105.0)/ (a*1260.0);
		t1 -= y;
		r = rt2pi*rta*std::exp(t1);
		
		//
		// 40
		//
		if (r == 0.0) {
		    // 420
		    if (x > a) {
			return {1.0, 0.0};
		    } else {
			return {0.0, 1.0};
		    }
		}
		
		if (x <= std::fmax(a, alog10)) {
		    // 50
		    // TAYLOR SERIES FOR P/R
		    //
		    apn = a + 1.;
		    t = x / apn;
		    wk[0] = t;
		    last_entry = 0;
		    for (n = 1; n < 20; n++) {
			apn += 1.0;
			t *= (x / apn);
			if (t <= 1e-3) { break; }
			wk[n] = t;
			last_entry = n;
		    }
		    ssum = t;
		    tol = 0.5 * acc;
		    while (1) {
			apn += 1.0;
			t *= x / apn;
			ssum += t;
			if (!(t > tol)) { break; }
		    }
		    
		    for (m = last_entry; m >= 0; m--)
			{
			    ssum += wk[m];
			}
		    ans = (r/a) * (1.0 + ssum);
		    qans = 0.5 + (0.5 - ans);
		    return {ans, qans};
		}
		
		if (x < x0) {
		    // 250
		    // CONTINUED FRACTION EXPANSION
		    //
		    tol = std::fmax(5.0*eps, acc);
		    a2nm1 = 1.0;
		    a2n = 1.0;
		    b2nm1 = x;
		    b2n = x + (1.0 - a);
		    c = 1.0;
		    
		    while (1) {
			a2nm1 = x*a2n + c*a2nm1;
			b2nm1 = x*b2n + c*b2nm1;
			am0 = a2nm1/b2nm1;
			c += 1.0;
			cma = c - a;
			a2n = a2nm1 + cma*a2n;
			b2n = b2nm1 + cma*b2n;
			an0 = a2n/b2n;
			if (!(std::fabs(an0-am0) >= tol*an0)) { break; }
		    }
		    
		    qans = r*an0;
		    ans = 0.5 + (0.5 - qans);
		    return {ans, qans};
		}
		// 100
		// ASYMPTOTIC EXPANSION
		//
		amn = a - 1.;
		t = amn / x;
		wk[0] = t;
		last_entry = 0;
		for (n = 1; n < 20; n++)
		    {
			amn -= 1.0;
			t *= amn / x;
			if (std::abs(t) <= 1e-3) { break; }
			wk[n] = t;
			// F77 code was using "n" later in the code
			last_entry = n;
		    }
		ssum = t;
		
		while (!(std::abs(t) <= acc)) {
		    amn -= 1.0;
		    t *= amn / x;
		    ssum += t;
		}
		for (m = last_entry; m >= 0; m--)
		    {
			ssum += wk[m];
		    }
		
		qans = (r/x) * (1. + ssum);
		ans = 0.5 + (0.5 - qans);
		return {ans, qans};
	    }
	    
	    twoa = a + a;
	    m = static_cast<int>(twoa);
	    
	    if (((a > x) || (x >= x0)) || (twoa != m)) {
		t1 = a*std::log(x) - x;
		r = std::exp(t1)/gamma(a);
		
		//
		// 40 Again - This time coming from 20
		//
		if (r == 0.0) {
		    if (x > a) {
			return {1.0, 0.0};
		    } else {
			return {0.0, 1.0};
		    }
		}
		
		if (x <= (a > alog10 ? a : alog10)) {
		    // 50
		    // TAYLOR SERIES FOR P/R
		    //
		    apn = a + 1.;
		    t = x / apn;
		    wk[0] = t;
		    last_entry = 0;
		    for (n =1; n < 20; n++) {
			apn += 1.0;
			t *= (x / apn);
			if (t <= 1e-3) { break; }
			wk[n] = t;
			last_entry = n;
		    }
		    ssum = t;
		    tol = 0.5 * acc;
		    while (!(t <= tol)) {
			apn += 1.0;
			t *= x / apn;
			ssum += t;
		    }
		    for (m = last_entry; m >= 0; m--)
			{
			    ssum += wk[m];
			}
		    ans = (r/a) * (1. + ssum);
		    qans = 0.5 + (0.5 - ans);
		    return {ans, qans};
		}
		
		if (x < x0) {
		    // 250
		    // CONTINUED FRACTION EXPANSION
		    //
		    tol = (5.0*eps > acc ? 5.0*eps : acc);
		    a2nm1 = 1.0;
		    a2n = 1.0;
		    b2nm1 = x;
		    b2n = x + (1.0 - a);
		    c = 1.0;
		    
		    while (1) {
			a2nm1 = x*a2n + c*a2nm1;
			b2nm1 = x*b2n + c*b2nm1;
			am0 = a2nm1/b2nm1;
			c += 1.0;
			cma = c - a;
			a2n = a2nm1 + cma*a2n;
			b2n = b2nm1 + cma*b2n;
			an0 = a2n/b2n;
			if (!(std::abs(an0-am0) >= tol*an0)) { break; }
		    }
		    
		    qans = r*an0;
		    ans = 0.5 + (0.5 - qans);
		    return {ans, qans};
		}
		// 100
		// ASYMPTOTIC EXPANSION
		//
		amn = a - 1.;
		t = amn / x;
		wk[0] = t;
		last_entry = 0;
		for (n = 1; n < 20; n++)
		    {
			amn -= 1.0;
			t *= amn / x;
			if (std::abs(t) <= 1e-3) { break; }
			wk[n] = t;
			last_entry = n;
		    }
		ssum = t;
		
		while (!(std::abs(t) <= acc)) {
		    amn -= 1.0;
		    t *= amn / x;
		    ssum += t;
		}
		for (m = last_entry; m >= 0; m--)
		    {
			ssum += wk[m];
		    }
		
		qans = (r/x) * (1. + ssum);
		ans = 0.5 + (0.5 - qans);
		return {ans, qans};
	    }
	    
	    i = m / 2;
	    if (a == i) {
		// 210
		ssum = std::exp(-x);
		t = ssum;
		n = 1;
		c = 0.;
	    } else {
		// 220
		rtx = std::sqrt(x);
		ssum = erfc1(0, rtx);
		t = std::exp(-x) / (rtpi*rtx);
		n = 0;
		c = -0.5;
	    }
	    
	    while (n < i) {
		// 230
		n += 1;
		c += 1.0;
		t *= x/c;
		ssum += t;
	    }
	    
	    qans = ssum;
	    ans = 0.5 + (0.5 - qans);
	    return {ans, qans};
	}
	
	if (a == 0.5) {
	    // 390
	    if (x >= 0.5) {
		qans = erfc1(0, std::sqrt(x));
		ans = 0.5 + (0.5 - qans);
	    } else {
		ans = erf(std::sqrt(x));
		qans = 0.5 + (0.5 - ans);
	    }
	    return {ans, qans};
	}

	if (x < 1.1) {
	    // 160
	    // TAYLOR SERIES FOR P(A,X)/X**A
	    //
	    an = 3.0;
	    c = x;
	    ssum = x / (a + 3.);
	    tol = 3.*acc / (a + 1.);
	    while (1) {
		an += 1.0;
		c *= -(x / an);
		t = c / (a + an);
		ssum += t;
		if (!(std::abs(t) > tol)) { break; }
	    }
	    j = a*x*((ssum / 6. - 0.5 / (a + 2.))*x + 1./(a + 1.));
	    z = a*std::log(x);
	    h = gam1(a);
	    g = 1. + h;
	    
	    if (((x < 0.25) && (z > -0.13394)) || (a < x/2.59)) {
		l = rexp(z);
		w = 0.5 + (0.5 + l);
		qans = (w*j - l)*g - h;
		if (qans < 0.) {
		    return {1.0, 0.0};
		}
		ans = 0.5 + (0.5 - qans);
	    } else {
		w = std::exp(z);
		ans = w*g*(0.5 + (0.5 - j));
		qans = 0.5 + (0.5 - ans);
	    }
	    return {ans, qans};
	}

	t1 = a*std::log(x) - x;
	u = a*std::exp(t1);
	if (u == 0.0) {
	    return {1.0, 0.0};
	}
	
	r = u * (1. + gam1(a));
	// 250
	// CONTINUED FRACTION EXPANSION
	//
	tol = (5.0*eps > acc ? 5.0*eps : acc);
	a2nm1 = 1.0;
	a2n = 1.0;
	b2nm1 = x;
	b2n = x + (1.0 - a);
	c = 1.0;
	
	while (1) {
	    a2nm1 = x*a2n + c*a2nm1;
	    b2nm1 = x*b2n + c*b2nm1;
	    am0 = a2nm1/b2nm1;
	    c += 1.0;
	    cma = c - a;
	    a2n = a2nm1 + cma*a2n;
	    b2n = b2nm1 + cma*b2n;
	    an0 = a2n/b2n;
	    if (!(std::abs(an0-am0) >= tol*an0)) { break; }
	}
	qans = r*an0;
	ans = 0.5 + (0.5 - qans);
	return {ans, qans};
    }

    XSF_HOST_DEVICE inline std::pair<double, int> bgrat(double a, double b, double x , double y, double w, double eps) {
	//    Asymptotic expansion for ix(a,b) when a is larger than b.
	//    The result of the expansion is added to w. It is assumed
	//    that a >= 15 And b <= 1.  Eps is the tolerance used.
	//    Ierr is a variable that reports the status of the results.
	double bp2n, cn, coef, dj, j, l, n2, q, r, s, ssum, t, t2, u, v;
	double c[30] = { 0.0 };
	double d[30] = { 0.0 };
	double bm1 = (b - 0.5) - 0.5;
	double nu = a + bm1*0.5;
	double lnx = (y > 0.375 ? std::log(x) : std::log1p(-y));
	double z = -nu*lnx;
	int i, n;

	if ((b*z) == 0.0) { return {w, 1}; }

	// COMPUTATION OF THE EXPANSION
	// SET R = EXP(-Z)*Z**B/GAMMA(B)
	r = b * (1. + gam1(b)) * std::exp(b*std::log(z));
	r *= std::exp(a*lnx) * std::exp(0.5*bm1*lnx);
	u = algdiv(b, a) + b*std::log(nu);
	u = r*std::exp(-u);
	if (u == 0.0) { return {w, 1}; }

	std::pair<double, double> ret = grat1(b, z, r, eps);
	q = ret.second;
	v = 0.25 * std::pow((1 / nu), 2);
	t2 = 0.25*lnx*lnx;
	l = w / u;
	j = q / r;
	ssum = j;
	t = 1.;
	cn = 1.;
	n2 = 0.;

	for (n = 1; n <= 30; n++) {
	    bp2n = b + n2;
	    j = (bp2n*(bp2n + 1.)*j + (z + bp2n + 1.)*t)*v;
	    n2 += 2.;
	    t *= t2;
	    cn *= 1/(n2*(n2 + 1.));
	    c[n-1] = cn;
	    s = 0.0;
	    if (n > 1) {
		coef = b - n;
		for (i = 1; i < n; i++) {
		    s += coef*c[i-1]*d[n-i-1];
		    coef += b;
		}
	    }
	    d[n-1] = bm1*cn + s/n;
	    dj = d[n-1]*j;
	    ssum += dj;
	    if (ssum <= 0.) { return {w, 1}; }
	    if (!(std::abs(dj) > eps*(ssum+l))) { break; }
	}
	return {w + u*ssum, 0};
    }

    XSF_HOST_DEVICE inline std::pair<double, double> cumgam(double x, double a) {
	//        Double precision cUMulative incomplete GAMma distribution
	//
	//
	//                            Function
	//
	//
	//    Computes   the  cumulative        of    the     incomplete   gamma
	//    distribution, i.e., the integral from 0 to X of
	//        (1/GAM(A))*EXP(-T)*T**(A-1) DT
	//    where GAM(A) is the complete gamma function of A, i.e.,
	//        GAM(A) = integral from 0 to infinity of
	//                EXP(-T)*T**(A-1) DT
	//
	//
	//                            Arguments
	//
	//
	//    X --> The upper limit of integration of the incomplete gamma.
	//                                            X is DOUBLE PRECISION
	//
	//    A --> The shape parameter of the incomplete gamma.
	//                                            A is DOUBLE PRECISION
	//
	//    CUM <-- Cumulative incomplete gamma distribution.
	//                                    CUM is DOUBLE PRECISION
	//
	//    CCUM <-- Compliment of Cumulative incomplete gamma distribution.
	//                                            CCUM is DOUBLE PRECISIO
	//
	//
	//                            Method
	//
	//
	//    Calls the routine GRATIO.
	if (x > 0.0) {
	    return gratio(a, x, 0);
	} else {
	    return {0.0, 1.0};
	}
    }
    
}

XSF_HOST_DEVICE inline double gdtria(double p, double b, double x) {
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

XSF_HOST_DEVICE inline double gdtrib(double a, double p, double x) {
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

XSF_HOST_DEVICE inline double gdtrix(double a, double b, double p) {
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
