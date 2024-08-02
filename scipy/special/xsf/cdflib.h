#include "config.h"
#include "error.h"

#include "tools.h"
#include "digamma.h"
#include "gamma.h"
#include "cephes/polevl.h"
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
	/* Included only for benchmarking, prefer xsf::gammaln */

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
	    return aug + log(x);
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
		return gamln(a) + algdiv(a,b);
	    } else {
		return gamln(a) + (gamln(b) - gamln(a+b));
	    }
	}
	
	if (a <= 2) {
	    if (b <= 2) {
		return gamln(a) + gamln(b) - gsumln(a, b);
	    }
	    if (b >= 8) {
		return gamln(a) + algdiv(a, b);
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
		    return w + gamln(a) + algdiv(a, b);
		}
	    } else {
		n = static_cast<int>(a - 1.);
		w = 1.0;
		for (i = 0; i < n; i++) {
		    a -= 1.0;
		    w *= a/(1. + (a/b));
		}
		return (std::log(w) - n*std::log(b)) + (gamln(a) + algdiv(a, b));
	    }
	}
	n = static_cast<int>(b - 1.);
	z = 1.0;
	for (i = 0; i < n; i++) {
	    b -= 1.0;
	    z *= b / (a + b);
	}
	return w + std::log(z) + (gamln(a) + gamln(b) - gsumln(a, b));
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
