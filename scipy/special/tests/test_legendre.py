import functools
import itertools
import operator
import platform
import sys

import numpy as np
from numpy import (array, isnan, r_, arange, finfo, pi, sin, cos, tan, exp,
        log, zeros, sqrt, asarray, inf, nan_to_num, real, arctan, double,
        array_equal)

import pytest
from pytest import raises as assert_raises
from numpy.testing import (assert_equal, assert_almost_equal,
        assert_array_equal, assert_array_almost_equal, assert_approx_equal,
        assert_, assert_allclose, assert_array_almost_equal_nulp,
        suppress_warnings)

from scipy import special

import math

# Base polynomials come from Abrahmowitz and Stegan
class TestLegendre:
    def test_legendre(self):
        leg0 = special.legendre(0)
        leg1 = special.legendre(1)
        leg2 = special.legendre(2)
        leg3 = special.legendre(3)
        leg4 = special.legendre(4)
        leg5 = special.legendre(5)
        assert_equal(leg0.c, [1])
        assert_equal(leg1.c, [1,0])
        assert_almost_equal(leg2.c, array([3,0,-1])/2.0, decimal=13)
        assert_almost_equal(leg3.c, array([5,0,-3,0])/2.0)
        assert_almost_equal(leg4.c, array([35,0,-30,0,3])/8.0)
        assert_almost_equal(leg5.c, array([63,0,-70,0,15,0])/8.0)

    @pytest.mark.parametrize('n', [1, 2, 3, 4, 5])
    @pytest.mark.parametrize('zr', [0.5241717, 12.80232, -9.699001,
                                    0.5122437, 0.1714377])
    @pytest.mark.parametrize('zi', [9.766818, 0.2999083, 8.24726, -22.84843,
                                    -0.8792666])
    def test_lpn_against_clpmn(self, n, zr, zi):
        reslpn = special.lpn(n, zr + zi*1j)
        resclpmn = special.clpmn(0, n, zr+zi*1j)
        assert_allclose(reslpn[0], resclpmn[0][0])
        assert_allclose(reslpn[1], resclpmn[1][0])


class TestLegendreFunctions:
    def test_clpmn(self):
        z = 0.5+0.3j
        clp = special.clpmn(2, 2, z, 3)
        assert_array_almost_equal(clp,
                   (array([[1.0000, z, 0.5*(3*z*z-1)],
                           [0.0000, sqrt(z*z-1), 3*z*sqrt(z*z-1)],
                           [0.0000, 0.0000, 3*(z*z-1)]]),
                    array([[0.0000, 1.0000, 3*z],
                           [0.0000, z/sqrt(z*z-1), 3*(2*z*z-1)/sqrt(z*z-1)],
                           [0.0000, 0.0000, 6*z]])),
                    7)

    def test_clpmn_close_to_real_2(self):
        eps = 1e-10
        m = 1
        n = 3
        x = 0.5
        clp_plus = special.clpmn(m, n, x+1j*eps, 2)[0][m, n]
        clp_minus = special.clpmn(m, n, x-1j*eps, 2)[0][m, n]
        assert_array_almost_equal(array([clp_plus, clp_minus]),
                                  array([special.lpmv(m, n, x),
                                         special.lpmv(m, n, x)]),
                                  7)

    def test_clpmn_close_to_real_3(self):
        eps = 1e-10
        m = 1
        n = 3
        x = 0.5
        clp_plus = special.clpmn(m, n, x+1j*eps, 3)[0][m, n]
        clp_minus = special.clpmn(m, n, x-1j*eps, 3)[0][m, n]
        assert_array_almost_equal(array([clp_plus, clp_minus]),
                                  array([special.lpmv(m, n, x)*np.exp(-0.5j*m*np.pi),
                                         special.lpmv(m, n, x)*np.exp(0.5j*m*np.pi)]),
                                  7)

    def test_clpmn_across_unit_circle(self):
        eps = 1e-7
        m = 1
        n = 1
        x = 1j
        for type in [2, 3]:
            assert_almost_equal(special.clpmn(m, n, x+1j*eps, type)[0][m, n],
                            special.clpmn(m, n, x-1j*eps, type)[0][m, n], 6)

    def test_inf(self):
        for z in (1, -1):
            for n in range(4):
                for m in range(1, n):
                    lp = special.clpmn(m, n, z)
                    assert_(np.isinf(lp[1][1,1:]).all())
                    lp = special.lpmn(m, n, z)
                    assert_(np.isinf(lp[1][1,1:]).all())

    def test_deriv_clpmn(self):
        # data inside and outside of the unit circle
        zvals = [0.5+0.5j, -0.5+0.5j, -0.5-0.5j, 0.5-0.5j,
                 1+1j, -1+1j, -1-1j, 1-1j]
        m = 2
        n = 3
        for type in [2, 3]:
            for z in zvals:
                for h in [1e-3, 1e-3j]:
                    approx_derivative = (special.clpmn(m, n, z+0.5*h, type)[0]
                                         - special.clpmn(m, n, z-0.5*h, type)[0])/h
                    assert_allclose(special.clpmn(m, n, z, type)[1],
                                    approx_derivative,
                                    rtol=1e-4)

    @pytest.mark.parametrize("shape", [(10,), (4, 9), (3, 5, 7, 10)])
    def test_lpmn(self, shape):
        rng = np.random.default_rng(1234)

        n = rng.integers(0, 10, shape)
        m = rng.integers(-10, 10, shape)
        x = rng.uniform(-0.99, 0.99, shape)

        p, p_jac, p_hess = special.lpmn(m, n, x, diff_n = 2, legacy = False)

        assert p.shape == shape
        assert p_jac.shape == p.shape
        assert p_hess.shape == p_jac.shape

        np.testing.assert_allclose((1 - x * x) * p_hess, 2 * x * p_jac - (n * (n + 1) - m * m / (1 - x * x)) * p, rtol = 1e-05, atol = 1e-08)

    @pytest.mark.parametrize("shape", [(10,), (4, 9), (3, 5, 7)])
    def test_lpmn_all(self, shape):
        rng = np.random.default_rng(1234)

        n_max = 20
        m_max = 20

        x = rng.uniform(-0.99, 0.99, shape)

        p, p_jac, p_hess = special.lpmn_all(m_max, n_max, x, diff_n = 2)

        m = np.concatenate([np.arange(m_max + 1), np.arange(-m_max, 0)])
        n = np.arange(n_max + 1)

        m = np.expand_dims(m, axis = tuple(range(1, x.ndim + 2)))
        n = np.expand_dims(n, axis = (0,) + tuple(range(2, x.ndim + 2)))
        np.testing.assert_allclose((1 - x * x) * p_hess, 2 * x * p_jac - (n * (n + 1) - m * m / (1 - x * x)) * p, rtol = 1e-05, atol = 1e-08)

    @pytest.mark.parametrize("shape", [(10,), (4, 9), (3, 5, 7)])
    def test_lpmn_all_exact(self, shape):
        rng = np.random.default_rng(1234)

        x = rng.uniform(-5, 5, shape)

        type = np.where(np.abs(x) <= 1, 2, 3)

        p, p_jac = special.lpmn_all(4, 4, x, diff_n = 1)

        np.testing.assert_allclose(p[0, 0], lpmn_ref(0, 0, x, type = type))

        np.testing.assert_allclose(p[0, 1], lpmn_ref(0, 1, x, type = type))
        np.testing.assert_allclose(p[1, 1], lpmn_ref(1, 1, x, type = type))
        np.testing.assert_allclose(p[-1, 1], lpmn_ref(-1, 1, x, type = type))

        np.testing.assert_allclose(p[0, 2], lpmn_ref(0, 2, x, type = type))
        np.testing.assert_allclose(p[1, 2], lpmn_ref(1, 2, x, type = type))
        np.testing.assert_allclose(p[2, 2], lpmn_ref(2, 2, x, type = type))
        np.testing.assert_allclose(p[-2, 2], lpmn_ref(-2, 2, x, type = type))
        np.testing.assert_allclose(p[-1, 2], lpmn_ref(-1, 2, x, type = type))

        np.testing.assert_allclose(p[0, 3], lpmn_ref(0, 3, x, type = type))
        np.testing.assert_allclose(p[1, 3], lpmn_ref(1, 3, x, type = type))
        np.testing.assert_allclose(p[2, 3], lpmn_ref(2, 3, x, type = type))
        np.testing.assert_allclose(p[3, 3], lpmn_ref(3, 3, x, type = type))
        np.testing.assert_allclose(p[-3, 3], lpmn_ref(-3, 3, x, type = type))
        np.testing.assert_allclose(p[-2, 3], lpmn_ref(-2, 3, x, type = type))
        np.testing.assert_allclose(p[-1, 3], lpmn_ref(-1, 3, x, type = type))

        np.testing.assert_allclose(p[0, 4], lpmn_ref(0, 4, x, type = type))
        np.testing.assert_allclose(p[1, 4], lpmn_ref(1, 4, x, type = type))
        np.testing.assert_allclose(p[2, 4], lpmn_ref(2, 4, x, type = type))
        np.testing.assert_allclose(p[3, 4], lpmn_ref(3, 4, x, type = type))
        np.testing.assert_allclose(p[4, 4], lpmn_ref(4, 4, x, type = type))
        np.testing.assert_allclose(p[-4, 4], lpmn_ref(-4, 4, x, type = type))
        np.testing.assert_allclose(p[-3, 4], lpmn_ref(-3, 4, x, type = type))
        np.testing.assert_allclose(p[-2, 4], lpmn_ref(-2, 4, x, type = type))
        np.testing.assert_allclose(p[-1, 4], lpmn_ref(-1, 4, x, type = type))

        np.testing.assert_allclose(p_jac[0, 1], lpmn_jac(0, 1, x, type = type))
        np.testing.assert_allclose(p_jac[1, 1], lpmn_jac(1, 1, x, type = type))
        np.testing.assert_allclose(p_jac[-1, 1], lpmn_jac(-1, 1, x, type = type))

        np.testing.assert_allclose(p_jac[0, 2], lpmn_jac(0, 2, x, type = type))
        np.testing.assert_allclose(p_jac[1, 2], lpmn_jac(1, 2, x, type = type))
        np.testing.assert_allclose(p_jac[2, 2], lpmn_jac(2, 2, x, type = type))
        np.testing.assert_allclose(p_jac[-2, 2], lpmn_jac(-2, 2, x, type = type))
        np.testing.assert_allclose(p_jac[-1, 2], lpmn_jac(-1, 2, x, type = type))

        np.testing.assert_allclose(p_jac[0, 3], lpmn_jac(0, 3, x, type = type))
        np.testing.assert_allclose(p_jac[1, 3], lpmn_jac(1, 3, x, type = type))
        np.testing.assert_allclose(p_jac[2, 3], lpmn_jac(2, 3, x, type = type))
        np.testing.assert_allclose(p_jac[3, 3], lpmn_jac(3, 3, x, type = type))
        np.testing.assert_allclose(p_jac[-3, 3], lpmn_jac(-3, 3, x, type = type))
        np.testing.assert_allclose(p_jac[-2, 3], lpmn_jac(-2, 3, x, type = type))
        np.testing.assert_allclose(p_jac[-1, 3], lpmn_jac(-1, 3, x, type = type))

    def test_lpmn_legacy(self):
        lp = special.lpmn(0, 2, .5)
        assert_array_almost_equal(lp,(array([[1.00000,
                                                      0.50000,
                                                      -0.12500]]),
                                      array([[0.00000,
                                                      1.00000,
                                                      1.50000]])), 4)

    @pytest.mark.parametrize("shape", [(1000,), (4, 9), (3, 5, 7)])
    @pytest.mark.parametrize("type", [2, 3])
    def test_clpmn_all_exact(self, shape, type):
        rng = np.random.default_rng(1234)

        z = rng.uniform(-100, 100, shape) + 1j * rng.uniform(-100, 100, shape)

#        p, p_jac = special.clpmn(4, 4, z, type = type, legacy = False)
        p, p_jac = special.clpmn_all(4, 4, type, z, n_diff = 1)

        np.testing.assert_allclose(p[0, 0], lpmn_ref(0, 0, z, type = type))

        np.testing.assert_allclose(p[0, 1], lpmn_ref(0, 1, z, type = type))
        np.testing.assert_allclose(p[1, 1], lpmn_ref(1, 1, z, type = type))
        np.testing.assert_allclose(p[-1, 1], lpmn_ref(-1, 1, z, type = type))

        np.testing.assert_allclose(p[0, 2], lpmn_ref(0, 2, z, type = type))
        np.testing.assert_allclose(p[1, 2], lpmn_ref(1, 2, z, type = type))
        np.testing.assert_allclose(p[2, 2], lpmn_ref(2, 2, z, type = type))
        np.testing.assert_allclose(p[-2, 2], lpmn_ref(-2, 2, z, type = type))
        np.testing.assert_allclose(p[-1, 2], lpmn_ref(-1, 2, z, type = type))
 
        np.testing.assert_allclose(p[0, 3], lpmn_ref(0, 3, z, type = type))
        np.testing.assert_allclose(p[1, 3], lpmn_ref(1, 3, z, type = type))
        np.testing.assert_allclose(p[2, 3], lpmn_ref(2, 3, z, type = type))
        np.testing.assert_allclose(p[3, 3], lpmn_ref(3, 3, z, type = type))
        np.testing.assert_allclose(p[-3, 3], lpmn_ref(-3, 3, z, type = type))
        np.testing.assert_allclose(p[-2, 3], lpmn_ref(-2, 3, z, type = type))
        np.testing.assert_allclose(p[-1, 3], lpmn_ref(-1, 3, z, type = type))

        np.testing.assert_allclose(p[0, 4], lpmn_ref(0, 4, z, type = type))
        np.testing.assert_allclose(p[1, 4], lpmn_ref(1, 4, z, type = type))
        np.testing.assert_allclose(p[2, 4], lpmn_ref(2, 4, z, type = type))
        np.testing.assert_allclose(p[3, 4], lpmn_ref(3, 4, z, type = type))
        np.testing.assert_allclose(p[4, 4], lpmn_ref(4, 4, z, type = type))
        np.testing.assert_allclose(p[-4, 4], lpmn_ref(-4, 4, z, type = type))
        np.testing.assert_allclose(p[-3, 4], lpmn_ref(-3, 4, z, type = type))
        np.testing.assert_allclose(p[-2, 4], lpmn_ref(-2, 4, z, type = type))
        np.testing.assert_allclose(p[-1, 4], lpmn_ref(-1, 4, z, type = type))

        np.testing.assert_allclose(p_jac[0, 0], lpmn_jac(0, 0, z, type = type))

        np.testing.assert_allclose(p_jac[0, 1], lpmn_jac(0, 1, z, type = type))
        np.testing.assert_allclose(p_jac[1, 1], lpmn_jac(1, 1, z, type = type))
        np.testing.assert_allclose(p_jac[-1, 1], lpmn_jac(-1, 1, z, type = type))

        np.testing.assert_allclose(p_jac[0, 2], lpmn_jac(0, 2, z, type = type))
        np.testing.assert_allclose(p_jac[1, 2], lpmn_jac(1, 2, z, type = type))
        np.testing.assert_allclose(p_jac[2, 2], lpmn_jac(2, 2, z, type = type))
        np.testing.assert_allclose(p_jac[-2, 2], lpmn_jac(-2, 2, z, type = type))
        np.testing.assert_allclose(p_jac[-1, 2], lpmn_jac(-1, 2, z, type = type))

        np.testing.assert_allclose(p_jac[0, 3], lpmn_jac(0, 3, z, type = type))
        np.testing.assert_allclose(p_jac[1, 3], lpmn_jac(1, 3, z, type = type))
        np.testing.assert_allclose(p_jac[2, 3], lpmn_jac(2, 3, z, type = type))
        np.testing.assert_allclose(p_jac[3, 3], lpmn_jac(3, 3, z, type = type))
        np.testing.assert_allclose(p_jac[-3, 3], lpmn_jac(-3, 3, z, type = type))
        np.testing.assert_allclose(p_jac[-2, 3], lpmn_jac(-2, 3, z, type = type))
        np.testing.assert_allclose(p_jac[-1, 3], lpmn_jac(-1, 3, z, type = type))

    @pytest.mark.parametrize("shape", [(10,), (4, 9), (3, 5, 7)])
    def test_lpn(self, shape):
        rng = np.random.default_rng(1234)

        n = rng.integers(0, 100, shape)
        x = rng.uniform(-1, 1, shape)

        p, p_jac, p_hess = special.lpn(n, x, diff_n = 2, legacy = False)

        assert p.shape == shape
        assert p_jac.shape == p.shape
        assert p_hess.shape == p_jac.shape

        err = (1 - x * x) * p_hess - 2 * x * p_jac + n * (n + 1) * p
        np.testing.assert_allclose(err, 0, atol = 1e-10)

    @pytest.mark.parametrize("n_max", [1, 2, 4, 8, 16, 32])
    @pytest.mark.parametrize("x_shape", [(10,), (4, 9), (3, 5, 7)])
    def test_lpn_all(self, n_max, x_shape):
        rng = np.random.default_rng(1234)

        x = rng.uniform(-1, 1, x_shape)
        p, p_jac, p_hess = special.lpn_all(n_max, x, diff_n = 2)

        n = np.arange(n_max + 1)
        n = np.expand_dims(n, axis = tuple(range(1, x.ndim + 1)))

        assert p.shape == (len(n),) + x.shape
        assert p_jac.shape == p.shape
        assert p_hess.shape == p_jac.shape

        err = (1 - x * x) * p_hess - 2 * x * p_jac + n * (n + 1) * p
        np.testing.assert_allclose(err, 0, atol = 1e-10)

    def test_lpn_legacy(self):
        p, pd = special.lpn(2, 0.5)
        assert_array_almost_equal(p, [1.00000, 0.50000, -0.12500], 4)
        assert_array_almost_equal(pd, [0.00000, 1.00000, 1.50000], 4)

    def test_lpmv(self):
        lp = special.lpmv(0,2,.5)
        assert_almost_equal(lp,-0.125,7)
        lp = special.lpmv(0,40,.001)
        assert_almost_equal(lp,0.1252678976534484,7)

        # XXX: this is outside the domain of the current implementation,
        #      so ensure it returns a NaN rather than a wrong answer.
        with np.errstate(all='ignore'):
            lp = special.lpmv(-1,-1,.001)
        assert_(lp != 0 or np.isnan(lp))

    def test_lqmn(self):
        lqmnf = special.lqmn(0,2,.5)
        lqf = special.lqn(2,.5)
        assert_array_almost_equal(lqmnf[0][0],lqf[0],4)
        assert_array_almost_equal(lqmnf[1][0],lqf[1],4)

    def test_lqmn_gt1(self):
        """algorithm for real arguments changes at 1.0001
           test against analytical result for m=2, n=1
        """
        x0 = 1.0001
        delta = 0.00002
        for x in (x0-delta, x0+delta):
            lq = special.lqmn(2, 1, x)[0][-1, -1]
            expected = 2/(x*x-1)
            assert_almost_equal(lq, expected)

    def test_lqmn_shape(self):
        a, b = special.lqmn(4, 4, 1.1)
        assert_equal(a.shape, (5, 5))
        assert_equal(b.shape, (5, 5))

        a, b = special.lqmn(4, 0, 1.1)
        assert_equal(a.shape, (5, 1))
        assert_equal(b.shape, (5, 1))

    def test_lqn(self):
        lqf = special.lqn(2,.5)
        assert_array_almost_equal(lqf,(array([0.5493, -0.7253, -0.8187]),
                                       array([1.3333, 1.216, -0.8427])),4)

    @pytest.mark.parametrize("function", [special.lpn, special.lqn])
    @pytest.mark.parametrize("n", [1, 2, 4, 8, 16, 32])
    @pytest.mark.parametrize("z_complex", [False, True])
    @pytest.mark.parametrize("z_inexact", [False, True])
    @pytest.mark.parametrize(
        "input_shape",
        [
            (), (1, ), (2, ), (2, 1), (1, 2), (2, 2), (2, 2, 1), (2, 2, 2)
        ]
    )
    def test_array_inputs_lxn(self, function, n, z_complex, z_inexact, input_shape):
        """Tests for correct output shapes."""
        rng = np.random.default_rng(1234)
        if z_inexact:
            z = rng.integers(-3, 3, size=input_shape)
        else:
            z = rng.uniform(-1, 1, size=input_shape)

        if z_complex:
            z = 1j * z + 0.5j * z

        P_z, P_d_z = function(n, z)
        assert P_z.shape == (n + 1, ) + input_shape
        assert P_d_z.shape == (n + 1, ) + input_shape

    @pytest.mark.parametrize("function", [special.lqmn])
    @pytest.mark.parametrize(
        "m,n",
        [(0, 1), (1, 2), (1, 4), (3, 8), (11, 16), (19, 32)]
    )
    @pytest.mark.parametrize("z_inexact", [False, True])
    @pytest.mark.parametrize(
        "input_shape", [
            (), (1, ), (2, ), (2, 1), (1, 2), (2, 2), (2, 2, 1)
        ]
    )
    def test_array_inputs_lxmn(self, function, m, n, z_inexact, input_shape):
        """Tests for correct output shapes and dtypes."""
        rng = np.random.default_rng(1234)
        if z_inexact:
            z = rng.integers(-3, 3, size=input_shape)
        else:
            z = rng.uniform(-1, 1, size=input_shape)

        P_z, P_d_z = function(m, n, z)
        assert P_z.shape == (m + 1, n + 1) + input_shape
        assert P_d_z.shape == (m + 1, n + 1) + input_shape

    @pytest.mark.parametrize("function", [special.clpmn, special.lqmn])
    @pytest.mark.parametrize(
        "m,n",
        [(0, 1), (1, 2), (1, 4), (3, 8), (11, 16), (19, 32)]
    )
    @pytest.mark.parametrize(
        "input_shape", [
            (), (1, ), (2, ), (2, 1), (1, 2), (2, 2), (2, 2, 1)
        ]
    )
    def test_array_inputs_clxmn(self, function, m, n, input_shape):
        """Tests for correct output shapes and dtypes."""
        rng = np.random.default_rng(1234)
        z = rng.uniform(-1, 1, size=input_shape)
        z = 1j * z + 0.5j * z

        P_z, P_d_z = function(m, n, z)
        assert P_z.shape == (m + 1, n + 1) + input_shape
        assert P_d_z.shape == (m + 1, n + 1) + input_shape

def lpmn_ref(m, n, z, *, type = 2):
    ls = np.where(type == 3, -1, 1)
    qs = np.where(type == 3, np.sign(np.real(z)), 1)

    if (n == 0):
        if (m == 0):
            return np.ones_like(z)

    if (n == 1):
        if (m == 0):
            return z

        if (m == 1):
            return qs * -ls * np.sqrt(ls * (1 - z * z))

        if (m == -1):
            return qs * np.sqrt(ls * (1 - z * z)) / 2

    if (n == 2):
        if (m == 0):
            return (3 * z * z - 1) / 2

        if (m == 1):
            return -3 * ls * qs * z * np.sqrt(ls * (1 - z * z))

        if (m == 2):
            return 3 * ls * (1 - z * z)

        if (m == -2):
            return ls * (1 - z * z) / 8

        if (m == -1):
            return qs * z * np.sqrt(ls * (1 - z * z)) / 2

    if (n == 3):
        if (m == 0):
            return (5 * z * z - 3) * z / 2

        if (m == 1):
            return 3 * (1 - 5 * z * z) * ls * qs * np.sqrt(ls * (1 - z * z)) / 2

        if (m == 2):
            return 15 * z * ls * (1 - z * z)

        if (m == 3):
            return -15 * qs * (1 - z * z) * np.sqrt(ls * (1 - z * z))

        if (m == -3):
            return (1 - z * z) * qs * ls * np.sqrt(ls * (1 - z * z)) / 48

        if (m == -2):
            return z * ls * (1 - z * z) / 8

        if (m == -1):
            return -qs * (1 - 5 * z * z) * np.sqrt(ls * (1 - z * z)) / 8

    if (n == 4):
        if (m == 0):
            return ((35 * z * z - 30) * z * z + 3) / 8

        if (m == 1):
            return 5 * ls * qs * z * np.sqrt(ls * (1 - z * z)) * (3 - 7 * z * z) / 2

        if (m == 2):
            return -15 * ls * (1 - (8 - 7 * z * z) * z * z) / 2

        if (m == 3):
            return -105 * z * qs * (1 - z * z) * np.sqrt(ls * (1 - z * z))

        if (m == 4):
            return 105 * np.square(z * z - 1)

        if (m == -1):
            return -z * qs * np.sqrt(ls * (1 - z * z)) * (3 - 7 * z * z) / 8

        if (m == -2):
            return -ls * (1 - (8 - 7 * z * z) * z * z) / 48

        if (m == -3):
            return z * ls * qs * (1 - z * z) * np.sqrt(ls * (1 - z * z)) / 48

        if (m == -4):
            return np.square(z * z - 1) / 384

    raise NotImplementedError

def lpmn_jac(m, n, z, type):
    ls = np.where(type == 3, -1, 1)
    qs = np.where(type == 3, np.sign(np.real(z)), 1)

    if (n == 0):
        if (m == 0):
            return np.zeros_like(z)

    if (n == 1):
        if (m == 0):
            return np.ones_like(z)

        if (m == 1):
            return qs * z / np.sqrt(ls * (1 - z * z))

        if (m == -1):
            return -qs * ls * z / (2 * np.sqrt(ls * (1 - z * z)))

    if (n == 2):
        if (m == 0):
            return 3 * z

        if (m == 1):
            return -qs * 3 * (1 - 2 * z * z) / np.sqrt(ls * (1 - z * z))

        if (m == 2):
            return ls * -6 * z

        if (m == -2):
            return ls * -z / 4

        if (m == -1):
            return ls * qs * (1 - 2 * z * z) / (2 * np.sqrt(ls * (1 - z * z)))

    if (n == 3):
        if (m == 0):
            return 3 * (5 * z * z - 1) / 2

        if (m == 1):
            return -qs * 3 * (11 - 15 * z * z) * z / (2 * np.sqrt(ls * (1 - z * z)))

        if (m == 2):
            return 15 * ls * (1 - 3 * z * z)

        if (m == 3):
            return 45 * z * qs * np.sqrt(ls * (1 - z * z))

        if (m == -3):
            return -z * qs * ls * np.sqrt(ls * (1 - z * z)) / 16

        if (m == -2):
            return ls * (1 - 3 * z * z) / 8

        if (m == -1):
            return qs * ls * (11 - 15 * z * z) * z / (8 * np.sqrt(ls * (1 - z * z)))

    raise NotImplementedError