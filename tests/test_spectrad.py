from __future__ import division

from itertools import product

from numpy import random, eye, array, empty, diag, allclose, log, sqrt
from numpy.testing import assert_allclose
from scipy.linalg import expm, det

import lib
import spectrad


def gen_Q(n, pi):
    Q = random.random((n,n))
    ix = eye(n, dtype=bool)
    Q[ix] = -Q.sum(1) + Q[ix]
    Q /= -Q.diagonal().dot(pi)
    return Q

def gen_J(n=4, t=0.25):
    pi = random.rand(n)
    pi /= pi.sum()
    order = range(1,4)
    random.shuffle(order)
    Ps = map(expm, [t*i*gen_Q(n, pi) for i in order])
    r = range(n)
    J = empty((n,)*3)
    for i, j, k in product(*[r]*3):
        J[i, j, k] = pi.dot(Ps[0][:,i]*Ps[1][:,j]*Ps[2][:,k])
    return J, pi, Ps

def wtf_paralinear(J):
    return -log(det(J) / sqrt((J.sum(0)*J.sum(1)).prod())) / J.shape[0]

def test_fit_paralinear():
    'fit_paralinear should accurately recover edgewise paralinear distances'
    for n in (4, 61):
        J, pi, Ps = gen_J(n)
        in_ds = [spectrad.calc_paralinear(diag(pi).dot(P)) for P in Ps]
        out_ds = spectrad.fit_paralinear(J)
        assert_allclose(in_ds, out_ds)

def test_triad():
    for n in (4, 61):
        J, pi, Ps = gen_J(n)
        fpi, fPs = spectrad.triad(J)
        assert_allclose(fpi, pi, rtol=1., atol=1e-4)
        for P, fP in zip(Ps, fPs):
            assert_allclose(P, fP, rtol=1., atol=1e-4)
