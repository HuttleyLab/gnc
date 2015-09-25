from __future__ import division

import sys
import os
from itertools import product

from numpy import array, isnan, log, empty, diag
from numpy.linalg import slogdet, inv, eig
from numpy.random import standard_normal
from numpy.testing import assert_allclose
from scipy.optimize import nnls

def bafz(pi):
    'bound away from zero'
    pi = array(pi)
    if (pi == 0.).any():
        pi[pi == 0.] = 1e-9
        pi = pi/pi.sum()
    return pi

def calc_paralinear(J):
    # -log(det(J) / sqrt((J.sum(0)*J.sum(1)).prod())) / J.shape[0]
    summands = [0.5*log(p) for p in bafz(J.sum(1))]
    summands += [0.5*log(p) for p in bafz(J.sum(0))]
    summands.append(-slogdet(J)[1])
    pos = sum(sorted(filter(lambda x: x >= 0., summands)))
    neg = sum(sorted(filter(lambda x: x < 0., summands), reverse=True))
    p = (pos + neg)/J.shape[0]
    assert not isnan(p), 'Paralinear failed\n' + repr(J)
    return p

def fit_paralinear(J):
    b = array([calc_paralinear(J.sum(c)) for c in (2, 0, 1)])
    A = array([[1.,1.,0.],[0.,1.,1.],[1.,0.,1.]])
    return nnls(A, b)[0]

def funny_P(J, a, b, c, gamma):
    p = bafz(J.sum(axis=(b, c)))
    ix = [None]*3
    ix[c] = gamma
    P = empty((J.shape[0],)*2)
    for i, j in product(*[range(J.shape[0])]*2):
        ix[a], ix[b] = i, j
        P[i,j] = J[tuple(ix)]/p[i]
    return P

def DLC_permutation(P):
    R = array([r == r.max() for r in P.T])
    if not (R.sum(0) == 1).all():
        print '**** FAILED DLC ****'
        print repr(P)
        print '**** FAILED DLC ****'
    assert (R.sum(0) == 1).all(), 'no DLC permutation of P exists' 
    return R

def pairwise_J(J, a, b, c):
    Jab = J.sum(c)
    if b < a:
        return Jab.T
    return Jab

def pairwise_P(J, a, b, c):
    Jab = pairwise_J(J, a, b, c)
    return diag(1./bafz(Jab.sum(1))).dot(Jab)

def triad(J):
    # figure out the paralinear branch lengths
    a, c, b = fit_paralinear(J).argsort()

    # use spectral method to get Par for shortest branch to minimise error
    U = standard_normal(J.shape[0])
    #U = array([1., 1., 1., 1.])
    Pabgammas = [funny_P(J, a, b, c, gamma) for gamma in range(J.shape[0])]
    PabU = sum(u*P for u, P in zip(U, Pabgammas))
    Jab = pairwise_J(J, a, b, c)
    pia = bafz(Jab.sum(1))
    Pab = diag(1/pia).dot(Jab)
    w, X = eig(PabU.dot(inv(Pab)))
    eta = inv(X).sum(1) # some mossel and roch magic here
    Par = X.dot(diag(eta))
    Par[Par < 0.] = 0.
    Par = diag(1./bafz(Par.sum(1))).dot(Par)

    # use Par to calculate Prb and Prc
    Prb = inv(Par).dot(Pab)
    Pac = pairwise_P(J, a, c, b)
    Prc = inv(Par).dot(Pac)

    assert_allclose(Par.dot(diag(Prc.dot(U))).dot(Prb), PabU)

    # calculate pir and Pra from Par and pia
    Jar = diag(pia).dot(Par)
    pir = bafz(Jar.sum(0))
    Pra = diag(1/pir).dot(Jar.T)

    # now reconstruct the rows
    R = DLC_permutation(Pra)
    Ps = [R.dot(P) for P in array((Pra, Prb, Prc))[array([a, b, c]).argsort()]]
    pir = R.dot(pir)

    for P in Ps:
        assert (P >= 0.).all(), 'P not stochastic\n' + repr(P)
        assert_allclose(P.sum(1), 1.), 'P not stochastic\n' + repr(P)

    return pir, Ps

def main():
    pass

if __name__ == '__main__':
    sys.exit(main())
