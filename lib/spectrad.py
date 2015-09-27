from __future__ import division

import sys
import os
from itertools import product

from numpy import array, isnan, log, empty, diag, isinf, logical_not,\
        real_if_close, eye, allclose, sign, ones
from numpy.linalg import slogdet, inv, eig
from numpy.random import standard_normal
from numpy.testing import assert_allclose
from scipy.optimize import nnls
from cvxopt import matrix, solvers
solvers.options['show_progress'] = False

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2015, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPLv3 or any later version'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Development'
__version__ = '0.0.1-dev'

def bafz(pi):
    'bound away from zero'
    pi = array(pi)
    if (pi == 0.).any():
        pi[pi == 0.] = 1e-9
        pi = pi/pi.sum()
    return pi

def calc_paralinear(J):
    # -log(det(J) / sqrt((J.sum(0)*J.sum(1)).prod())) / J.shape[0]
    summands = [0.5*log(p) for p in J.sum(1)]
    summands += [0.5*log(p) for p in J.sum(0)]
    summands.append(-slogdet(J)[1])
    pos = sum(sorted(filter(lambda x: x >= 0., summands)))
    neg = sum(sorted(filter(lambda x: x < 0., summands), reverse=True))
    p = (pos + neg)/J.shape[0]
    assert not isnan(p), 'Paralinear failed\n' + repr(J)
    return p

def fix_J(J):
    J = array(J)
    pc = 0.5*J[J != 0.].min()
    k = len(J.shape)
#    sums = array([J.sum(tuple(range(i)+range(i+1,k))) for i in range(k)]).T
#    for i in range(J.shape[0]):
#        if (sums[i] == 0.).any():
#            J[(i,)*k] = pc
    for i in range(J.shape[0]):
        if J[(i,)*k] == 0.:
            J[(i,)*k] = pc
    J /= J.sum()
    return J

def fit_paralinear(J):
    b = array([calc_paralinear(J.sum(c)) for c in (2, 0, 1)])
    if sum(isinf(b)) == 3:
        return array([1., 1., 1.])
    else:
        b[isinf(b)] = 2.*b[logical_not(isinf(b))].max()
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

def DLCesque_permutation(P):
    m = P.shape[0]
    not_done_rows = range(m)
    not_done_cols = range(m)
    column_map = [None]*m
    while not_done_rows:
        candidates = [[] for i in range(m)]
        for i in not_done_cols:
            col = P[not_done_rows,i]
            preference = col.argmax()
            colmax = col[preference]
            preference = not_done_rows[preference]
            if len(col[col < colmax]) > 0.:
                weight = colmax - col[col < colmax].max()
            else:
                weight = 0.
            candidates[preference].append((weight, i))
        
        for row, candidate in enumerate(candidates):
            if candidate:
                column_map[row] = max(candidate)[1]
                not_done_rows.remove(row)
                not_done_cols.remove(column_map[row])

    return eye(m)[:,column_map]

def pairwise_J(J, a, b, c):
    Jab = J.sum(c)
    if b < a:
        return Jab.T
    return Jab

def pairwise_P(J, a, b, c):
    Jab = pairwise_J(J, a, b, c)
    return diag(1./bafz(Jab.sum(1))).dot(Jab)

def is_stochastic(P):
    return (P >= 0.).all() and allclose(P.sum(1), 1.)

def _calcQq(Ps):
    Ps = array(Ps)
    m = Ps.shape[2]**Ps.shape[0]
    n = Ps.shape[1]
    A = empty((m,n))
    b = empty(m)
    dx = tuple(range(Ps.shape[0]))
    for i, ix in enumerate(product(*[range(Ps.shape[2])]*Ps.shape[0])):
        A[i] = Ps[dx,:,ix].prod(0)
        b[i] = J[ix]
    Q = A.T.dot(A)
    q = -b.T.dot(A)
    return Q, q

def best_pi(J, Ps):
    def _calcQq(Ps):
        Ps = array(Ps)
        m = Ps.shape[2]**Ps.shape[0]
        n = Ps.shape[1]
        A = empty((m,n))
        b = empty(m)
        dx = tuple(range(Ps.shape[0]))
        for i, ix in enumerate(product(*[range(Ps.shape[2])]*Ps.shape[0])):
            A[i] = Ps[dx,:,ix].prod(0)
            b[i] = J[ix]
        Q = A.T.dot(A)
        q = -b.T.dot(A)
        return Q, q

    Q, q = map(matrix, _calcQq(Ps))
    n = Ps[0].shape[0]
    G = matrix(-eye(n))
    h = matrix(array([-1e-9]*n))
    A = matrix(ones(n).reshape((1,n)))
    b = matrix([1.])
    initvals = {'x' : matrix([1./n]*n)}
    sol = solvers.qp(Q, q, G, h, A, b, initvals=initvals)
    assert sol['status'] == 'optimal'
    return array(sol['x'].T)[0]

def triad(J):
    assert len(J.shape) == 3, 'triad only works for three taxa'
    assert_allclose(J.shape[0], J.shape[1:], err_msg='sets of states differ')

    # groom J a little
    J = fix_J(J)

    # figure out the paralinear branch lengths
    b, a, c = fit_paralinear(J).argsort()

    # use spectral method to get Par for shortest branch to minimise error
    U = standard_normal(J.shape[0])
    Pabgammas = [funny_P(J, a, b, c, gamma) for gamma in range(J.shape[0])]
    PabU = sum(u*P for u, P in zip(U, Pabgammas))
    Jab = pairwise_J(J, a, b, c)
    pia = bafz(Jab.sum(1))
    Pab = diag(1/pia).dot(Jab)

    invPmb = eig(inv(Pab).dot(PabU))[1]
    Pmb = inv(invPmb).real
    Pmb = diag(sign(Pmb.sum(1))).dot(Pmb)
    Pmb[Pmb < 1e-9] = 1e-9
    Pmb = diag(1./bafz(Pmb.sum(1))).dot(Pmb)
    invPmb = inv(Pmb)
    assert is_stochastic(Pmb), 'Pmb not stochastic'

    Jma = (Jab.dot(invPmb)).T
    Jma[Jma < 1e-9] = 1e-9
    Jma /= Jma.sum()
    pim = bafz(Jma.sum(1))
    Pma = diag(1./pim).dot(Jma)
    assert is_stochastic(Pma), 'Pma not stochastic'

    Jcb = pairwise_J(J, c, b, a)
    Jmc = (Jcb.dot(invPmb)).T
    Jmc[Jmc < 1e-9] = 1e-9
    Jmc /= Jmc.sum()
    pim = Jmc.sum(1)
    Pmc = diag(1./pim).dot(Jmc)
    assert is_stochastic(Pmc), 'Pmc not stochastic'

    # now reconstruct the rows
    R = DLCesque_permutation(Pma)
    Ps = [R.dot(P) for P in array((Pma, Pmb, Pmc))[array([a, b, c]).argsort()]]
    pim = best_pi(J, Ps)

    return pim, Ps

    X = eig(PabU.dot(inv(Pab)))[1].real # making it real is a bit dodgey
    eta = inv(X).sum(1) # some mossel and roch magic here
    Par = X.dot(diag(eta))
    Par[Par < 0.] = 0.
    Par = diag(1./bafz(Par.sum(1))).dot(Par)

    # calculate pir and Pra from Par and pia
    Jar = diag(pia).dot(Par)
    pir = bafz(Jar.sum(0))
    Pra = diag(1/pir).dot(Jar.T)
    assert (Pra >= 0.).all(), 'Pra not stochastic'
    for row in Pra:
        if not allclose(row.sum(), 1.):
            print row
    assert allclose(Pra.sum(1), 1.), 'Pra not stochastic'

    # use Par to calculate Prb and Prc
    Prb = inv(Par).dot(Pab)
    Pac = pairwise_P(J, a, c, b)
    Prc = inv(Par).dot(Pac)
    assert (Prb >= 0.).all(), 'Prb not stochastic'
    assert (Prc >= 0.).all(), 'Prc not stochastic'
    assert allclose(Prb.sum(1), 1.), 'Prb not stochastic'
    assert allclose(Prc.sum(1), 1.), 'Prc not stochastic'

    # the following is not true for noisy data
    # assert_allclose(Par.dot(diag(Prc.dot(U))).dot(Prb), PabU)

    # now reconstruct the rows
    Pra = real_if_close(Pra)
    R = DLCesque_permutation(Pra)
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
