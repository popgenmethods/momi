from __future__ import division
import numpy
import scipy.misc
import operator
import math
from util import memoize_instance
import warnings
from size_history import ConstantTruncatedSizeHistory
import numpy as np
from convolution_momi import convolve_chen

math_mod = math
myint,myfloat = int,float

## UNCOMMENT FOR HIGHER PRECISION
# import gmpy2
# math_mod = gmpy2
# gmpy2.get_context().precision=100
# myint,myfloat = gmpy2.mpz, gmpy2.mpfr

'''
Formulas from Hua Chen 2012, Theoretical Population Biology
Note that for all formulas from that paper, N = diploid population size
'''

class _SumProduct_Chen(object):
    '''
    compute sfs of data via Hua Chen's sum-product algorithm
    '''
    def __init__(self, demography):
        self.G = demography
        attach_Chen(self.G)

    def p(self):
        '''Return the likelihood for the data'''
        return self.partial_likelihood_bottom(self.G.root)[1]

    def partial_likelihood_top(self, node):
        lik, sfs = self.partial_likelihood_bottom(node)
        return self.G.chen[node].apply_transition(lik), sfs

    def partial_likelihood_bottom(self, node):
        sfs = 0.0
        if self.G.is_leaf(node):
            n = self.G.n_lineages_subtended_by[node]
            n_der = self.G._n_derived_subtended_by[node]
            lik = np.zeros((len(n_der), n+1, n+1))
            lik[range(len(n_der)), n-n_der, n_der] = 1.0
        else:
            children = tuple(self.G[node])
            ch_liks, ch_sfs = zip(*[self.partial_likelihood_top(ch)
                                    for ch in children])
            ch_liks = [l * self.combinatorial_factors(ch)
                       for l,ch in zip(ch_liks, children)]
            lik = convolve_chen(*ch_liks) / self.combinatorial_factors(node)

            for ch1, ch2 in ((0,1),(1,0)):
                sfs += ch_sfs[ch1] * (self.G._n_derived_subtended_by[children[ch2]] == 0)
        sfs += (lik * self.truncated_sfs(node)).sum(axis=(1,2))
        return lik,sfs

    def combinatorial_factors(self, node):
        n_node = self.G.n_lineages_subtended_by[node]

        n_der = np.outer(np.ones(n_node+1), np.arange(n_node+1))
        n_anc = np.outer(np.arange(n_node+1), np.ones(n_node+1))

        return scipy.misc.comb(n_der + n_anc, n_der)

    def truncated_sfs(self, node):
        n_node = self.G.n_lineages_subtended_by[node]
        sfs = np.zeros((n_node+1,n_node+1))
        for n_der in range(1,n_node+1):
            for n_anc in range(n_node-n_der+1):
                sfs[n_anc,n_der] = self.G.chen[node].freq(n_der, n_der+n_anc)
        return sfs


def attach_Chen(tree):
    '''Attach Hua Chen equations to each node of tree.
    Does nothing if these formulas have already been added.'''
    if not hasattr(tree, "chen"):
        tree.chen = {}
        for node in tree:
            size_model = tree._node_data[node]['model']
            if type(size_model) is not ConstantTruncatedSizeHistory:
                raise NotImplementedError("Hua Chen's equations only implemented for constant population size along each branch")
            tree.chen[node] = SFS_Chen(size_model.N / 2.0, size_model.tau, tree.n_lineages_subtended_by[node])

class SFS_Chen(object):
    def __init__(self, N_diploid, timeLen, max_n):
        self.timeLen = timeLen
        self.N_diploid = N_diploid
        self.max_n = max_n
        # precompute
        for n in range(1,max_n+1):
            for i in range(1,n+1):
                self.freq(i,n)

            max_m = n
            if timeLen == float('inf'):
                max_m = 1
            for m in range(1,max_m+1):
                self.g(n,m)

    @memoize_instance
    def g(self, n, m):
        return g(n, m, self.N_diploid, self.timeLen)

    @memoize_instance
    def ET(self, i, n, m):
        try:
            return ET(i, n, m, self.N_diploid, self.timeLen)
        except ZeroDivisionError:
            warnings.warn("divide by zero in hua chen formula")
            return 0.0

    @memoize_instance
    def ES_i(self, i, n, m):
        '''TPB equation 4'''
        assert n >= m
        return math.fsum([p_n_k(i, n, k) * k * self.ET(k, n, m) for k in range(m, n + 1)])

    @memoize_instance
    def freq(self, i, n):
        max_m = n-i+1
        if self.timeLen == float('inf'):
            max_m = 1

        ret = 0.0
        for m in range(1,max_m+1):
            ret += self.ES_i(i, n, m)
        return ret

    def apply_transition(self, likelihoods):
        ## einsum should be faster, but causes memory overflow for our machines/tests
        # return np.einsum('ijkl,mij->mkl',
        #                  self.transition_tensor(),
        #                  likelihoods)
        ## just use for loops
        n = likelihoods.shape[-1]-1
        assert likelihoods.shape[1:] == (n+1,n+1)
        ret = np.zeros(likelihoods.shape)
        for n_top in range(1,n+1):
            for n_bottom in range(n_top,n+1):
                for n_derived_bottom in range(n_bottom+1):
                    for n_derived_top in range(n_derived_bottom+1):
                        n_ancestral_bottom = n_bottom - n_derived_bottom
                        n_ancestral_top = n_top - n_derived_top

                        tmp = np.array(likelihoods[:,n_ancestral_bottom,n_derived_bottom])
                        tmp *= self.g(n_bottom,
                                      n_top) * math.exp(log_urn_prob(n_derived_top,
                                                                     n_ancestral_top,
                                                                     n_derived_bottom,
                                                                     n_ancestral_bottom))

                        ret[:,n_ancestral_top,n_derived_top] += tmp
        return ret

    @memoize_instance
    def transition_tensor(self):
        n = self.max_n

        ret = np.zeros((n+1,n+1,n+1,n+1))
        for n_top in range(1,n+1):
            for n_bottom in range(n_top,n+1):
                for n_derived_bottom in range(n_bottom+1):
                    for n_derived_top in range(n_derived_bottom+1):
                        n_ancestral_bottom = n_bottom - n_derived_bottom
                        n_ancestral_top = n_top - n_derived_top
                        ret[n_ancestral_bottom,
                            n_derived_bottom,
                            n_ancestral_top,
                            n_derived_top] = self.g(n_bottom,
                                                    n_top) * math.exp(log_urn_prob(n_derived_top,
                                                                                   n_ancestral_top,
                                                                                   n_derived_bottom,
                                                                                   n_ancestral_bottom))
        return ret


def log_factorial(n):
    return math_mod.lgamma(n+1)

def log_rising(n,k):
    return log_factorial(n+k-1) - log_factorial(n-1)

def log_falling(n,k):
    return log_factorial(n) - log_factorial(n-k)

def gcoef(k, n, m, N_diploid, tau):
    k, n, m = map(myint, [k, n, m])
    N_diploid = myfloat(N_diploid)
    tau = myfloat(tau)
    return (2*k - 1) * (-1)**(k - m) * math_mod.exp(log_rising(m, k-1) + log_falling(n, k) - log_factorial(m) - log_factorial(k - m) - log_rising(n, k))
    #return (2*k - 1) * (-1)**(k - m) * rising(m, k-1) * falling(n, k) / math_mod.factorial(m) / math_mod.factorial(k - m) / rising(n, k)


def g_sum(n, m, N_diploid, tau):
    if tau == float("inf"):
        if m == 1:
            return 1.0
        return 0.0
    tau = myfloat(tau)
    return float(sum([gcoef(k, n, m, N_diploid, tau) * math_mod.exp(-k * (k - 1) * tau / 4 / N_diploid) for k in range(m, n + 1)]))


g = g_sum

def formula1(n, m, N_diploid, tau):
    def expC2(k):
        return math_mod.exp(-k * (k - 1) / 4 / N_diploid * tau)
    r = sum(gcoef(k, n, m, N_diploid, tau) *
            ((expC2(m) - expC2(k)) / (k - m) / (k + m - 1) - (tau / 4 / N_diploid * expC2(m)))
            for k in range(m + 1, n + 1))
    #q = 4 * N_diploid / g(n, m, N_diploid, tau)
    q = 4 * N_diploid
    return float(r * q)


def formula3(j, n, m, N_diploid, tau):
    # Switch argument to j here to stay consistent with the paper.
    j, n, m = map(myint, [j, n, m])
    tau, N_diploid = map(myfloat, [tau, N_diploid])
    def expC2(kk):
        return math_mod.exp(-kk * (kk - 1) / 4 / N_diploid * tau)
    r = sum(gcoef(k, n, j, N_diploid, tau) * # was gcoef(k, n, j + 1, N_diploid, tau) *
            sum(gcoef(ell, j, m, N_diploid, tau) * ( # was gcoef(ell, j - 1, m, N_diploid, tau) * (
                    (
                        expC2(j) * (tau / 4 / N_diploid - ((k - j) * (k + j - 1) + (ell - j)*(ell + j - 1)) / # tau / 4 / N_diploid was 1 in this
                             (k - j) / (k + j- 1) / (ell - j) / (ell + j - 1))
                    )
                    +
                    (
                        expC2(k) * (ell - j) * (ell + j - 1) / (k - j) / (k + j - 1) / (ell - k) / (ell + k - 1)
                    )
                    -
                    (
                        expC2(ell) * (k - j) * (k + j - 1) / (ell - k) / (ell + k - 1) / (ell - j) / (ell + j - 1)
                    )
                )
                for ell in range(m, j)
                )
            for k in range(j + 1, n + 1)
            )
    #q = 4 * N_diploid / myfloat(g(n, m, N_diploid, tau))
    q = 4 * N_diploid
    return float(q * r)

def formula2(n, m, N_diploid, tau):
    def expC2(k):
        return math_mod.exp(-k * (k - 1) / 4 / N_diploid * tau)
    r = sum(gcoef(k, n, m, N_diploid, tau) *
            ((expC2(k) - expC2(n)) / (n - k) / (n + k - 1) - (tau / 4 / N_diploid * expC2(n)))
            for k in range(m, n))
    #q = 4 * N_diploid / g(n, m, N_diploid, tau)
    q = 4 * N_diploid
    return float(r * q)

def ET(i, n, m, N_diploid, tau):
    '''Starting with n lineages in a population of size N_diploid,
    expected time when there are i lineages conditional on there
    being m lineages at time tau in the past.'''
    if tau == float("inf"):
        if m != 1 or i == 1:
            return 0.0
        return 2 * N_diploid / float(nChoose2(i)) * g(n, m, N_diploid, tau)
    if n == m:
        return tau * (i == n) * g(n, m, N_diploid, tau)
    if m == i:
        return formula1(n, m, N_diploid, tau)
    elif n == i:
        return formula2(n, m, N_diploid, tau)
    else:
        return formula3(i, n, m, N_diploid, tau)

def p_n_k(i, n, k):
    if k == 1:
        return int(i == n)
    else:
        #return scipy.misc.comb(n-i-1,k-2) / scipy.misc.comb(n-1,k-1)
        return math.exp(log_binom(n - i - 1, k - 2) - log_binom(n - 1, k - 1))

def nChoose2(n):
    return (n * (n-1)) / 2

def log_binom(n, k):
    if k < 0 or k > n:
        return -float('inf')
    return log_factorial(n) - log_factorial(n - k) - log_factorial(k)

def log_urn_prob(n_parent_derived, n_parent_ancestral, n_child_derived, n_child_ancestral):
    n_parent = n_parent_derived + n_parent_ancestral
    n_child = n_child_derived + n_child_ancestral
    if n_child_derived >= n_parent_derived and n_parent_derived > 0 and n_child_ancestral >= n_parent_ancestral and n_parent_ancestral > 0:
        return log_binom(n_child_derived - 1, n_parent_derived - 1) + log_binom(n_child_ancestral - 1, n_parent_ancestral - 1) - log_binom(n_child-1, n_parent-1)
    elif n_child_derived == n_parent_derived == 0 or n_child_ancestral == n_parent_ancestral == 0:
        return 0.0
    else:
        return float("-inf")
