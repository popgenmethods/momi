from __future__ import print_function, division
import scipy.misc, scipy.signal
import math
import numpy as np
from convolution_momi import convolve_momi

class _SumProduct(object):
    '''
    compute joint SFS entry via a sum-product like algorithm,
    using Moran model to do transitions up the tree,
    and using Polanski Kimmel to compute frequency of mutations at each node in tree
    '''
    def __init__(self, demography):
        self.G = demography

    def p(self):
        '''Return joint SFS entry for the demography'''
        return self.partial_likelihood_bottom(self.G.root)[1]

    def leaf_likelihood_bottom(self, leaf):
        n_der = self.G._n_derived_subtended_by[leaf]
        n_node = self.G._node_data[leaf]['lineages']

        ret = np.zeros((len(n_der), n_node + 1))
        ret[zip(*enumerate(n_der))] = 1.0
        return ret

    def combinatorial_factors(self, node):
        n_node = self.G.n_lineages_subtended_by[node]
        return scipy.misc.comb(n_node, np.arange(n_node + 1))

    def truncated_sfs(self, node):
        n_node = self.G.n_lineages_subtended_by[node]
        sfs = np.array([self.G._node_data[node]['model'].freq(n_derived, n_node) for n_derived in range(n_node + 1)])
        sfs[sfs == float("inf")] = 0.0
        return sfs

    def partial_likelihood_top(self, node):
        bottom_likelihood, sfs = self.partial_likelihood_bottom(node)
        return np.einsum('ij,kj->ik',
                         bottom_likelihood,
                         self.G._node_data[node]['model'].transition_matrix), sfs

    def partial_likelihood_bottom(self, node):
        sfs = 0.0
        if self.G.is_leaf(node):
            lik = self.leaf_likelihood_bottom(node)
        else:
            children = tuple(self.G[node])
            ch_liks, ch_sfs = zip(*[self.partial_likelihood_top(ch) for ch in children])

            ch_liks = [l * self.combinatorial_factors(ch)
                       for l,ch in zip(ch_liks, children)]

            lik = convolve_momi(*ch_liks) / self.combinatorial_factors(node)

            for ch1,ch2 in ((0,1),(1,0)):
                sfs += ch_sfs[ch1] * (self.G._n_derived_subtended_by[children[ch2]] == 0)

        sfs += (lik * self.truncated_sfs(node)).sum(axis=1)
        return lik,sfs
