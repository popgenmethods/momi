import networkx as nx
from Bio import Phylo
from cStringIO import StringIO
from util import cached_property
from size_history import ConstantTruncatedSizeHistory, ExponentialTruncatedSizeHistory, PiecewiseHistory
import numpy as np
from sum_product import _SumProduct
from huachen_eqs import _SumProduct_Chen

class Demography(nx.DiGraph):
    @classmethod
    def from_newick(cls, newick, default_lineages=None, default_N=1.0):
        t = cls(_newick_to_nx(newick, default_lineages))
        # add models to edges
        for v0, v1, d in t.edges(data=True):
            n_sub = t.n_lineages_subtended_by[v1]
            nd = t._node_data[v1]
            _set_model(nd, n_sub, d['branch_length'], default_N=default_N)
        _set_model(t._node_data[t.root], t.n_lineages_subtended_by[t.root], float('inf'), default_N=default_N)
        return t

    def __init__(self, *args, **kwargs):
        super(Demography, self).__init__(*args, **kwargs)
        nd = self._node_data
        if not all('lineages' in nd[k] for k in self.leaves):
            raise Exception("'lineages' attribute must be set for each leaf node.")

    def to_newick(self):
        return _to_newick(self, self.root)

    def compute_sfs(self, config_list, use_chen_eqs = False):
        if not all([set(c.keys()) == self.leaves for c in config_list]):
            raise Exception("Derived allele counts should be specified at all leafs, and no other vertices")

        config_list = {l: np.array([c[l] for c in config_list], dtype=np.int) for l in self.leaves}

        if any([np.any(d * (self.n_lineages_subtended_by[l] - d) < 0) for l,d in config_list.iteritems()]):
            raise Exception("Derived allele counts should all be between 0 and n")

        sum_d = np.sum(list(config_list.values()), axis=0)
        sum_n = sum([self.n_lineages_subtended_by[l]
                     for l in self.leaves])
        if np.any(sum_d * (sum_n - sum_d) == 0):
            raise Exception("Sites with all ancestral or all derived alleles are not allowed")

        self._n_derived_subtended_by = {v: np.sum([config_list[l] for l in self.leaves_subtended_by[v]], axis=0) for v in self}

        sp = self._sum_product(use_chen_eqs)
        return sp.p()

    def sfs(self, state, use_chen_eqs = False):
        """Deprecated, use compute_sfs instead"""
        if any([v['derived'] + v['ancestral'] != self.n_lineages_subtended_by[k] for k,v in state.iteritems()]):
            raise Exception("Sum of ancestral, derived alleles must be n")
        config_list = [{k: v['derived'] for k,v in state.iteritems()}]
        return self.compute_sfs(config_list, use_chen_eqs=use_chen_eqs)

    @cached_property
    def root(self):
        nds = [node for node, deg in dict(self.in_degree()).items() if deg == 0]
        assert len(nds) == 1
        return nds[0]

    @cached_property
    def leaves(self):
        return set([k for k, v in dict(self.out_degree()).items() if v == 0])

    def is_leaf(self, node):
        return node in self.leaves

    @cached_property
    def leaves_subtended_by(self):
        return {v: self.leaves & set(nx.dfs_preorder_nodes(self, v)) for v in self}

    @cached_property
    def n_lineages_subtended_by(self):
        nd = self._node_data
        return {v: sum(nd[l]['lineages'] for l in self.leaves_subtended_by[v]) for v in self}

    ### "private" methods

    @cached_property
    def _node_data(self):
        return dict(self.nodes(data=True))

    def _sum_product(self, use_chen_eqs):
        if use_chen_eqs:
            return _SumProduct_Chen(self)
        else:
            return _SumProduct(self)


_field_factories = {
    "N": float, "lineages": int, "ancestral": int,
    "derived": int, "model": str, "n": int,
    "tau": float,
    }
def _extract_momi_fields(comment):
    for field in comment.split("&&"):
        if field.startswith("momi:"):
            attrs = field.split(":")
            assert attrs[0] == "momi"
            attrs = [a.split("=") for a in attrs[1:]]
            attrdict = dict((a, _field_factories[a.split('_')[0]](b)) for a, b in attrs)
            return attrdict
    return {}

def _newick_to_nx(newick, default_lineages=None):
    newick = StringIO(newick)
    phy = Phylo.read(newick, "newick")
    phy.rooted = True
    edges = []
    nodes = []
    node_data = {}
    clades = [phy.root]
    phy.root.name = phy.root.name or "root"
    i = 0
    while clades:
        clade = clades.pop()
        nd = _extract_momi_fields(clade.comment or "")
        if 'lineages' not in nd and default_lineages is not None:
            nd['lineages'] = default_lineages
        nodes.append((clade.name, nd))
        for c_clade in clade.clades:
            clades += clade.clades
            if c_clade.name is None:
                c_clade.name = "node%d" % i
                i += 1
            ed = {'branch_length': c_clade.branch_length}
            edges.append((clade.name, (c_clade.name), ed))
    t = nx.DiGraph()
    t.add_edges_from(edges)
    t.add_nodes_from(nodes)
    tn = dict(t.nodes(data=True))
    for node in node_data:
        tn[node].update(node_data[node])
    return t

def _to_newick(G, root):
    parent = list(G.predecessors(root))
    try:
        edge_length = str(G[parent[0]][root]['branch_length'])
    except IndexError:
        edge_length = None
    if not G[root]:
        assert edge_length is not None
        return root + ":" + edge_length
    else:
        children = list(G[root])
        assert len(children) == 2
        dta = [(_to_newick(G, child),) for child in children]
        ret = "(%s,%s)" % (dta[0] + dta[1])
        if edge_length:
            ret += ":" + edge_length
        return ret

def _set_model(nd, n_sub, branch_length, default_N):
    if 'model' not in nd or nd['model'] == "constant":
        nd['model'] = ConstantTruncatedSizeHistory(
                N=nd.get('N', default_N),
                tau=branch_length,
                n_max=n_sub)
    elif nd['model'] == "exponential":
        nd['model'] = ExponentialTruncatedSizeHistory(n_max=n_sub,
                                                      tau=branch_length,
                                                      N_top=nd['N_top'],
                                                      N_bottom=nd['N_bottom'])
    elif nd['model'] == 'piecewise':
        i = -1
        epochs = []
        while True:
            i += 1
            try:
                model = nd['model_'+str(i)]
            except KeyError:
                break
            if  model == 'constant':
                epochs.append(ConstantTruncatedSizeHistory(n_max=n_sub,
                                                           tau=nd['tau_'+str(i)],
                                                           N=nd['N_'+str(i)]))
            elif model == 'exponential':
                epochs.append(ExponentialTruncatedSizeHistory(n_max=n_sub,
                                                              tau=nd['tau_'+str(i)],
                                                              N_top=nd['N_top_'+str(i)],
                                                              N_bottom=nd['N_bottom_'+str(i)]))
            else:
                raise Exception("Unsupported model type: %s" % model)
        nd['model'] = PiecewiseHistory(epochs)
        if not np.allclose(nd['model'].tau , branch_length):
            raise Exception("Sum of piecewise epoch lengths does not equal branch length")
    else:
        raise Exception("Unsupported model type: %s" % nd['model'])

