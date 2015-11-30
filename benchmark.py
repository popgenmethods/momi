## usage: python benchmark.py /path/to/ms [--threads num_threads] [--reset]

import argparse, time, subprocess, os, sqlite3, random, itertools
import networkx as nx
import multiprocessing as mp
import numpy as np
from collections import Counter, defaultdict

from momi import Demography

## comand line arguments
parser = argparse.ArgumentParser()
parser.add_argument("ms", type=str, help="Path to ms/scrm")
parser.add_argument("--threads", type=int, default=1, help="Number of parallel threads")
parser.add_argument("--reset", action="store_true", help="Reset results database")

## the sql database to store results
conn = sqlite3.connect('benchmark_results.db')
cur = conn.cursor()

## number of repetitions for each sample size and number of demes
reps = 20

def main():
    args = parser.parse_args()

    ## arguments for time_runs() jobs
    args_list = []

    ## benchmarking for n=2,4,8,...,256 samples
    for n_total_log2 in range(1,9):
        ## samples come from D=1,2,4,...,n demes
        for n_demes_log2 in range(1, n_total_log2+1):

            n_demes = 2**n_demes_log2
            n_per_deme = 2**(n_total_log2 - n_demes_log2)

            ## don't use coalescent for n>=256 and D>=32, due to long running time
            moran_only = (n_total_log2 >= 8) and (n_demes >= 32)

            args_list += [(args.ms, n_demes, n_per_deme, moran_only)] * reps

    ## run the jobs
    results_list = mp.Pool(processes=args.threads).map(time_runs, args_list)

    ## store the results in sqlite
    store_results_list(results_list, args.reset)


def time_runs(args):
    ms_path, n_taxa, lineages_per_taxon, moranOnly = args

    print "# Starting job with n_per_deme=%d, n_demes=%d" % (lineages_per_taxon, n_taxa)

    # Get a random phylogeny
    newick_str = random_binary_tree(n_taxa)
    demo = Demography.from_newick(newick_str, lineages_per_taxon)

    n = n_taxa * lineages_per_taxon
    results = []
    snp_list = run_simulation(ms_path, demo, 100, lineages_per_taxon)
    for snp,state in enumerate(snp_list):
        #print(state)
        state_tuple = tuple([v['derived'] for k,v in sorted(state.iteritems())])
        rid = random.getrandbits(32)

        method_list = [("moran",False)]
        if not moranOnly:
            method_list += [("chen",True)]

        for name,use_chen_eqs in method_list:
            with Timer() as t:
                ret = demo.sfs(state, use_chen_eqs)
            results.append((name,n_taxa, lineages_per_taxon, snp, t.interval, ret, rid, str(state_tuple), newick_str))

    print "# Finished job with n_per_deme=%d, n_demes=%d" % (lineages_per_taxon, n_taxa)
    return results


class Timer:
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start
        #print('Call took %.03f sec.' % self.interval)

## functions for generating a random binary tree demography

def random_binary_tree(n):
    g = nx.DiGraph()
    nodes = ["l%i" % i for i in range(n)]
    i = 0

    lengths = {v: 0.0 for v in nodes}
    for _ in range(n - 1):
        coal = random.sample(nodes, 2)
        int_node = "i%i" % i

        t = random.expovariate(1) / float(n) * 2.0
        lengths = {k: v + t for k, v in lengths.items()}

        g.add_edges_from([(int_node, c, {'edge_length': lengths[c]}) for c in coal])

        lengths[int_node] = 0.0

        i += 1
        nodes = list((set(nodes) - set(coal)) | set([int_node]))
    assert len(nodes) == 1
    return newick_helper(g, nodes[0])

def newick_helper(g, node):
    children = g.successors(node)
    el = None
    try:
        parent = g.predecessors(node)[0]
        el = g[parent][node]['edge_length']
    except IndexError:
        parent = None
    if children:
        ret = "(%s,%s)" % tuple([newick_helper(g, c) for c in children])
        if parent is not None:
            ret += ":%f" % el
        return ret
    else:
        ret = "%s" % node
        if el is not None:
            ret += ":%f" % el
        return ret

## functions for generating dataset with ms/scrm

def build_command_line(ms_path, demo, L, lineages_per_taxon):
    '''Given a tree, build a ms command line which will simulate from it.'''
    ejopts = []
    Iopts = []

    tfac = 0.5
    theta = 1.

    lineages = []
    lineage_map = {}
    for i, leaf_node in list(enumerate(sorted(demo.leaves), 1)):
        nsamp = lineages_per_taxon
        Iopts.append(nsamp)
        lineage_map[leaf_node] = i
        lineages += [leaf_node] * nsamp
        age = demo._node_data[leaf_node]['model'].tau * tfac

        p, = demo.predecessors(leaf_node)
        while True:
            if p not in lineage_map:
                lineage_map[p] = i
                tau = demo._node_data[p]['model'].tau
                #if p.edge_length == float("inf"):
                if tau == float('inf'):
                    break
                age += tau * tfac
                old_p = p
                p, = demo.predecessors(p)
            else:
                # We have a join-on time
                ejopts.append((age, i, lineage_map[p]))
                break

    cmdline = ["-I %d %s" % (len(Iopts), " ".join(map(str, Iopts)))]
    for ej in ejopts:
        cmdline.append("-ej %g %d %d" % ej)
    cmdline = ["%s %d 1 -t %g" % (ms_path, sum(Iopts), theta)] + cmdline
    #print(cmdline)
    return lineages, " ".join(cmdline)

def run_simulation(ms_path, tree, L, lineages_per_taxon):
    lineages, cmd = build_command_line(ms_path, tree, L, lineages_per_taxon)
    species = list(set(lineages))
    n_lineages = Counter(lineages)

    try:
        lines = subprocess.check_output(cmd.split())
    except subprocess.CalledProcessError, e:
        ## ms gives really weird error codes, so ignore them
        lines = e.output
    lines = lines.split("\n")

    #print(cmd)
    output = [l.strip() for l in lines]
    def f(x):
        if x == "//":
            f.i += 1
        return f.i
    f.i = 0
    for k, lines in itertools.groupby(output, f):
        if k == 0:
            continue
        # Skip preamble
        next(lines)
        # segsites
        segsites = int(next(lines).split(" ")[1])
        # positions
        next(lines)
        # at haplotypes
        lin_counts = defaultdict(lambda: np.zeros(segsites, dtype=int))
        for hap, lin in zip(lines, lineages):
            hap = list(map(int, hap))
            lin_counts[lin] += hap
    return [{lin: {'derived':lin_counts[lin][i],
                   'ancestral':n_lineages[lin] - lin_counts[lin][i]}
             for lin in lineages}
            for i in range(segsites)]

### functions for manipulating the database of results (benchmark_results.db)

def store_results_list(results_list, reset_db):
    if reset_db:
        cur.execute("drop table results")

    try:
        create_table()
    except sqlite3.OperationalError:
        pass

    for results in results_list:
        for entry in results:
            store_result(*entry)
    conn.commit()
    conn.close()

def create_table():
    cur.execute("""create table results (model varchar, n integer, """
                 """lineages integer, site integer, time real, result real, run_id integer, state varchar, tree varchar)""")

def store_result(name, n, l, i, t, res, rid, state, tree):
    cur.execute("insert into results values (?, ?, ?, ?, ?, ?, ?, ?, ?)", (name, n, l, i, t, res, rid, state, tree))


if __name__=="__main__":
    main()
