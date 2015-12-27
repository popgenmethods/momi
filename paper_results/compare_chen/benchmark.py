## usage: python benchmark.py [--help] [--ms_path /path/to/ms] [--threads num_threads]

description = """
Generate and plot the figures in the paper, comparing accuracy and speed of momi vs. algorithm of Chen 2012.
"""

import argparse, time, subprocess, os, random, itertools, re
import networkx as nx
import multiprocessing as mp
import numpy as np
from collections import Counter, defaultdict
import cPickle as pickle

from momi import Demography

results_file = 'benchmark.results.pickle'
r_file = 'benchmark.plot.R'

## comand line arguments
parser = argparse.ArgumentParser(description=description)
parser.add_argument("--ms_path", type=str, default=None,
                    help="Path to ms/scrm. If provided, will simulate new results and overwrite the results file (%s)" % results_file)
parser.add_argument("--threads", type=int, default=1,
                    help="Number of parallel threads")

## number of repetitions for each sample size and number of demes
reps = 20
loci_per_rep = 100
theta_per_locus = 1.0

def main():
    args = parser.parse_args()
    if args.ms_path is None:
        print "--ms_path not provided, using existing results from %s" % results_file
    else:
        print "Overwriting results in %s" % results_file
        run_benchmarks(args.ms_path, args.threads)

    print "Plotting results with %s" % r_file
    with file(results_file,'r') as f:
        results_list = pickle.load(f)
    plot_accuracy(results_list)
    plot_times(results_list)

def run_benchmarks(ms_path, n_threads):
    ## arguments for time_runs() jobs
    jobs_list = []
    ## benchmarking for n=2,4,8,...,256 samples
    for n_total_log2 in range(1,9):
        ## samples come from D=1,2,4,...,n demes
        for n_demes_log2 in range(1, n_total_log2+1):

            n_demes = 2**n_demes_log2
            n_per_deme = 2**(n_total_log2 - n_demes_log2)

            ## don't use coalescent for n>=256 and D>=32, due to long running time
            moran_only = (n_total_log2 >= 8) and (n_demes >= 32)

            jobs_list += [(ms_path, n_demes, n_per_deme, moran_only, random.randint(0,999999999))] * reps

    ## run the jobs
    results_list = mp.Pool(processes=n_threads).map(do_job, jobs_list, chunksize=1)

    ## store the results
    with file(results_file,'w') as f:
        pickle.dump(results_list, f)


def do_job((ms_path, n_demes, n_per_pop, moranOnly, seed)):
    print "# Starting job with n_per_deme=%d, n_demes=%d" % (n_per_pop, n_demes)

    random.seed(seed)

    # Get a random phylogeny
    newick_str = random_binary_tree(n_demes)
    demo = Demography.from_newick(newick_str, n_per_pop)

    n = n_demes * n_per_pop
    snp_list = run_simulation(ms_path, demo, loci_per_rep, n_per_pop)
    snp_list = [{pop: counts['derived'] for pop,counts in snp.iteritems()}
                for snp in snp_list]

    # only keep unique snps
    snp_list = set([snp.iteritems() for snp in snp_list])
    snp_list = map(dict, snp_list)

    method_list = [("momi",False)]
    if not moranOnly:
        method_list += [("Chen",True)]

    results = {}
    for name,use_chen_eqs in method_list:
        with Timer() as t:
            demo.compute_sfs([snp_list[0]], use_chen_eqs)
        precompute_t = t.interval

        with Timer() as t:
            results[name] = {'sfs': demo.compute_sfs(snp_list, use_chen_eqs)}
        results[name].update({'timing': {"'Per SNP'" : t.interval / float(len(snp_list)), "'Precomputation'" : precompute_t}})

    print "# Finished job with n_per_deme=%d, n_demes=%d" % (n_per_pop, n_demes)
    return {'n': n_demes * n_per_pop, 'n_pops' : n_demes, 'results': results, 'seed' : seed}

class Timer:
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start
        #print('Call took %.03f sec.' % self.interval)


## functions for plotting

def plot_times(results_list):
    dataframe = ""
    dataframe += "\t".join(["method","n","D","component","time"])
    for job_result in results_list:
        for method, res in job_result['results'].iteritems():
            for component, time in res['timing'].iteritems():
                dataframe += "\n" + "\t".join(map(str,[method, job_result['n'], job_result['n_pops'], component, time]))

    p = subprocess.Popen(["Rscript",r_file,'timing'], stdin=subprocess.PIPE)
    p.communicate(dataframe)

def plot_accuracy(results_list):
    dataframe = ""
    dataframe += "\t".join(["n","D","momi","Chen"])
    for res in results_list:
        n,D = res['n'],res['n_pops']
        res = res['results']
        try:
            momi, chen = res['momi']['sfs'], res['Chen']['sfs']
        except KeyError:
            continue
        assert len(momi) == len(chen)
        for m,c in zip(momi,chen):
            dataframe += "\n" + "\t".join(map(str, [n,D,m,c]))
    p = subprocess.Popen(["Rscript",r_file,'accuracy'], stdin=subprocess.PIPE)
    p.communicate(dataframe)


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

def build_command_line(ms_path, demo, L, n_per_pop):
    '''Given a tree, build a ms command line which will simulate from it.'''
    ejopts = []
    Iopts = []

    tfac = 0.5

    haps2pops = []
    lineage_map = {}
    for i, leaf_node in list(enumerate(sorted(demo.leaves), 1)):
        nsamp = n_per_pop
        Iopts.append(nsamp)
        lineage_map[leaf_node] = i
        haps2pops += [leaf_node] * nsamp
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
    cmdline = ["%s %d %d -t %g" % (ms_path, sum(Iopts), L, theta_per_locus)] + cmdline
    #print(cmdline)
    cmdline += ['-seeds'] + [str(random.randint(0,999999999)) for _ in range(3)]
    return haps2pops, " ".join(cmdline)

def run_simulation(ms_path, tree, L, n_per_pop):
    haps2pops, cmd = build_command_line(ms_path, tree, L, n_per_pop)

    try:
        lines = subprocess.check_output(cmd.split())
    except subprocess.CalledProcessError, e:
        ## ms gives really weird error codes, so ignore them
        lines = e.output

    lines_by_rep = lines.split("\n//\n")[1:]

    ret = []
    for lines in lines_by_rep:
        lines = [l.strip() for l in lines.split("\n") if l.strip() != ""]
        segsites = int(re.match("segsites: (\d+)", lines[0]).group(1))

        assert segsites == 0 or len(lines) == 2 + n_per_pop * len(tree.leaves)

        pop_counts = defaultdict(lambda: np.zeros(segsites, dtype=int))
        for hap, pop in zip(lines[2:], haps2pops):
            hap = list(map(int, hap))
            pop_counts[pop] += hap

        for i in range(segsites):
            ret += [{pop: {'derived':pop_counts[pop][i],
                           'ancestral':n_per_pop - pop_counts[pop][i]}
                     for pop in tree.leaves}]
    return ret

### functions for manipulating the database of results (benchmark_results.db)


if __name__=="__main__":
    main()
