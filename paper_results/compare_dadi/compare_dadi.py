from __future__ import division
import dadi, momi, time, itertools, subprocess, sys, pandas, random
import numpy as np
import cPickle as pickle
import pandas as pd
from StringIO import StringIO

#n_per_deme_list = [16,32]
#n_grid_dadi_list = [16,32,64,128]

n_per_deme_list = [16,32,64,128]
n_grid_dadi_list = [16,32,64,128,256,512,1024]

results_file = 'figures/compare_dadi.pickle'

def main():
    try:
        with file(results_file,'r') as f:
            results = pickle.load(f)
    except IOError:
        results = {}

    if sys.argv[-1] == "--reset":
        results = {}
        
    jobs = []
    old_results = dict(results)
    for n_per_deme in n_per_deme_list:
        jobs += [{'method': 'momi',
                  'n_per_deme': n_per_deme}]
        for n_pts in n_grid_dadi_list:
            if n_pts < n_per_deme:
                continue
            jobs += [{'method' : 'dadi',
                      'n_pts': n_pts,
                      'n_per_deme' : n_per_deme}]
    dadi_phis = {}
    for job in jobs:
        key = tuple(sorted(job.items()))
        if key in results:
            continue

        print job
        n_per_deme = job['n_per_deme']

        if job['method'] == 'momi':
            with Timer() as t:
                results[key] = {'sfs': compute_momi(n_per_deme)}
            results[key]['time'] = t.interval
        elif job['method'] == 'dadi':
            phi_key = job['n_pts']

            if phi_key not in dadi_phis:
                with Timer() as t:
                    phi,xx = compute_dadi(job['n_pts'])
                dadi_phis[phi_key] = (phi,xx,t.interval)

            phi,xx,phi_t = dadi_phis[phi_key]
            with Timer() as t:
                results[key] = {'sfs': 2.0 * np.array(dadi.Spectrum.from_phi(phi, (n_per_deme,)*3, (xx,)*3))}
            results[key]['time'] = phi_t + t.interval

    # save results
    if results != old_results:
        with file(results_file,'w') as f:
            pickle.dump(results,f,-1)

    # plot results
    summarize_results(results)
    plot_timing(results)
    plot_accuracy(results)



# ancestral size
N_a = 7300

# scaled parameters
N_af = 12300 / N_a
N_b = 2100 / N_a
N_eu0 = 1000 / N_a
N_as0 = 510 / N_a

r_eu = .004 * 4 * N_a
r_as = .0055 * 4 * N_a

T_af = 220e3 / 25 / N_a / 2
T_b = 140e3 / 25 / N_a / 2
T_eu_as = 21.2e3 / 25 / N_a / 2

N_eu_f = N_eu0 * np.exp(r_eu/2 * T_eu_as)
N_as_f = N_as0 * np.exp(r_as/2 * T_eu_as)


## compute SFS with dadi

def compute_dadi(pts):
    xx = dadi.Numerics.default_grid(pts)
    # at ancestor
    phi = dadi.PhiManip.phi_1D(xx)
    #  out of africa
    phi = dadi.Integration.one_pop(phi, xx, T_af - T_b, nu = N_af)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # integrate two pops
    phi = dadi.Integration.two_pops(phi, xx, T_b - T_eu_as, nu1=N_af, nu2=N_b)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)

    # integrate three pops
    phi = dadi.Integration.three_pops(phi, xx, T_eu_as,
                                      nu1=N_af,
                                      nu2=lambda t:N_as0*np.exp(r_as/2*t),
                                      nu3=lambda t:N_eu0*np.exp(r_eu/2*t))
    return phi,xx
    #return dadi.Spectrum.from_phi(phi, (n_per_deme,)*3, (xx,)*3)



### compute SFS with momi


def compute_momi(n_per_deme):
    demo_newick_str = """
    ((
    eu:%f[&&momi:lineages=%d:model=exponential:N_top=%f:N_bottom=%f],
    as:%f[&&momi:lineages=%d:model=exponential:N_top=%f:N_bottom=%f]
    ):%f[&&momi:N=%f],
    af:%f[&&momi:lineages=%d:N=%f]
    )[&&momi:model=piecewise:model_0=constant:tau_0=%f:N_0=%f:model_1=constant:N_1=1:tau_1=inf]
    """ % (T_eu_as, n_per_deme, N_eu0, N_eu_f,
           T_eu_as, n_per_deme, N_as0, N_as_f,
           T_b - T_eu_as, N_b,
           T_b, n_per_deme, N_af,
           T_af-T_b, N_af)

    demo = momi.Demography.from_newick(demo_newick_str)

    idxs = [(i,j,k)
            for i,j,k in itertools.product(range(n_per_deme+1),repeat=3)
            if not (i == j and i ==k and (i==0 or i == n_per_deme))]
    config_list = [{'af':i,'as':j,'eu':k}
                   for (i,j,k) in idxs]

    ret = np.zeros((n_per_deme+1,)*3)
    ret[zip(*idxs)] = demo.compute_sfs(config_list)
    return ret


### plotting functions

r_file = 'plot.r'

def plot_timing(results):
    df = "\t".join(["G","n","seconds"])
    for k,v in results.items():
        k = dict(k)
        if k['method'] == 'momi':
            G = 0
        elif k['method'] == 'dadi':
            G = k['n_pts']
        df += "\n" + "\t".join(map(str,[G,k['n_per_deme'],v['time']]))

    p = subprocess.Popen(["Rscript",r_file,"timing"], stdin=subprocess.PIPE,cwd='figures')
    p.communicate(df)

def summarize_results(results):
    df = "\t".join(["n","n.pts","min","neg_entries","total_entries"]) + "\n"
    for k,v in sorted(results.items()):
        k = dict(k)
        if k['method'] == 'momi':
            n_pts = 'momi'
        else:
            n_pts = k['n_pts']
        n = k['n_per_deme']

        v = np.array(v['sfs'])
        v[0,0,0] = v[-1,-1,-1] = float('inf')
        v = pd.Series(np.reshape(v,-1))

        df += "\t".join(map(str,
                            [n,n_pts,
                             v.min(),
                             np.sum(v < 0),
                             v.size])) + "\n"
    with file('figures/neg_entries.txt','w') as f:
        f.write(df)

def plot_accuracy(results):
    results_list = []

    n_list = set([dict(k)['n_per_deme'] for k in results.keys()])
    idxs = {n: [(i,j,k)
                for i,j,k in itertools.product(range(n+1), repeat=3)
                if not (i==j and i==k and (i==0 or i==n))]
            for n in n_list}
    # for n,idx in idxs.iteritems():
    #     try:
    #         # sample down number of points to plot
    #         idxs[n] = random.sample(idx, 10000)
    #     except ValueError:
    #         pass

    for k,v in results.items():
        k = dict(k)
        if k['method'] == 'dadi':
            n_pts = k['n_pts']
        else:
            n_pts = 'momi'
        n = k['n_per_deme']
        v = v['sfs']
        assert v.shape == (n+1,)*3
        for i,j,k in idxs[n]:
            results_list += [(n_pts,n,i,j,k,v[i,j,k])]

    momi_results = {x[1:-1]: x[-1] for x in results_list if x[0] == 'momi'}
    results_list = [x for x in results_list if x[0] != 'momi']

    df = "\t".join(['n','n.pts','val','momi'])
    for x in results_list:
        df += "\n" + "\t".join(map(str, ["\n" , x[1], x[0], x[-1], momi_results[x[1:-1]]]))

    p = subprocess.Popen(["Rscript",r_file,"accuracy"], stdin=subprocess.PIPE,cwd='figures')
    p.communicate(df)


class Timer:
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start

if __name__=="__main__":
    main()
