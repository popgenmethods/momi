from momi import Demography, ConstantTruncatedSizeHistory, ExponentialTruncatedSizeHistory, PiecewiseHistory

## First, let's compute the SFS for a single population with a three-epoch history

n = 20 ## sample size

## specify the epochs backwards in time from the present
epochs = [ExponentialTruncatedSizeHistory(n_max=n, tau=0.25, N_top=1.0, N_bottom=10.0), # an exponential growth period for 0.25 units of time. "top"=past, "bottom"=present
          ConstantTruncatedSizeHistory(n_max = n, tau=0.1, N=0.1), # a bottleneck for 0.1 units of time
          ConstantTruncatedSizeHistory(n_max = n, tau=float('inf'), N=1.0), # ancestral population size = 1
          ]
## turn the list of epochs into a single population history
demo1 = PiecewiseHistory(epochs)

## print the SFS entries
print "Printing SFS entries for three epoch history"
print [demo1.freq(i, n) for i in range(1,n)]


## Next, print out the "truncated SFS" for just the middle epoch
## i.e., the frequency spectrum for mutations that occur within the middle epoch
print "\nPrinting truncated SFS for middle epoch"
print [epochs[1].freq(i,n) for i in range(1,n)]


## Finally, let's compute SFS entries for a multipopulation demography
# For our demography, we'll have three leaf populations: 'a','b','c', with 10,5,8 sampled alleles respectively

## specify demography via a newick string
## to specify additional population parameters, follow branch length with [&&momi:...]
## where ... contains:
##    lineages= # alleles (if population is leaf)
##    model= population history (default=constant)
##    N, N_top, N_bottom= parameters for constant/exponential size history
##    model_i= model for i-th epoch of piecewise history (either constant or exponential)
##    N_i, N_top_i, N_bottom_i, tau_i= parameters for i-th epoch of piecewise history
newick_str = """
((
a:.25[&&momi:lineages=10:model=constant:N=10.0],
b:.3[&&momi:lineages=5:model=exponential:N_top=1.0:N_bottom=10.0]
):.1[&&momi:model=constant:N=1.5],
c:.3[&&momi:lineages=8:model=piecewise:model_0=exponential:tau_0=.2:N_top_0=.1:N_bottom_0=1.0:model_1=constant:tau_1=.1:N_1=.3]
)[&&momi:N=3.0]
"""

demo3 = Demography.from_newick(newick_str)

# construct a [list] of SNPs to compute SFS for
# each SNP is represented as a {dict} giving the derived allele counts

derived_counts_list = [{'a': 0, 'b': 0, 'c':1}, # singleton SNP in population 'c'
                       {'a': 2, 'b':0, 'c':0}, # doubleton SNP in population 'a'
                       {'a': 10,'b':5, 'c':0}, # SNP with derived allele fixed in 'a','b', but not present in 'c'
                       ]

# compute the SFS
sfs_list = demo3.compute_sfs(derived_counts_list)

print "Printing a few SFS entries for a 3 population demography"

print "Derived_Counts","\t","SFS_Value"
for derived_count,sfs_val in zip(derived_counts_list, sfs_list):
    print derived_count, "\t", sfs_val


### For benchmarking in paper (see paper_results/compare_chen/benchmark.py), we also implemented Hua Chen's formulas,
### but ONLY for the special case of constant population size along each branch.
try:
    # set use_chen_eqs=True to use Chen's formulas
    print demo3.compute_sfs(derived_counts_list, use_chen_eqs=True)
except NotImplementedError:
    # Chen's formulas not implemented for demo3, due to changing size along branches
    pass


