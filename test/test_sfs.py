from momi import Demography
import pytest
from test_demography import demo
import cPickle as pickle
import numpy as np
import itertools

results_file = 'test_sfs.pickle'

try:
    with file(results_file,'r') as f:
        results = pickle.load(f)
except:
    pass

@pytest.mark.parametrize("use_chen_eqs", (True,False))
def test_compute_sfs(demo, use_chen_eqs):
    assert np.allclose(demo.compute_sfs(results['derived_counts'], use_chen_eqs=use_chen_eqs),
                       results['sfs_vals'])

@pytest.mark.parametrize("use_chen_eqs", (True,False))    
def test_old_sfs(demo, use_chen_eqs):
    for i,config in enumerate(results['derived_counts']):
        config = {l : {'derived' : d,
                       'ancestral' : demo.n_lineages_subtended_by[l] - d} for l,d in config.iteritems()}
        assert np.isclose(demo.sfs(config, use_chen_eqs=use_chen_eqs), results['sfs_vals'][i])
       
if __name__=="__main__":
    derived_counts = [{'a': i, 'b': j}
                      for i,j in itertools.product(range(11), range(9))
                      if not ((i == 0 and j==0) or (i==10 and j==8))]

    sfs_vals = demo().compute_sfs(derived_counts)
    with file(results_file,'w') as f:
        pickle.dump({'derived_counts': derived_counts,
                     'sfs_vals' : sfs_vals},
                    f)
    
