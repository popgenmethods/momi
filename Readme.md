# momi

This repository implements the algorithm in [this paper](http://www.tandfonline.com/doi/abs/10.1080/10618600.2016.1159212) ([preprint here](http://arxiv.org/abs/1503.01133)).

In particular, this repository computes the expected joint site frequency spectrum (SFS) for a tree-shaped demography without migration,
via a multipopulation Moran model.

It also computes the "truncated site frequency spectrum" for a single population, i.e. the frequency
spectrum for mutations arising after a certain point in time. This can be used in both Moran and coalescent
approaches to computing the multipopulation SFS.

## Installation

Prerequisites:
* Scientific distribution of Python 2.7, e.g. [Anaconda](http://continuum.io/downloads), [Enthought Canopy](https://www.enthought.com/products/canopy/)
  * Alternatively, custom installation of pip, and the SciPy stack
* gcc

To install, simply type
```
pip install momi
```

## Getting started

See [tutorial.py](tutorial.py) to get started.

See [this repository](https://github.com/jackkamm/momi1_paper_results) to reproduce the results in the paper.
