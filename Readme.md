# momi

This repository implements the algorithm and results in this [paper](http://arxiv.org/abs/1503.01133).
In particular, it computes the expected joint site frequency spectrum (SFS) for a tree-shaped demography without migration,
via a multipopulation Moran model.

It also computes the "truncated site frequency spectrum" for a single population, i.e. the frequency
spectrum for mutations arising after a certain point in time. This can be used in both Moran and coalescent
approaches to computing the multipopulation SFS.

## Installation

To install, type
```
pip install momi
```
Prerequisites:
* Scientific distribution of Python 2.7, e.g. [Anaconda](http://continuum.io/downloads), [Enthought Canopy](https://www.enthought.com/products/canopy/)
  * Alternatively, custom installation of pip, and the SciPy stack
* gcc

## Getting started

See [tutorial.py](tutorial.py) to get started.

The results from the paper are implemented in [benchmark.py](benchmark.py), and plotted with [plot.R](plot.R). To run [benchmark.py](benchmark.py) you will also need [ms](http://home.uchicago.edu/rhudson1/source/mksamples.html) or [scrm](https://scrm.github.io/), or similar program. Then to redo the results in the paper, type
```
python benchmark.py /path/to/ms/or/scrm [--threads num_threads]
Rscript plot.R
```

## Upcoming features

The current version of momi is 1.1. Upcoming features for the next version (2.0) of momi are:
* Pulse migration/admixture
* Parameter inference via gradient descent, automatic differentiation
* Improvements to user interface and computational efficiency

Version 2.0 will be released in early 2016, to be accompanied with a second paper.