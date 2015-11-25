# momi

This repository implements the algorithm and results in this [paper](http://arxiv.org/abs/1503.01133).
In particular, it computes the expected joint site frequency spectrum (SFS) for a tree-shaped demography without migration,
via a multipopulation Moran model.

It also computes the "truncated site frequency spectrum" for a single population, i.e. the frequency
spectrum for mutations arising after a certain point in time. This can be used in both Moran and coalescent
approaches to computing the multipopulation SFS.

## Installation

Prerequisites:
* Scientific distribution of Python 2.7, e.g. [Anaconda](http://continuum.io/downloads), [Enthought Canopy](https://www.enthought.com/products/canopy/)
  * Alternatively, custom installation of pip, cython, the SciPy stack
* gcc
  * May be required for installation of some dependencies

To install, in the top-level directory of momi (where "setup.py" lives), type
```
pip install .
```

## Getting started

See [example.py](example.py) for an example of how to construct a population history and compute its SFS entries.

The results from the paper are implemented in [benchmark.py](benchmark.py) and [benchmark_run.py](benchmark_run.py).


## Upcoming features

The current version of momi is 1.0. Upcoming features for the next version (1.0) of momi are:
* Pulse migration/admixture
* Parameter inference via gradient descent, automatic differentiation
* Improvements to user interface and computational efficiency

Version 2.0 will be released late 2015 or early 2016, to be accompanied with a second paper.
Stay tuned!