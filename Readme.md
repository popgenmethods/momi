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
  * Alternatively, custom installation of pip, and the SciPy stack
* gcc

To install, simply type
```
pip install momi
```

## Getting started

See [tutorial.py](tutorial.py) to get started.

The results from the paper are implemented in [paper_results/](paper_results/).

## Upcoming features

The current version of momi is 1.2.3. Upcoming features for the next version (2.0.0) of momi are:
* Pulse migration/admixture
* Parameter inference via gradient descent, automatic differentiation
* Improvements to user interface

Version 2.0 will be released in early 2016, to be accompanied with a second paper.