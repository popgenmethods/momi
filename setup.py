#!/usr/bin/env python

from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy

extensions = [Extension("convolution_momi",
                        sources=["momi/convolution_momi.pyx"],
                        include_dirs=[numpy.get_include()])]

setup(name='momi',
      version='1.2.4',
      description='Compute the site frequency spectrum (SFS) of population genetics',
      author='Jack Kamm, Jonathan Terhorst, Yun S. Song',
      author_email='jkamm@stat.berkeley.edu, terhorst@stat.berkeley.edu, yss@eecs.berkeley.edu',
      packages=['momi'],
      install_requires=[
            'biopython',
            'numpy',
            'scipy',
            'networkx<2.0'],
      url='https://github.com/popgenmethods/momi',
      keywords=['population genetics', 'statistics', 'site frequency spectrum', 'coalescent'],
      ext_modules=cythonize(extensions),
      )
