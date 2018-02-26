#!/usr/bin/env python

from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

cmdclass = {}
ext_modules = []

if use_cython:
    ext_modules.append(
        Extension("convolution_momi",
                  sources=["momi/convolution_momi.pyx"],
                  include_dirs=[numpy.get_include()]))
    cmdclass['build_ext'] = build_ext
else:
    ext_modules.append(
        Extension("convolution_momi",
                  sources=["momi/convolution_momi.c"],
                  include_dirs=[numpy.get_include()]))


# see https://stackoverflow.com/a/4515279/3718509 and https://packaging.python.org/tutorials/distributing-packages/ for uploading to pypi

setup(name='momi',
      version='1.2.5',
      description='Compute the site frequency spectrum (SFS) of population genetics',
      author='Jack Kamm, Jonathan Terhorst, Yun S. Song',
      author_email='jkamm@stat.berkeley.edu, terhorst@stat.berkeley.edu, yss@eecs.berkeley.edu',
      packages=['momi'],
      install_requires=[
            'biopython',
            'numpy',
            'scipy',
            'networkx'],
      url='https://github.com/popgenmethods/momi',
      keywords=['population genetics', 'statistics', 'site frequency spectrum', 'coalescent'],
      ext_modules=ext_modules,
      cmdclass=cmdclass)
