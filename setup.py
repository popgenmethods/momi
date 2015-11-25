#!/usr/bin/env python

from distutils.core import setup

setup(name='momi',
      version='0.1',
      description='Moran model for inference',
      author='Jack Kamm, Jonathan Terhorst, Yun S. Song',
      author_email='terhorst@stat.berkeley.edu',
      install_requires=['numpy','scipy','biopython','networkx'],
      url='https://github.com/terhorst/momi',
      )
