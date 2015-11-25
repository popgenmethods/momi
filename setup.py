#!/usr/bin/env python

from distutils.core import setup

setup(name='momi',
      version='0.1',
      description='Moran model for inference',
      author='Jack Kamm, Jonathan Terhorst, Yun S. Song',
      author_email='jkamm@stat.berkeley.edu, terhorst@stat.berkeley.edu, yss@eecs.berkeley.edu',
      packages=['momi'],      
      install_requires=['biopython','numpy','scipy','networkx'],
      url='https://github.com/jackkamm/momi',
      download_url='https://github.com/jackkamm/momi/tarball/1.0',
      keywords=['population genetics','statistics','site frequency spectrum','coalescent'],      
      )
