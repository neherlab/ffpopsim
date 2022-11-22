#!/usr/bin/env python
# vim: fdm=indent

import sys

# Python 3.6 compatibility: Python 3.6 has no compatible version of setuptools that supports
# PEP-517 (the `pyproject.toml`-based config), so we use a shim package `ppsetuptools` instead.
if sys.version_info >= (3, 7):
  from setuptools import setup, Extension
else:
  from ppsetuptools import setup, Extension

import numpy as np

setup_args = dict(
  py_modules=["FFPopSim"],
  package_dir={ '': 'src/python' },
  ext_modules = [
    Extension(
      '_FFPopSim',
      [
        'src/python/FFPopSim_wrap.cpp',
        'src/haploid_highd.cpp',
        'src/haploid_lowd.cpp',
        'src/hivpopulation.cpp',
        'src/hivgene.cpp',
        'src/rootedTree.cpp',
        'src/multiLocusGenealogy.cpp',
        'src/hypercube_lowd.cpp',
        'src/hypercube_highd.cpp',
      ],
      include_dirs=[np.get_include()],
      libraries=['gsl', 'gslcblas'],
      py_limited_api = True
    )
  ]
)
setup(**setup_args)
