#!/usr/bin/env python
# vim: fdm=indent

import sys
import subprocess
import os

# Python 3.6 compatibility: Python 3.6 has no compatible version of setuptools that supports
# PEP-517 (the `pyproject.toml`-based config), so we use a shim package `ppsetuptools` instead.
if sys.version_info >= (3, 7):
  from setuptools import setup, Extension
else:
  from ppsetuptools import setup, Extension

import numpy as np

# Get GSL include paths
def get_gsl_include():
    try:
        output = subprocess.check_output(['gsl-config', '--cflags'], stderr=subprocess.DEVNULL)
        cflags = output.decode().strip().split()
        return [flag[2:] for flag in cflags if flag.startswith('-I')]
    except:
        return []

# Get Boost include paths
def get_boost_include():
    for path in ['/opt/homebrew/opt/boost/include', '/usr/local/opt/boost/include', '/opt/local/include']:
        if os.path.isdir(path):
            return [path]
    return []

# Get GSL library paths
def get_gsl_lib_dirs():
    try:
        output = subprocess.check_output(['gsl-config', '--libs'], stderr=subprocess.DEVNULL)
        ldflags = output.decode().strip().split()
        return [flag[2:] for flag in ldflags if flag.startswith('-L')]
    except:
        return []

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
      include_dirs=[np.get_include()] + get_gsl_include() + get_boost_include(),
      library_dirs=get_gsl_lib_dirs(),
      libraries=['gsl', 'gslcblas'],
      extra_compile_args=['-std=c++11'],
      py_limited_api = True
    )
  ]
)
setup(**setup_args)
