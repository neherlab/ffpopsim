#!/usr/bin/env python2

"""
setup.py file for SWIG python bindings
"""

from distutils.core import setup, Extension


hivpython_module = Extension('_hivpython',
                           sources=['hivpython_wrap.cpp',
                                    'hivpython.cpp'],
                           #TODO: this flag should be generated dynamically
                           # (using a C++ interface file? use Makefile +
                           # command-line args?)
                           include_dirs=['/usr/lib/python2.7/site-packages/numpy/core/include',
                                         '/home/fabio/university/phd/artificial_evolution/libraries/PopGenLib/src'],
                           library_dirs=['/home/fabio/university/phd/artificial_evolution/libraries/PopGenLib/src',
                                         '/home/fabio/university/phd/artificial_evolution/libraries/HandyTools/src'],
                           libraries=['HandyTools', 'PopGenLib'],
                           )

setup (name = 'popgenlib',
       author      = "Richard Neher, Boris Shraiman, Fabio Zanini",
       description = """PopGenLib library""",
       ext_modules = [hivpython_module],
       py_modules = ["hivpython"],
       )
