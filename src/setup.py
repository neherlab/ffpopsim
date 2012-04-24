#!/usr/bin/env python2

"""
setup.py file for SWIG python bindings
"""

from distutils.core import setup, Extension


hivpopulation_module = Extension('_hivpopulation',
                           sources=['hivpopulation_test_wrap.cpp',
                                    'hivpopulation_test.cpp'],
                           #TODO: this flag should be generated dynamically
                           # (using a C++ interface file?)
                           include_dirs=['/usr/lib/python2.7/site-packages/numpy/core/include'],
                           library_dirs=['/home/fabio/university/phd/artificial_evolution/libraries/PopGenLib/src',
                                         '/home/fabio/university/phd/artificial_evolution/libraries/HandyTools/src'],
                           libraries=['HandyTools', 'PopGenLib'],
                           )

setup (name = 'popgenlib',
       author      = "Richard Neher, Boris Shraiman, Fabio Zanini",
       description = """PopGenLib library""",
       ext_modules = [hivpopulation_module],
       py_modules = ["hivpopulation"],
       )
