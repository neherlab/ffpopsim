#!/usr/bin/env python2

"""
setup.py file for SWIG python bindings
"""

from distutils.core import setup, Extension


popgenlib_module = Extension('_PopGenLib', sources=['PopGenLib_wrap.cpp'])

setup (name = 'popgenlib',
       author      = "Richard Neher, Fabio Zanini",
       description = """PopGenLib library""",
       ext_modules = [popgenlib_module],
       py_modules = ["PopGenLib"],
       )
