#!/usr/bin/env python
#
from distutils.core import setup
from Cython.Build import cythonize
#from distutils.extension import Extension
#from Cython.Distutils import build_ext
#import sympy
#import timeit
#
#setup(
#      cmdclass = {'build_ext': build_ext},
#      ext_modules = [Extension("testsympy", ["testsympy.pyx"])]
#     )
setup(
      name = 'test4',
      ext_modules = cythonize("test4.pyx"),
     )
#
