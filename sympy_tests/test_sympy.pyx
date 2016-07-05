#!/usr/bin/env python
#
# test.pyx
# cython: profile=True
#import cython
from sympy import *
#import timeit
#@cython.boundscheck(False)
def main_routine():
  a= symbols('a')
  s= sin (a)
  c= cos(a)
  pprint( simplify( s**2+c**2 ))
