#!/usr/bin/env python
#
#import csympy
import sys
sys.path.append( '/home/nixuser/Documents/code/sympy_tests/common/' )
from utils import *
from sympy import *
#
#csympy.__version__
#vf, d, a, vi, t = S('vf d a vi t'.split())
#xi, tau = symbols ( 'xi, tau' )
x, y, s, z, a, w = symbols ( 'x, y, s, z, a, w' )
#U = Function( 'u' )
#u = U( xi, tau )
#F = Function( 'f' )
#u = F( (Abs(u))**2 )
#x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = S('x0 x1 x2 x3 x4 x5 x6 x7 x8 x9 x10'.split())
#a, b, c, x = symbols(  'a, b, c, x'  )
#equations = ([
# Eq(vf, vi+a*t),
# Eq(d, vi*t + a*t**2/2),
# Eq(a, 10),
# Eq(d, 60),
# Eq(vi, 5)
#])
#equations = Eq(0, a*x**2 + b*x + c)
#equations = ([
# Eq( -44.000000,   6.000000*x0 - 124.000000*x1 - 80.000000*x2 + 52.000000*x3 ),
# Eq(-101.000000, -12.000000*x0 +  12.000000*x1 - 69.000000*x2 - 60.000000*x3 ),
# Eq(  70.000000, -95.000000*x0 +   1.000000*x1 - 76.000000*x2 - 84.000000*x3 ),
# Eq(  81.000000,  -5.000000*x0 - 109.000000*x1 + 15.000000*x2 + 69.000000*x3 )
#])
#pprint(solve(equations, x))
#pprint(solve(equations, [x0, x1, x2, x3]))
#genform = diff ( x * log ( 2 * ( sqrt ( x**2 + y**2 + z**2 ) + y ) ) + y * log ( 2 * ( sqrt ( x**2 + y**2 + z**2 ) + x ) ) - z * atanh ( x * y / ( z * sqrt ( x**2 + y**2 + z**2 ) ) ) + z * atanh ( x / z ) - x, x )
genform = x+y+z+2*x*sin(a)+s**3 + exp(w)
#
#pprint(genform)
result_stream = open( 'test2.sympy.tex', 'w' )
result_stream.write( '\\documentclass{article}\n\\usepackage{amsmath}\n' )
result_stream.write( '\\usepackage{mathtools}\n\\begin{document}\n' )
genformsimpl = simplifySingleExpressioninMaxima ( genform )
printTexExpressionInline ( result_stream, genformsimpl )
#print ( repr ( genformsimpl ) )
result_stream.write( '\\end{document}' )
result_stream.close()
#
