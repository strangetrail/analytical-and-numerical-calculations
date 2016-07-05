#!/bin/python
#
###############################################################################
#                                                                             #
#                                                                             #
#                                                                             #
#                                   test6.py                                  #
#                                                                             #
#                                                                             #
#                                                                             #
###############################################################################
#
# Solving common integrals for square cube anisotropic potential source.
# MatLab and Wolfram use constant optimization - should be aware of it.
#
###############################################################################
#                                   IMPORTS                                   #
###############################################################################
#
import sys
sys.path.append( '/home/nixuser/Documents/code/sympy_tests/common/' )
from utils import *
from sympy import *
#
###############################################################################
#^^^^^                              IMPORTS                              ^^^^^#
###############################################################################
#
beginSimpleTiming()
#
###############################################################################
#                                   SETTINGS                                  #
###############################################################################
#
bIntegrateSimpleSqrtEquationU1 = False
bComputeSeries                 = False
bSolveLogArgumentEquation      = False
bIntegrate1stBaseEquationX     = False
bIntegrate1stBaseEquationY     = False
bIntegrate2ndBaseEquationXY    = False
bIntegrate2ndBaseEquationXYZ   = False
bIntegrateSimpleAsinhEquationX = False
bSimplifyAsinExpandArgument    =  True
#
###############################################################################
#^^^^^                              SETTINGS                             ^^^^^#
###############################################################################
#
#
###############################################################################
#                                   SYMBOLS                                   #
###############################################################################
#
x, y, z, a, b, u_1, u_2, u_3, u, X, Y, Z, u =                         (
 symbols( 'x, y, z, a, b, u_1, u_2, u_3, u X, Y, Z, u', real = True )
                                                                      )
#
l_0 = symbols( 'l_0', positive = True, real = True )
n, m = symbols( 'n, m', positive = True, integer = True )
#
###############################################################################
#^^^^^                              SYMBOLS                              ^^^^^#
###############################################################################
#
if ( bSimplifyAsinExpandArgument ):
  genform_ =                           (
             (
              l_0 * (2*m - 1) +
              sqrt
              (
               l_0**2 * (2*m - 1)**2 +
               l_0**2 * (2*n + 1)**2 +
               4 * z**2
              )
             )**2 *
             (
              l_0**2 * (2*m + 1)**2 +
              4 * z**2
             ) -
             (
              l_0 * (2*m + 1) *
              (
               l_0 * (2*m - 1) +
               sqrt
               (
                l_0**2 * (2*m - 1) +
                l_0**2 * (2*n + 1) +
                4 * z**2
               )
              )
              +
              4 * z**2
             )**2
                                       )
  genform_ = simplify( powdenest( genform_ , force = True ) )
pass# IF ( bSimplifyAsinExpandArgument )
#
if bIntegrate1stBaseEquationX:
  genform = simplify( integrate( x**2 / sqrt( x**2 + y**2 + z**2 ), x ) )
pass# IF
#
if bIntegrate1stBaseEquationY:
  genform = simplify( integrate( genform, y ) )
pass# IF
#
if bIntegrate2ndBaseEquationXY:
  genform = simplify( integrate( integrate( (x*y) / sqrt( x**2 + y**2 + z**2 ), x ), y ) )
pass# IF
#
if bIntegrate2ndBaseEquationXYZ:
  genform = simplify(
             integrate( integrate( integrate( (x*y) / sqrt( x**2 + y**2 + z**2 ), x ), y ), z )
            )
pass# IF
#
if bIntegrateSimpleAsinhEquationX:
  genform = integrate( x**2 * asinh( a / sqrt( x**2 + b ) ), x )
pass# IF
#
if bIntegrateSimpleSqrtEquationU1:
  genform = integrate( (( sqrt( ( u_1 - x )**2 - x**2 - z**2 ) )**3) / u_1, u_1 )
pass# IF
#
#eq_to_differentiate = log( Abs( 2 * sqrt( u**2 - 2*x*u - y**2 ) + 2*u - y**2 ) )
#genform = simplify( diff( eq_to_differentiate, u ) )
#eq_to_simplify = (-u_3) - 2 * I * ( x * ( x + sqrt( x**2 + z**2 + Y ) ) + Y ) + sqrt( ( x + sqrt( x**2 + z**2 + Y ) )**2 * ( 4 * x**2 + 4 * Y ) - ( 2 * x * ( x + sqrt( x**2 + z**2 + Y ) ) + 2 * Y )**2 )
#eq_to_simplify_LHS = u_3 + 2 * I * Y + 2 * I * x**2
#eq_to_simplify_RHS = 2 * sqrt(Y) * Abs(z) - 2 * I * x * sqrt( Y + x**2 + z**2 )
#eq_to_simplify_LHS = (u_3)**2 + 4 * I * u_3 * x**2 + 4 * x**2 * z**2 - 4 * Y**2 + 4 * I * Y * u_3 - 4 * Y * x**2 - 4 * Y * z**2
#eq_to_simplify_RHS = -8 * I * sqrt(Y) * x * Abs(z) * sqrt( Y + x**2 + z**2 )
#genform = simplify( diff( genform, x ) )
#roots = solve( Eq( genform, 0 ) , Y )
#pprint( roots )
#
#genform = simplify( diff( x/(2*x*sqrt( x**2 + y**2 + z**2 ))*( -sqrt( x**2 + y**2 + z**2 ) *(y**2+z**2) * asinh(x/(sqrt(y**2+z**2))) + x*(x**2+y**2+z**2)   ) ,  x) )
#
if bSolveLogArgumentEquation:
  l_x = l_0
  l_y = l_0
pass# IF
#
#x_1 = l_x * ( 2*n + 1 ) / 2
#x_2 = l_x * ( 2*n - 1 ) / 2
#y_1 = l_y * ( 2*m + 1 ) / 2
#y_2 = l_y * ( 2*m - 1 ) / 2
#
if bSolveLogArgumentEquation:
  x_1 = l_x * ( 2*n + 1 )
  x_2 = l_x * ( 2*n - 1 )
  y_1 = l_y * ( 2*m + 1 )
  y_2 = l_y * ( 2*m - 1 )
pass# IF
#
if bSolveLogArgumentEquation:
  Arg2 = y_2 + sqrt( y_2**2 + x_1**2 + 4 * z**2 )
  Arg1 = y_2 * Arg2 + 4 * z**2
  Arg3 = y_1**2 + 4 * z**2
  #Argument = Arg1 / ( Arg2 * Arg3 )
  #Argument = sqrt( Arg2**2 * Arg3 - Arg1**2 ) - I * Arg1
  Argument = Arg2**2 * Arg3 - Arg1**2
  #genform_ =                                                                                (
  #          simplify
  #          (
  #           powdenest
  #           (
  #            ( Arg2 * Arg3 ) *
  #            simplify( powdenest( I * Argument + sqrt( 1 - Argument**2 ), force = True ) )
  #            ,
  #            force = True
  #           )
  #          )
  #                                                                                          )
  genform_ = simplify( powdenest( Argument, force = True ) )
  #
  if ( bComputeSeries ):
    log_common_series = 0
    for i in range(1,7):
      log_common_series = (-1)**(i+1) * (Argument**i / i)
    pass# FOR i in range(0,3)
  pass# IF ( bComputeSeries )
  #
  #genform = Arg3**2 * Arg2**2 - ( Arg1 )**2
  #genform0 = simplify( genform )
  #genform0_ = Arg1
  #genform_0 = simplify( Arg1 )
  #gf = l_0 * 2 * z * (2*n - 1) + a * (Arg1) - u
  #roots = solve( Eq( gf, 0 ), z )
  #pprint( roots )
pass# IF
#
#Arg11 = (-2) * ( x_2 * ( x_2 + sqrt( x_2**2 + y_1**2 + z**2 ) ) + z**2 )
#Arg12 = x_2 + sqrt( x_2**2 + y_1**2 + z**2 )
#Arg13 = 2 * sqrt( x_2**2 + z**2 )
#Argument1 = Arg11 / ( Arg12 * Arg13 )
#genform1 = I * Argument1 + sqrt( 1 - Argument1**2 )
#Arg21 = (-2) * ( x_2 * ( x_2 + sqrt( x_2**2 + y_2**2 + z**2 ) ) + z**2 )
#Arg22 = x_2 + sqrt( x_2**2 + y_2**2 + z**2 )
#Arg23 = 2 * sqrt( x_2**2 + z**2 )
#Argument2 = Arg21 / ( Arg22 * Arg23 )
#genform2 = I * Argument2 + sqrt( 1 - Argument2**2 )
#genform = simplify( genform1 / genform2 )
#
#Arg11 =                                                                                      (
#        (-2) *
#        ( x_2 * ( x_2 + powdenest( sqrt( x_2**2 + y_1**2 + z**2 ), force = True ) ) + z**2 )
#                                                                                             )
#Arg12 = x_2 + powdenest( sqrt( x_2**2 + y_1**2 + z**2 ), force = True )
#Arg13 = 2 * powdenest( sqrt( x_2**2 + z**2 ), force = True )
#Argument1 = Arg11 / ( Arg12 * Arg13 )
#genform1 = powdenest( simplify( powdenest( I * Argument1 + powdenest( sqrt( 1 - Argument1**2 ), force = True ), force=True ) ), force=True )
#Arg21 = (-2) * ( x_2 * ( x_2 + powdenest( sqrt( x_2**2 + y_2**2 + z**2 ), force = True ) ) + z**2 )
#Arg22 = x_2 + powdenest( sqrt( x_2**2 + y_2**2 + z**2 ), force = True )
#Arg23 = 2 * powdenest( sqrt( x_2**2 + z**2 ), force = True )
#Argument2 = Arg21 / ( Arg22 * Arg23 )
#genform2 = powdenest( simplify( powdenest( I * Argument2 + powdenest( sqrt( 1 - Argument2**2 ), force = True ), force=True ) ), force=True )
#genform = powdenest( simplify( powdenest( simplify( powdenest( genform1 / genform2, force=True ) ), force=True ) ), force = True )
#
print "Done."
f = open( 'test6.sympy.tex', 'w' )
f.write( '\\documentclass{article}\n\\usepackage{amsmath}\n' )
f.write( '\\usepackage{mathtools}\n\\begin{document}\n' )
f.write( '\\begin{equation}\n\\begin{split}\n' )
f.write( latex( genform_, mode='inline' ) )
f.write( '\n\\end{split}\n\\end{equation}\n' )
#f.write( '\\begin{equation}\n\\begin{split}\n' )
#f.write( latex( Argument, mode='inline' ) )
#f.write( '\n\\end{split}\n\\end{equation}\n' )
#f.write( '\\begin{equation}\n\\begin{split}\n' )
#f.write( latex( genform0, mode='inline' ) )
#f.write( '\n\\end{split}\n\\end{equation}\n' )
#f.write( '\\begin{equation}\n\\begin{split}\n' )
#f.write( latex( genform0_, mode='inline' ) )
#f.write( '\n\\end{split}\n\\end{equation}\n' )
#f.write( '\\begin{equation}\n\\begin{split}\n' )
#f.write( latex( genform_0, mode='inline' ) )
#f.write( '\n\\end{split}\n\\end{equation}\n' )
#for r in roots:
  #f.write( '\\begin{equation}\n\\begin{split}\n' )
  #f.write( latex( r, mode='inline' ) )
  #f.write( '\n\\end{split}\n\\end{equation}\n' )
  #gen = simplify( genform.subs( z, r ) )
  #z = r
  #f.write( '\\begin{equation}\n\\begin{split}\n' )
  #f.write( latex( genform, mode='inline' ) )
  #f.write( '\n\\end{split}\n\\end{equation}\n' )
#pass# FOR r in roots
f.write( '\\end{document}' )
f.close()
#
