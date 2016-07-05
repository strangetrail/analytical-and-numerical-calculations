#!/usr/bin/env python
#
###############################################################################
#                                                                             #
#                                                                             #
#                                                                             #
#                                   test5.py                                  #
#                                                                             #
#                                                                             #
#                                                                             #
###############################################################################
#
import sys
sys.path.append( '/home/nixuser/Documents/code/sympy_tests/common/' )
from sympy import *
from utils import *
#
###############################################################################
#                                    TYPES                                    #
###############################################################################
#
enumExprType = enum( SELF = 0, OTHER = 1 )
#
###############################################################################
#^^^^^                               TYPES                               ^^^^^#
###############################################################################
#
#
###############################################################################
#                                   SETTINGS                                  #
###############################################################################
#
bIntegrate     = False
bDifferentiate =  True
bSolve         = False
#
iExprType = enumExprType.OTHER
#
iIntIdx = -1
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
y, a, b, c      = symbols ( 'y, a, b, c',    real = True, )
if ( bSolve ) :
  x, d, e = symbols ( 'x, a, c, d, e', real = True )
else:
  x, b_p        = symbols ( 'x, b_p',        real = True, positive = True )
  b_n           = symbols ( 'b_n',           real = True, negative = True )
pass# IF ( bSolve )
#
###############################################################################
#^^^^^                              SYMBOLS                              ^^^^^#
###############################################################################
#
if ( bSolve ) :
  R = a*x**4 + b*x**3 + c*x**2 + d*x + e
else:
  if ( iExprType == enumExprType.SELF ) :
    R_p   = sqrt ( x**2 + b_p*x - 1 )
    R_n   = sqrt ( x**2 + b_n*x - 1 )
  else:
    if   ( iIntIdx == 0 ) :
      R   = sqrt ( x**2 + b_n*x - 1 )
    elif ( iIntIdx == 1 ) :
      R   = sqrt ( ( x - b_p ) * ( x - b_p ) )
    else:
      R   = sqrt ( a*y**2 + b*y + c )
    pass# IFCASE ( iIntIdx )
  pass# IF ( iExprType == enumExprType.SELF )
pass# IF ( bSolve )
#
if ( bDifferentiate ) :
  if ( iExprType == enumExprType.SELF ) :
    genform_p = diff ( log ( Abs( 2*R_p + 2*x + b_p ) ), x )
    genform_n = diff ( log ( Abs( 2*R_n + 2*x + b_n ) ), x )
  else:
    if   ( iIntIdx == 0 ) :
      genform = diff ( log ( Abs( 2*R + 2*x + b_n ) ), y )
    elif ( iIntIdx == 1 ) :
      genform = diff ( log ( Abs( 2*R + 2*x + b_n ) ), y )
    else:
      #genform = simplify( diff ( log ( Abs( 2*R + 2*y + b ) ), y ) )
      genform = simplify( diff ( log ( 2 * sqrt( a ) * R + 2*a*y + b ) / sqrt ( a ), y ) )
    pass# IFCASE ( iIntIdx )
  pass# IF ( iExprType == enumExprType.SELF )
pass# IF ( bDifferentiate )
#
if ( bIntegrate ) :
  if ( iExprType == enumExprType.SELF ) :
    genform_p   = integrate ( 1 / R_p, x )
    genform_n   = integrate ( 1 / R_n, x )
  else:
    if   ( iIntIdx == 0 ) :
      genform   = integrate ( 1 / R, x )
    elif ( iIntIdx == 1 ) :
      genform   = integrate ( R, x )
    else:
      genform   = simplify( integrate ( 1 / R, y ) )
    pass# IFCASE ( iIntIdx )
  pass# IF ( iExprType == enumExprType.SELF )
pass# IF ( bIntegrate )
#
if ( bSolve ) :
  geneq = Eq( 0, R )
  roots = solve( geneq, x )
pass# IF ( bSolve )
#
###############################################################################
#                            PRINTING OUTPUT TO FILE                          #
###############################################################################
#
OutputFile = open( 'test5.sympy.tex', 'w' )
OutputFile.write( '\\documentclass{article}\n' )
OutputFile.write( '\\usepackage{amsmath}\n' )
OutputFile.write( '\\begin{document}\n' )
if ( bSolve ) :
  if ( roots ) :
    for item in roots:
      printTexExpressionInline( OutputFile, item )
    pass# FOR item in roots
  pass# IF ( items )
  print '\n'
  pprint ( geneq )
  print '\n'
  printTexExpressionInline ( OutputFile, 0 )
else:
  if ( iExprType == enumExprType.SELF ) :
    printTexExpressionInline ( OutputFile, simplify( genform_p ) )
    printTexExpressionInline ( OutputFile, simplify( genform_n ) )
  else:
    printTexExpressionInline ( OutputFile, simplify( genform ) )
    pprint ( genform )
  pass# IF ( iExprType == enumExprType.SELF )
pass# IF ( bSolve )
OutputFile.write( '\\end{document}' )
OutputFile.close()
#
###############################################################################
#^^^^^                       PRINTING OUTPUT TO FILE                     ^^^^^#
###############################################################################
#
