#!/usr/bin/env python
#
###############################################################################
#                                                                             #
#                                                                             #
#                                                                             #
#                                   test3.py                                  #
#                                                                             #
#                                                                             #
#                                                                             #
###############################################################################
#
import sys
sys.path.append( '/home/nixuser/Documents/code/sympy_tests/common/' )
from utils import *
from sympy import *
#
beginSimpleTiming()
#
enumArgTypeAB = enum( BOTH = 0, ATYPE = 1, BTYPE = 2 )
enumArgType12 = enum( BOTH = 0, TYPE1 = 1, TYPE2 = 2 )
#
###############################################################################
#                                   SETTINGS                                  #
###############################################################################
#
bDoNotSimplify             = False
bPrintArcSines             =  True
bTry2IntegrateArcsines     = False
bExpandArcsines            =  True
bDecomposeLogExpansion2Sum =  True
#
iIsLimitsDependsOnZ = 1
#
eArcsineArgTypeAB = enumArgTypeAB.ATYPE
eArcsineArgType12 = enumArgType12.TYPE1
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
x, y, z, q, u = symbols( 'x, y, z, q, u' )
l_0 = symbols( 'l_0', positive = True, real = True )
n, m = symbols( 'n, m', positive = True, integer = True )
#
###############################################################################
#^^^^^                              SYMBOLS                              ^^^^^#
###############################################################################
#
#
###############################################################################
#                                   GLOBALS                                   #
###############################################################################
#
l_1 = l_0
l_2 = l_0
#
if( iIsLimitsDependsOnZ ):
  lmA = ( l_2 * (2*m - 1) ) / (2*z)
  lmB = ( l_2 * (2*m + 1) ) / (2*z)
  lnA = ( l_1 * (2*n - 1) ) / (2*z)
  lnB = ( l_1 * (2*n + 1) ) / (2*z)
else:
  lmA = ( l_2 * (2*m - 1) ) / 2
  lmB = ( l_2 * (2*m + 1) ) / 2
  lnA = ( l_1 * (2*n - 1) ) / 2
  lnB = ( l_1 * (2*n + 1) ) / 2
pass# IF ( iIsLimitsDependsOnZ )
#
if ( iIsLimitsDependsOnZ ):
  u1 = lmA + sqrt( lnB**2 + 1 + lmA**2 )
  u2 = lmA + sqrt( lnA**2 + 1 + lmA**2 )
  u3 = lmB + sqrt( lnB**2 + 1 + lmB**2 )
  u4 = lmB + sqrt( lnA**2 + 1 + lmB**2 )
else:
  u1 = lmA + sqrt( lnB**2 + z**2 + lmA**2 )
  u2 = lmA + sqrt( lnA**2 + z**2 + lmA**2 )
  u3 = lmB + sqrt( lnB**2 + z**2 + lmB**2 )
  u4 = lmB + sqrt( lnA**2 + z**2 + lmB**2 )
pass# IF ( iIsLimitsDependsOnZ )
#
Aa = 1
if ( iIsLimitsDependsOnZ ):
  Ba = (-2) * lmA
  Ca = -1
else:
  Ba = (-2) * lmA
  Ca = -(z**2)
pass# IF ( iIsLimitsDependsOnZ )
#
Ab = 1
if ( iIsLimitsDependsOnZ ):
  Bb = (-2) * lmB
  Cb = -1
else:
  Bb = (-2) * lmB
  Cb = -(z**2)
pass# IF ( iIsLimitsDependsOnZ )
#
#m = 0
#n = 0
#l_1 = 1
#l_2 = 1
#
###############################################################################
#^^^^^                              GLOBALS                              ^^^^^#
###############################################################################
#
getSimpleTimingData( 'Start timer' )
#
def Ra(U):
  return powdenest( sqrt( Aa*(U**2) + Ba*U + Ca ) )
pass# DEF Ra ( U )
#
def Rb(U):
  return powdenest( sqrt( Ab*(U**2) + Bb*U + Cb ) )
pass# DEF Rb ( U )
#
#def arsinh_symbolic(X):
#  return log( X + sqrt( 1 + X**2 ) )
#
def asinh_arg_A(U):
  return asinh( simplify( powdenest( ( 2*U + Ba ) / sqrt( (-4) - (Ba)**2 ), force=True ) ) )
pass# DEF asinh_arg_A ( U )
#
def asinh_arg_B(U):
  return asinh( simplify( powdenest( ( 2*U + Bb ) / sqrt( (-4) - (Bb)**2 ), force=True ) ) )
pass# DEF asinh_arg_B(U)
#
def asinh_argonly_A(U):
  if ( iIsLimitsDependsOnZ ):
    if ( bDoNotSimplify ):
      return (-I) * ( 2*U + Ba ) / sqrt( (4) + (Ba)**2 )
    else:
      return simplify( (-I) * ( 2*U + Ba ) / sqrt( (4) + (Ba)**2 ) )
  else:
    if ( bDoNotSimplify ):
      return ( 2*U + Ba ) / sqrt( ( 4*Ca ) - (Ba)**2 )
    else:
      return simplify( ( 2*U + Ba ) / sqrt( ( 4*Ca ) - (Ba)**2 ) )
pass# DEF asinh_argonly_A ( U )
#
def asinh_argonly_B(U):
  if ( iIsLimitsDependsOnZ ):
    if ( bDoNotSimplify ):
      return (-I) * ( 2*U + Bb ) / sqrt( (4) + (Bb)**2 )
    else:
      return simplify( (-I) * ( 2*U + Bb ) / sqrt( (4) + (Bb)**2 ) )
  else:
    if ( bDoNotSimplify ):
      return ( 2*U + Bb ) / sqrt( ( 4*Cb ) - (Bb)**2 )
    else:
      return simplify( ( 2*U + Bb ) / sqrt( ( 4*Cb ) - (Bb)**2 ) )
pass# DEF asinh_argonly_B ( U )
#
def asin_argonly_A( U ):
  if ( iIsLimitsDependsOnZ ):
    return                                                                      (
           simplify
           (
            powdenest( ( Ba*U - 2 ) / ( U * sqrt( Ba**2 + 4 ) ), force = True )
           )
                                                                                )
  else:
    return                                                                                 (
           simplify
           (
            powdenest( ( Ba*U + 2*Ca ) / ( Abs(U) * sqrt( (Ba)**2 - 4*Ca ) ), force=True )
           )
                                                                                           )
  pass# IF ( iIsLimitsDependsOnZ )
pass# DEF asin_argonly_A ( U )
#
def asin_argonly_B( U ):
  if ( iIsLimitsDependsOnZ ):
    return                                                                      (
           simplify
           (
            powdenest( ( Bb*U - 2 ) / ( U * sqrt( Bb**2 + 4 ) ), force = True )
           )
                                                                                )
  else:
    return                                                                                 (
           simplify
           (
            powdenest( ( Bb*U + 2*Cb ) / ( Abs(U) * sqrt( (Bb)**2 - 4*Cb ) ), force=True )
           )
                                                                                           )
  pass# IF ( iIsLimitsDependsOnZ )
pass# DEF asin_argonly_B ( U )
#
def asin_arg_A( U ):
  return asin( asin_argonly_A( U ) )
pass# DEF asinh_arg_A ( U )
#
def asin_arg_B(U):
  return asin( asin_argonly_B( U ) )
pass# DEF asin_arg_B ( U )
#
#
###############################################################################
#                                     TEST                                    #
###############################################################################
#
if(False):
  V = (
       z * (
            (
                  Ba * asinh_arg_A( u1 ) - Ba * asinh_arg_A( u2 )
            )/2 +
            (
                  Bb * asinh_arg_B( u4 ) - Bb * asinh_arg_B( u3 )
            )/2
           )
      )
#
#
###############################################################################
#^^^^^                                TEST                               ^^^^^#
###############################################################################
#
def re_collect( arg ):
  collection_ReIm = collect( arg, I, evaluate = False, exact = False )
  if ( S.One in collection_ReIm ):
    return collection_ReIm[S.One]
  else:
    print 'Something is wrong - there is no zero-order term in collection with respect to I.\n'
    return 0
  pass# IF ( S.One in collection_ReIm )
pass# DEF re_collect ( arg )
#
def im_collect( arg ):
  collection_ReIm = collect( arg, I, evaluate = False, exact = False )
  if ( I in collection_ReIm ):
    return collection_ReIm[I]
  else:
    print 'Something is wrong - there is no first-order term in collection with respect to I.\n'
    return 0
  pass# IF ( S.One in collection_ReIm )
pass# DEF im_collect ( arg )
#
def asin_expand( z ):
  return (-I) * log( I*z + sqrt( 1 - z**2 ) )
pass# DEF asin_expand ( z )
#
###############################################################################
#                                     TEST                                    #
###############################################################################
#
def arg_explicit_inequality():
  if ( iIsLimitsDependsOnZ ):
    lmA_raw = lmA * 2*z
    lnB_raw = lnB * 2*z
    lmB_raw = lmB * 2*z
    lnA_raw = lnA * 2*z
  else:
    lmA_raw = lmA * 2
    lnB_raw = lnB * 2
    lmB_raw = lmB * 2
    lnA_raw = lnA * 2
  pass# IF ( iIsLimitsDependsOnZ )
  arg0 = lmB_raw + sqrt( lmB_raw**2 + lnA_raw**2 + 4*z**2 )
  arg1 = lmB_raw**2 + 4*z**2
  arg2 = lmB_raw * arg0 + 4*z**2
  arg_expl = arg0**2 * arg1 - arg2**2
  return arg_expl
pass# DEF arg_explicit_inequality (  )
#
def arg_nonintegrable():
  if ( iIsLimitsDependsOnZ ):
    lmA_raw = lmA * 2*z
    lnB_raw = lnB * 2*z
  else:
    lmA_raw = lmA * 2
    lnB_raw = lnB * 2
  pass# IF ( iIsLimitsDependsOnZ )
  arg0 = lmA_raw + sqrt( lmA_raw**2 + lnB_raw**2 + 4*z**2 )
  arg2 = lmA_raw * arg0 + 4*z**2
  arg3 = 2*l_0*z*(2*n + 1) - arg2
  return arg3
pass# DEF arg_nonintegrable (  )
#
def explicit_arg_diff_squares():
  if ( iIsLimitsDependsOnZ ):
    lmA_raw = lmA * 2*z
    lnB_raw = lnB * 2*z
  else:
    lmA_raw = lmA * 2
    lnB_raw = lnB * 2
  pass# IF ( iIsLimitsDependsOnZ )
  arg_0 = sqrt( lmA_raw**2 + lnB_raw**2 + 4*z**2 )
  arg0 = lmA_raw + arg_0
  arg2 = lmA_raw * arg0 + 4*z**2
  arg3 = ((2*l_0*z*(2*n + 1))**2 - arg2**2)
  return arg3
pass# DEF arg_nonintegrable (  )
#
def arg_nonintegrable_var_subst():
  argNonintegrable = arg_nonintegrable()
  print 'Solving ...\n'
  getSimpleTimingData( 'Solve timer begin' )
  return solve( Eq( argNonintegrable, u ), z )
pass# DEF arg_nonintegrable_var_subst (  )
#
###############################################################################
#^^^^^                                TEST                               ^^^^^#
###############################################################################
#
def taylor_series_log_expansion( arg, idx_last ):
  log_expansion = 0
  for i in range ( 1, idx_last ):
    if ( i % 2 ):
      one_multiplier = 1
    else:
      one_multiplier = -1;
    pass# IF ( i % 2 )
    log_expansion += one_multiplier * (arg - 1)**i / i
  pass# FOR i in range( 0, 3 )
  return log_expansion
pass# DEF taylor_series_log_expansion (  )
#
def atan_expand( arg ):
  return ( I * ( log( 1 - I * arg ) - log( 1 + I * arg ) ) ) / 2
pass# DEF atan_expand ( arg )
#
###############################################################################
#                                     TEST                                    #
###############################################################################
#
def atan_expand_test( arg ):
  #return diff( ( log( 1 - arg ) - log( 1 + arg ) ), z )
  return I * ( log( 1 - I * arg ) - log( 1 + I * arg ) ) / 2
  #return 1 / arg
pass# DEF atan_expand_test ( arg )
#
def log_complex_expand_test( arg1, arg2 ):
  return                                                                        (
         log( sqrt( ( re_collect( arg1 ) )**2 + ( im_collect( arg1 ) )**2 ) ) +
         I * atan_expand_test( im_collect( arg2 ) / re_collect( arg2 ) )
                                                                                )
pass# DEF log_complex_expand_test ( arg )
#
def asin_expand_deep_test( arg ):
  return (-I) * log_complex_expand_test( I*arg + sqrt( 1 - arg**2 ), I*arg + sqrt( -1 + arg**2 ) )
pass# DEF asin_expand_test ( arg )
#
def asinh_expand_deep_test ( arg ):
  return log_complex_expand_test( arg + sqrt( arg**2 + 1 ), arg + sqrt( arg**2 + 1 ) )
pass# DEF asinh_expand_deep_test ( arg )
#
###############################################################################
#^^^^^                                TEST                               ^^^^^#
###############################################################################
#
def log_complex_expand( arg ):
  return                                                                  (
         log( sqrt( (re_collect( arg ))**2 + (im_collect( arg ))**2 ) ) +
         I * atan_expand( im_collect( arg ) / re_collect( arg ) )
                                                                          )
pass# DEF log_complex_expand ( arg )
#
def asin_expand_deep( arg ):
  return (-I) * log_complex_expand( I*arg + sqrt( 1 - arg**2 ) )
pass# DEF asin_expand ( arg )
#
def asinh_expand(z):
  return log( z + sqrt( z**2 + 1 ) )
pass# DEF asinh_expand ( z )
#
def asinh_expand_only(z):
  return ( z + sqrt( z**2 + 1 ) )
pass# DEF asinh_expand_only ( z )
#
def asin_sum_argonly ( alpha, beta ) :
  return alpha * sqrt ( 1 - beta**2 ) + beta * sqrt ( 1 - alpha**2 )
pass# DEF asin_sum_argonly ( arg1, arg2 )
#
def log_argonly ( U ) :
  return lnA * ( log ( u2 ) - log ( u4 ) ) + lnB * ( log ( u3 ) - log ( u1 ) )
pass# DEF log_argonly ( U )
#
# Sets of arguments:
#
if ( (eArcsineArgTypeAB == enumArgTypeAB.BOTH) and (eArcsineArgType12 == enumArgType12.BOTH) ):
  iA1EnaFlag = 1
  iA2EnaFlag = 1
  iB1EnaFlag = 1
  iB2EnaFlag = 1
else:
  iA1EnaFlag =                                                 (
               1
               if
               (
                (eArcsineArgTypeAB == enumArgTypeAB.ATYPE) and
                (eArcsineArgType12 == enumArgType12.TYPE1)
               )
               else
               0
                                                               )
  iA2EnaFlag =                                                 (
               1
               if
               (
                (eArcsineArgTypeAB == enumArgTypeAB.ATYPE) and
                (eArcsineArgType12 == enumArgType12.TYPE2)
               )
               else
               0
                                                               )
  iB1EnaFlag =                                                 (
               1
               if
               (
                (eArcsineArgTypeAB == enumArgTypeAB.BTYPE) and
                (eArcsineArgType12 == enumArgType12.TYPE1)
               )
               else
               0
                                                               )
  iB2EnaFlag =                                                 (
               1
               if
               (
                (eArcsineArgTypeAB == enumArgTypeAB.BTYPE) and
                (eArcsineArgType12 == enumArgType12.TYPE2)
               )
               else
               0
                                                               )
pass# IF
pass# (
pass#   (eArcsineArgTypeAB == enumArgTypeAB.BOTH) and
pass#   (eArcsineArgType12 == enumArgType12.BOTH)
pass# )
#
# Arguments:
#
if ( bExpandArcsines ):
  if ( bDecomposeLogExpansion2Sum ):
    #
    ### TEMP ###asinArgonlyAu1 = asin_expand_deep( asin_argonly_A( u1 ) )### TEST ###
    ### TEMP ###asinArgonlyAu2 = asin_expand_deep( asin_argonly_A( u2 ) )### TEST ###
    ### TEMP ###asinArgonlyBu3 = asin_expand_deep( asin_argonly_B( u3 ) )### TEST ###
    ### TEMP ###asinArgonlyBu4 = asin_expand_deep( asin_argonly_B( u4 ) )### TEST ###
    #
    ###########################################################################
    #                                 TEST                                    #
    ###########################################################################
    #
    asinArgonlyAu1 = 0
    asinArgonlyAu2 = 0
    asinArgonlyBu3 = 0
    asinArgonlyBu4 = 0
    #
    #asinArgonlyAu1 = 1
    #
    #getSimpleTimingData( 'Solve start timer' )
    #asinArgonlyAu1 = solve( Eq( explicit_arg_diff_squares(), 0 ), z )
    #getSimpleTimingData( 'Solve finish timer' )
    #
    #asinArgonlyAu1 = asinh_expand_deep_test ( asinh_argonly_A( u1 ) )
    #asinArgonlyAu2 = asinh_expand_deep_test ( asinh_argonly_A( u2 ) )
    #asinArgonlyBu3 = asinh_expand_deep_test ( asinh_argonly_B( u3 ) )
    #asinArgonlyBu4 = asinh_expand_deep_test ( asinh_argonly_B( u4 ) )
    #
    asinArgonlyAu1 = asin_argonly_A( u1 )
    asinArgonlyAu2 = asin_argonly_A( u2 )
    asinArgonlyBu3 = asin_argonly_B( u3 )
    asinArgonlyBu4 = asin_argonly_B( u4 )
    #
    #getSimpleTimingData( 'Calculating nonintegrable begin timer' )
    #asinArgonlyAu1 = arg_nonintegrable()
    #getSimpleTimingData( 'Calculating nonintegrable end timer' )
    #getSimpleTimingData( 'Calculating series begin timer' )
    #asinArgonlyAu1 = taylor_series_log_expansion( asinArgonlyAu1, 3 )
    #getSimpleTimingData( 'Calculating series end timer' )
    #getSimpleTimingData( 'Expanding series begin timer' )
    #asinArgonlyAu1 = expand( asinArgonlyAu1 )
    #getSimpleTimingData( 'Expanding series end timer' )
    #getSimpleTimingData( 'Integrate series begin timer' )
    #asinArgonlyAu1 = integrate( asinArgonlyAu1, z )
    #getSimpleTimingData( 'Integrate series end timer' )
    #
    #asinArgonlyRoots = arg_nonintegrable_var_subst()
    #getSimpleTimingData( 'Solve timer end' )
    #
    ###########################################################################
    #^^^^^                            TEST                               ^^^^^#
    ###########################################################################
    #
  else:
    asinArgonlyAu1 = asin_expand( asin_argonly_A( u1 ) )
    asinArgonlyAu2 = asin_expand( asin_argonly_A( u2 ) )
    asinArgonlyBu3 = asin_expand( asin_argonly_B( u3 ) )
    asinArgonlyBu4 = asin_expand( asin_argonly_B( u4 ) )
  pass# IF ( bDecomposeLogExpansion2Sum )
else:
  asinArgonlyAu1 = asin_arg_A( u1 )
  asinArgonlyAu2 = asin_arg_A( u2 )
  asinArgonlyBu3 = asin_arg_B( u3 )
  asinArgonlyBu4 = asin_arg_B( u4 )
pass# IF ( bExpandArcsines )
#
V1 =                                  (
     (iIsLimitsDependsOnZ * 1) *
     (
       (
        iA2EnaFlag * asinArgonlyAu2 -
        iA1EnaFlag * asinArgonlyAu1
       ) +
       (
        iB1EnaFlag * asinArgonlyBu3 -
        iB2EnaFlag * asinArgonlyBu4
       )
     )
                                      )
#
Vu1 = asinh_expand_only( asinh_argonly_A( u1 ) )
Vu2 = asinh_expand_only( asinh_argonly_A( u2 ) )
Vu3 = asinh_expand_only( asinh_argonly_B( u3 ) )
Vu4 = asinh_expand_only( asinh_argonly_B( u4 ) )
#
if(False):
  if ( iIsLimitsDependsOnZ ):
    if ( bDoNotSimplify ):
      genform1 = sympify( Vu1 )
      genform2 = sympify( Vu2 )
      genform3 = sympify( Vu3 )
      genform4 = sympify( Vu4 )
    else:
      genform1 = simplify( sympify( Vu1 ) )
      # + simplify( powdenest( sympify( re(V3) ), force=True ) )
      genform2 = simplify( sympify( Vu2 ) )
      # + simplify( powdenest( sympify( re(V3) ), force=True ) )
      genform3 = simplify( sympify( Vu3 ) )
      # + simplify( powdenest( sympify( re(V3) ), force=True ) )
      genform4 = simplify( sympify( Vu4 ) )
      # + simplify( powdenest( sympify( re(V3) ), force=True ) )
    pass# IF ( bDoNotSimplify )
  else:
    if ( bDoNotSimplify ):
      genform1 = sympify( Vu1 )
      genform2 = sympify( Vu2 )
      genform3 = sympify( Vu3 )
      genform4 = sympify( Vu4 )
    else:
      genform1 = simplify( sympify( Vu1 ) )
      genform2 = simplify( sympify( Vu2 ) )
      genform3 = simplify( sympify( Vu3 ) )
      genform4 = simplify( sympify( Vu4 ) )
    pass# IF ( bDoNotSimplify )
  pass# IF ( iIsLimitsDependsOnZ )
pass# IF ( True )
#
### TEMP ###genform = simplify( powdenest( sympify( V1 ), force = True ) )### TEST ###
#
###############################################################################
#                                     TEST                                    #
###############################################################################
#
genform =                                                           (
          simplify
          (
           powdenest
           (
            simplify
            (
             asin_sum_argonly
             (
              asin_sum_argonly ( asinArgonlyAu2, -asinArgonlyAu1 ),
              asin_sum_argonly ( asinArgonlyBu3, -asinArgonlyBu4 )
             )
            )
            ,
            force = True
           )
          )
                                                                    )
#
#genform = simplify( powdenest( sympify( V1 ), force = True ) )
#
#genform1 = simplify( powdenest( sympify( asinArgonlyAu1 ), force = True ) )
#genform2 = simplify( powdenest( sympify( asinArgonlyAu2 ), force = True ) )
#genform3 = simplify( powdenest( sympify( asinArgonlyBu3 ), force = True ) )
#genform4 = simplify( powdenest( sympify( asinArgonlyBu4 ), force = True ) )
#
###############################################################################
#^^^^^                                TEST                               ^^^^^#
###############################################################################
#
OutputFile = open( 'test3.sympy.tex', 'w' )
OutputFile.write( '\\documentclass{article}\n' )
OutputFile.write( '\\usepackage{amsmath}\n' )
OutputFile.write( '\\begin{document}\n' )
if ( bPrintArcSines ):
  #
  ### TEMP ###printTexExpressionInline( OutputFile, genform )### TEST ###
  #
  #############################################################################
  #                                   TEST                                    #
  #############################################################################
  #
  printTexExpressionInline( OutputFile, genform )
  #
  #printTexExpressionInline( OutputFile, genform1 )
  #printTexExpressionInline( OutputFile, genform2 )
  #printTexExpressionInline( OutputFile, genform3 )
  #printTexExpressionInline( OutputFile, genform4 )
  #
  #for item in asinArgonlyRoots:
  #  printTexExpressionInline( OutputFile, item )
  #pass# FOR item in asinArgonlyRoots
  #
  #############################################################################
  #^^^^^                              TEST                               ^^^^^#
  #############################################################################
  #
  if ( bTry2IntegrateArcsines and iIsLimitsDependsOnZ ):
    print '\nAttempting to integrate arcsines ...\n'
    printTexExpressionInline     (
     OutputFile,
     integrate
     (
      simplify
      (
       powdenest
       (
        z**2 * asin_arg_A( u2 ),
        force = True
       )
      ),
      z
     )
                                 )
    printTexExpressionInline     (
     OutputFile,
     integrate
     (
      simplify
      (
       powdenest
       (
        z**2 * asin_arg_B( u2 ),
        force = True
       )
      ),
      z
     )
                                 )
  pass# IF ( bTry2IntegrateArcsines )
else:
  printTexExpressionInline( OutputFile, genform1 )
  printTexExpressionInline( OutputFile, genform2 )
  printTexExpressionInline( OutputFile, genform3 )
  printTexExpressionInline( OutputFile, genform4 )
pass# IF ( bPrintArcSines )
OutputFile.write( '\\end{document}' )
OutputFile.close()
getSimpleTimingData( 'Finish timer' )
#
