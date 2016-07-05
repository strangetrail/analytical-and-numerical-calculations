#!/usr/bin/env python
#
import sys
sys.path.append( '/home/thisuser/Documents/code/sympy_tests/common/' )
from utils import *
from mathutils import *
import math as pythmath
import random as pythrand
from sympy import *
#
###############################################################################
#                                    TYPES                                    #
###############################################################################
#
eTransProfiles = enum(TEM=1, TE=2, TM=3)
#
###############################################################################
#^^^^^                               TYPES                               ^^^^^#
###############################################################################
#
#
###############################################################################
#                                   SYMBOLS                                   #
###############################################################################
#
t, w, tau, infSymb, VARNONE, t_0, k, w_0, x, y, z =   (
 symbols
 (
  't, w, tau, infSymb, VARNONE, t_0, k, w_0, x, y, z'
 )
                                             )
chi_1221, chi_1122, chi_1212 =   (
 symbols
 (
  'chi_1221, chi_1122, chi_1212'
 )
                                 )
r_p, a_p, t_p, A_1, R_1, A_2, R_2 =   (
 symbols
 (
  'r_p, a_p, t_p, A_1, R_1, A_2, R_2'
 )
                                      )
M11, M12, M13, M21, M22, M23, M31, M32, M33 =   (
 symbols
 (
  'M11, M12, M13, M21, M22, M23, M31, M32, M33'
 )
                                                )
batch_max = 0
#
###############################################################################
#^^^^^                              SYMBOLS                              ^^^^^#
###############################################################################
#
#
###############################################################################
#                                   SETTINGS                                  #
###############################################################################
#
bSwitchExec2LibOnly     =  True
bQuasiSteadyState       =  True
bLaguerreSymbolic       = False
bDoNotEvaluateDiffs     = False
bDoNotEvaluateRatios    =  True
bUsePump2               = False
bPump1QWP               =  True# TODO # HERE
bPump2QWP               =  True# TODO
bProbeQWP               =  True# TODO
bIsNumericConstants     = False# Does w_0, k, chi_1221, chi_1122, chi_1212 an \
#                              # arbitrary numeric constants?                  
bSimplifyEachEquation   = False
bSendSimplifying2Maxima = False
bDEBUGdoubleSimplify    = False
#
iRotatePump1Pump2 = 0# 90 180 270 # TODO # HERE
iRotatePump1Probe = 0# 90 180 270 # TODO
#
rot180x =            (
          [
           [1,0, 0],
           [0,0,-1],
           [0,1, 0]
          ]
                     )
rot90z =            (
         [
          [0,-1,0],
          [1, 0,0],
          [0, 0,1]
         ]
                    )
rot180z =             (
          [
           [-1, 0,0],
           [ 0,-1,0],
           [ 0, 0,1]
          ]
                      )
rot270z =            (
          [
           [ 0,1,0],
           [-1,0,0],
           [ 0,0,1]
          ]
                     )
#
def conditional_rotate ( HE_beam, cond ) :
  if   ( cond == 90 ) :
    rot_matrix = rot90z
  elif ( cond == 180 ) :
    rot_matrix = rot180z
  elif ( cond == 270 ) :
    rot_matrix = rot270z
  else :
    rot_matrix = Ident_M
  pass# IFCASE ( cond )
  #
  if ( cond != 0 ) :
    return                                         (
           dot_prod_l2r_tens_vect
           (
            rot_matrix,
            HE_beam,
            bForceSimplify = bSimplifyEachEquation
           )
                                                   )
  else :
    return HE_beam
  pass# IF ( cond != 0 )
pass# DEF conditional_rotate ( HE_beam )
#
fPump2ZDisplacement = 0.0# TODO # HERE
fProbeZDisplacement = 0.0# TODO
#
eSwitch_mode_TE2TM = eTransProfiles.TEM
#
###############################################################################
#^^^^^                              SETTINGS                             ^^^^^#
###############################################################################
#
#
###############################################################################
#                                  CONSTANTS                                  #
###############################################################################
#
wavelength = 630e-9
#
if ( bIsNumericConstants ) :
  k = 2 * pythmath.pi / wavelength
  if ( bSimplifyEachEquation ) :
    k = simplify ( k )
  pass# IF ( bSimplifyEachEquation )
  w_0 = 1e-3
  chi_1chloronaphtallene = ( 1.632 * 3.04e+14 * 488e-9 ) / ( 4 * pythmath.pi )
  chi_1chloronaphtallene_diff = chi_1chloronaphtallene - chi_1chloronaphtallene * 0.05
  chi_1221 = chi_1chloronaphtallene_diff + chi_1chloronaphtallene * 0.1 * pythrand.random()
  chi_1122 = chi_1chloronaphtallene_diff + chi_1chloronaphtallene * 0.1 * pythrand.random()
  chi_1212 = chi_1chloronaphtallene_diff + chi_1chloronaphtallene * 0.1 * pythrand.random()
  # HERE!!!02.01.15-22:10!!! Assign correct values for 1-cloronaphtallene or anethole # DONE !!!02.10.15-23:10!!!
pass# IF ( bIsNumericConstants )
#
###############################################################################
#^^^^^                             CONSTANTS                             ^^^^^#
###############################################################################
#
#x = rho * cos ( phi )
#y = rho * sin ( phi )
#
z_r = k*(w_0**2)/2
if ( bSimplifyEachEquation ) :
  z_r = simplify ( z_r )
pass# IF ( bSimplifyEachEquation )
#
rho = sqrt( x**2 + y**2 )
if ( bSimplifyEachEquation ) :
  rho = simplify ( rho )
pass# IF ( bSimplifyEachEquation )
#
def r_vector( z_displ ):
  r_vector_res = sqrt( rho**2 + ( z + z_displ )**2 )
  if ( bSimplifyEachEquation ) :
    r_vector_res = simplify( r_vector_res )
  pass# IF ( bSimplifyEachEquation )
  #
  return r_vector_res
pass# DEF r_vector ( z_displ )
#
f = 1 / ( k * w_0 )
if ( bSimplifyEachEquation ) :
  f = simplify ( f )
pass# IF ( bSimplifyEachEquation )
#
phi = atan( y / x )
if ( bSimplifyEachEquation ) :
  phi = simplify ( phi )
pass# IF ( bSimplifyEachEquation )
#
L_ra = Function( 'L_ra' )
L = L_ra( r_p, a_p, t_p )
#
def E_TETM_comon ( radial, azimuth, z_displ ) :
  if ( bLaguerreSymbolic ) :
    E_TETM_comon_res =                                          (
             -(
               ( I**( 2*radial - azimuth + 1 ) * z_r ) /
               ( ( r_vector( z_displ ) )**2 * rho**2 )
              ) *
             (rho / ( sqrt(2) * r_vector(z_displ) * f ))**azimuth *
             L *
             exp
             (
              I*k*( r_vector( z_displ ) ) -
              rho**2 / ( 4* f**2 * ( r_vector( z_displ ) )**2 )
             ) *
             cos( azimuth * phi )
                                                                )
  else:
    E_TETM_comon_res =                                                  (
           -(
             ( I**( 2*radial - azimuth + 1 ) * z_r ) /
             ( ( r_vector( z_displ ) )**2 * rho**2 )
            ) *
           (rho / ( sqrt(2) * ( r_vector( z_displ ) ) * f ))**azimuth *
           assoc_laguerre
           (
            radial,
            azimuth,
            rho**2 / ( 2 * f**2 * ( r_vector( z_displ ) )**2 )
           ) *
           exp
           (
            I*k*( r_vector( z_displ ) ) -
            rho**2 / ( 4* f**2 * ( r_vector( z_displ ) )**2 )
           ) *
           cos( azimuth * phi )
                                                                        )
  pass# IF ( bLaguerreSymbolic )
  #
  if ( bSimplifyEachEquation ) :
    E_TETM_comon_res = simplify ( E_TETM_comon_res )
  pass# IF ( bSimplifyEachEquation )
  #
  return E_TETM_comon_res
pass# DEF E_TETM_comon ( radial, azimuth, z_displ )
#
def E_TM_comon_xyz ( radial, azimuth, z_displ ) :
  E_TM_comon_xyz_res = E_TETM_comon( radial, azimuth, z_displ ) * y * ( z + z_displ )
  if ( bSimplifyEachEquation ) :
    E_TM_comon_xyz_res = simplify ( E_TM_comon_xyz_res )
  pass# IF ( bSimplifyEachEquation )
  #
  return E_TM_comon_xyz_res
pass# DEF E_TM_comon_xyz ( radial, azimuth, z_displ )
#
def E_TE_common_xyz ( radial, azimuth, z_displ ) :
  E_TE_common_xyz_res = E_TETM_comon( radial, azimuth, z_displ ) * x
  if ( bSimplifyEachEquation ) :
    E_TE_common_xyz_res = simplify ( E_TE_common_xyz_res )
  pass# IF ( bSimplifyEachEquation )
  #
  return E_TE_common_xyz_res
pass# DEF E_TE_common_xyz ( radial, azimuth, z_displ )
#
def E_TE_comon ( radial, azimuth, z_displ ) :
  E_TE_common = [0 for i in xrange(3)]
  E_TE_common[0] = E_TE_common_xyz( radial, azimuth, z_displ ) * y
  E_TE_common[1] = -E_TE_common_xyz( radial, azimuth, z_displ ) * x
  E_TE_common[2] = 0
  if ( bSimplifyEachEquation ) :
    for i in range(0,3):
      E_TE_common[i] = simplify ( E_TE_common[i] )
    pass# FOR i in range(0,3):
  pass# IF ( bSimplifyEachEquation )
  #
  return E_TE_common
pass# DEF E_TE_comon ( radial, azimuth, z_displ )
#
def E_TM_comon ( radial, azimuth, z_displ ) :
  E_TM_common = [0 for i in xrange(3)]
  E_TM_common[0] = E_TE_common_xyz( radial, azimuth, z_displ ) * x * ( z + z_displ )
  E_TM_common[1] = E_TE_common_xyz( radial, azimuth, z_displ ) * y * ( z + z_displ )
  E_TM_common[2] = -E_TE_common_xyz( radial, azimuth, z_displ ) * rho**2
  if ( bSimplifyEachEquation ) :
    for i in range(0,3):
      E_TM_common[i] = simplify ( E_TM_common[i] )
    pass# FOR i in range(0,3):
  pass# IF ( bSimplifyEachEquation )
  #
  return E_TM_common
pass# DEF E_TM_comon ( radial, azimuth, z_displ )
#
def E_TEM_common( radial, azimuth, z_displ ) :
  E_TM = E_TM_comon( radial, azimuth, z_displ )
  E_TE = E_TE_comon( radial, azimuth, z_displ )
  #
  E_TEM = [0 for i in xrange(3)]
  for i in range(0,3):
    E_TEM[i] = E_TM[i] + E_TE[i]
    if ( bSimplifyEachEquation ) :
      E_TEM[i] = simplify ( E_TEM[i] )
    pass# IF ( bSimplifyEachEquation )
  pass# FOR i in range(0,3)
  #
  return E_TEM
pass# DEF E_TEM_common
#
if ( bIsNumericConstants ) :
  exp_QWP = exp( I*pythmath.pi/2 )
else:
  exp_QWP = exp( I*pi/2 )
pass# IF ( bIsNumericConstants )
if ( bSimplifyEachEquation ) :
  exp_QWP = simplify ( exp_QWP )
pass# IF ( bSimplifyEachEquation )
#
QWP =                  (
      [
       [exp_QWP, 0],
       [0, -I*exp_QWP]
      ]
                       )
#
if ( bSimplifyEachEquation ) :
  for i in range(0,2):
    for j in range(0,2):
      QWP[i][j] = simplify( QWP[i][j] )
    pass# FOR j in range(0,2):
  pass# FOR i in range(0,2)
pass# IF ( bSimplifyEachEquation )
#
def jones_calc( M, V ) :
  U = [0 for i in xrange(3)]
  for i in range(0,2):
    for j in range(0,2):
      U[i] = U[i] + V[j] * M[i][j]
      if ( bSimplifyEachEquation ) :
        U[i] = simplify ( U[i] )
      pass# IF ( bSimplifyEachEquation )
    pass# FOR j in range(0,2)
  pass# FOR i in range(0,2)
  U[2] = V[2]
  if ( bSimplifyEachEquation ) :
    U[2] = simplify( U[2] )
  pass# IF ( bSimplifyEachEquation )
  #
  return U
pass# DEF jones_calc( M, V )
#
idx_azimuthal_pump = 1
idx_radial_pump = 0
idx_azimuthal_pump_back = -1
idx_radial_pump_back = 0
idx_azimuthal_probe = 2
idx_radial_probe = 0
#
if ( bUsePump2 ) :
  if ( eSwitch_mode_TE2TM == eTransProfiles.TEM ) :
    E_pump_backward =                             (
                      jones_calc
                      (
                       QWP
                       ,
                       E_TEM_common
                       (
                        idx_radial_pump_back,
                        idx_azimuthal_pump_back,
                        fPump2ZDisplacement
                       )
                      )
                                                  )
    print 'I\'m here, where are TEM modes are calculated.'
  elif ( eSwitch_mode_TE2TM == eTransProfiles.TM ):
    E_pump_backward =                             (
                      jones_calc
                      (
                       QWP,
                       E_TM_comon
                       (
                        idx_radial_pump_back,
                        idx_azimuthal_pump_back,
                        fPump2ZDisplacement
                       )
                      )
                                                  )
    print 'I\'m here, where are TM modes are calculated.'
  elif ( eSwitch_mode_TE2TM == eTransProfiles.TE ):
    E_pump_backward =                             (
                      jones_calc
                      (
                       QWP,
                       E_TE_comon
                       (
                        idx_radial_pump_back,
                        idx_azimuthal_pump_back,
                        fPump2ZDisplacement
                       )
                      )
                                                  )
    print 'I\'m here, where are TE modes are calculated.'
  pass# IFCASE ( eSwitch_mode_TE2TM )
  #
  E_pump_backward =                                         (
                    dot_prod_l2r_tens_vect
                    (
                     rot180x,
                     E_pump_backward,
                     bForceSimplify = bSimplifyEachEquation
                    )
                                                            )
  #
  E_pump_backward = conditional_rotate ( E_pump_backward, iRotatePump1Pump2 )
pass# IF ( bUsePump2 )
#
if ( eSwitch_mode_TE2TM == eTransProfiles.TEM ):
  E_pump_forward = jones_calc( QWP, E_TEM_common( idx_radial_pump, idx_azimuthal_pump, 0.0 ) )
  print 'I\'m here, where are TEM modes are calculated.'
elif ( eSwitch_mode_TE2TM == eTransProfiles.TM ):
  E_pump_forward = jones_calc( QWP, E_TM_comon( idx_radial_pump, idx_azimuthal_pump, 0.0 ) )
  print 'I\'m here, where are TM modes are calculated.'
elif ( eSwitch_mode_TE2TM == eTransProfiles.TE ):
  E_pump_forward = jones_calc( QWP, E_TE_comon( idx_radial_pump, idx_azimuthal_pump, 0.0 ) )
  print 'I\'m here, where are TE modes are calculated.'
pass# IFCASE ( eSwitch_mode_TE2TM )
#
if ( bUsePump2 ):
  E_pump = [a+b for a,b in zip( E_pump_backward, E_pump_forward )]
  if ( bSimplifyEachEquation ) :
    for i in range(0,3):
      E_pump[i] = simplify ( E_pump[i] )
    pass# FOR i in range(0,2):
  pass# IF ( bSimplifyEachEquation )
else:
  #E_pump = [a for a in zip( E_pump_forward )]
  E_pump = E_pump_forward
pass# IF ( bUsePump2 )
#
if ( eSwitch_mode_TE2TM == eTransProfiles.TEM ):
  E_probe =                                                                          (
           jones_calc
           (
            QWP,
            E_TEM_common( idx_radial_probe, idx_azimuthal_probe, fProbeZDisplacement )
           )
                                                                                     )
  print 'I\'m here, where are TEM modes are calculated.'
elif ( eSwitch_mode_TE2TM == eTransProfiles.TM ):
  E_probe =                                                                        (
           jones_calc
           (
            QWP,
            E_TM_comon( idx_radial_probe, idx_azimuthal_probe, fProbeZDisplacement )
           )
                                                                                   )
  print 'I\'m here, where are TM modes are calculated.'
elif ( eSwitch_mode_TE2TM == eTransProfiles.TE ):
  E_probe =                                                                        (
           jones_calc
           (
            QWP,
            E_TE_comon( idx_radial_probe, idx_azimuthal_probe, fProbeZDisplacement )
           )
                                                                                   )
  print 'I\'m here, where are TE modes are calculated.'
pass# IFCASE ( eSwitch_mode_TE2TM )
#
E_probe = conditional_rotate ( E_probe, iRotatePump1Probe )
#
# Optical Kerr quasi steady-state limit:
chi_1111 = chi_1221 + chi_1122 + chi_1212
#
def x_3rd_order (  ):
  x_3rd = [[[[0 for i in xrange(3)] for j in xrange(3)] for k in xrange(3)] for p in xrange(3)]
  for i in range(0,3):
    for j in range(0,3):
      for k in range(0,3):
        for p in range(0,3):
          if ( (i == j) and (j == k) and (k == p) ):
            x_3rd[i][j][k][p] = chi_1111
          else:
            if ( (i == j) and (k == p) ):
              x_3rd[i][j][k][p] = chi_1122
            else:
              if ( (i == k) and (j == p) ):
                x_3rd[i][j][k][p] = chi_1212
              else:
                if ( (i == p) and (j == k) ):
                  x_3rd[i][j][k][p] = chi_1221
                else:
                  x_3rd[i][j][k][p] = 0
                pass# IF ( (i == p) and (j == k) )
              pass# IF ( (i == k) and (j == p) )
            pass# IF ( (i == j) and (k == p) )
          pass# IF ( (i == j) and (j == k) and (k == p) )
        pass# FOR p in range(0,3)
      pass# FOR k in range(0,3)
    pass# FOR j in range(0,3)
  pass# FOR i in range(0,3)
  #
  return x_3rd
pass# DEF x_3rd_order()
#
x_3rd = x_3rd_order()
#
def x_2nd_order_rel ():
  x_2nd_ord = [[0 for i in xrange(3)] for j in xrange(3)]
  for i in range(0,3):
    P_3nd_ord_i = 0;
    for j in range(0,3):
      for k in range(0,3):
        for p in range(0,3):
          x_2nd_ord[i][j] += x_3rd[i][j][k][p] * E_pump[k] * conjugate( E_pump[p] )
          if ( bSimplifyEachEquation ) :
            x_2nd_ord[i][j] = simplify( x_2nd_ord[i][j] )
          pass# IF ( bSimplifyEachEquation )
        pass# FOR
      pass# FOR
    pass# FOR
  pass# FOR
  #
  return x_2nd_ord
pass# DEF
#
#def x_2nd_order_rel ():
#  x_2nd_ord = [[0 for i in xrange(3)] for j in xrange(3)]
#  for i in range(0,3):
#    P_3nd_ord_i = 0;
#    for j in range(0,3):
#      for k in range(0,3):
#        for p in range(0,3):
#          P_3nd_ord_i += x_3rd[i][j][k][p] * E_probe[j] * E_pump[k] * conjugate( E_pump[p] )
#          if ( bSimplifyEachEquation ) :
#            P_3nd_ord_i = simplify( P_3nd_ord_i )
#          pass# IF ( bSimplifyEachEquation )
#        pass# FOR
#      pass# FOR
#    pass# FOR
#    for j in range(0,3):
#      x_2nd_ord[i][j] = P_3nd_ord_i / E_probe[j];
#      if ( bSimplifyEachEquation ) :
#        x_2nd_ord[i][j] = simplify( x_2nd_ord[i][j] )
#      pass# IF ( bSimplifyEachEquation )
#    pass# FOR
#  pass# FOR
#  #
#  return x_2nd_ord
#pass# DEF
#
def get_x_2nd():
  x_2nd = x_2nd_order_rel()
  for i in range(0,3):
    for j in range(0,3):
      if ( bQuasiSteadyState ):
        print '\nCalculating x_2nd['+repr(i)+']['+repr(j)+']:'
        getSimpleTimingData ( "start timer" )
        #
        if ( bSendSimplifying2Maxima ) :
          printMaximaRedirectMessage ()
          if ( not bDEBUGdoubleSimplify ) :
            ###x_2nd[i][j] = x_2nd[i][j]
            x_2nd[i][j] = simplifySingleExpressioninMaxima( x_2nd[i][j] )
          else :
            x_2nd[i][j] = simplify(simplifySingleExpressioninMaxima( x_2nd[i][j] ) - simplify( powdenest(x_2nd[i][j], force=True) ))### DEBUG
          pass# IF ( bSendSimplifying2Maxima )
        else :
          x_2nd[i][j] = simplify( powdenest(x_2nd[i][j], force=True) )
          #x_2nd[i][j] = x_2nd[i][j]
        pass# IF ( bSendSimplifying2Maxima )
        #
        getSimpleTimingData ( "finish timer" )
        print 'Done calculating x_2nd['+repr(i)+']['+repr(j)+'].'
      else:
        chi_integral = integrate                                                   (
                        x_2nd[i][j] * exp( -I * w * t_0 ) * exp( -(t - t_0)/tau ),
                        (t_0, -oo, t)
                                                                                   )
        print '\nCalculating x_2nd['+repr(i)+']['+repr(j)+'] non quasi steady-state case:'
        getSimpleTimingData ( "start timer" )
        #
        if ( bSendSimplifying2Maxima ) :
          printMaximaRedirectMessage ()
          if ( not bDEBUGdoubleSimplify ) :
            ###x_2nd[i][j] = chi_integral
            x_2nd[i][j] = simplifySingleExpressioninMaxima( chi_integral )
          else :
            x_2nd[i][j] = simplify( powdenest(simplifySingleExpressioninMaxima( chi_integral ), force=True) )### DEBUG
          pass# IF ( not bDEBUGdoubleSimplify )
        else :
          x_2nd[i][j] = simplify( powdenest( chi_integral, force=True ) )
        pass# IF ( bSendSimplifying2Maxima )
        #
        getSimpleTimingData ( "finish timer" )
        print 'Done calculating x_2nd['+repr(i)+']['+repr(j)+'] non quasi steady-state case.'
      pass# IF ( bQuasiSteadyState )
    pass# FOR j in range(0,3)
  pass# FOR i in range(0,3)
  return x_2nd
pass# DEF get_x_2nd
#
#@cython.boundscheck(False)
def do_main():
  #
  M = [[M11, M12, M13],[M21, M22, M23],[M31, M32, M33]]
  # AB consists of:
  # 1st row with 3 pairs are off-diagonal terms with respect to main diagonal \
  #     elements;
  # 2nd row with 3 pairs are main diagonal elements;
  # 3rd row with 4 pairs are term 00 compared with off-diagonal terms with    \
  #     respect to main diagonal elements;
  # 4th row with 4 pairs are off-diagonal terms with respect to main diagonal \
  #     elements compared to each other;
  # 5th row with 3 pairs are off-diagonal terms with respect to auxilary      \
  #     diagonal elements;
  # 6th row with 4 pairs are off-diagonal terms with respect to auxilary      \
  #     diagonal elements compared to each other;
  # 7th row with 4 pairs are off-diagonal terms from bottom half-matrix with  \
  #     respect to main diagonal elements compared to each other;
  # 8th row with 4 pairs are off-diagonal terms from bottom half-matrix with  \
  #     respect to auxilary diagonal elements compared to each other;
  #AB_count = 35
  #AB =                                                  (
  #     [
  #      [
  #       [0, 1], [0, 2], [1, 2],
  #       [0, 0], [1, 1], [0, 0],
  #       [0, 0], [0, 0], [0, 0], [0, 0],
  #       [0, 1], [0, 2], [0, 1], [0, 1],
  #       [0, 1], [0, 0], [1, 0],
  #       [0, 1], [0, 0], [0, 1], [0, 1],
  #       [1, 0], [2, 0], [1, 0], [1, 0],
  #       [1, 2], [2, 1], [1, 2], [1, 2],
  #       [0, 1], [0, 2], [1, 0], [2, 0], [1, 2], [2, 1]
  #      ],
  #      [
  #       [1, 0], [2, 0], [2, 1],
  #       [1, 1], [2, 2], [2, 2],
  #       [0, 1], [0, 2], [1, 2], [1, 1],
  #       [0, 2], [1, 2], [1, 2], [0, 2],
  #       [1, 2], [2, 2], [2, 1],
  #       [0, 0], [1, 0], [1, 0], [0, 0],
  #       [2, 0], [2, 1], [2, 1], [2, 0],
  #       [2, 1], [2, 2], [2, 2], [2, 1],
  #       [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]
  #      ]
  #     ]
  #                                                      )
  #
  chi_2nd_differences = [VARNONE for i in xrange(3*3*3*3)]
  chi_2nd_ratios      = [VARNONE for i in xrange(3*3*3*3)]
  #
  x_2nd_result = get_x_2nd ()
  #
  for i in range ( 0, 3 ) :
    for j in range ( 0, 3 ) :
      if ( ( i < 2 ) or ( j < 2 ) ) :
        if ( ( j == 2 ) and ( i < 2 ) ) :
          k_min = i + 1
        else:
          k_min = i
        pass# IF ( ( j == 2 ) and ( i < 2 ) )
        #
        for k in range ( k_min, 3 ) :
          if ( ( j < 2 ) and ( k == i ) ) :
            l_min = j + 1
          else :
            l_min = 0
          pass# IF ( ( j < 2 ) and ( k == i ) )
          #
          for l in range ( l_min, 3 ) :
            if ( not bDoNotEvaluateRatios ) :
              print ( '\nCalculating ratio\n' )
              print ( 'CHI['+repr(i)+']['+repr(j)+']   over   CHI['+repr(k)+']['+repr(l)+']:\n' )
              getSimpleTimingData( "Ratio evaluation start" )
              #
              if                                                              (
                 ( i == 2 ) and ( j == 1 ) and ( k == 2 ) and ( l == 2  ) and
                 ( not bSendSimplifying2Maxima )
                                                                              ) :
                print ( "=====SIMPLIFICATION SKIPPED=====" )
                chi_2nd_ratios[27*i+9*j+3*k+l] =                       (
                                                  x_2nd_result[i][j] /
                                                  x_2nd_result[k][l]
                                                                       )
              else :
                if ( bSendSimplifying2Maxima ) :
                  printMaximaRedirectMessage ()
                  if ( not bDEBUGdoubleSimplify ) :
                    chi_2nd_ratios[27*i+9*j+3*k+l] =                                  (
                                                     simplifySingleExpressioninMaxima
                                                     (
                                                       x_2nd_result[i][j] /
                                                       x_2nd_result[k][l]
                                                      )
                                                                                      )
                  else :
                    chi_2nd_ratios[27*i+9*j+3*k+l] =                          simplify(powdenest(
                                                     simplifySingleExpressioninMaxima
                                                     (
                                                      x_2nd_result[i][j] /
                                                      x_2nd_result[k][l]
                                                     )
                                                                                      ), force=True)### DEBUG
                  pass# IF ( not bDEBUGdoubleSimplify )
                else :
                  chi_2nd_ratios[27*i+9*j+3*k+l] =                       (
                                                   simplify
                                                   (
                                                    x_2nd_result[i][j] /
                                                    x_2nd_result[k][l]
                                                   )
                                                                         )
                pass# IF ( bSendSimplifying2Maxima )
              pass# IF ( ( i == 2 ) and ( j >= 0 ) and ( k == 2 ) and ( l == 2  ) )
              #
              getSimpleTimingData( "Ratio evaluation finish" )
              print ( 'Done calculating ratio' )
              print ( 'chi['+repr(i)+']['+repr(j)+'] over chi['+repr(k)+']['+repr(l)+'].\n' )
            pass# IF ( not bDoNotEvaluateRatios )
            #
            if ( not bDoNotEvaluateDiffs ):
              print ( '\nCalculating difference between\n' )
              print ( 'CHI['+repr(i)+']['+repr(j)+']   and   CHI['+repr(k)+']['+repr(l)+']:\n' )
              getSimpleTimingData( "Difference evaluation start" )
              #
              if ( bSendSimplifying2Maxima ) :
                printMaximaRedirectMessage ()
                if ( not bDEBUGdoubleSimplify ) :
                  chi_2nd_differences[27*i+9*j+3*k+l] =                                  (
                                                        simplifySingleExpressioninMaxima
                                                        (
                                                         x_2nd_result[i][j] -
                                                         x_2nd_result[k][l]
                                                        )
                                                                                         )
                else :
                  chi_2nd_differences[27*i+9*j+3*k+l] =               simplify(powdenest(
                                                        simplifySingleExpressioninMaxima
                                                        (
                                                         x_2nd_result[i][j] -
                                                         x_2nd_result[k][l]
                                                        )
                                                                                    ), force=True)### DEBUG
                pass# IF ( not bDEBUGdoubleSimplify )
              else :
                chi_2nd_differences[27*i+9*j+3*k+l] =                       (
                                                      simplify
                                                      (
                                                       x_2nd_result[i][j] -
                                                       x_2nd_result[k][l]
                                                      )
                                                                            )
              pass# IF ( bSendSimplifying2Maxima )
              #
              getSimpleTimingData( "Difference evaluation finish" )
              print ( 'Done calculating difference between' )
              print ( 'chi['+repr(i)+']['+repr(j)+'] and chi['+repr(k)+']['+repr(l)+'].\n' )
            pass# IF ( not bDoNotEvaluateDiffs )
          pass#FOR l in range(0,3)
        pass#FOR k in range(0,3)
      pass# IF ( ( i < 2 ) or ( j < 2 ) )
    pass# FOR j in range(0,3)
  pass# FOR i in range(0,3)
  f = open( 'test4.sympy.tex', 'w' )
  f.write( '\\documentclass{article}\n\\usepackage{amsmath}\n' )
  f.write( '\\usepackage{mathtools}\n\\begin{document}\n' )
  x_2nd_all = 0;
  for i in range(0,3):
    for j in range(0,3):
      f.write ( '\\text{Chi ['+repr(i)+'] ['+repr(j)+'] :}\\newline\n' )
      printTexExpressionInline ( f, x_2nd_result[i][j] )
      #x_2nd_all += x_2nd_result[i][j] * M[i][j]
    pass# FOR
  pass# FOR
  #print ( '\n Calculating sum of all together: ' )
  #x_2nd_all_simplified = simplify( x_2nd_all )
  #print ( '\n Done calculating sum of all together. ' )
  #f.write( '\\text{All together: }\\newline\n' )
  #f.write( '\\begin{equation}\n\\begin{gathered}\n\\begin{multlined}\n' )
  #f.write( latex( x_2nd_all_simplified, mode = 'inline' ) )
  #f.write( '\n\\end{multlined}\n\\end{gathered}\n\\end{equation}\n' )
  #
  f.write ( '\\text{E Pump:}\\newline\n' )
  printTexExpressionInline ( f, E_pump )
  #
  f.write ( '\\text{E Probe:}\\newline\n' )
  printTexExpressionInline ( f, E_probe )
  #
  for i in range(0,3):
    for j in range(0,3):
      for k in range(0,3):
        for l in range(0,3):
          if ( not bDoNotEvaluateDiffs ):
            f.write                                           (
               '\\text{Chi['+repr(i)+']['+repr(j)+'] -' +
               ' Chi['+repr(k)+']['+repr(l)+']: }\\newline\n'
                                                              )
            printTexExpressionInline( f, chi_2nd_differences[27*i+9*j+3*k+l] )
          pass# IF ( not bDoNotEvaluateDiffs )
          f.write                                           (
             '\\text{Chi['+repr(i)+']['+repr(j)+'] /' +
             ' Chi['+repr(k)+']['+repr(l)+']: }\\newline\n'
                                                            )
          printTexExpressionInline( f, chi_2nd_ratios[27*i+9*j+3*k+l] )
        pass#FOR l in range(0,3)
      pass#FOR k in range(0,3)
    pass# FOR j in range(0,3)
  pass# FOR i in range(0,3)
  f.write( '\\end{document}' )
  f.close()
  pprint( x_3rd )
#
if ( not bSwitchExec2LibOnly ):
  beginSimpleTiming()
  do_main()
  getSimpleTimingData( 'Calculation timer' )
  getTotalTime ()
pass# IF
#
