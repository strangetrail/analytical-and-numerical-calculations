#!/usr/bin/env python
#
# test.pyx
# cython: profile=True
#
import timeit
#import sys
#sys.path.append( '/home/nixuser/Documents/code/sympy_tests/common/' )
#from utils import *
#
# utils INLINE
#
def enum(**enums):
  return type('Enum', (), enums)
#
bgn_timestamp = timeit.default_timer()
end_timestamp = timeit.default_timer()
#
def getSimpleTimingData( timingMessageText ):
  global end_timestamp
  global bgn_timestamp
  end_timestamp = timeit.default_timer()
  print timingMessageText + ':  '
  print end_timestamp - bgn_timestamp
  bgn_timestamp = timeit.default_timer()
pass# DEF
#
def beginSimpleTiming():
  getSimpleTimingData( 'Simple test timer' )
pass# DEF
#
def printTexExpressionInline( texfile, expr ):
  texfile.write( '\\begin{equation}\n\\begin{split}\n' )
  texfile.write( latex( expr, mode = 'inline' ) )
  texfile.write( '\n' )
  texfile.write( '\\end{split}\n\\end{equation}\n' )
pass# DEF
#
#from mathutils import *
#
# mathutils INLINE
#
#
# Dot prod.:
#
def dot_product(U,V):
  ##return simplify( sympify( U[0] * V[0] + U[1] * V[1] + U[2] * V[2] ) )
  return U[0] * V[0] + U[1] * V[1] + U[2] * V[2]
pass# DEF dot_product ( U, V )
#
# Dot prod. from l2r of vect. with 2nd. ord. tens.:
#
def dot_prod_l2r_tens_vect( V, M ):
  R = [0 for i in xrange(3)]
  U = [0 for i in xrange(3)]
  for i in range(0,3):
    for j in range(0,3):
      R[j] = M[i][j]
    pass# FOR j in range(0,3)
    U[i] = dot_product( R, V )
  pass# FOR i in range(0,3)
  return U
#
import math as pythmath
import random as pythrand
from sympy import *
#
def the_main_routine():
  #
  eTransProfiles = enum(TEM=1, TE=2, TM=3)
  #
  #############################################################################
  #                                 SYMBOLS                                   #
  #############################################################################
  #
  t, w, tau, infSymb, t_0, k, w_0, x, y, z =   (
   symbols
   (
    't, w, tau, infSymb, t_0, k, w_0, x, y, z'
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
  bSwitchExec2LibOnly   = False
  bQuasiSteadyState     =  True
  bLaguerreSymbolic     = False
  bSubstituteValues     =  True
  bCheckEquals          = False
  bDoNotEvaluateDiffs   =  True
  bUsePump2             = False
  bPump1QWP             =  True# TODO # HERE
  bPump2QWP             =  True# TODO
  bProbeQWP             =  True# TODO
  bIsNumericConstants   =  True# Is w_0, k, chi_1221, chi_1122, chi_1212 an     \
  #                            # arbitrary numeric constants?
  bSimplifyEachEquation = False#
  bDoNotDenestPowers    =  True
  #
  iRotatePump1Pump2 = 0# 90 180 270 # TODO # HERE
  iRotatePump1Probe = 0# 90 180 270 # TODO
  #
  fPump2ZDisplacement = 0.0# TODO # HERE
  fProbeZDisplacement = 0.0# TODO
  #
  eSwitch_probe_TE2TM = eTransProfiles.TEM
  eSwitch_pump_TE2TM  = eTransProfiles.TEM
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
               (rho / ( sqrt(2) * r * f ))**azimuth *
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
  idx_azimuthal_probe = 3
  idx_radial_probe = 0
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
  if ( bUsePump2 ) :
    if ( eSwitch_probe_TE2TM == eTransProfiles.TEM ) :
      E_pump_backward =                             (
                        jones_calc
                        (
                         QWP
                         ,
                         E_TEM_common
                         (
                          idx_radial_probe_back,
                          idx_azimuthal_probe_back,
                          fPump2ZDisplacement
                         )
                        )
                                                    )
      print 'I\'m here, where are TEM modes are calculated.'
    elif ( eSwitch_probe_TE2TM == eTransProfiles.TM ):
      E_pump_backward =                             (
                        jones_calc
                        (
                         QWP,
                         E_TM_comon
                         (
                          idx_radial_probe_back,
                          idx_azimuthal_probe_back,
                          fPump2ZDisplacement
                         )
                        )
                                                    )
      print 'I\'m here, where are TM modes are calculated.'
    elif ( eSwitch_probe_TE2TM == eTransProfiles.TE ):
      E_pump_backward =                             (
                        jones_calc
                        (
                         QWP,
                         E_TE_comon
                         (
                          idx_radial_probe_back,
                          idx_azimuthal_probe_back,
                          fPump2ZDisplacement
                         )
                        )
                                                    )
      print 'I\'m here, where are TE modes are calculated.'
    pass# IFCASE ( eSwitch_probe_TE2TM )
    #
    E_pump_backward =                                         (
                      dot_prod_l2r_tens_vect
                      (
                       rot180x,
                       E_pump_backward,
                       bForceSimplify = bSimplifyEachEquation
                      )
                                                              )
  pass# IF ( bUsePump2 )
  #
  if ( eSwitch_probe_TE2TM == eTransProfiles.TEM ):
    E_pump_forward = jones_calc( QWP, E_TEM_common( idx_radial_probe, idx_azimuthal_probe, 0.0 ) )
    print 'I\'m here, where are TEM modes are calculated.'
  elif ( eSwitch_probe_TE2TM == eTransProfiles.TM ):
    E_pump_forward = jones_calc( QWP, E_TM_comon( idx_radial_probe, idx_azimuthal_probe, 0.0 ) )
    print 'I\'m here, where are TM modes are calculated.'
  elif ( eSwitch_probe_TE2TM == eTransProfiles.TE ):
    E_pump_forward = jones_calc( QWP, E_TE_comon( idx_radial_probe, idx_azimuthal_probe, 0.0 ) )
    print 'I\'m here, where are TE modes are calculated.'
  pass# IFCASE ( eSwitch_probe_TE2TM )
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
  if ( eSwitch_pump_TE2TM == eTransProfiles.TEM ):
    E_probe =                                                                          (
             jones_calc
             (
              QWP,
              E_TEM_common( idx_radial_pump, idx_azimuthal_pump, fProbeZDisplacement )
             )
                                                                                       )
    print 'I\'m here, where are TEM modes are calculated.'
  elif ( eSwitch_pump_TE2TM == eTransProfiles.TM ):
    E_probe =                                                                        (
             jones_calc
             (
              QWP,
              E_TM_comon( idx_radial_pump, idx_azimuthal_pump, fProbeZDisplacement )
             )
                                                                                     )
    print 'I\'m here, where are TM modes are calculated.'
  elif ( eSwitch_pump_TE2TM == eTransProfiles.TE ):
    E_probe =                                                                        (
             jones_calc
             (
              QWP,
              E_TE_comon( idx_radial_pump, idx_azimuthal_pump, fProbeZDisplacement )
             )
                                                                                     )
    print 'I\'m here, where are TE modes are calculated.'
  pass# IFCASE ( eSwitch_pump_TE2TM )
  #
  # Optical Kerr quasi steady-state linit:
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
            P_3nd_ord_i += x_3rd[i][j][k][p] * E_probe[j] * E_pump[k] * conjugate( E_pump[p] )
            if ( bSimplifyEachEquation ) :
              P_3nd_ord_i = simplify( P_3nd_ord_i )
            pass# IF ( bSimplifyEachEquation )
          pass# FOR
        pass# FOR
      pass# FOR
      for j in range(0,3):
        x_2nd_ord[i][j] = P_3nd_ord_i / E_probe[j];
        if ( bSimplifyEachEquation ) :
          x_2nd_ord[i][j] = simplify( x_2nd_ord[i][j] )
        pass# IF ( bSimplifyEachEquation )
      pass# FOR
    pass# FOR
    #
    return x_2nd_ord
  pass# DEF
  #
  def get_x_2nd():
    x_2nd = x_2nd_order_rel()
    for i in range(0,3):
      for j in range(0,3):
        if ( bQuasiSteadyState ):
          print '\n Calculating x_2nd['+repr(i)+']['+repr(j)+']: \n'
          if ( bDoNotDenestPowers ) :
            x_2nd[i][j] = simplify( x_2nd[i][j] )
          else:
            x_2nd[i][j] = simplify( powdenest( x_2nd[i][j], force=True ) )
          pass# IF ( bDoNotDenestPowers )
          print '\n Done calculating x_2nd['+repr(i)+']['+repr(j)+']. \n'
        else:
          chi_integral = integrate                                                   (
                          x_2nd[i][j] * exp( -I * w * t_0 ) * exp( -(t - t_0)/tau ),
                          (t_0, -oo, t)
                                                                                     )
          print '\n Calculating x_2nd['+repr(i)+']['+repr(j)+'] non quasi steady-state case: \n'
          if ( bDoNotDenestPowers ) :
            x_2nd[i][j] = simplify( powdenest( chi_integral, force=True ) )
          else:
            x_2nd[i][j] = simplify( chi_integral )
          print '\n Done calculating x_2nd['+repr(i)+']['+repr(j)+'] non quasi steady-state case. \n'
        pass# IF ( bQuasiSteadyState )
      pass# FOR j in range(0,3)
    pass# FOR i in range(0,3)
    return x_2nd
  pass# DEF get_x_2nd
  #
  #@cython.boundscheck(False)
  def do_main():
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
    chi_2nd_differences = [0 for i in xrange(3*3*3*3)]
    chi_2nd_ratios = [0 for i in xrange(3*3*3*3)]
    x_2nd_result = get_x_2nd()
    for i in range(0,3):
      for j in range(0,3):
        if ( ( not ( j < 2 ) ) and ( i < 2 ) ) :
          k_min = i + 1
        else:
          k_min = i
        pass# IF ( not ( j < 2 ) ) and ( i < 2 )
        #
        for k in range(k_min,3):
          if ( j == 2 ) :
            print "\nWierd...\n"
          pass# IF ( j < 2 )
          #
          if ( ( i < 2 ) and ( j < 2 ) ) :
            if ( ( k == i ) and ( j < 2 ) ) :
              l_min = j + 1
            else:
              l_min = 0
            pass# IF ( ( k == i ) and ( j < 2 ) )
          pass# IF ( ( i < 2 ) and ( j < 2 ) )
          #
          for l in range(l_min,3):
            print                                                                              (
                  '\n Calculating ' +
                  'ratio chi['+repr(i)+']['+repr(j)+'] over chi['+repr(k)+']['+repr(l)+']: \n'
                                                                                               )
            getSimpleTimingData( "Difference evaluation start" )
            chi_2nd_ratios[27*i+9*j+3*k+l] =                       (
                                             simplify
                                             (
                                              x_2nd_result[i][j] /
                                              x_2nd_result[k][l]
                                             )
                                                                   )
            getSimpleTimingData( "Difference evaluation finish" )
            print                                                                              (
                  '\n Done calculating ' +
                  'ratio chi['+repr(i)+']['+repr(j)+'] over chi['+repr(k)+']['+repr(l)+']. \n'
                                                                                               )
            if ( not bDoNotEvaluateDiffs ):
              print                                                                       (
                    '\n Calculating difference between ' +
                    'chi['+repr(i)+']['+repr(j)+'] and chi['+repr(k)+']['+repr(l)+']: \n'
                                                                                          )
              getSimpleTimingData( "Ratio evaluation start" )
              chi_2nd_differences[27*i+9*j+3*k+l] =                       (
                                                    simplify
                                                    (
                                                     x_2nd_result[i][j] -
                                                     x_2nd_result[k][l]
                                                    )
                                                                          )
              getSimpleTimingData( "Ratio evaluation finish" )
              print                                                                       (
                    '\n Done calculating difference between ' +
                    'chi['+repr(i)+']['+repr(j)+'] and chi['+repr(k)+']['+repr(l)+']. \n'
                                                                                          )
            pass# IF ( not bDoNotEvaluateDiffs )
          pass#FOR l in range(0,3)
        pass#FOR k in range(0,3)
      pass# FOR j in range(0,3)
    pass# FOR i in range(0,3)
    f = open( 'test4.sympy.tex', 'w' )
    f.write( '\\documentclass{article}\n\\usepackage{amsmath}\n' )
    f.write( '\\usepackage{mathtools}\n\\begin{document}\n' )
    x_2nd_all = 0;
    for i in range(0,3):
      for j in range(0,3):
        f.write( '\\begin{equation}\n\\begin{gathered}\n\\begin{multlined}\n' )
        f.write( latex( x_2nd_result[i][j], mode='inline' ) )
        f.write( '\n\\end{multlined}\n\\end{gathered}\n\\end{equation}\n' )
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
    f.write( '\\begin{equation}\n\\begin{gathered}\n\\begin{multlined}\n' )
    f.write( latex( E_pump, mode='inline' ) )
    f.write( '\n\\end{multlined}\n\\end{gathered}\n\\end{equation}\n' )
    #
    f.write( '\\begin{equation}\n\\begin{gathered}\n\\begin{multlined}\n' )
    f.write( latex( E_probe, mode='inline' ) )
    f.write( '\n\\end{multlined}\n\\end{gathered}\n\\end{equation}\n' )
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
    if (bCheckEquals):
      for i in range(0,3):
        for j in range(0,3):
          for k in range(0,3):
            for l in range(0,3):
              if ( (i != k) and (j != l) ):
                print                                                                       (
                      '\n (Search for equals) Calculating difference between ' +
                      'chi['+repr(i)+']['+repr(j)+'] and chi['+repr(k)+']['+repr(l)+']: \n'
                                                                                            )
                difference = simplify( Abs( x_2nd[i][j] ) - Abs( x_2nd[k][l] ) )
                print                                                                       (
                      '\n (Search for equals) Done calculating difference between ' +
                      'chi['+repr(i)+']['+repr(j)+'] and chi['+repr(k)+']['+repr(l)+']. \n'
                                                                                            )
                if ( difference == 0 ):
                  print( 'x_2nd[' + repr(i) + '][' + repr(j) + '] =' )
                  print( '= x_2nd[' + repr(k) + '][' + repr(l) + ']' )
                else:
                  print( '\n' )
                  print 'i=', i, 'j=', j, 'k=', k, 'l=', l
                  print( '\n' )
                  pprint( difference )
                  print( '\n' )
                  print( '\n' )
                pass
              pass
            pass
          pass
        pass
      pass
    pass
  #
  if ( not bSwitchExec2LibOnly ):
    beginSimpleTiming()
    do_main()
    getSimpleTimingData( 'Calculation timer' )
  pass# IF
  #
pass# DEF the_main_routine()
#
