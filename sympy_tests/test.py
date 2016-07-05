#!/usr/bin/env python
#
from sympy import *
import timeit
#
t, w, tau, infSymb, t_0, k, x, y, z, w_0 = symbols( 't, w, tau, infSymb, t_0, k, x, y, z, w_0' )
chi_1221, chi_1122, chi_1212 = symbols( 'chi_1221, chi_1122, chi_1212' )
r_p, a_p, t_p, A_1, R_1, A_2, R_2 = symbols( 'r_p, a_p, t_p, A_1, R_1, A_2, R_2' )
#
z_r = k*(w_0**2)/2
rho = sqrt( x**2 + y**2 )
r = sqrt( x**2 + y**2 + z**2 )
f = 1 / ( k * w_0 )
phi = atan( y /x )
#
L_ra = Function( 'L_ra' )
L = L_ra( r_p, a_p, t_p )
#
bQuasiSteadyState = True
bLaguerreSymbolic = False
bSubstituteValues = True
bCheckEquals = False
#
def det2x2(M1,M2,M3,M4):
  ##return simplify( sympify( M1*M4 - M2*M3 ) )
  return M1*M4 - M2*M3
def det3x3(M):
  ##return simplify( sympify( M[0][0]*det2x2(M[1][1],M[1][2],M[2][1],M[2][2])-M[0][1]*det2x2(M[1][0],M[1][2],M[2][0],M[2][2])+M[0][2]*det2x2(M[1][0],M[1][1],M[2][0],M[2][1]) ) )
  return M[0][0]*det2x2(M[1][1],M[1][2],M[2][1],M[2][2])-M[0][1]*det2x2(M[1][0],M[1][2],M[2][0],M[2][2])+M[0][2]*det2x2(M[1][0],M[1][1],M[2][0],M[2][1])
#
def E_TETM_comon( radial, azimuth ):
  if ( bLaguerreSymbolic ):
    return(
                -( ( I**( 2*radial - azimuth + 1 ) * z_r ) /
                ( r**2 * rho**2 ) ) *
                (rho / ( sqrt(2) * r * f ))**azimuth *
                L *
                exp( I*k*r - rho**2 / ( 4* f**2 * r**2 ) ) * cos( azimuth * phi )
                )
  else:
    return(
                -( ( I**( 2*radial - azimuth + 1 ) * z_r ) /
                ( r**2 * rho**2 ) ) *
                (rho / ( sqrt(2) * r * f ))**azimuth *
                assoc_laguerre( radial, azimuth,
                 rho**2 / ( 2 * f**2 * r**2 )
                ) *
                exp( I*k*r - rho**2 / ( 4* f**2 * r**2 ) ) * cos( azimuth * phi )
                )
#
def E_TM_comon_xyz( radial, azimuth ):
  return E_TETM_comon( radial, azimuth ) * y * z
#
def E_TE_common_xyz ( radial, azimuth ):
  return E_TETM_comon( radial, azimuth ) * x
#
def E_TE_comon( radial, azimuth ):
  E_common = [0 for i in xrange(3)]
  E_common[0] = E_TE_common_xyz( radial, azimuth ) * y
  E_common[1] = -E_TE_common_xyz( radial, azimuth ) * x
  E_common[2] = 0
  return E_common
#
def E_TM_comon( radial, azimuth ):
  E_common = [0 for i in xrange(3)]
  E_common[0] = E_TE_common_xyz( radial, azimuth ) * x * z
  E_common[1] = E_TE_common_xyz( radial, azimuth ) * y * z
  E_common[2] = -E_TE_common_xyz( radial, azimuth ) * rho**2
  return E_common
#
#
#
exp_QWP = exp( I*pi/2 )
QWP =                  (
      [
       [exp_QWP, 0],
       [0, -I*exp_QWP]
      ]
                       )
#
def jones_calc( M, V ):
  U = [0 for i in xrange(3)]
  for i in range(0,2):
    for j in range(0,2):
      U[i] = U[i] + V[j] * M[i][j]
  U[2] = V[2]
  return U
#
if ( bSubstituteValues ):
  E_pump = jones_calc( QWP, E_TM_comon( 0, 1 ) )
  E_probe = jones_calc( QWP, E_TM_comon( 0, 1 ) )
else:
  E_pump = jones_calc( QWP, E_TM_comon( R_1, A_1 ) )
  E_probe = jones_calc( QWP, E_TM_comon( R_2, A_2 ) )
#
# optical kerr quasi steady-state linit:
chi_1111 = chi_1221 + chi_1122 + chi_1212
def x_3rd_order (  ):
  x_3rd = [[[[0 for i in xrange(3)] for j in xrange(3)] for k in xrange(3)] for p in xrange(3)]
  for i in range(0,3):
    for j in range(0,3):
      for k in range(0,3):
        for p in range(0,3):
          if ( (i == j) and (j == k) and (k==p) ):
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
  return x_3rd
x_3rd = x_3rd_order()
def x_2nd_order_rel ():
  x_2nd_ord = [[0 for i in xrange(3)] for j in xrange(3)]
  for i in range(0,3):
    for j in range(0,3):
      for k in range(0,3):
        for p in range(0,3):
          x_2nd_ord[i][j] += x_3rd[i][j][k][p] * E_pump[k] * conjugate( E_pump[p] )
  return x_2nd_ord
#
def get_x_2nd():
  x_2nd = x_2nd_order_rel()
  for i in range(0,3):
    for j in range(0,3):
      if ( bQuasiSteadyState ):
        x_2nd[i][j] = simplify( powdenest( x_2nd[i][j], force=True ) )
      else:
        chi_integral = integrate                                                   (
                        x_2nd[i][j] * exp( -I * w * t_0 ) * exp( -(t - t_0)/tau ),
                        (t_0, -oo, t)
                                                                                   )
        x_2nd[i][j] = simplify( powdenest( chi_integral, force=True ) )
  return x_2nd
def do_main():
  x_2nd_result = get_x_2nd()
  f = open( 'test4.sympy.tex', 'w' )
  f.write( '\\documentclass{article}\n\\usepackage{amsmath}\n' )
  f.write( '\\usepackage{mathtools}\n\\begin{document}\n' )
  f.write( '\\begin{equation}\n\\begin{gathered}\n\\begin{multlined}\n' )
  for i in range(0,3):
    for j in range(0,3):
      f.write( latex( x_2nd_result[i][j], mode='inline' ) )
      f.write( '\n\\end{multlined}\n\\end{gathered}\n\\end{equation}\n' )
      f.write( '\\begin{equation}\n\\begin{gathered}\n\\begin{multlined}\n' )
  f.write( latex( E_pump, mode='inline' ) )
  f.write( '\n\\end{multlined}\n\\end{gathered}\n\\end{equation}\n' )
  f.write( '\\begin{equation}\n\\begin{gathered}\n\\begin{multlined}\n' )
  f.write( latex( E_probe, mode='inline' ) )
  f.write( '\n\\end{multlined}\n\\end{gathered}\n\\end{equation}\n\\end{document}' )
  f.close()
  pprint( x_3rd )
  if (bCheckEquals):
    for i in range(0,3):
      for j in range(0,3):
        for k in range(0,3):
          for l in range(0,3):
            if ( (i != k) and (j != l) ):
              difference = simplify( Abs( x_2nd[i][j] ) - Abs( x_2nd[k][l] ) )
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
do_main()
