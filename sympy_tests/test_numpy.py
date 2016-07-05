#!/usr/bin/env python
# Note that : (sic!) N2 is multiplied on lapl. from l2r. in the main formulae with 18 derivatives
#
#
##
###from sympy.solvers import pdsolve
###from sympy.abc import x, y, a, b, c, rho, phi, z
###from sympy import Function, pprint
###from sympy.mpmath import *
###a, b, c, d = symbols( 'a, b, c, d' )
###f_Vxyz = Function( 'Vxyz' )
###f = Function( 'f' )
###G = Function( 'G' )
###Vxyz = f_Vxyz( a, b, c )
###u = f( x, y )
###def try_to_change_list_contents(the_list):
###    print 'got', the_list
###    the_list[0] = 'minus_one'
###    print 'changed to', the_list
###outer_list = ['one', 'two', 'three']
###print 'before, outer_list =', outer_list
###try_to_change_list_contents(outer_list)
###print 'after, outer_list =', outer_list
###aaaaa = (z**2 * diff(Hrho,z) + diff(Hrho,phi) + sin(z) * diff(Hrho,z)) / diff(diff(Hrho,z),z)
###aaaaa = numer(aaaaa)
###aaaaa = aaaaa.subs(diff(Hrho,z),dHrho_dz)
###aaaaa = aaaaa.subs(diff(Hrho,phi),dHrho_dphi)
###gf_c1 = collect(aaaaa,dHrho_dz,evaluate=False)
###gf_c2 = collect(aaaaa,dHrho_dphi,evaluate=False)
###genform = a*u + b*ux + c*uy - G(x,y)
###genform = dVxyz_dc + diff(Vxyz*b**2,b)
###ux = u.diff(x)
###uy = u.diff(y)
###import timeit
###start = timeit.default_timer()
###stop = timeit.default_timer()
###print stop - start
##
#
#
from numpy import *
from sympy import *
import timeit
#
bgn_timestamp = timeit.default_timer()
end_timestamp = timeit.default_timer()
print 'test timer:'
print end_timestamp - bgn_timestamp
print '\n\n'
#
e_rho, e_phi, e_z, alpha, beta, epsilon, rho, phi, x, y, z, k, detM, dHrho_dz, dHrho_dphi, dHrho_drho, d2Hrho_dphidz, d2Hrho_drhodz, d2Hrho_drhodphi, d2Hrho_dz2, d2Hrho_dphi2, d2Hrho_drho2, C_r, n_0, n_e = symbols('e_rho, e_phi, e_z, alpha, beta, epsilon, rho, phi, x, y, z, k, detM, dHrho_dz, dHrho_dphi, dHrho_drho, d2Hrho_dphidz, d2Hrho_drhodz, d2Hrho_drhodphi, d2Hrho_dz2, d2Hrho_dphi2, d2Hrho_drho2, C_r, n_0, n_e')
#
# Tens. components:
f_M11 = Function('M11')
f_M12 = Function('M12')
f_M13 = Function('M13')
f_M21 = Function('M21')
f_M22 = Function('M22')
f_M23 = Function('M23')
f_M31 = Function('M31')
f_M32 = Function('M32')
f_M33 = Function('M33')
#
# Dummy tens. components:
f_MD11 = Function('MD11')
f_MD12 = Function('MD12')
f_MD13 = Function('MD13')
f_MD21 = Function('MD21')
f_MD22 = Function('MD22')
f_MD23 = Function('MD23')
f_MD31 = Function('MD31')
f_MD32 = Function('MD32')
f_MD33 = Function('MD33')
#
# Dummy tens. components:
f_NN11 = Function('NN11')
f_NN12 = Function('NN12')
f_NN13 = Function('NN13')
f_NN21 = Function('NN21')
f_NN22 = Function('NN22')
f_NN23 = Function('NN23')
f_NN31 = Function('NN31')
f_NN32 = Function('NN32')
f_NN33 = Function('NN33')
#
# Dummy tens. components:
f_N11 = Function('N11')
f_N12 = Function('N12')
f_N13 = Function('N13')
f_N21 = Function('N21')
f_N22 = Function('N22')
f_N23 = Function('N23')
f_N31 = Function('N31')
f_N32 = Function('N32')
f_N33 = Function('N33')
#
# Basic vect./tens. field components:
f_Hrho = Function('Hrho')
f_Hphi = Function('Hphi')
f_Hx = Function('Hx')
f_Hy = Function('Hy')
f_Hz = Function('H_z')# check this part # replace all with SCRIPT_subscript notations - for tex compatibility
#
Hrho = f_Hrho(rho,phi,z)
Hphi = f_Hphi(rho,phi,z)
Hz_ = f_Hz(rho,phi,z)
#
Hx = f_Hx(x,y,z)
Hy = f_Hy(x,y,z)
Hz = f_Hz(x,y,z)
#
#M11 = f_M11(x,y,z)
#M12 = f_M12(x,y,z)
#M13 = f_M13(x,y,z)
#M21 = f_M21(x,y,z)
#M22 = f_M22(x,y,z)
#M23 = f_M23(x,y,z)
#M31 = f_M31(x,y,z)
#M32 = f_M32(x,y,z)
#M33 = f_M33(x,y,z)
#
#M11 = f_M11(rho,phi,z)
#M12 = f_M12(rho,phi,z)
#M13 = f_M13(rho,phi,z)
#M21 = f_M21(rho,phi,z)
#M22 = f_M22(rho,phi,z)
#M23 = f_M23(rho,phi,z)
#M31 = f_M31(rho,phi,z)
#M32 = f_M32(rho,phi,z)
#M33 = f_M33(rho,phi,z)
#
n_d = n_e - n_0
#
M11 = n_0 + n_d * (cos( C_r * rho ))**2 * (cos( phi ))**2
M12 = -n_d * (cos( C_r * rho ))**2 * cos( phi ) * sin( phi )
M13 = -n_d * cos( C_r * rho ) * sin( C_r * rho ) * cos( phi )
M21 = -n_d * (cos( C_r * rho ))**2 * cos( phi ) * sin( phi )
M22 = n_0 + n_d * (cos( C_r * rho ))**2 * (sin( phi ))**2
M23 = n_d * cos( C_r * rho ) * sin( C_r * rho ) * sin( phi )
M31 = -n_d * cos( C_r * rho ) * sin( C_r * rho ) * cos( phi )
M32 = n_d * cos( C_r * rho ) * sin( C_r * rho ) * sin( phi )
M33 = n_0 + n_d * (sin( C_r * rho ))**2
#
# Dummies:
#
MD11 = f_MD11(rho,phi,z)
MD12 = f_MD12(rho,phi,z)
MD13 = f_MD13(rho,phi,z)
MD21 = f_MD21(rho,phi,z)
MD22 = f_MD22(rho,phi,z)
MD23 = f_MD23(rho,phi,z)
MD31 = f_MD31(rho,phi,z)
MD32 = f_MD32(rho,phi,z)
MD33 = f_MD33(rho,phi,z)
#
NN11 = f_NN11(rho,phi,z)
NN12 = f_NN12(rho,phi,z)
NN13 = f_NN13(rho,phi,z)
NN21 = f_NN21(rho,phi,z)
NN22 = f_NN22(rho,phi,z)
NN23 = f_NN23(rho,phi,z)
NN31 = f_NN31(rho,phi,z)
NN32 = f_NN32(rho,phi,z)
NN33 = f_NN33(rho,phi,z)
#
N11 = f_N11(rho,phi,z)
N12 = f_N12(rho,phi,z)
N13 = f_N13(rho,phi,z)
N21 = f_N21(rho,phi,z)
N22 = f_N22(rho,phi,z)
N23 = f_N23(rho,phi,z)
N31 = f_N31(rho,phi,z)
N32 = f_N32(rho,phi,z)
N33 = f_N33(rho,phi,z)
#
# Vect. field:
Hexp_part = exp( I * beta * z )
H = array([( Hexp_part * Hrho ), ( Hexp_part * Hphi ), ( Hexp_part * Hz_ )])
#H = [0 for i in xrange(3)]
#H[0] = exp( I * beta * z ) * Hrho
#H[1] = exp( I * beta * z ) * Hphi
#H[2] = exp( I * beta * z ) * Hz_
#
# Vect. field in cartes. coord. expressed in terms of cyl. coord.:
Hxyz = array([( Hexp_part * Hrho * cos( Hphi ) ), ( Hexp_part * Hrho * sin( Hphi ) ), ( Hexp_part * Hrho * sin( Hz_ ) )])
#Hxyz = [0 for i in xrange(3)]
#Hxyz[0] = exp( I * beta * z ) * Hrho * cos( Hphi )
#Hxyz[1] = exp( I * beta * z ) * Hrho * sin( Hphi )
#Hxyz[2] = exp( I * beta * z ) * Hz_
#
# Vect. field in cartes. coord.:
Hcartes = array([Hx, Hy, Hz])
#Hcartes = [0 for i in xrange(3)]
#Hcartes[0] = Hx
#Hcartes[1] = Hy
#Hcartes[2] = Hz
#
# Cyl. coord.:
v = x * sin( alpha * z ) + y * cos( alpha * z )
u = x * cos( alpha * z ) - y * sin( alpha * z )
#
# Basic 2nd. ord. tens. for inverse squared refr. index:
M = array([[M11, M12, M13],
           [M21, M22, M23],
           [M31, M32, M33]])
#M = [[0 for i in xrange(3)] for i in xrange(3)]
#M = [
#     [M11, M12, M13],
#     [M21, M22, M23],
#     [M31, M32, M33]
#    ]
M_dummy = array([[MD11, MD12, MD13],
                 [MD21, MD22, MD23],
                 [MD31, MD32, MD33]])
#M_dummy = [
#           [MD11, MD12, MD13],
#           [MD21, MD22, MD23],
#           [MD31, MD32, MD33]
#          ]
#
rot_left = array([[(cos(phi)), (sin(phi)), 0],
                  [(-sin(phi)), (cos(phi)), 0],
                  [0, 0, 1]])
#rot_left = [
#            [ cos( phi ), sin( phi ), 0],
#            [-sin( phi ), cos( phi ), 0],
#            [   0,          0,        1]
#           ]
rot_right = array([[(cos(phi)), -sin(phi), 0],
                   [sin(phi), cos(phi), 0],
                   [0, 0, 1]])
#rot_right = [
#             [cos( phi ), -sin( phi ), 0],
#             [sin( phi ),  cos( phi ), 0],
#             [  0,           0,        1]
#            ]
#
# Cyl. coord. covar. metric:
g_covar = array([1, 0, 0],
                [0, (rho**2), 0],
                [0, 0, 1])
#g_covar = [
#           [1,   0,    0],
#           [0, rho**2, 0],
#           [0,   0,    1]
#          ]
g_contravar = array([[1, 0, 0],
                     [0, (1 / rho**2), 0],
                     [0, 0, 1]])
#g_contravar = [
#               [1,   0,        0],
#               [0, 1 / rho**2, 0],
#               [0,   0,        1]
#              ]
#
# Basic tens. as metric for rotating frame:
#M[0][0] = epsilon + alpha * v**2
#M[0][1] = ( (-1) * alpha**2 ) * u * v
#M[0][2] = -alpha * v
#M[1][0] = ( (-1) * alpha**2 ) * u * v
#M[1][1] = epsilon + alpha * u**2
#M[1][2] = -alpha * u
#M[2][0] = -alpha * v
#M[2][1] = -alpha * u
#M[2][2] = epsilon
#
# More precisely for inverse matrix:
#M[0][0] = ( epsilon + alpha * u**2 * ( 1 - alpha ) ) / detM
#M[0][1] = 0
#M[0][2] = ( alpha * v * ( 1 + alpha * u**2 * ( 1 - alpha ) ) ) / detM
#M[1][0] = 0
#M[1][1] = ( epsilon + alpha * v**2 * ( 1 - alpha ) ) / detM
#M[1][2] = ( alpha * u * ( 1 + alpha * v**2 * ( 1 - alpha ) ) ) / detM
#M[2][0] = ( alpha * v * ( 1 + alpha * u**2 * ( 1 - alpha ) ) ) / detM
#M[2][1] = ( alpha * u * ( 1 + alpha * u**2 * ( 1 - alpha ) ) ) / detM
#M[2][2] = ( epsilon + alpha * ( u**2 * v**2 * alpha * ( 1 - alpha**2 ) + alpha * ( u**2 + v**2 ) ) ) / detM
#
# Additional 2nd ord. tens. for auxilary expression:
N1 = array([[M11, M21, M31],
            [M12, M22, M32],
            [M13, M23, M33]])
#N1 = [
#     [M11, M21, M31],
#     [M12, M22, M32],
#     [M13, M23, M33]
#    ]
N1_dummy = array([[N11, N12, N13],
                  [N21, N22, N23],
                  [N31, N32, N33]])
#N1_dummy = [
#            [N11, N12, N13],
#            [N21, N22, N23],
#            [N31, N32, N33]
#           ]
#
# Additional 2nd ord. tens. for lapl. expression:
N2 = array([[(M33 + M22), (-M21), (-M31)],
            [(-M12), (M33 + M11), (-M32)],
            [(-M13), (-M23), (M11 + M22)]])
#N2 = [[0 for i in xrange(3)] for i in xrange(3)]
N2_dummy = array([[NN11, NN12, NN13],
                  [NN21, NN22, NN23],
                  [NN31, NN32, NN33]])
#N2_dummy = [
#            [NN11, NN12, NN13],
#            [NN21, NN22, NN23],
#            [NN31, NN32, NN33]
#           ]
#
# and it's components in general form:
#N2[0][0] = M22 + M33 + M33 + M11 + M22 - M11 - M22 - M33
#N2[0][1] = M32 + M13 - M32 - M13 - M21
#N2[0][2] = M12 + M23 - M31 - M12 - M23
#N2[1][0] = M23 + M31 - M23 - M31 - M12
#N2[1][1] = M33 + M11 + M22 + M33 + M11 - M22 - M33 - M11
#N2[1][2] = M13 + M21 - M21 - M32 - M13
#N2[2][0] = M21 + M32 - M13 - M21 - M32
#N2[2][1] = M31 + M12 - M12 - M23 - M31
#N2[2][2] = M11 + M22 + M11 + M22 + M33 - M33 - M11 - M22
#
N2[0][0] = M33 + M22
N2[0][1] = -M21
N2[0][2] = -M31
N2[1][0] = -M12
N2[1][1] = M33 + M11
N2[1][2] = -M32
N2[2][0] = -M13
N2[2][1] = -M23
N2[2][2] = M11 + M22
#
# Covar. basis (identical to basis {rho, phi, z}) coord. function names (for differentiating metric tens. in Chris. coeff. calculations):
cov_b = array([rho, phi, z])
#cov_b = [rho, phi, z]
#
# Levi-Civita non-zero pairs relative to 3rd index:
LC_nzPairs = array([[1, 2, 2, 1],
                    [2, 0, 0, 2],
                    [0, 1, 1, 0]])
#LC_nzPairs = [# i  j  i  j
#              [ 1, 2, 2, 1 ],# k = 1
#              [ 2, 0, 0, 2 ],# k = 2
#              [ 0, 1, 1, 0 ]#  k = 3
#             ]# +  +  -  -
#
# Chris. symbols of 2nd kind (I've been too lasy enough to derive them in paper) for cyl. coord.:
def Chris_2nd_cyl( i_con, j_cov, k_cov ):
  ChGamma = 0
  for m in range( 0, 3 ):
    ChGamma += g_contravar[i_con, m] * ( diff( g_covar[m, j_cov], cov_b[k_cov] ) + diff( g_covar[m, k_cov], cov_b[j_cov] ) - diff( g_covar[j_cov, k_cov], cov_b[m] ) )
  ##return simplify( sympify( ChGamma / 2 ) )
  return ChGamma / 2
#
# Covariant derivatives (I've been too lasy enough to derive them in paper) used to calculate curvilinear cyl. curl of 2nd. ord. covar. tens.:
# def cov_diff_2nd_ord_tens_single_comp( M, i_cov, j_con, m_diff ):# M[i][j] - mixed: covariant i, contravariant j
  # Mdiff = diff( M[i_cov][j_con], cov_b[m_diff] )
  # for k in range(0,3):
    # Mdiff += Chris_2nd_cyl( j_con, k, m_diff ) * M[i_cov][k] - Chris_2nd_cyl( k, i_cov, m_diff ) * M[k][j_con]
  # return simplify( Mdiff )
#
def cov_diff_vect_single___tens( Stens, i_cov, j_con, m_diff ):
  Mdiff = diff( Stens[i_cov, j_con], cov_b[m_diff] )
  for k in range(0,3):
    Mdiff -= Chris_2nd_cyl( k, m_diff, i_cov ) * Stens[k, j_con]
  ##return simplify( sympify( Mdiff ) )
  return Mdiff
#
# Matrix multiplication:
def mult1_2nd_ord_tens( M1, M2 ):#RENAME
 MM = zeros((3,3), dtype=numpy.int)
 MM.fill(S.Zero)
 #MM = [[0 for i in xrange(3)] for i in xrange(3)]
 for i in range(0,3):
  for j in range(0,3):
   for k in range(0,3):
    MM[i, j] += M1[i, k] * M2[k, j]
   ##MM[i][j] = simplify( sympify( MM[i][j] ) )
 return MM
#
# Convert vect. from physical representation in normal covar. basis {e_rho, e_phi, e_z} to general representation in nonnormalized contravar. basis {b_rho, b_phi, b_z} for further calculations:
def conv_vect_eee2bbb_cyl(Veee):
  U = zeros((3,), dtype=numpy.int)
  U.fill(S.Zero)
  #U = [0 for i in xrange(3)]
  for i in range(0,3):
    ##U[i] = simplify( sympify( powdenest ( Veee[i] / sqrt( g_covar[i][i] ) ) ) )
    U[i] = powdenest( Veee[i] / sqrt( g_covar[i, i] ), force=True )
  return U
#
# Convert tens. from physical representation in normal covar. basis {e_rho, e_phi, e_z} to general representation in nonnormalized contravar. basis {b_rho, b_phi, b_z} for further calculations:
def conv_tens_eee2bbb_cyl(Meee):
  Mbbb = zeros((3,3), dtype=numpy.int)
  Mbbb.fill(S.Zero)
  #Mbbb = [[0 for i in xrange(3)] for i in xrange(3)]
  for i in range(0,3):
    for j in range(0,3):
      ##Mbbb[i][j] = simplify( sympify( powdenest( Meee[i][j] / sqrt( g_covar[i][i] * g_covar[j][j] ), force=True ) ) )
      Mbbb[i, j] = powdenest( Meee[i, j] / sqrt( g_covar[i, i] * g_covar[j, j] ), force=True )
  return Mbbb
#
# Convert vect. back to phys. representation in norm. covar. basis {e_rho, e_phi, e_z} from general representation in nonnorm. contravar. basis {b_rho, b_phi, b_z} for further calculations:
def conv_vect_bbb2eee_cyl( Vbbb ):
  U = zeros((3,), dtype=numpy.int)
  U.fill(S.Zero)
  #U = [0 for i in xrange(3)]
  for i in range(0,3):
    ##U[i] = simplify( sympify( powdenest ( Veee[i] * sqrt( g_covar[i][i] ) ) ) )
    U[i] = powdenest( Vbbb[i] * sqrt( g_covar[i, i] ), force=True )
  return U
#
# Convert tens. back to phys. representation in norm. covar. basis {e_rho, e_phi, e_z} from general representation in nonnorm. contravar. basis {b_rho, b_phi, b_z} for further calculations:
def conv_tens_bbb2eee_cyl( Mbbb ):
  Meee = zeros((3,3), dtype=numpy.int)
  Meee.fill(S.Zero)
  #Meee = [[0 for i in xrange(3)] for i in xrange(3)]
  for i in range(0,3):
    for j in range(0,3):
      ##Meee[i][j] = simplify( sympify( powdenest( Mbbb[i][j] * sqrt( g_covar[i][i] * g_covar[j][j] ), force=True ) ) )
      Meee[i, j] = powdenest( Mbbb[i, j] * sqrt( g_covar[i, i] * g_covar[j, j] ), force=True )
  return Meee
#
def convert_2nd_order_tens( Mxyz ):
  Mrpz = zeros((3,3), dtype=numpy.int)
  Mrpz.fill(S.Zero)
  #Mrpz = [[0 for i in xrange(3)] for i in xrange(3)]
  Mrpz = mult1_2nd_ord_tens( rot_left, Mxyz )
  Mrpz = mult1_2nd_ord_tens( Mrpz, rot_right )
  return Mrpz
def prod_scal2vect( s, V ):
  U = zeros((3,), dtype=numpy.int)
  Mrpz.fill(S.Zero)
  #U = [0 for i in xrange(3)]
  ##U[0] = simplify( sympify( s * V[0] ) )
  U[0] = s * V[0]
  ##U[1] = simplify( sympify( s * V[1] ) )
  U[1] = s * V[1]
  ##U[2] = simplify( sympify( s * V[2] ) )
  U[2] = s * V[2]
  return U
# Part. deriv. by x in cyl. coord.:
def diffx_rho_phi( A ):
  ##return simplify( sympify( diff(A,rho)*cos(phi)-1/rho*diff(A,phi)*sin(phi) ) )
  return diff( A, rho ) * cos( phi ) - (1 / rho) * diff( A, phi ) * sin( phi )
# Part. deriv. by y in cyl. coord.:
def diffy_rho_phi( A ):
  ##return simplify( sympify( diff(A,rho)*sin(phi)+1/rho*diff(A,phi)*cos(phi) ) )
  return diff(A,rho)*sin(phi)+1/rho*diff(A,phi)*cos(phi)
# Div. in cyl. coord.:
def div_cyl( V ):
  ##return simplify( sympify( (1/rho)*diff(rho*V[0],rho)+(1/rho)*diff(V[1],phi)+diff(V[2],z) ) )
  return (1/rho)*diff(rho*V[0],rho)+(1/rho)*diff(V[1],phi)+diff(V[2],z)
# Div. in cyl. coord. expressed in cartes. coord.:
def div_cartes_2_cyl( V ):
  ##return simplify( sympify( diffx_rho_phi(V[0])+diffy_rho_phi(V[1])+diff(V[2],z) ) )
  return diffx_rho_phi(V[0])+diffy_rho_phi(V[1])+diff(V[2],z)
#
# Div. of norm. ( eee ) 2nd. ord. tens. in cyl. coord.:
def div_tens_cyl( M ):
  Udiv = zeros((3,), dtype=numpy.int)
  Udiv.fill(S.Zero)
  #Udiv = [0 for i in xrange(3)]
  ##Udiv[0] = simplify( sympify( diff( M[0][0], rho ) + (M[0][0] / rho) + ((1 / rho) * diff( M[1][0], phi )) + diff( M[2][0], z ) + diff( M[1][1], rho ) ) )
  Udiv[0] = diff( M[0, 0], rho ) + (M[0, 0] / rho) + ((1 / rho) * diff( M[1, 0], phi )) + diff( M[2, 0], z ) + diff( M[1, 1], rho )
  ##Udiv[1] = simplify( sympify( ((1 / rho) * diff( M[1][1], phi )) + diff( M[0][1] , rho ) + (M[0][1] / rho) + (M[1][0] / rho) + diff( M[2][1], z ) ) )
  Udiv[1] = ((1 / rho) * diff( M[1, 1], phi )) + diff( M[0, 1] , rho ) + (M[0, 1] / rho) + (M[1, 0] / rho) + diff( M[2, 1], z )
  ##Udiv[2] = simplify( sympify( diff( M[2][2], z ) + diff( M[0][2], rho ) + (M[0][2] / rho) + ((1 / rho) * diff( M[1][2], phi )) ) )
  Udiv[2] = diff( M[2, 2], z ) + diff( M[0, 2], rho ) + (M[0, 2] / rho) + ((1 / rho) * diff( M[1, 2], phi ))
  return Udiv
#
# Cross prod. in cartes. coord.:
def cross_prod_cartes( U, V ):
  T = zeros((3,), dtype=numpy.int)
  T.fill(S.Zero)
  #T = [0 for i in xrange(3)]
  ##T[0] = simplify( sympify( U[1] * V[2] - U[2] * V[1] ) )
  T[0] = U[1] * V[2] - U[2] * V[1]
  ##T[1] = simplify( sympify( U[2] * V[0] - U[0] * V[2] ) )
  T[1] = U[2] * V[0] - U[0] * V[2]
  ##T[2] = simplify( sympify( U[0] * V[1] - U[1] * V[0] ) )
  T[2] = U[0] * V[1] - U[1] * V[0]
  return T
# Cross prod. in cyl. coord.:
X1 = rho*cos(phi)
X2 = rho*sin(phi)
X3 = z
# Jacob_cyl_cartes = [[0 for i in xrange(3)] for i in xrange(3)]
Jacob_cyl_cartes = [
                    [(diff(X1,rho)), (diff(X1,phi)), (diff(X1,z))]
                    [(diff(X2,rho)), (diff(X2,phi)), (diff(X2,z))]
                    [(diff(X3,rho)), (diff(X3,phi)), (diff(X3,z))]
                   ]
#Jacob_cyl_cartes = [
#                    [diff(X1,rho),diff(X1,phi),diff(X1,z)],
#                    [diff(X2,rho),diff(X2,phi),diff(X2,z)],
#                    [diff(X3,rho),diff(X3,phi),diff(X3,z)]
#                   ]
def det2x2(M1,M2,M3,M4):
  ##return simplify( sympify( M1*M4 - M2*M3 ) )
  return M1*M4 - M2*M3
def det3x3(M):
  ##return simplify( sympify( M[0][0]*det2x2(M[1][1],M[1][2],M[2][1],M[2][2])-M[0][1]*det2x2(M[1][0],M[1][2],M[2][0],M[2][2])+M[0][2]*det2x2(M[1][0],M[1][1],M[2][0],M[2][1]) ) )
  return M[0, 0]*det2x2(M[1, 1],M[1, 2],M[2, 1],M[2, 2])-M[0, 1]*det2x2(M[1, 0],M[1, 2],M[2, 0],M[2, 2])+M[0, 2]*det2x2(M[1, 0],M[1, 1],M[2, 0],M[2, 1])
Jacob_det = det3x3(Jacob_cyl_cartes)
# Cross prod. of two covar. vect. (i.g. their components are covar. and have lower indices, but they are given in contrav. nonnormalized(sic!) basis {b_rho; b_phi; b_z} with roof indices rho, phi, z):
def cross_prod_cyl( U, V ):
  T = zeros((3,), dtype=numpy.int)
  T.fill(S.Zero)
  #T = [0 for i in xrange(3)] # (i,j,k)
  ##T[0] = simplify( sympify( (1 / Jacob_det) * ( U[1] * V[2] - U[2] * V[1] ) ) )# (2,3,1) - (3,2,1) # rho
  T[0] = (1 / Jacob_det) * ( U[1] * V[2] - U[2] * V[1] )# (2,3,1) - (3,2,1) # rho
  ##T[1] = simplify( sympify( (1 / Jacob_det) * ( U[2] * V[0] - U[0] * V[2] ) ) )# (3,1,2) - (1,3,2) # phi
  T[1] = (1 / Jacob_det) * ( U[2] * V[0] - U[0] * V[2] )# (3,1,2) - (1,3,2) # phi
  ##T[2] = simplify( sympify( (1 / Jacob_det) * ( U[0] * V[1] - U[1] * V[0] ) ) )# (1,2,3) - (2,1,3) # z
  T[2] = (1 / Jacob_det) * ( U[0] * V[1] - U[1] * V[0] )# (1,2,3) - (2,1,3) # z
  return T# Contravar. vect. with components having roof indices, but given in covar. basis {b_rho; b_phi; b_z} with lower indices rho, phi, z
#
# Cross prod. spec. case - dot prod. of nabla and 2nd rank tens. on the left side:
def cross_prod_spec_case( M, V ):
  U = zeros((3,), dtype=numpy.int)
  U.fill(S.Zero)
  #U = [0 for i in xrange(3)]
  ##U[0] = simplify( sympify( (M[1][0]*diff(V[2],x)+M[1][1]*diff(V[2],y)+M[1][2]*diff(V[2],z)-M[2][0]*diff(V[1],x)+M[2][1]*diff(V[1],y)+M[2][2]*diff(V[1],z)) ) )
  U[0] = (M[1, 0]*diff(V[2],x)+M[1, 1]*diff(V[2],y)+M[1, 2]*diff(V[2],z)-M[2, 0]*diff(V[1],x)+M[2, 1]*diff(V[1],y)+M[2, 2]*diff(V[1],z))
  ##U[1] = simplify( sympify( (M[2][0]*diff(V[0],x)+M[2][1]*diff(V[0],y)+M[2][2]*diff(V[0],z)-M[0][0]*diff(V[2],x)+M[0][1]*diff(V[2],y)+M[0][2]*diff(V[2],z)) ) )
  U[1] = (M[2, 0]*diff(V[0],x)+M[2, 1]*diff(V[0],y)+M[2, 2]*diff(V[0],z)-M[0, 0]*diff(V[2],x)+M[0, 1]*diff(V[2],y)+M[0, 2]*diff(V[2],z))
  ##U[2] = simplify( sympify( (M[0][0]*diff(V[1],x)+M[0][1]*diff(V[1],y)+M[0][2]*diff(V[1],z)-M[1][0]*diff(V[0],x)+M[1][1]*diff(V[0],y)+M[1][2]*diff(V[0],z)) ) )
  U[2] = (M[0, 0]*diff(V[1],x)+M[0, 1]*diff(V[1],y)+M[0, 2]*diff(V[1],z)-M[1, 0]*diff(V[0],x)+M[1, 1]*diff(V[0],y)+M[1, 2]*diff(V[0],z))
  return U
#
# Cross prod. from l2r. of mixed basis `b_cov (X) b_con' tens. and vect. in cyl. coord.:
def cross_prod_tens_l2r_vect_cyl( M, V ):# M have roof 1st. ind. and floor 2nd. ind.
  T = [0 for i in xrange(3)]
### T_cov = index_lower_both_2ndord_tens_cyl(M)
### T_mix = index_rise_2ndord_tens_cyl(T_cov)
### V_cov = index_lower_vect_cyl(V)
  T[0] = cross_prod_cyl( M[0], V )
  T[1] = cross_prod_cyl( M[1], V )
  T[2] = cross_prod_cyl( M[2], V )
  return T
#
# Curl in cartes. cord.:
def curl_simple(X,Y,Z):
  U = [0 for i in xrange(3)]
  ##U[0] = simplify( sympify( diff(Z,y)-diff(Y,z) ) )
  U[0] = diff(Z,y)-diff(Y,z)
  ##U[1] = simplify( sympify( diff(X,z)-diff(Z,x) ) )
  U[1] = diff(X,z)-diff(Z,x)
  ##U[2] = simplify( sympify( diff(Y,x)-diff(X,y) ) )
  U[2] = diff(Y,x)-diff(X,y)
  return U
# Curl applied to 2nd rank tensor from left side in cartes. coord.:
def curl_parallel_l2r_tensor_2nd(V):
  U = [[0 for i in xrange(3)] for i in xrange(3)]
  T = curl_simple(V[0][0],V[1][0],V[2][0])
  U[0][0] = T[0]
  U[1][0] = T[1]
  U[2][0] = T[2]
  T = curl_simple(V[0][1],V[1][1],V[2][1])
  U[0][1] = T[0]
  U[1][1] = T[1]
  U[2][1] = T[2]
  T = curl_simple(V[0][2],V[1][2],V[2][2])
  U[0][2] = T[0]
  U[1][2] = T[1]
  U[2][2] = T[2]
  return U
# Curl spec. case - cartes. coord. derivatives expressed with cyl. counterparts, except z:
def curl_parallel_2_cyl(V):
  U = [0 for i in xrange(3)]
  ##U[0] = simplify( sympify( diffy_rho_phi(V[2])-diff(V[1],z) ) )
  U[0] = diffy_rho_phi(V[2])-diff(V[1],z)
  ##U[1] = simplify( sympify( diff(V[0],z)-diffx_rho_phi(V[2]) ) )
  U[1] = diff(V[0],z)-diffx_rho_phi(V[2])
  ##U[2] = simplify( sympify( diffx_rho_phi(V[1])-diffy_rho_phi(V[0]) ) )
  U[2] = diffx_rho_phi(V[1])-diffy_rho_phi(V[0])
  return U
# same, but only rho component:
# def curl_parallel_2_cyl_rho(V):
#   return diffy_rho_phi(V[2])-diff(V[1],z)
# Curl in cyl. coordinates:
def curl(V):
  U = [0 for i in xrange(3)]
# U = [
  ##U[0] = simplify( sympify( diff(V[2],phi)/rho-diff(V[1],z) ) )
  U[0] = diff(V[2],phi)/rho-diff(V[1],z)
  ##U[1] = simplify( sympify( diff(V[0],z)-diff(V[2],rho) ) )
  U[1] = diff(V[0],z)-diff(V[2],rho)
  ##U[2] = simplify( sympify( (diff(rho*V[1],rho)-diff(V[0],phi))/rho ) )
  U[2] = (diff(rho*V[1],rho)-diff(V[0],phi))/rho
#     ]
  return U
#
# Curvil. cyl. curl of 2nd. ord. covar. tens.:
# def curl_cyl_2ndord_tens(M):
  # Mcurl = [[0 for i in xrange(3)] for i in xrange(3)]
  # for j in range(0,3):# columns
    # for i in range(0,3):# rows
      # Mcurl[i][j] = (1 / rho) * ( cov_diff_2nd_ord_tens_single_comp( M, LC_nzPairs[i][1], j, LC_nzPairs[i][0] ) - cov_diff_2nd_ord_tens_single_comp( M, LC_nzPairs[i][3], j, LC_nzPairs[i][2] ) )
  # return simplify( Mcurl )
#
def curl_cyl_2ndord_tens___2( S ):
  Mcurl = [[0 for i in xrange(3)] for i in xrange(3)]
  for j in range(0,3):# columns
    for i in range(0,3):# rows
      ##Mcurl[i][j] = simplify( sympify( (1 / rho) * ( cov_diff_vect_single___tens( S, LC_nzPairs[i][1], j, LC_nzPairs[i][0] ) - cov_diff_vect_single___tens( S, LC_nzPairs[i][3], j, LC_nzPairs[i][2] ) ) ) )
      Mcurl[i][j] = (1 / rho) * ( cov_diff_vect_single___tens( S, LC_nzPairs[i][1], j, LC_nzPairs[i][0] ) - cov_diff_vect_single___tens( S, LC_nzPairs[i][3], j, LC_nzPairs[i][2] ) )
  return Mcurl
#
# Dot prod.:
def dot_product(U,V):
  ##return simplify( sympify( U[0] * V[0] + U[1] * V[1] + U[2] * V[2] ) )
  return U[0] * V[0] + U[1] * V[1] + U[2] * V[2]
#
# def dot_prod_cyl( U, V ):
  # DP = 0
  # for i in range(0,3):
    # for j in range(0,3):
      # DP = DP + g_covar[i][j] * U[i] * V[j]
  # return DP
#
# Dot prod. from l2r of 2nd. ord. tens. with vect.:
def dot_prod_l2r_tens_vect( M, V ):
  U = [0 for i in xrange(3)]
  U[0] = dot_product( M[0], V )
  U[1] = dot_product( M[1], V )
  U[2] = dot_product( M[2], V )
  return U
#
# def dot_prod_l2r_tens_vect_cyl( M, V ):
  # U = [0 for i in xrange(3)]
  # U[0] = dot_prod_cyl( M[0], V )
  # U[1] = dot_prod_cyl( M[1], V )
  # U[2] = dot_prod_cyl( M[2], V )
  # return U
#
# Index lowering and rising used to evaluste cross prod. of 2nd. ord. covar. tens. l2r with vect. in cyl. coord.:
def index_lower_vect_cyl( V ):
  return dot_prod_l2r_tens_vect( g_covar, V )
#
# def index_rise_2ndord_tens_cyl( M ):
  # return mult1_2nd_ord_tens( g_contravar, M )
#
def index_lower_1st_2ndord_tens_cyl( M ):
  return mult1_2nd_ord_tens( g_covar, M )
#
def index_rise_1st_2ndord_tens_cyl( M ):
  return mult1_2nd_ord_tens( g_contravar, M )
#
def index_lower_both_2ndord_tens_cyl( M ):
 T = [[0 for i in xrange(3)] for i in xrange(3)]
 for i in range(0,3):
  for j in range(0,3):
   for k in range(0,3):
    for l in range(0,3):
     T[i][j] += g_covar[i][k] * g_covar[l][j] * M[k][l]
   ##T[i][j] = simplify( sympify( T[i][j] ) )
 return T
#
# Cartes. coord. to cyl. coord. STRICT transformation:
def cartes_2_cyl( V ):
  U = [0 for i in xrange(3)]
# U = [
  ##U[0] = simplify( sympify( sqrt( (V[0])**2 + (V[1])**2 ) ) )
  U[0] = sqrt( (V[0])**2 + (V[1])**2 )
  ##U[1] = simplify( sympify( atan(V[1]/V[0]) ) )
  U[1] = atan(V[1]/V[0])
  ##U[2] = simplify( sympify( V[2] ) )
  U[2] = V[2]
#     ]
  return U
def cyl_2_cartes(V):
  U = [0 for i in xrange(3)]
  ##U[0] = simplify( sympify( V[0]*cos(V[1]) ) )
  U[0] = V[0]*cos(V[1])
  ##U[1] = simplify( sympify( V[1]*sin(V[1]) ) )
  U[1] = V[1]*sin(V[1])
  ##U[2] = simplify( sympify( V[2] ) )
  U[2] = V[2]
  return U
# Vect. laplacian in cartes. coord:
def lapl(V):
  U = [0 for i in xrange(3)]
  for i in range(0,3):
    ##U[i] = simplify( sympify( diff(diff(V[i],x),x)+diff(diff(V[i],y),y)+diff(diff(V[i],z),z) ) )
    U[i] = diff(diff(V[i],x),x)+diff(diff(V[i],y),y)+diff(diff(V[i],z),z)
  return U
def lapl_cyl(V):
  U = [0 for i in xrange(3)]
  ##U[0] = simplify( sympify( diff(diff(V[0],rho),rho)+(1/rho**2)*diff(diff(V[0],phi),phi)+diff(diff(V[0],z),z)+(1/rho)*diff(V[0],rho)-(2/rho**2)*diff(V[1],phi)-V[0]/rho**2 ) )
  U[0] = diff(diff(V[0],rho),rho)+(1/rho**2)*diff(diff(V[0],phi),phi)+diff(diff(V[0],z),z)+(1/rho)*diff(V[0],rho)-(2/rho**2)*diff(V[1],phi)-V[0]/rho**2
  ##U[1] = simplify( sympify( diff(diff(V[1],rho),rho)+(1/rho**2)*diff(diff(V[1],phi),phi)+diff(diff(V[1],z),z)+(1/rho)*diff(V[1],rho)+(2/rho**2)*diff(V[0],phi)-V[1]/rho**2 ) )
  U[1] = diff(diff(V[1],rho),rho)+(1/rho**2)*diff(diff(V[1],phi),phi)+diff(diff(V[1],z),z)+(1/rho)*diff(V[1],rho)+(2/rho**2)*diff(V[0],phi)-V[1]/rho**2
  ##U[2] = simplify( sympify( diff(diff(V[2],rho),rho)+(1/rho**2)*diff(diff(V[2],phi),phi)+diff(diff(V[2],z),z)+(1/rho)*diff(V[2],rho) ) )
  U[2] = diff(diff(V[2],rho),rho)+(1/rho**2)*diff(diff(V[2],phi),phi)+diff(diff(V[2],z),z)+(1/rho)*diff(V[2],rho)
  return U
# Vect lapl. cyl. to cartes. mapping:
def lapl_cartes2cyl(V):
  U = [0 for i in xrange(3)]
  for i in range(0,3):
    ##U[i] = simplify( sympify( diffx_rho_phi(diffx_rho_phi(V[i]))+diffy_rho_phi(diffy_rho_phi(V[i]))+diff(diff(V[i],z),z) ) )
    U[i] = diffx_rho_phi(diffx_rho_phi(V[i]))+diffy_rho_phi(diffy_rho_phi(V[i]))+diff(diff(V[i],z),z)
  return U
# Matrix transponse:
def matr_trans(M):
  MT = [[0 for i in xrange(3)] for i in xrange(3)]
  for i in range(0,3):
    for j in range(0,3):
      MT[i][j] = M[j][i]
  return MT
#
##
#def mult1(a, b):
#  return a*b
##
#
#dmult1_drho = mult1.diff(rho)
#
##
#dVxyz_dc = Vxyz.diff(c)
##
#
#
#
#
#
#
end_timestamp = timeit.default_timer()
print 'precalculation timer:'
print end_timestamp - bgn_timestamp
print '\n\n'
bgn_timestamp = timeit.default_timer()
test = False
if test:
  Mtest = conv_tens_bbb2eee_cyl( curl_cyl_2ndord_tens___2( index_lower_1st_2ndord_tens_cyl( conv_tens_eee2bbb_cyl( M ) ) ) )
  # Mtest = curl_cyl_2ndord_tens___2( index_lower_1st_2ndord_tens_cyl( conv_tens_eee2bbb_cyl( M ) ) )
  # Mtest = index_lower_1st_2ndord_tens_cyl( conv_tens_eee2bbb_cyl( M ) )
  # Mtest = conv_tens_eee2bbb_cyl( M )
  # pprint(simplify( cov_diff_2nd_ord_tens_single_comp( M, 1, 0, 0 ) - cov_diff_2nd_ord_tens_single_comp( M, 0, 0, 1 )), wrap_line=True )
  pprint( Mtest, wrap_line=True )
  # print( latex( Mtest, mode='equation' ) )
else:
  # Main calculations.
  #
  # 1) k**2 * H[0..2] = curl( M[0..2][0..2] * (curl(H[0..2])) )
  #
  # 1.1) Right part
  #
  # 1.1.a) Cartesian coordinates
  #
  # M[0..2][0..2] * (curl(H[0..2])):
  ##M_curlH__ = dot_prod_l2r_tens_vect( M, curl_simple( Hcartes[0], Hcartes[1], Hcartes[2] ) )
  #
  # RHS:
  ##curl_rhs__ = curl_simple( M_curlH__[0], M_curlH__[1], M_curlH__[2] )
  #
  #
  ##### 1.1.b) Cyl. coord.
  #####
  ##### Vect. field curl:
  ####curlH = curl( H )
  ##### Convert 2nd. ord. tens. from cartes. to cyl. row-wise:
  ####Mcyl = [0 for i in xrange(3)]
  ####Mcyl[0] = cartes_2_cyl( M[0] )
  ####Mcyl[1] = cartes_2_cyl( M[1] )
  ####Mcyl[2] = cartes_2_cyl( M[2] )
  ##### Dot prod. of 2nd. ord. tens. and vect. field curl:
  ####M_curlH = dot_prod_l2r_tens_vect( Mcyl, curlH )
  #####
  ##### Convert above expression from cartes. to cyl. coord.:
  ####f_MH = cartes_2_cyl( M_curlH )
  #####
  ##### RHS:
  ####curl_rhs = curl( f_MH )
  #
  #
  # 1.1.c) Cart. to cyl. coord.
  #
  # Dot prod. of 2nd order tens. and a curl of vect. field, but with spec. case of a curl:
  ##MH_ = dot_prod_l2r_tens_vect(M,curl_parallel_2_cyl(Hxyz))
  #
  # RHS:
  ##curl_rhs_ = curl_parallel_2_cyl(MH_)
  #
  #
  # 1.1.d) Cyl. coord. tens. transform.:
  #
  # curl(H):
  _curlH = curl( H )
  #print "9\n\n"
  #pprint( _curlH, wrap_line=True )
  #
  # Convert M[0..2][0..2] from cartes. to cyl. coord.:
  _Mcyl0 = convert_2nd_order_tens( M )
  _Mcyl = M_dummy
  #
  # Convert M[0..2][0..2] from norm. covar. basis {e_rho; e_phi; e_z} to nonnorm. covar. basis {b_rho; b_phi; b_z}:
  _Mbbb = conv_tens_eee2bbb_cyl( _Mcyl )
  #
  # Convert of curl(H[0..2]) from norm. covar. basis {e_rho; e_phi; e_z} to nonnorm. covar. basis {b_rho; b_phi; b_z}:
  _curlH_bbb = conv_vect_eee2bbb_cyl( _curlH )
  #print "8\n\n"
  #pprint( _curlH_bbb, wrap_line=True )
  #
  # Convert curl(H[0..2]) from covar. to contravar. basis:
  _curlH_cov = index_lower_vect_cyl( _curlH_bbb )
  #print "7\n\n"
  #pprint( _curlH_cov, wrap_line=True )
  #
  # M[0..2][0..2] * curl(H[0..2]):
  _M_curlH = dot_prod_l2r_tens_vect( _Mbbb, _curlH )
  #
  # Convert `M[0..2][0..2] * curl(H[0..2])' back to norm. covar. basis {e_rho; e_phi; e_z} from nonnorm. covar. basis {b_rho; b_phi; b_z}:
  _McurlH_eee = conv_vect_bbb2eee_cyl( _M_curlH )
  #
  # RHS `curl( M[0..2][0..2] * curl(H[0..2]) )':
  _curl_rhs = curl( _McurlH_eee )
  #
  #
  #
  # 1.2) Left part
  #
  # 1.2.a) Cartesian coordinates
  #
  # LHS:
  ##curl_lhs__ = prod_scal2vect( k**2, Hcartes )
  #
  #
  ##### 1.2.b) Cyl. coord.
  #####
  ##### LHS:
  ####curl_lhs = prod_scal2vect( k**2, H )
  #####
  #####
  # 1.2.c) Cart. to cyl. coord.
  #
  # LHS:
  ##curl_lhs_ = prod_scal2vect(k**2,Hxyz)
  #
  #
  # 1.2.d) Cyl. coord. tens. transform.:
  #
  # LHS:
  _curl_lhs = prod_scal2vect( k**2, H )
  #
  #
  #
  #
  # 2) ( N2[0..2][0..2] * nabla**2 + k**2 ) * H[0..2] = (curl(M[0..2][0..2])) * (curl(H[0..2])) + ( N1[0..2][0..2] * nabla ) [X] (curl(H[0..2]))
  #
  # 2.1) Right part
  #
  # 2.1.a) Cartesian coordinates
  #
  # (curl(M[0..2][0..2])) * (curl(H[0..2])):
  ##curlM_curlH__ = dot_prod_l2r_tens_vect( curl_parallel_l2r_tensor_2nd( M ), curl_simple( Hcartes[0], Hcartes[1], Hcartes[2] ) )
  #
  # ( N1[0..2][0..2] * nabla ) [X] (curl(H[0..2])):
  ##N1nabla_X_curlH__ = cross_prod_spec_case( N1, curl_simple( Hcartes[0], Hcartes[1], Hcartes[2] ) )
  #
  # RHS:
  ##lapl_rhs__ = [a+b for a,b in zip(N1nabla_X_curlH__,curlM_curlH__)]
  #
  #
  ##### 2.1.b) Cylindrical coordinates
  #####
  ##### Cartes. mapping of curl(M[0..2][0..2]) in cyl. coord.:
  ####curlMV1 = cyl_2_cartes(curl(cartes_2_cyl([row[0] for row in M])))
  ####curlMV2 = cyl_2_cartes(curl(cartes_2_cyl([row[1] for row in M])))
  ####curlMV3 = cyl_2_cartes(curl(cartes_2_cyl([row[2] for row in M])))
  ####MV = [[0 for i in xrange(3)] for i in xrange(3)]
  ####for i in range(0,3):
  ######MV[i][0] = curlMV1[i]
  ######MV[i][1] = curlMV2[i]
  ######MV[i][2] = curlMV3[i]
  #####
  ##### (curl(M[0..2][0..2])) * (curl(H[0..2])):
  ####curlM_curlH = cartes_2_cyl(dot_prod_l2r_tens_vect(MV,cyl_2_cartes(curl(H))))
  #####
  ##### Cartesian mapping of nabla * N1[0..2][0..2] in polar coordinates:
  ####divN1 = [0 for i in xrange(3)]
  ####divN1[0] = div_cyl(cartes_2_cyl([row[0] for row in N1]))
  ####divN1[1] = div_cyl(cartes_2_cyl([row[1] for row in N1]))
  ####divN1[2] = div_cyl(cartes_2_cyl([row[2] for row in N1]))
  #####
  ##### (nabla * N1[0..2][0..2]) [X] (curl(H[0..2])):
  ####divN1_x_curlH = cross_prod_cyl(cartes_2_cyl(divN1),curl(H))
  #####
  ##### N1[0..2][0..2] in cyl. coord. row-wise:
  ####N1cyl = [0 for i in xrange(3)]
  ####N1cyl[0] = cartes_2_cyl(N1[0])
  ####N1cyl[1] = cartes_2_cyl(N1[1])
  ####N1cyl[2] = cartes_2_cyl(N1[2])
  #####
  ##### Cartes. mapping of `N1[0..2][0..2] [X] curl(H[0..2])' in cyl. coord. row-wise:
  ####mapN1curlH = [0 for i in xrange(3)]
  ####mapN1curlH[0] = cyl_2_cartes(cross_prod_cyl(N1cyl[0],curl(H)))
  ####mapN1curlH[1] = cyl_2_cartes(cross_prod_cyl(N1cyl[1],curl(H)))
  ####mapN1curlH[2] = cyl_2_cartes(cross_prod_cyl(N1cyl[2],curl(H)))
  #####
  ##### nabla * (N1[0..2][0..2] [X] curl(H[0..2])):
  ####div_N1curlH = [0 for i in xrange(3)]
  ####div_N1curlH[0] = div_cyl(cartes_2_cyl([row[0] for row in mapN1curlH]))
  ####div_N1curlH[1] = div_cyl(cartes_2_cyl([row[1] for row in mapN1curlH]))
  ####div_N1curlH[2] = div_cyl(cartes_2_cyl([row[2] for row in mapN1curlH]))
  #####
  ##### RHS:
  ####lapl_rhs = [a+b-c for a,b,c in zip(curlM_curlH,divN1_x_curlH,div_N1curlH)]
  #####
  #####
  # 2.1.c) Cart. to cyl. coord.
  #
  # curl(M[0..2][0..2]):
  ##curlM_ = [0 for i in xrange(3)]
  ##curlM_[0] = curl_parallel_2_cyl([row[0] for row in M])
  ##curlM_[1] = curl_parallel_2_cyl([row[1] for row in M])
  ##curlM_[2] = curl_parallel_2_cyl([row[2] for row in M])
  #
  # (curl(M[0..2][0..2])) * (curl(H[0..2])):
  ##curlM_curlH_ = dot_prod_l2r_tens_vect(curlM_,curl_parallel_2_cyl(Hxyz))
  #
  # nabla * N1[0..2][0..2]:
  ##divN1_ = [0 for i in xrange(3)]
  ##divN1_[0] = div_cartes_2_cyl([row[0] for row in N1])
  ##divN1_[1] = div_cartes_2_cyl([row[1] for row in N1])
  ##divN1_[2] = div_cartes_2_cyl([row[2] for row in N1])
  #
  # (nabla * N1[0..2][0..2]) [X] (curl(H[0..2])):
  ##divN1_x_curlH_ = cross_prod_cartes(divN1_,curl_parallel_2_cyl(Hxyz))
  #
  # N1[0..2][0..2] [X] curl(H[0..2]):
  ##curlH_ = curl_parallel_2_cyl(Hxyz)
  ##N1_x_curlH_ = [0 for i in xrange(3)]
  ##N1_x_curlH_[0] = cross_prod_cartes(N1[0],curlH_)
  ##N1_x_curlH_[1] = cross_prod_cartes(N1[1],curlH_)
  ##N1_x_curlH_[2] = cross_prod_cartes(N1[2],curlH_)
  #
  # nabla * (N1[0..2][0..2] [X] curl(H[0..2])):
  ##div_N1curlH_ = [0 for i in xrange(3)]
  ##div_N1curlH_[0] = div_cartes_2_cyl([row[0] for row in N1_x_curlH_])
  ##div_N1curlH_[1] = div_cartes_2_cyl([row[1] for row in N1_x_curlH_])
  ##div_N1curlH_[2] = div_cartes_2_cyl([row[2] for row in N1_x_curlH_])
  #
  # RHS:
  ##lapl_rhs_ = [a+b-c for a,b,c in zip(curlM_curlH_,divN1_x_curlH_,div_N1curlH_)]
  #
  #
  # 2.1.d) Cyl. coord. tens. transform.:
  ##########
  # M[0..2][0..2] in cyl. coord.:
  # _Mcyl = convert_2nd_order_tens(M)# DUPLICATE
  #
  # Convertion of M[0..2][0..2] from normalized covariant {e_rho; e_phi; e_z} to nonnormalized covariant basis {b_rho; b_phi; b_z}:
  # _Mbbb = conv_tens_eee2bbb_cyl(_Mcyl)# DUPLICATE
  #
  # M[0..2][0..2] in dyadic representation `b_covar (X) b_contravar':
  _Mcovcon = index_lower_1st_2ndord_tens_cyl( _Mbbb )
  #
  # curl(H[0..2]):
  # _curlH = curl(H)# DUPLICATE
  #
  # Convert of curl(H[0..2]) from norm. covar. basis {e_rho; e_phi; e_z} to nonnorm. covar. basis {b_rho; b_phi; b_z}:
  # _curlH_bbb = conv_vect_eee2bbb_cyl(_curlH)# DUPLICATE
  #
  # Convert curl(H[0..2]) from covar. to contravar. basis:
  # _curlH_cov = index_lower_vect_cyl(_curlH_bbb)# DUPLICATE
  #
  # (curl(M[0..2][0..2])) * (curl(H[0..2])):
  _curlM_curlH = dot_prod_l2r_tens_vect( curl_cyl_2ndord_tens___2( _Mcyl ), _curlH_cov )
  #
  _curlMcurlH_eee = conv_vect_bbb2eee_cyl( _curlM_curlH )
  #
  ##########
  # N1[0..2][0..2] in cyl. coord.:
  _N1cyl0 = convert_2nd_order_tens( N1 )
  _N1cyl = N1_dummy
  #print "7\n\n"
  #pprint( _N1cyl, wrap_line=True )
  #
  # div(N1[0..2][0..2]):
  _divN1 = div_tens_cyl( _N1cyl )
  #
  # Convert `div(N1[0..2][0..2])' from norm. covar. basis {e_rho; e_phi; e_z} to nonnorm. covar. basis {b_rho; b_phi; b_z}:
  _divN1_bbb = conv_vect_eee2bbb_cyl( _divN1 )
  #
  # Convert `div(N1[0..2][0..2])' from covar. to contravar. basis:
  _divN1_cov = index_lower_vect_cyl( _divN1_bbb )
  #
  # div(N1[0..2][0..2]) [X] curl(H[0..2]):
  _divN1_x_curlH = cross_prod_cyl( _divN1_cov, _curlH_cov )
  #
  _divN1xcurlH_eee = conv_vect_bbb2eee_cyl( _divN1_x_curlH )
  #
  ##########
  # Convert `N1[0..2][0..2]' from norm. covar. basis {e_rho; e_phi; e_z} to nonnorm. covar. basis {b_rho; b_phi; b_z}:
  _N1bbb = conv_tens_eee2bbb_cyl( _N1cyl )
  #print "6\n\n"
  #pprint( _N1bbb, wrap_line=True )
  #
  # Convert `N1[0..2][0..2]' from dyad. represent. of type `b_covar (X) b_covar' to dyad. represent. of type `b_covar (X) b_contravar':
  _N1concov = index_rise_1st_2ndord_tens_cyl( index_lower_both_2ndord_tens_cyl( _N1bbb ) )
  #print "5\n\n"
  #pprint( index_lower_both_2ndord_tens_cyl( _N1bbb ), wrap_line=True )
  #print "4\n\n"
  #pprint( _N1concov, wrap_line=True )
  #
  # N1[0..2][0..2] [X] curl(H[0..2]):
  _N1_x_curlH = cross_prod_tens_l2r_vect_cyl( _N1concov, _curlH_cov )
  #print "3\n\n"
  #pprint( _N1_x_curlH, wrap_line=True )
  _N1xcurlH_eee = conv_tens_bbb2eee_cyl( _N1_x_curlH )
  #print "2\n\n"
  #pprint( _N1xcurlH_eee, wrap_line=True )
  _div_N1xcurlH = div_tens_cyl( _N1xcurlH_eee )
  #print "1\n\n"
  #pprint( _div_N1xcurlH, wrap_line=True )
  ##########
  _lapl_rhs = [a+b-c for a,b,c in zip( _curlMcurlH_eee, _divN1xcurlH_eee, _div_N1xcurlH )]
  #
  #
  # 2.2) Left part
  #
  # 2.2.a) Cartesian coordinates
  #
  # N2[0..2][0..2] * nabla**2 * H[0..2]:
  ##N2_nabla2H__ = dot_prod_l2r_tens_vect( N2, lapl( Hcartes ) )
  #
  # k**2 * H[0..2]:
  ##k2_H__ = prod_scal2vect( k**2, Hcartes )
  # LHS:
  ##lapl_lhs__ = [a+b for a,b in zip( N2_nabla2H__, k2_H__ )]
  #
  ##### 2.2.b) Cyl. coord.
  #####
  ##### N2[0..2][0..2] in cyl. coord.:
  ####N2cyl = [0 for i in xrange(3)]
  ####N2cyl[0] = cartes_2_cyl( N2[0] )
  ####N2cyl[1] = cartes_2_cyl( N2[1] )
  ####N2cyl[2] = cartes_2_cyl( N2[2] )
  ##### N2[0..2][0..2] * nabla**2 * H[0..2]:
  ####N2_nabla2H = dot_prod_l2r_tens_vect( N2cyl, lapl_cyl( H ) )
  ##### N2[0..2][0..2] * nabla**2 * H[0..2] in cyl. coord.:
  ####N2nabla2Hcyl = cartes_2_cyl( N2_nabla2H )
  ##### k**2 * H[0..2]:
  ####k2_H = prod_scal2vect( k**2, H )
  ##### Sum:
  ####lapl_lhs = [a+b for a,b in zip( N2nabla2Hcyl, k2_H )]
  #####
  # 2.2.c) Cart. to cyl. coord.
  #
  # N2[0..2][0..2] * nabla**2 * H[0..2]:
  ##N2_nabla2H_ = dot_prod_l2r_tens_vect(N2,lapl_cartes2cyl(Hxyz))
  # k**2 * H[0..2] in cyl. coord.:
  ##k2_H_ = prod_scal2vect(k**2,Hxyz)
  # Sum:
  ##lapl_lhs_ = [a+b for a,b in zip( N2_nabla2H_, k2_H_ )]
  #
  #
  # 2.2.d) Cyl. coord. tens. transform.
  #
  # lapl(H[0..2]):
  _laplH = lapl_cyl( H )
  #
  # Convert lapl(H[0..2]) from norm. covar. basis {e_rho; e_phi; e_z} to nonnorm. covar. basis {b_rho; b_phi; b_z}:
  _laplH_bbb = conv_vect_eee2bbb_cyl( _laplH )
  #
  # Convert lapl(H[0..2]) from covar. basis to contravar. basis:
  _laplH_cov = index_lower_vect_cyl( _laplH_bbb )
  #
  #
  _N2cyl0 = convert_2nd_order_tens( N2 )
  _N2cyl = N2_dummy
  #
  # Convert N2 from norm. covar. basis {e_rho; e_phi; e_z} to nonnorm. covar. basis {b_rho; b_phi; b_z}:
  _N2bbb = conv_tens_eee2bbb_cyl( _N2cyl )
  #
  # N2[0..2][0..2] * lapl(H[0..2]):
  _N2_laplH = dot_prod_l2r_tens_vect( _N2bbb, _laplH_cov )
  # Convert N2[0..2][0..2] * lapl(H[0..2]) back to norm. covar basis {e_rho; e_phi; e_z} from nonnorm. covar. basis {b_rho; b_phi; b_z}:
  _N2laplH_eee = conv_vect_bbb2eee_cyl( _N2_laplH )
  ##########
  # k**2 * H[0..2]:
  _k2_H = prod_scal2vect( k**2, H )
  ##########
  _lapl_lhs = [a+b for a,b in zip( _N2laplH_eee, _k2_H )]
  #
  #
  #
  #
  #
  end_timestamp = timeit.default_timer()
  print 'main calculation timer:'
  print end_timestamp - bgn_timestamp
  print '\n\n'
  bgn_timestamp = timeit.default_timer()
  # Result:
  # genform = numer(simplify(curl_lhs_[0]-curl_rhs_[0]-lapl_lhs_[0]+lapl_rhs_[0]))
  genform0 = simplify( sympify( _lapl_lhs[0] - _lapl_rhs[0] ) )
  genform1 = simplify( sympify( _lapl_lhs[1] - _lapl_rhs[1] ) )
  genform2 = simplify( sympify( _lapl_lhs[2] - _lapl_rhs[2] ) )
  #
  genformA = simplify( sympify( _curl_lhs[0] - _curl_rhs[0] ) )
  genformB = simplify( sympify( _curl_lhs[1] - _curl_rhs[1] ) )
  genformC = simplify( sympify( _curl_lhs[2] - _curl_rhs[2] ) )
  #
  #genform = genform.subs(diff(diff(Hrho,rho),rho),d2Hrho_drho2)
  #genform = genform.subs(diff(diff(Hrho,phi),phi),d2Hrho_dphi2)
  #genform = genform.subs(diff(diff(Hrho,z),z),d2Hrho_dz2)
  #
  #genform = genform.subs(diff(diff(Hrho,rho),phi),d2Hrho_drhodphi)
  #genform = genform.subs(diff(diff(Hrho,phi),rho),d2Hrho_drhodphi)
  #
  #genform = genform.subs(diff(diff(Hrho,rho),z),d2Hrho_drhodz)
  #genform = genform.subs(diff(diff(Hrho,phi),z),d2Hrho_dphidz)
  #
  #genform = genform.subs(diff(diff(Hrho,z),rho),d2Hrho_drhodz)
  #genform = genform.subs(diff(diff(Hrho,z),phi),d2Hrho_dphidz)
  #
  #genform = genform.subs(diff(Hrho,rho),dHrho_drho)
  #genform = genform.subs(diff(Hrho,phi),dHrho_dphi)
  #genform = genform.subs(diff(Hrho,z),dHrho_dz)
  #
  end_timestamp = timeit.default_timer()
  print 'vector simplification timer:'
  print end_timestamp - bgn_timestamp
  print '\n\n'
  bgn_timestamp = timeit.default_timer()
  #
  genforms = [genform0, genform1, genform2, genformA, genformB, genformC]
  rhophiz = [rho, phi, z]
  Hrhophiz = [Hrho, Hphi, Hz]
  rhophiz_text = ["rho", "phi", "z"]
  Hrhophiz_text = ["Hrho", "Hphi", "Hz"]
  #
  # Fourier method for variable separation:
  #
  #f_Hrho_R = Function('Hrho_R')
  #f_Hrho_P = Function('Hrho_P')
  #f_Hrho_T = Function('Hrho_T')
  #f_Hrho_Z = Function('Hrho_Z')
  #f_Hphi_R = Function('Hphi_R')
  #f_Hphi_P = Function('Hphi_P')
  #f_Hphi_T = Function('Hphi_T')
  #f_Hphi_Z = Function('Hphi_Z')
  #f_Hrho_R = Function('Hz_R')
  #f_Hrho_P = Function('Hz_P')
  #f_Hrho_T = Function('Hz_T')
  #f_Hrho_Z = Function('Hz_Z')
  #
  #Hrho_R = Hrho_R( rho )
  #Hrho_P = Hrho_P( phi )
  #Hrho_T = Hrho_T( rho, phi )
  #Hrho_Z = Hrho_Z( z )
  #Hphi_R = Hphi_R( rho )
  #Hphi_P = Hphi_P( phi )
  #Hphi_T = Hphi_T( rho, phi )
  #Hphi_Z = Hphi_Z( z )
  #Hrho_R = Hz_R( rho )
  #Hrho_P = Hz_P( phi )
  #Hrho_T = Hz_T( rho, phi )
  #Hrho_Z = Hz_Z( z )
  #
  #Hrhophiz_SV1 = [Hrho_T, Hrho_Z, Hphi_T, Hphi_Z, Hz_T, Hz_Z]
  #Hrhophiz_SV2 = [Hrho_R, Hrho_P, Hphi_R, Hphi_P, Hz_R, Hz_P]
  #
  #for i in range(0,3):
  #  for j in range(0,3):
  #    genforms[i] = genforms[i].subs( Hrhophiz[j], Hrhophiz_SV1[j] * Hrhophiz_SV1[j+1] )
  #  genforms[i] = simplify( genforms[i] )
  #
  #for i in range(0,3):
  #  for j in range(0,3):
  #    genforms[i] = genforms[i].subs( Hrhophiz[j], Hrhophiz_SV2[j] * Hrhophiz_SV2[j+1] )
  #  genforms[i] = simplify( genforms[i] )
  #
  for i in range(0,3):
    for j in range(0,3):
      for k in range(0,3):
        genforms[i] = genforms[i].subs( M_dummy[j][k], _Mcyl0[j][k] )
        genforms[i] = genforms[i].subs( N1_dummy[j][k], _N1cyl0[j][k] )
        genforms[i] = genforms[i].subs( N2_dummy[j][k], _N2cyl0[j][k] )
    genforms[i] = simplify( genforms[i] )
  #
  end_timestamp = timeit.default_timer()
  print 'substitution timer:'
  print end_timestamp - bgn_timestamp
  print '\n\n'
  bgn_timestamp = timeit.default_timer()
  #
  #
  # Output:
  #genform_collect_Hrho = collect(genform,Hrho,evaluate=False)
  #genform_collect_Hrhodz = collect(genform,dHrho_dz,evaluate=False)
  #genform_collect_Hrhodphi = collect(genform,dHrho_dphi,evaluate=False)
  #genform_collect_Hrhodrho = collect(genform,dHrho_drho,evaluate=False)
  #genform_collect_Hrhodphidz = collect(genform,d2Hrho_dphidz,evaluate=False)
  #genform_collect_Hrhodrhodz = collect(genform,d2Hrho_drhodz,evaluate=False)
  #genform_collect_Hrhodrhodphi = collect(genform,d2Hrho_drhodphi,evaluate=False)
  #genform_collect_Hrhodz2 = collect(genform,d2Hrho_dz2,evaluate=False)
  #genform_collect_Hrhodphi2 = collect(genform,d2Hrho_dphi2,evaluate=False)
  #genform_collect_Hrhodrho2 = collect(genform,d2Hrho_drho2,evaluate=False)
  #
  #pprint( genform0, wrap_line=True )
  #print( "\n\n" )
  #pprint( genform1, wrap_line=True )
  #print( "\n\n" )
  #pprint( genform2, wrap_line=True )
  #print( "\n\n" )
  #
  #genform_s = simplify( sympify( genform_collect_Hrhodphi2[S.One] ) )
  for l in range(0,3):
    print "\n"
    print "\n"
    print "\n"
    print "\n"
    print "\n"
    for i in range(0,3):
      for j in range(0,3):
        print "%s / %s :: " % (Hrhophiz_text[i], rhophiz_text[j])
        pprint( genforms[l].has( diff( Hrhophiz[i], rhophiz[j] ) ) )
        pprint( genforms[l+3].has( diff( Hrhophiz[i], rhophiz[j] ) ) )
        print "\n"
    for i in range(0,3):
      for j in range(0,3):
        for k in range(0,3):
          print "%s / %s %s :: " % (Hrhophiz_text[i], rhophiz_text[j], rhophiz_text[k])
          pprint( genforms[l].has( diff( diff( Hrhophiz[i], rhophiz[j] ), rhophiz[k] ) ) )
          pprint( genforms[l+3].has( diff( diff( Hrhophiz[i], rhophiz[j] ), rhophiz[k] ) ) )
          print "\n"
  #
  end_timestamp = timeit.default_timer()
  print 'search timer:'
  print end_timestamp - bgn_timestamp
  print '\n\n'
  bgn_timestamp = timeit.default_timer()
  #
  f = open( 'test.sympy.tex', 'w' )
  f.write( '\\documentclass{article}\n\\usepackage{amsmath}\n\\begin{document}\n\\begin{equation*}\n\\begin{split*}\n&' )
  f.write( latex( genform0, mode='inline') )# add \documentclass{article}\usepackage{amsmath}\begin{document}<...>\end{document}
  f.write( ' \\\n&' )
  f.write( latex( genform1, mode='inline') )# add \documentclass{article}\usepackage{amsmath}\begin{document}<...>\end{document}
  f.write( ' \\\n&' )
  f.write( latex( genform2, mode='inline') )# add \documentclass{article}\usepackage{amsmath}\begin{document}<...>\end{document}
  f.write( ' \\\n&' )
  f.write( ' \\\n&' )
  f.write( ' \\\n&' )
  f.write( latex( genformA, mode='inline') )# add \documentclass{article}\usepackage{amsmath}\begin{document}<...>\end{document}
  f.write( ' \\\n&' )
  f.write( latex( genformB, mode='inline') )# add \documentclass{article}\usepackage{amsmath}\begin{document}<...>\end{document}
  f.write( ' \\\n&' )
  f.write( latex( genformC, mode='inline') )# add \documentclass{article}\usepackage{amsmath}\begin{document}<...>\end{document}
  f.write( '\n\\end{split*}\n\\end{equation*}\n\\end{document}' )
  f.close()
  #
  end_timestamp = timeit.default_timer()
  print 'print timer:'
  print end_timestamp - bgn_timestamp
  print '\n\n'
  bgn_timestamp = timeit.default_timer()
print 'v1.1.build61'
