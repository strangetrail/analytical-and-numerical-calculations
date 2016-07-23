#!/usr/bin/env python
#
# Note that : (sic!) N2 is multiplied on lapl. from l2r. in the main formulae \
# with 18 derivatives
#
import sys
sys.path.append( '/home/thisuser/Documents/code/sympy_tests/common/' )
from functools import reduce
from sympy import *
from utils import *
from mathutils import *
from test4 import get_x_2nd
#
beginSimpleTiming()
#
eTensors = enum(DIFF=1, LAPL=2, CURL=3)
#
###############################################################################
#                                   SETTINGS                                  #
###############################################################################
#
x, y, z = [S.Zero for idx in xrange(3)]
bCheckAutoCalculated     =  True
bSymmetric               =  True
bMTens_XYZ               =  True
bKerr_2ndOrd_tens        = False
bSolveEigenValue         = False
bDoNotHandleDeterminant  = False# It's better to rather not simplify it, then \
#                               # skip expanding.
bForceExpandSimplified   = False# Whether to expand generic forms in tex file
bRemovePlainWaveExponent =  True
bPrintAll                =  True
bDoNotCountDerivatives   =  True
bCartesian               =  True
bCurvilinear             = False
bFourierDomain           =  True
bFourPlaneWaveOnly       = False
bVarSeparate             = False
bDoNotSubstitute         =  True# DO NOT Replace the function                 \
#                               #  Gamma ( alpha_1, alpha_2 ) with variable   \
#                               #  Gamma_s.                                    
bSkipInitialTests        =  True
#
eTensor4Eigenvalues      = eTensors.DIFF
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
i,j,k,l,n,m,t = 0,0,0,0,0,0,0
SPACE, gamma_TEST, gamma, u_i, u_j, u_k, alpha_1, alpha_2, e_rho, e_phi =   (
 symbols
 (
  'SPACE, gamma_TEST, gamma, u_i, u_j, u_k, alpha_1, alpha_2, e_rho, e_phi'
 )
                                                                            )
#
e_z, alpha, beta, epsilon, rho, phi, theta_1, theta_2, theta_3, x, y, z, k_0, isum, jsum, ksum, nsum, msum, lsum, tsum, detM, dHrho_dz = (
 symbols
 (
  'e_z, alpha, beta, epsilon, rho, phi, theta_1, theta_2, theta_3, x, y, z, k_0, i_s, j_s, k_s, n_s, m_s, l_s, t_s, detM, dHrho_dz'
 )
                                                                            )
#
dHrho_dphi, dHrho_drho, d2Hrho_dphidz, d2Hrho_drhodz, d2Hrho_drhodphi =   (
 symbols
 (
  'dHrho_dphi, dHrho_drho, d2Hrho_dphidz, d2Hrho_drhodz, d2Hrho_drhodphi'
 )
                                                                          )
#
d2Hrho_dz2, d2Hrho_dphi2, d2Hrho_drho2, C_r, n_0, n_e, Gamma_s, mu =   (
 symbols
 (
  'd2Hrho_dz2, d2Hrho_dphi2, d2Hrho_drho2, C_r, n_0, n_e, Gamma_s, mu'
 )
                                                                       )
#
###############################################################################
#^^^^^                              SYMBOLS                              ^^^^^#
###############################################################################
#
alphas = [alpha_1, alpha_2]
Gamma_f = Function( 'Gamma_f' )
Gamma = Gamma_f( alpha_1, alpha_2 )
#Gamma = Gamma_s
plain_wave_exp_curv = exp ( I*k_0 * ( alpha_1*theta_1 + alpha_2*theta_2 + Gamma*theta_3 ) )
plain_wave_exp_cyl = exp ( I*k_0 * ( alpha_1*rho + alpha_2*phi + Gamma*z ) )
plain_wave_exp = exp ( I*k_0 * ( alpha_1*x + alpha_2*y + Gamma*z ) )
print "Substitution test\n"
pprint( (6*plain_wave_exp + tan(phi/plain_wave_exp)).subs(plain_wave_exp, S.One) )
print "\n"
FDvars = [alpha_1, alpha_2, Gamma]
def removePWE ( expr ):
  #print "\nBEFORE:\n"
  #pprint (expr)
  expr2 = (expand(simplify(sympify(expr)))).subs( plain_wave_exp, S.One )
  expr3 = (expr2).subs( plain_wave_exp_cyl, S.One )
  expr4 = (expr3).subs( plain_wave_exp_curv, S.One )
  #print "\nAFTER:\n"
  #pprint (expr2)
  return expr4
  #return sum(collect( expr, plain_wave_exp, evaluate = False, exact = False ).values())
pass# DEF removePWE ( expr )
#
###############################################################################
#                                   CLASSES                                   #
###############################################################################
#
###class FunctionFD ( Function ):
###  """METADATA"""
###  def __new__(cls, *args, **options):
###    result = super(Function, cls).__new__(cls, *args, **options)
###    return result
###  #
###  # Fourier domain derivative:
###  #
###  def fdiff ( self, argindex=1 ):
###    if ( not (1 <= argindex <= len(self.args)) ):
###      raise ArgumentIndexError(self, argindex)
###    return I * FDvars[argindex] * Function ( self.__class__.__name__ )(*alphas)
#
###############################################################################
#^^^^^                              CLASSES                              ^^^^^#
###############################################################################
#
if ( not bSkipInitialTests ):
  collect_test = collect( alpha**8 + alpha**3 + beta, alpha**3, exact=True, evaluate=False )
  pprint( collect_test )
  collect_test = collect( alpha**8 + alpha**3 + beta, alpha**3, exact=True, evaluate=True )
  pprint( collect_test )
  collect_test = collect( alpha**8 + alpha**3 + beta, alpha**3, exact=False, evaluate=False )
  pprint( collect_test )
  collect_test = collect( alpha**8 + alpha**3 + beta, alpha**3, exact=False, evaluate=True )
  pprint( collect_test )
  collect_test = collect( alpha**8 + alpha**3 + beta, alpha**8, exact=True, evaluate=False )
  pprint( collect_test )
  collect_test = collect( alpha**8 + alpha**3 + beta, alpha**8, exact=True, evaluate=True )
  pprint( collect_test )
  collect_test = collect( alpha**8 + alpha**3 + beta, alpha**8, exact=False, evaluate=False )
  pprint( collect_test )
  collect_test = collect( alpha**8 + alpha**3 + beta, alpha**8, exact=False, evaluate=True )
  pprint( collect_test )
#
#
# Tens. components:
f_M11 = Function( 'M_11' )
f_M12 = Function( 'M_12' )
f_Maa = Function( 'M_a' )
f_M13 = Function( 'M_13' )
f_Mbb = Function( 'M_b' )
f_M21 = Function( 'M_21' )
f_M22 = Function( 'M_22' )
f_M23 = Function( 'M_23' )
f_Mcc = Function( 'M_c' )
f_M31 = Function( 'M_31' )
f_M32 = Function( 'M_32' )
f_M33 = Function( 'M_33' )
f_mu11 = Function( 'mu_11' )
f_mu12 = Function( 'mu_12' )
f_mu13 = Function( 'mu_13' )
f_mu21 = Function( 'mu_21' )
f_mu22 = Function( 'mu_22' )
f_mu23 = Function( 'mu_23' )
f_mu31 = Function( 'mu_31' )
f_mu32 = Function( 'mu_32' )
f_mu33 = Function( 'mu_33' )
f_eps11 = Function( 'eps_11' )
f_eps12 = Function( 'eps_12' )
f_eps13 = Function( 'eps_13' )
f_eps21 = Function( 'eps_21' )
f_eps22 = Function( 'eps_22' )
f_eps23 = Function( 'eps_23' )
f_eps31 = Function( 'eps_31' )
f_eps32 = Function( 'eps_32' )
f_eps33 = Function( 'eps_33' )
#
# Dummy tens. components:
f_MD11 = Function('MD_11')
f_MD12 = Function('MD_12')
f_MD13 = Function('MD_13')
f_MD21 = Function('MD_21')
f_MD22 = Function('MD_22')
f_MD23 = Function('MD_23')
f_MD31 = Function('MD_31')
f_MD32 = Function('MD_32')
f_MD33 = Function('MD_33')
#
# Dummy tens. components:
f_NN11 = Function('NN_11')
f_NN12 = Function('NN_12')
f_NN13 = Function('NN_13')
f_NN21 = Function('NN_21')
f_NN22 = Function('NN_22')
f_NN23 = Function('NN_23')
f_NN31 = Function('NN_31')
f_NN32 = Function('NN_32')
f_NN33 = Function('NN_33')
#
# Dummy tens. components:
f_N11 = Function('N_11')
f_N12 = Function('N_12')
f_N13 = Function('N_13')
f_N21 = Function('N_21')
f_N22 = Function('N_22')
f_N23 = Function('N_23')
f_N31 = Function('N_31')
f_N32 = Function('N_32')
f_N33 = Function('N_33')
#
# Basic vect./tens. field components:
f_Hrho = Function('H_rho')
f_Hphi = Function('H_phi')
f_Hx = Function('H_x')
f_Hy = Function('H_y')
f_Hz = Function('H_z')# Check this part.
f_Ex = Function('E_x')
f_Ey = Function('E_y')
f_Ez = Function('E_z')
#
# Replace all with SCRIPT_subscript notations - for tex compatibility
#
xyz = [x, y, z]
#
Hrho = f_Hrho(rho,phi,z)
Hphi = f_Hphi(rho,phi,z)
Hz_ = f_Hz(rho,phi,z)
#
Hx = f_Hx(x,y,z)
Hy = f_Hy(x,y,z)
Hz = f_Hz(x,y,z)
Ex = f_Ex(x,y,z)
Ey = f_Ey(x,y,z)
Ez = f_Ez(x,y,z)
#==================================================================================================
#
#
#
# Index of refraction tens.:
if ( bKerr_2ndOrd_tens ):
  #
  #
  #  Tens. components in quasi steady-state optical Kerr eff. case:
  Chi_2nd = get_x_2nd()# Iteration instead of loop or explicit matrix initialization is used here \
  #                    # because the function for propagators, given below, was derived manually  \
  #                    # without automation and assistance of SymPy objects.
  M11 = Chi_2nd[0][0]
  M12 = Chi_2nd[0][1]
  M13 = Chi_2nd[0][2]
  M21 = Chi_2nd[1][0]
  M22 = Chi_2nd[1][1]
  M23 = Chi_2nd[1][2]
  M31 = Chi_2nd[2][0]
  M32 = Chi_2nd[2][1]
  M33 = Chi_2nd[2][2]
else:
  #
  #
  #  Tens. components in common case:
  if ( bMTens_XYZ ):
    if ( not bSymmetric ):
      #
      # X, Y, Z basis:
      M11 = f_M11(x,y,z)
      M12 = f_M12(x,y,z)
      M13 = f_M13(x,y,z)
      M21 = f_M21(x,y,z)
      M22 = f_M22(x,y,z)
      M23 = f_M23(x,y,z)
      M31 = f_M31(x,y,z)
      M32 = f_M32(x,y,z)
      M33 = f_M33(x,y,z)
    else:
      #
      # Consider symmetric 2nd ord. Kerr tens.
      #
      M11 = f_M11(x,y,z)
      M12 = f_Maa(x,y,z)
      M13 = f_Mbb(x,y,z)
      M21 = f_Maa(x,y,z)# Minus comes from 2D tensor components in     \
      #                        #  principal basis.                             
      M22 = f_M22(x,y,z)
      M23 = f_Mcc(x,y,z)
      M31 = f_Mbb(x,y,z)
      M32 = f_Mcc(x,y,z)
      M33 = f_M33(x,y,z)
    pass# IF
  else:
    #
    # Rho, Phi, Z basis:
    M11 = f_M11(rho,phi,z)
    M12 = f_M12(rho,phi,z)
    M13 = f_M13(rho,phi,z)
    M21 = f_M21(rho,phi,z)
    M22 = f_M22(rho,phi,z)
    M23 = f_M23(rho,phi,z)
    M31 = f_M31(rho,phi,z)
    M32 = f_M32(rho,phi,z)
    M33 = f_M33(rho,phi,z)
  pass# IF
pass# IF
#==================================================================================================
#
n_d = n_e - n_0
#
#M11 = n_0 + n_d * (cos( C_r * rho ))**2 * (cos( phi ))**2
#M12 = -n_d * (cos( C_r * rho ))**2 * cos( phi ) * sin( phi )
#M13 = -n_d * cos( C_r * rho ) * sin( C_r * rho ) * cos( phi )
#M21 = -n_d * (cos( C_r * rho ))**2 * cos( phi ) * sin( phi )
#M22 = n_0 + n_d * (cos( C_r * rho ))**2 * (sin( phi ))**2
#M23 = n_d * cos( C_r * rho ) * sin( C_r * rho ) * sin( phi )
#M31 = -n_d * cos( C_r * rho ) * sin( C_r * rho ) * cos( phi )
#M32 = n_d * cos( C_r * rho ) * sin( C_r * rho ) * sin( phi )
#M33 = n_0 + n_d * (sin( C_r * rho ))**2
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
H = [0 for i in xrange(3)]
H[0] = exp( I * beta * z ) * Hrho
H[1] = exp( I * beta * z ) * Hphi
H[2] = exp( I * beta * z ) * Hz_
#
# Vect. field in cartes. coord. expressed in terms of cyl. coord.:
Hxyz = [0 for i in xrange(3)]
Hxyz[0] = exp( I * beta * z ) * Hrho * cos( Hphi )
Hxyz[1] = exp( I * beta * z ) * Hrho * sin( Hphi )
Hxyz[2] = exp( I * beta * z ) * Hz_
#
# Vect. field in cartes. coord.:
Hcartes = [0 for i in xrange(3)]
Hcartes[0] = Hx
Hcartes[1] = Hy
Hcartes[2] = Hz
#
#
#
# Fourier domain representation for plain wave decomposition:
#
F = [0 for i in xrange(3)]
Hfour = [0 for i in xrange(3)]
__Hfour = [0 for i in xrange(3)]
Efour = [0 for i in xrange(3)]
#
Fx_f = Function( 'F_x' )
F_x = Fx_f( alpha_1, alpha_2 )
F[0] = F_x
Fy_f = Function( 'F_y' )
F_y = Fy_f( alpha_1, alpha_2 )
F[1] = F_y
Fz_f = Function( 'F_z' )
F_z = Fz_f( alpha_1, alpha_2 )
F[2] = F_z
#
#
#
Hdummy = [0 for i in xrange(3)]
Hdummy[0] = Hx
Hdummy[1] = Hy
Hdummy[2] = 0
Edummy = [0 for i in xrange(3)]
Edummy[0] = Ex
Edummy[1] = Ey
Edummy[2] = 0
#
#
Hfour[2] = F[2] * plain_wave_exp
Hfour[1] = F[1] * plain_wave_exp
Hfour[0] = F[0] * plain_wave_exp
__Hfour[2] = F[2] * plain_wave_exp_curv
__Hfour[1] = F[1] * plain_wave_exp_curv
__Hfour[0] = F[0] * plain_wave_exp_curv
Efour[2] = F[2] * plain_wave_exp_cyl
Efour[1] = F[1] * plain_wave_exp_cyl
Efour[0] = F[0] * plain_wave_exp_cyl
#
#
# Test values for symmetric tensor M.
#
EigenFTLaplPropagator =                   (
 [
  [
   -(alpha_1**2) * k_0**2 * M22
   - alpha_1**2 * k_0**2 * M33
   - alpha_1 * alpha_2 * k_0**2 * M12
   - 2 * alpha_2**2 * k_0**2 * M22
   - alpha_2**2 * k_0**2 * M33
   - I * alpha_2 * k_0 * diff( M23, z )
   + I * alpha_2 * k_0 * diff( M33, y )
   + k_0**2,
   2 * alpha_1**2 * k_0**2 * M12
   + alpha_1 * alpha_2 * k_0**2 * M22
   + I * alpha_1 * k_0 * diff( M23, z )
   - I * alpha_1 * k_0 * diff( M33, y )
   + alpha_2**2 * k_0**2 * M12,
   2 * alpha_1**2 * k_0**2 * M13
   - alpha_1 * alpha_2 * k_0**2 * M23
   - I * alpha_1 * k_0 * diff( M22, z )
   + I * alpha_1 * k_0 * diff( M23, y )
   + alpha_2**2 * k_0**2 * M13
   + I * alpha_2 * k_0 * diff( M12, z )
   - I * alpha_2 * k_0 * diff( M13, y )
  ],
  [
   alpha_1**2 * k_0**2 * M12
   + alpha_1 * alpha_2 * k_0**2 * M11
   + I * alpha_2 * k_0 * diff( M13, z )
   - I * alpha_2 * k_0 * diff( M33, x ),
   -2 * alpha_1**2 * k_0**2 * M11
   - alpha_1**2 * k_0**2 * M33
   + alpha_1 * alpha_2 * k_0**2 * M12
   - I * alpha_1 * k_0 * diff( M13, z )
   + I * alpha_1 * k_0 * diff( M33, x )
   - alpha_2**2 * k_0**2 * M11
   - alpha_2**2 * k_0**2 * M33
   + k_0**2,
   alpha_1**2 * k_0**2 * M23
   + alpha_1 * alpha_2 * k_0**2 * M13
   + I * alpha_1 * k_0 * diff( M12, z )
   - I * alpha_1 * k_0 * diff( M23, x )
   + 2 * alpha_2**2 * k_0**2 * M23
   - I * alpha_2 * k_0 * diff( M11, z )
   + I * alpha_2 * k_0 * diff( M13, x )
  ],
  [
   alpha_1**2 * k_0**2 * M13
   + alpha_2**2 * k_0**2 * M13
   - I * alpha_2 * k_0 * diff( M13, y )
   + I * alpha_2 * k_0 * diff( M23, x ),
   alpha_1**2 * k_0**2 * M23
   + I * alpha_1 * k_0 * diff( M13, y )
   - I * alpha_1 * k_0 * diff( M23, x )
   + alpha_2**2 * k_0**2 * M23,
   -2 * alpha_1**2 * k_0**2 * M11
   - alpha_1**2 * k_0**2 * M22
   -2 * alpha_1 * alpha_2 * k_0**2 * M12
   - I * alpha_1 * k_0 * diff( M12, y )
   + I * alpha_1 * k_0 * diff( M22, x )
   - alpha_2**2 * k_0**2 * M11
   + I * alpha_2 * k_0 * diff( M11, y )
   - I * alpha_2 * k_0 * diff( M12, x )
   + k_0**2
  ]
 ]
                                      )
EigenFTLapltens1 =(
 [
  [
   -alpha_1 * k_0**2 * M13
   + I * k_0 * diff( M22, z )
   - I * k_0 * diff( M23, y )
   ,
   alpha_1 * k_0**2 * M23
   - I * k_0 * diff(M12, z)
   + I * k_0 * diff( M13, y )
   ,
   - alpha_1 * k_0**2 * M33
   ],
 [
  - alpha_2 * k_0**2 * M13
  - I * k_0 * diff( M12, z )
  + I * k_0 * diff( M23, x ),
  -alpha_2 * k_0**2 * M23
  + I * k_0 * diff( M11, z )
  - I * k_0 * diff( M13, x ),
  alpha_2 * k_0**2 * M33
 ],
 [
  alpha_2 * k_0**2 * M12
  + I * k_0 * diff( M12, y )
  - I * k_0 * diff( M22, x )
  + alpha_1 * k_0**2 * M11
  ,
  alpha_1 * k_0**2 * M12
  - alpha_2 * k_0**2 * M22
  - I * k_0 * diff( M11, y )
  + I * k_0 * diff( M12, x )
  ,
  -alpha_1 * k_0**2 * M13
  + alpha_2 * k_0**2 * M23
 ]
 ])
EigenFTLapltens2 = (
 [
  [
   -k_0**2 * M22,
   k_0**2 * M12,
   k_0**2 * M13
  ],
  [
   k_0**2 * M12,
   -k_0**2 * M11 - 2 * k_0**2 * M33,
   k_0**2 * M23
  ],
  [
   2 * k_0**2 * M13
   ,
   0
   ,
   -k_0**2 * M11 - k_0**2 * M22
  ]
 ]
               )
#
#
#
#
differenceEqGamma1Tens = (
[
 [
  alpha_1 * k_0**2 * M13
  + 2 * alpha_2 * k_0**2 * M23
  ,
  -2 * alpha_1 * k_0**2 * M23
  - alpha_2 * k_0**2 * M13
  ,
  alpha_1 * k_0**2 * M22
  + alpha_1 * k_0**2 * M33
  - alpha_2 * k_0**2 * M12
 ],
 [
  -alpha_1 * k_0**2 * M23
  ,
  2 * alpha_1 * k_0**2 * M13
  + alpha_2 * k_0**2 * M23
  ,
  -alpha_1 * k_0**2 * M12
  + alpha_2 * k_0**2 * M11
  - alpha_2 * k_0**2 * M33
 ],
 [
  -alpha_1 * k_0**2 * M11
  + alpha_1 * k_0**2 * M22
  - 2 * alpha_2 * k_0**2 * M12
  ,
  -2 * alpha_1 * k_0**2 * M12
  + alpha_2 * k_0**2 * M11
  + alpha_2 * k_0**2 * M22
  ,
  alpha_1 * k_0**2 * M13
  - alpha_2 * k_0**2 * M23
 ]
]
)
differenceEqGamma2Tens = (
[
 [
  0
  ,
  0
  ,
  -k_0**2 * M13
 ],
 [
  0
  ,
  2 * k_0**2 * M33
  ,
  -k_0**2 * M23
 ],
 [
  -2 * k_0**2 * M13
  ,
  0
  ,
  k_0**2 * M11
  + k_0**2 * M22
 ]
]
)
differenceEqTens = (
[
 [
  alpha_1**2 * k_0**2 * M22
  + alpha_1**2 * k_0**2 * M33
  + alpha_1 * alpha_2 * k_0**2 * M12
  + 2 * alpha_2**2 * k_0**2 * M22
  ,
  -2 * alpha_1**2 * k_0**2 * M12
  - alpha_1 * alpha_2 * k_0**2 * M22
  + alpha_1 * alpha_2 * k_0**2 * M33
  - alpha_2**2 * k_0**2 * M12
  ,
  -2 * alpha_1**2 * k_0**2 * M13
 ],
 [
  -alpha_1**2 * k_0**2 * M12
  - alpha_1 * alpha_2 * k_0**2 * M11
  + alpha_1 * alpha_2 * k_0**2 * M33
  ,
  2 * alpha_1**2 * k_0**2 * M11
  - alpha_1 * alpha_2 * k_0**2 * M12
  + alpha_2**2 * k_0**2 * M11
  + alpha_2**2 * k_0**2 * M33
  ,
  -2 * alpha_1 * alpha_2 * k_0**2 * M13
  - 2 * alpha_2**2 * k_0**2 * M23
 ],
 [
  -alpha_1**2 * k_0**2 * M13
  - alpha_1 * alpha_2 * k_0**2 * M23
  ,
  -alpha_1 * alpha_2 * k_0**2 * M13
  - alpha_2**2 * k_0**2 * M23
  ,
  2 * alpha_1**2 * k_0**2 * M11
  + 4 * alpha_1 * alpha_2 * k_0**2 * M12
 ]
]
)
#
auto_differenceEqGamma2Tens = [[0 for i in xrange(3)] for i in xrange(3)]
auto_differenceEqGamma1Tens = [[0 for i in xrange(3)] for i in xrange(3)]
auto_differenceEqTens       = [[0 for i in xrange(3)] for i in xrange(3)]
#
auto_EigenFTLapltens2       = [[0 for i in xrange(3)] for i in xrange(3)]
auto_EigenFTLapltens1       = [[0 for i in xrange(3)] for i in xrange(3)]
auto_EigenFTLaplPropagator  = [[0 for i in xrange(3)] for i in xrange(3)]
#
auto_EigenFTCurltens2       = [[0 for i in xrange(3)] for i in xrange(3)]
auto_EigenFTCurltens1       = [[0 for i in xrange(3)] for i in xrange(3)]
auto_EigenFTCurlPropagator  = [[0 for i in xrange(3)] for i in xrange(3)]
#
auto_greensfunction_E       = [[0 for i in xrange(3)] for i in xrange(3)]
#
#muT = [
#       [f_mu11(x, y, z), f_mu12(x, y, z), f_mu13(x, y, z)],
#       [f_mu21(x, y, z), f_mu22(x, y, z), f_mu23(x, y, z)],
#       [f_mu31(x, y, z), f_mu32(x, y, z), f_mu33(x, y, z)],
#      ]
muT = [
       [mu, 0, 0],
       [0, mu, 0],
       [0, 0, mu],
      ]
epsilonT = [
            [f_eps11(x, y, z), f_eps12(x, y, z), f_eps13(x, y, z)],
            [f_eps21(x, y, z), f_eps22(x, y, z), f_eps23(x, y, z)],
            [f_eps31(x, y, z), f_eps32(x, y, z), f_eps33(x, y, z)],
           ]
#
#
def summ2( T, V, idx ):
  summ2_result = 0;
  for i in range(0,2):
    summ2_result += T[idx][i] * V[i]
  return summ2_result
def summ3( T, V, idx ):
  summ3_result = 0;
  for i in range(0,3):
    summ3_result += T[idx][i] * V[i]
  return summ3_result
#
#
# Helicoidal (rotating frame) coord.:
v = x * sin( alpha * z ) + y * cos( alpha * z )
u = x * cos( alpha * z ) - y * sin( alpha * z )
#
# Basic 2nd. ord. tens. for inverse squared refr. index:
M = [[0 for i in xrange(3)] for i in xrange(3)]
#M =           (
#    [
#     [0,0,0],
#     [0,0,0],
#     [0,0,0]
#    ]         )
#
# Symmetric variant of tensor M for test values given above:
#
M = [
     [M11, M12, M13],
     [M21, M22, M23],
     [M31, M32, M33]
    ]
epsilon_1 = 1/epsilon
rhophiz = [rho, phi, z]
theta = [theta_1, theta_2, theta_3]
lambda_11 = Function('lambda_11')(*alphas) * plain_wave_exp
__M = [
       [epsilon_1, epsilon_1,         0],
       [epsilon_1, epsilon_1,         0],
       [        0,         0, lambda_11]
      ]
#M = [
#     [ M11, M12, M13],
#     [-M12, M22, M23],
#     [-M13, M32, M33]
#    ]
#M = [
#     [M11, M12, M13],
#     [M21, M22, M23],
#     [M31, M32, M22]
#    ]
#M =        (
#    [
#     [M11, M12, M13],
#     [M21, M22, M23],
#     [M31, M32, M33]
#    ]
#           )
M_dummy = [
           [MD11, MD12, MD13],
           [MD21, MD22, MD23],
           [MD31, MD32, MD33]
          ]
#
rot_left = [
            [ cos( phi ), sin( phi ), 0],
            [-sin( phi ), cos( phi ), 0],
            [   0,          0,        1]
           ]
rot_right = [
             [cos( phi ), -sin( phi ), 0],
             [sin( phi ),  cos( phi ), 0],
             [  0,           0,        1]
            ]
#
Identity_M = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
#
# Cyl. coord. metrics and Lame coefficients:
g_covar = [
           [1,   0,    0],
           [0, rho**2, 0],
           [0,   0,    1]
          ]
g_contravar = [
               [1,   0,        0],
               [0, 1 / rho**2, 0],
               [0,   0,        1]
              ]
g_scalar = prod ( [g_covar[i][i] for i in xrange(3)] )
#
# Curv. coord. metric and scale factors:
#
f_MD21 = Function('g_11')
gCurv_cov =                                         (
            [
             [
              (
               Function('g_'+repr(idx1)+repr(idx2))
              )(*alphas)
              * plain_wave_exp
              for idx2 in xrange(3)
             ]
             for idx1 in xrange(3)
            ]                                       )
#
gCurv_con =                                         (
            [
             [
              (
               Function('g__'+repr(idx1)+repr(idx2))
              )(*alphas)
              * plain_wave_exp
              for idx2 in xrange(3)
             ]
             for idx1 in xrange(3)
            ]                                       )
#
# Basis vectors triple product:
g_sqrt = (Function( 'g_sqrt' ))(*alphas) * plain_wave_exp
g_sqrt1 = (Function( 'g_sqrt1' ))(*alphas) * plain_wave_exp
gCurv = g_sqrt
#
# Orthogonal basis scale factors:
#
h = [gCurv_cov[idx0][idx0] for idx0 in xrange(3)]
h1 = [Function('g1_'+repr(idx0)+repr(idx0))(*alphas) * plain_wave_exp for idx0 in xrange(3)]
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
N1 =                                           (
     [
      [M[1-1][1-1], M[2-1][1-1], M[3-1][1-1]],
      [M[1-1][2-1], M[2-1][2-1], M[3-1][2-1]],
      [M[1-1][3-1], M[2-1][3-1], M[3-1][3-1]]
     ]                                         )
__N1 = ([
        [__M[1-1][1-1], __M[2-1][1-1], __M[3-1][1-1]],
        [__M[1-1][2-1], __M[2-1][2-1], __M[3-1][2-1]],
        [__M[1-1][3-1], __M[2-1][3-1], __M[3-1][3-1]]
       ])
N1_dummy = [
            [N11, N12, N13],
            [N21, N22, N23],
            [N31, N32, N33]
           ]
#
# Additional 2nd ord. tens. for lapl. expression ... :
N2 = [[0 for i in xrange(3)] for i in xrange(3)]
__N2 = [[0 for i in xrange(3)] for i in xrange(3)]
N2_dummy = [
            [NN11, NN12, NN13],
            [NN12, NN22, NN23],
            [NN31, NN32, NN33]
           ]
#
# ... and it's components in general form:
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
# Same as above, but simplified (bevel basis test in comments):
N2[0][0] =  M[3-1][3-1] + M[2-1][2-1]# 33 + 22 = 11 # V
N2[0][1] = -M[2-1][1-1]#                  - 21 = 12 # V
N2[0][2] = -M[3-1][1-1]#                  -  0 = 13 # V
N2[1][0] = -M[1-1][2-1]#                  - 12 = 21 # V
N2[1][1] =  M[3-1][3-1] + M[1-1][1-1]# 33 + 11 = 22 # V
N2[1][2] = -M[3-1][2-1]#                  -  0 = 23 # V
N2[2][0] = -M[1-1][3-1]#                  -  0 = 31 # V
N2[2][1] = -M[2-1][3-1]#                  -  0 = 32 # V
N2[2][2] =  M[1-1][1-1] + M[2-1][2-1]# 11 + 22 = 33 # V
#
__N2[0][0] =  __M[3-1][3-1] + __M[2-1][2-1]# 33 + 22 = 11 # V
__N2[0][1] = -__M[2-1][1-1]#                    - 21 = 12 # V
__N2[0][2] = -__M[3-1][1-1]#                    -  0 = 13 # V
__N2[1][0] = -__M[1-1][2-1]#                    - 12 = 21 # V
__N2[1][1] =  __M[3-1][3-1] + __M[1-1][1-1]# 33 + 11 = 22 # V
__N2[1][2] = -__M[3-1][2-1]#                    -  0 = 23 # V
__N2[2][0] = -__M[1-1][3-1]#                    -  0 = 31 # V
__N2[2][1] = -__M[2-1][3-1]#                    -  0 = 32 # V
__N2[2][2] =  __M[1-1][1-1] + __M[2-1][2-1]# 11 + 22 = 33 # V
#
# Covar. basis (identical to basis {rho, phi, z}) coord. function names (for differentiating metric tens. in Chris. coeff. calculations):
cov_b = [rho, phi, z]
#
# Levi-Civita non-zero pairs relative to 3rd index:
LC_nzPairs = [# i  j  i  j
              [ 1, 2, 2, 1 ],# k = 1
              [ 2, 0, 0, 2 ],# k = 2
              [ 0, 1, 1, 0 ]#  k = 3
             ]# +  +  -  -
#
def prod ( iterable ) :
  return reduce ( operator.mul, iterable, 1 )
pass# DEF prod ( iterable )
#
# Chris. symbols of 2nd kind (I've been too lasy enough to derive them in paper) for cyl. coord.:
#
def Chris_2nd_cyl( i_con, j_cov, k_cov ):
  ChGamma = 0
  for m in range( 0, 3 ):
    ChGamma += g_contravar[i_con][m] * ( diff( g_covar[m][j_cov], cov_b[k_cov] ) + diff( g_covar[m][k_cov], cov_b[j_cov] ) - diff( g_covar[j_cov][k_cov], cov_b[m] ) )
  ##return simplify( sympify( ChGamma / 2 ) )
  return ChGamma / 2
pass# DEF Chris_2nd_cyl ( i_con, j_cov, k_cov )
#
# Chris. symbols of 2nd kind for general curv. coord.:
#
def Chris_2nd_curv( i_con, j_cov, k_cov, g_cov = gCurv_cov, g_con = gCurv_con, g_vars = theta ):
  ChGamma = 0
  for m in range( 0, 3 ):
    ChGamma += (g_con[i_con][m] * ( diff( g_cov[j_cov][m], g_vars[k_cov] )
               + diff( g_cov[k_cov][m], g_vars[j_cov] )
               - diff( g_cov[j_cov][k_cov], g_vars[m] ) ))
  return ChGamma / 2
pass# DEF Chris_2nd_cyl ( i_con, j_cov, k_cov, g_con, g_cov )
#
# Covariant derivatives (I've been too lasy enough to derive them in paper) used to calculate curvilinear cyl. curl of 2nd. ord. covar. tens.:
# def cov_diff_2nd_ord_tens_single_comp( M, i_cov, j_con, m_diff ):# M[i][j] - mixed: covariant i, contravariant j
  # Mdiff = diff( M[i_cov][j_con], cov_b[m_diff] )
  # for k in range(0,3):
    # Mdiff += Chris_2nd_cyl( j_con, k, m_diff ) * M[i_cov][k] - Chris_2nd_cyl( k, i_cov, m_diff ) * M[k][j_con]
  # return simplify( Mdiff )
#
def cov_diff_vect_single___tens( Stens, i_cov, j_con, m_diff ):
  Mdiff = diff( Stens[i_cov][j_con], cov_b[m_diff] )
  for k in range(0,3):
    Mdiff -= Chris_2nd_cyl( k, m_diff, i_cov ) * Stens[k][j_con]
  ##return simplify( sympify( Mdiff ) )
  return Mdiff
#
# Matrix multiplication:
#
def mult1_2nd_ord_tens( M1, M2 ):#RENAME
 MM = [[0 for i in xrange(3)] for i in xrange(3)]
 for i in range(0,3):
  for j in range(0,3):
   for k in range(0,3):
    MM[i][j] += M1[i][k] * M2[k][j]
   pass# FOR
   ##MM[i][j] = simplify( sympify( MM[i][j] ) )
  pass# FOR
 pass# FOR
 return MM
pass# DEF
#
# Convert vect. from physical representation in normal covar. basis {e_rho, e_phi, e_z} to general representation in nonnormalized contravar. basis {b_rho, b_phi, b_z} for further calculations:
def conv_vect_eee2bbb_cyl(Veee):
  U = [0 for i in xrange(3)]
  for i in range(0,3):
    ##U[i] = simplify( sympify( powdenest ( Veee[i] / sqrt( g_covar[i][i] ) ) ) )
    U[i] = powdenest( Veee[i] / sqrt( g_covar[i][i] ), force=True )
  return U
#
# Convert tens. from physical representation in normal covar. basis {e_rho, e_phi, e_z} to general representation in nonnormalized contravar. basis {b_rho, b_phi, b_z} for further calculations:
def conv_tens_eee2bbb_cyl(Meee):
  Mbbb = [[0 for i in xrange(3)] for i in xrange(3)]
  for i in range(0,3):
    for j in range(0,3):
      ##Mbbb[i][j] = simplify( sympify( powdenest( Meee[i][j] / sqrt( g_covar[i][i] * g_covar[j][j] ), force=True ) ) )
      Mbbb[i][j] = powdenest( Meee[i][j] / sqrt( g_covar[i][i] * g_covar[j][j] ), force=True )
  return Mbbb
#
# Convert vect. back to phys. representation in norm. covar. basis {e_rho, e_phi, e_z} from general representation in nonnorm. contravar. basis {b_rho, b_phi, b_z} for further calculations:
def conv_vect_bbb2eee_cyl( Vbbb ):
  U = [0 for i in xrange(3)]
  for i in range(0,3):
    ##U[i] = simplify( sympify( powdenest ( Veee[i] * sqrt( g_covar[i][i] ) ) ) )
    U[i] = powdenest( Vbbb[i] * sqrt( g_covar[i][i] ), force=True )
  return U
#
# Convert tens. back to phys. representation in norm. covar. basis {e_rho, e_phi, e_z} from general representation in nonnorm. contravar. basis {b_rho, b_phi, b_z} for further calculations:
#
def conv_tens_bbb2eee_cyl( Mbbb ):
  Meee = [[0 for i in xrange(3)] for i in xrange(3)]
  for i in range(0,3):
    for j in range(0,3):
      ##Meee[i][j] = simplify( sympify( powdenest( Mbbb[i][j] * sqrt( g_covar[i][i] * g_covar[j][j] ), force=True ) ) )
      Meee[i][j] = powdenest( Mbbb[i][j] * sqrt( g_covar[i][i] * g_covar[j][j] ), force=True )
  return Meee
#
def convert_2nd_order_tens( Mxyz ):
  Mrpz = [[0 for i in xrange(3)] for i in xrange(3)]
  Mrpz = mult1_2nd_ord_tens( rot_left, Mxyz )
  Mrpz = mult1_2nd_ord_tens( Mrpz, rot_right )
  return Mrpz
def prod_scal2matr( s, M ):
  Tresult = [[0 for i in xrange(3)] for i in xrange(3)]
  for i in range(0,3):
    for j in range(0,3):
      Tresult[i][j] = s * M[i][j]
    pass# IF
  pass# IF
  return Tresult
pass# DEF
def summ_tens( T1, T2 ):
  T0 = [[0 for i in xrange(3)] for i in xrange(3)]
  for i in range(0,3):
    for j in range(0,3):
      T0[i][j] = T1[i][j] + T2[i][j]
    pass# DEF
  pass# DEF
  return T0
pass# DEF
def prod_scal2vect( s, V ):
  U = [0 for i in xrange(3)]
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
#
# Grad in general curv. coord.:
#
def grad_curv_con( S, g_cov = gCurv_cov, g_con = gCurv_con, g_vars = theta ):
  U =                                                        (
      [
       g_con[i][i] * diff ( S, g_vars[i] )
       for i in xrange(3)
      ]                                                      )
  return U
pass# DEF grad_curv_con ( V )
#
# Div. in cyl. coord.:
#
def div_cyl( V ):
  ##return simplify( sympify( (1/rho)*diff(rho*V[0],rho)+(1/rho)*diff(V[1],phi)+diff(V[2],z) ) )
  return (1/rho)*diff(rho*V[0],rho)+(1/rho)*diff(V[1],phi)+diff(V[2],z)
# Div. in cyl. coord. expressed in cartes. coord.:
def div_cartes_2_cyl( V ):
  ##return simplify( sympify( diffx_rho_phi(V[0])+diffy_rho_phi(V[1])+diff(V[2],z) ) )
  return diffx_rho_phi(V[0])+diffy_rho_phi(V[1])+diff(V[2],z)
#
# Div. of norm. ( eee ) 2nd. ord. tens. in cyl. coord.:
#
def div_tens_cyl( M ):
  Udiv = [0 for i in xrange(3)]
  ##Udiv[0] = simplify( sympify( diff( M[0][0], rho ) + (M[0][0] / rho) + ((1 / rho) * diff( M[1][0], phi )) + diff( M[2][0], z ) + diff( M[1][1], rho ) ) )
  Udiv[0] = diff( M[0][0], rho ) + (M[0][0] / rho) + ((1 / rho) * diff( M[1][0], phi )) + diff( M[2][0], z ) + diff( M[1][1], rho )
  ##Udiv[1] = simplify( sympify( ((1 / rho) * diff( M[1][1], phi )) + diff( M[0][1] , rho ) + (M[0][1] / rho) + (M[1][0] / rho) + diff( M[2][1], z ) ) )
  Udiv[1] = ((1 / rho) * diff( M[1][1], phi )) + diff( M[0][1] , rho ) + (M[0][1] / rho) + (M[1][0] / rho) + diff( M[2][1], z )
  ##Udiv[2] = simplify( sympify( diff( M[2][2], z ) + diff( M[0][2], rho ) + (M[0][2] / rho) + ((1 / rho) * diff( M[1][2], phi )) ) )
  Udiv[2] = diff( M[2][2], z ) + diff( M[0][2], rho ) + (M[0][2] / rho) + ((1 / rho) * diff( M[1][2], phi ))
  return Udiv
pass# DEF div_tens_cyl ( M )
#
# Div. in general curv. coord.:
#
def div_curv_ort( V_con ):
  P_h = prod( h )
  return P_h * sum( diff( V_con[tsum] * P_h, theta[tsum] ) for tsum in xrange(3) )
pass# DEF div_curv_ort ( V_con )
#
def div_curv_from_con ( V_con, g_cov = gCurv_cov, g_con = gCurv_con, g_vars = theta ) :
  S = sum ( diff ( V_con[ksum], g_vars[ksum] )
            + sum ( Chris_2nd_curv ( ksum, ksum, tsum, g_cov, g_con, g_vars )
                    * V_con[tsum]
                    for tsum in xrange(3)
                  )
            for ksum in xrange(3)
          )
  return S
pass# DEF div_curv_from_con ( V_con, g_con = gCurv_con, g_cov = gCurv_cov )
#
# `( nabla dot ( T dot u ) ) (X) u', e.g.                 \
# curv. column-wise div. from 2nd. ord. contravar. tens.:  
def div_tens_curv_con ( T_concon, g_cov = gCurv_cov, g_con = gCurv_con, g_vars = theta ) :
  U_RES = ( [ div_curv_from_con ( [ T_concon[k][n] for k in xrange(3) ],
                                  g_cov, g_con, g_vars
                                )
              for n in xrange(3)
            ]
          )
  return U_RES
pass# DEF div_tens_curv_con ( T_concon, g_cov = gCurv_cov, g_con = gCurv_con )
#
# Cross prod. in cartes. coord.:
def cross_prod_cartes( U, V ):
  T = [0 for i in xrange(3)]
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
                    [diff(X1,rho),diff(X1,phi),diff(X1,z)],
                    [diff(X2,rho),diff(X2,phi),diff(X2,z)],
                    [diff(X3,rho),diff(X3,phi),diff(X3,z)]
                   ]
#
def det2x2( M1, M2, M3, M4 ):
  ##return simplify( sympify( M1*M4 - M2*M3 ) )
  return M1*M4 - M2*M3
pass# DEF
#
def det3x3( M ):
  return                                                          (
         M[0][0] * det2x2( M[1][1], M[1][2], M[2][1], M[2][2] ) -
         M[0][1] * det2x2( M[1][0], M[1][2], M[2][0], M[2][2] ) +
         M[0][2] * det2x2( M[1][0], M[1][1], M[2][0], M[2][1] )
                                                                  )
pass# DEF
#
def Cmatrix3x3( M ):
  return                                                  (
         [
          [
            det2x2( M[1][1], M[1][2], M[2][1], M[2][2] ),
           -det2x2( M[1][0], M[1][2], M[2][0], M[2][2] ),
            det2x2( M[1][0], M[1][1], M[2][0], M[2][1] )
          ],
          [
           -det2x2( M[0][1], M[0][2], M[2][1], M[2][2] ),
            det2x2( M[0][0], M[0][2], M[2][0], M[2][2] ),
           -det2x2( M[0][0], M[0][1], M[2][0], M[2][1] )
          ],
          [
            det2x2( M[0][1], M[0][2], M[1][1], M[1][2] ),
           -det2x2( M[0][0], M[0][2], M[1][0], M[1][2] ),
            det2x2( M[0][0], M[0][1], M[1][0], M[1][1] )
          ]
         ]
                                                          )
pass# DEF
#
def inverse_matr( M ):
  M_inv = [[0 for i in xrange(3)] for i in xrange(3)]
  M_C = Cmatrix3x3( matr_trans( M ) )
  print '\nSimplifying determinant ...\n'
  M_DET = simplify( powdenest( det3x3( M ), force = True ) )
  pprint ( M_DET )
  for i in range(0,3):
    for j in range(0,3):
      M_inv[i][j] = M_C[i][j] / M_DET
    pass# FOR
  pass# FOR
  return M_inv
pass# DEF
#
Jacob_det = det3x3(Jacob_cyl_cartes)
#
# Cross prod. of two covar. vect. (i.g. their components are covar. and have lower indices, but they are given in contrav. nonnormalized(sic!) basis {b_rho; b_phi; b_z} with roof indices rho, phi, z):
#
def cross_prod_cyl( U, V ):
  T = [0 for i in xrange(3)] # (i,j,k)
  ##T[0] = simplify( sympify( (1 / Jacob_det) * ( U[1] * V[2] - U[2] * V[1] ) ) )# (2,3,1) - (3,2,1) # rho
  T[0] = (1 / Jacob_det) * ( U[1] * V[2] - U[2] * V[1] )# (2,3,1) - (3,2,1) # rho
  ##T[1] = simplify( sympify( (1 / Jacob_det) * ( U[2] * V[0] - U[0] * V[2] ) ) )# (3,1,2) - (1,3,2) # phi
  T[1] = (1 / Jacob_det) * ( U[2] * V[0] - U[0] * V[2] )# (3,1,2) - (1,3,2) # phi
  ##T[2] = simplify( sympify( (1 / Jacob_det) * ( U[0] * V[1] - U[1] * V[0] ) ) )# (1,2,3) - (2,1,3) # z
  T[2] = (1 / Jacob_det) * ( U[0] * V[1] - U[1] * V[0] )# (1,2,3) - (2,1,3) # z
  return T# Contravar. vect. with components having roof indices, but given in covar. basis {b_rho; b_phi; b_z} with lower indices rho, phi, z
pass# DEF cross_prod_cyl ( U, V )
#
# Cross prod. in gen. curv. coord. for covar. vectors:
def cross_prod_curv_con ( U_cov, V_cov, g_cov = gCurv_cov, g_con = gCurv_con,
                          g_scal = gCurv_cov
                        ) :
  T = ( [ sum ( sum ( LeviCivita ( isum+1, jsum+1, k+1 )
                      * U_cov[isum] * V_cov[jsum]
                      / sqrt ( g_scal )
                      for jsum in xrange(3)
                    )
                for isum in xrange(3)
              )
          for k in xrange(3)
        ]
      )
  return T
pass# DEF cross_prod_curv_con ( U_cov, V_cov, g_cov = gCurv_cov,      \
#   #                           g_con = gCurv_con, g_scal = gCurv_cov \
#   #                         )                                        
#
# Cross prod. in gen. curv. coord. for contravar. vectors:
def cross_prod_curv_cov ( U_con, V_con, g_cov = gCurv_cov, g_con = gCurv_con, g_scal = gCurv ) :
  T =                                                        (
      [
       sum(sum(
        LeviCivita ( isum+1, jsum+1, k+1 ) * U_con[isum] * V_con[jsum] * sqrt ( g_scal )
        for jsum in xrange(3)  ) for isum in xrange(3)  )
       for k in xrange(3)
      ]                                                      )
  return T
pass# DEF cross_prod_curv_con ( U_con, V_con, g_cov = gCurv_cov, g_con = gCurv_con )
#
#
def cross_prod_curv_concon_l2r_tens_vect ( T_concon, V_con, g_cov = gCurv_cov, g_con = gCurv_con, g_scal = gCurv ):
  M_concov = [ cross_prod_curv_cov ( T_concon[k], V_con, g_cov, g_con, g_scal )  for k in xrange(3)]
  M_RES = mult1_2nd_ord_tens ( M_concov, g_con )
  return M_RES
pass# DEF cross_prod_curv_con_l2r_tens_vect ( M, V )
#
# Cross prod. spec. case - dot prod. of nabla and 2nd rank tens. on the left side:
def cross_prod_spec_case( M, V ):
  U = [0 for i in xrange(3)]
  ##U[0] = simplify( sympify( (M[1][0]*diff(V[2],x)+M[1][1]*diff(V[2],y)+M[1][2]*diff(V[2],z)-M[2][0]*diff(V[1],x)+M[2][1]*diff(V[1],y)+M[2][2]*diff(V[1],z)) ) )
  U[0] = (M[1][0]*diff(V[2],x)+M[1][1]*diff(V[2],y)+M[1][2]*diff(V[2],z)-M[2][0]*diff(V[1],x)+M[2][1]*diff(V[1],y)+M[2][2]*diff(V[1],z))
  ##U[1] = simplify( sympify( (M[2][0]*diff(V[0],x)+M[2][1]*diff(V[0],y)+M[2][2]*diff(V[0],z)-M[0][0]*diff(V[2],x)+M[0][1]*diff(V[2],y)+M[0][2]*diff(V[2],z)) ) )
  U[1] = (M[2][0]*diff(V[0],x)+M[2][1]*diff(V[0],y)+M[2][2]*diff(V[0],z)-M[0][0]*diff(V[2],x)+M[0][1]*diff(V[2],y)+M[0][2]*diff(V[2],z))
  ##U[2] = simplify( sympify( (M[0][0]*diff(V[1],x)+M[0][1]*diff(V[1],y)+M[0][2]*diff(V[1],z)-M[1][0]*diff(V[0],x)+M[1][1]*diff(V[0],y)+M[1][2]*diff(V[0],z)) ) )
  U[2] = (M[0][0]*diff(V[1],x)+M[0][1]*diff(V[1],y)+M[0][2]*diff(V[1],z)-M[1][0]*diff(V[0],x)+M[1][1]*diff(V[0],y)+M[1][2]*diff(V[0],z))
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
def curl_simple ( X, Y, Z ):
  U = [0 for i in xrange(3)]
  ##U[0] = simplify( sympify( diff(Z,y)-diff(Y,z) ) )
  U[0] = diff(Z,y)-diff(Y,z)
  ##U[1] = simplify( sympify( diff(X,z)-diff(Z,x) ) )
  U[1] = diff(X,z)-diff(Z,x)
  ##U[2] = simplify( sympify( diff(Y,x)-diff(X,y) ) )
  U[2] = diff(Y,x)-diff(X,y)
  return U
pass# DEF curl_simple ( X, Y, Z )
#
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
pass# DEF curl_parallel_l2r_tensor_2nd ( V )
#
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
#
# Curl in cyl. coordinates:
def curl(V):
  U = [0 for i in xrange(3)]
  U[0] = diff(V[2],phi)/rho-diff(V[1],z)
  U[1] = diff(V[0],z)-diff(V[2],rho)
  U[2] = (diff(rho*V[1],rho)-diff(V[0],phi))/rho
  return U
#
# Curvil. cyl. curl of 2nd. ord. covar. tens.:
#def curl_cyl_2ndord_tens ( M ) :
#  Mcurl = [[0 for i in xrange(3)] for i in xrange(3)]
#  for j in range(0,3):# columns
#    for i in range(0,3):# rows
#      Mcurl[i][j] = (1 / rho) * ( cov_diff_2nd_ord_tens_single_comp( M, LC_nzPairs[i][1], j, LC_nzPairs[i][0] ) - cov_diff_2nd_ord_tens_single_comp( M, LC_nzPairs[i][3], j, LC_nzPairs[i][2] ) )
#  return simplify( Mcurl )
#pass# DEF curl_cyl_2ndord_tens ( M )
#
# Curvil. cyl. curl of 2nd. ord. covar. tens.:
def curl_cyl_2ndord_tens___2 ( S ) :
  Mcurl = [[0 for i in xrange(3)] for i in xrange(3)]
  for j in range(0,3):# columns
    for i in range(0,3):# rows
      ##Mcurl[i][j] = simplify( sympify( (1 / rho) * ( cov_diff_vect_single___tens( S, LC_nzPairs[i][1], j, LC_nzPairs[i][0] ) - cov_diff_vect_single___tens( S, LC_nzPairs[i][3], j, LC_nzPairs[i][2] ) ) ) )
      Mcurl[i][j] = (1 / rho) * ( cov_diff_vect_single___tens( S, LC_nzPairs[i][1], j, LC_nzPairs[i][0] ) - cov_diff_vect_single___tens( S, LC_nzPairs[i][3], j, LC_nzPairs[i][2] ) )
  return Mcurl
pass# DEF curl_cyl_2ndord_tens___2 ( S )
#
# Curl in general curv. ort. coord.:
def curl_curv_ort_con( V ):
  P_h = prod( h )
  U =                                     (
      [
       sum(sum(P_h * LeviCivita( isum+1, jsum+1, k+1 )
       * diff( V[jsum] * h[jsum]**2, theta[isum] ) for jsum in xrange(3)) for isum in xrange(3) )
       for k in xrange(3)
      ]
                                          )
  return U
pass# DEF curl_curv_ort_con ( V )
#
# Curl in general curv. coord.:
def curl_curv_con( V_cov, g_cov = gCurv_cov, g_con = gCurv_con, g_scal = gCurv, g_vars = theta ):
  U = (
       [
        sum
        (
         sum
         (
          LeviCivita ( isum+1, jsum+1, k+1 )
          * ( diff( V_cov[jsum], g_vars[isum] )
              - sum ( Chris_2nd_curv ( msum, isum, jsum, g_cov, g_con, g_vars )
                      * V_cov[msum]
                      for msum in xrange(3)
                    )
            )
          / sqrt ( g_scal )
          for jsum in xrange(3)
         )
         for isum in xrange(3)
        )
        for k in xrange(3)
       ]
      )
  return U
pass# DEF curl_curv_con ( V_cov, g_cov = gCurv_cov, g_con = gCurv_con )
#
# `( nabla X ( T dot u ) ) (X) u', e.g.                   \
# curv. column-wise curl from 2nd. ord. contravar. tens.:  
def curl_curv_concon ( T_concon, g_cov = gCurv_cov, g_con = gCurv_con, g_scal = gCurv, g_vars = theta ) :
  Tconcon_RES = [[0 for j in xrange(3)] for i in xrange(3)]
  for i in range ( 0, 3 ) :
    Ccon = curl_curv_con ( [T[j][i] for j in xrange(3)], g_cov, g_con, g_scal, g_vars )
    for j in range ( 0, 3 ) :
      Tconcon_RES[j][i] = Ccon[j]
    pass# FOR j in range ( 0, 3 )
  pass# FOR i in range ( 0, 3 )
  return Tconcon_RES
pass# DEF curl_curv_concon ( T, g_cov, g_con )
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
pass# DEf dot_prod_l2r_tens_vect( M, V )
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
pass# DEF cartes_2_cyl ( V )
#
def cyl_2_cartes(V):
  U = [0 for i in xrange(3)]
  ##U[0] = simplify( sympify( V[0]*cos(V[1]) ) )
  U[0] = V[0]*cos(V[1])
  ##U[1] = simplify( sympify( V[1]*sin(V[1]) ) )
  U[1] = V[1]*sin(V[1])
  ##U[2] = simplify( sympify( V[2] ) )
  U[2] = V[2]
  return U
pass# DEF cyl_2_cartes ( V )
#
# Vect. lapl. in cartes. coord:
#
def lapl(V):
  U = [0 for i in xrange(3)]
  for i in range(0,3):
    ##U[i] = simplify( sympify( diff(diff(V[i],x),x)+diff(diff(V[i],y),y)+diff(diff(V[i],z),z) ) )
    U[i] = diff(diff(V[i],x),x)+diff(diff(V[i],y),y)+diff(diff(V[i],z),z)
  return U
pass# DEF lapl ( V )
#
# Vect. laplacian in cyl. coord:
#
def lapl_cyl(V):
  U = [0 for i in xrange(3)]
  ##U[0] = simplify( sympify( diff(diff(V[0],rho),rho)+(1/rho**2)*diff(diff(V[0],phi),phi)+diff(diff(V[0],z),z)+(1/rho)*diff(V[0],rho)-(2/rho**2)*diff(V[1],phi)-V[0]/rho**2 ) )
  U[0] = diff(diff(V[0],rho),rho)+(1/rho**2)*diff(diff(V[0],phi),phi)+diff(diff(V[0],z),z)+(1/rho)*diff(V[0],rho)-(2/rho**2)*diff(V[1],phi)-V[0]/rho**2
  ##U[1] = simplify( sympify( diff(diff(V[1],rho),rho)+(1/rho**2)*diff(diff(V[1],phi),phi)+diff(diff(V[1],z),z)+(1/rho)*diff(V[1],rho)+(2/rho**2)*diff(V[0],phi)-V[1]/rho**2 ) )
  U[1] = diff(diff(V[1],rho),rho)+(1/rho**2)*diff(diff(V[1],phi),phi)+diff(diff(V[1],z),z)+(1/rho)*diff(V[1],rho)+(2/rho**2)*diff(V[0],phi)-V[1]/rho**2
  ##U[2] = simplify( sympify( diff(diff(V[2],rho),rho)+(1/rho**2)*diff(diff(V[2],phi),phi)+diff(diff(V[2],z),z)+(1/rho)*diff(V[2],rho) ) )
  U[2] = diff(diff(V[2],rho),rho)+(1/rho**2)*diff(diff(V[2],phi),phi)+diff(diff(V[2],z),z)+(1/rho)*diff(V[2],rho)
  return U
pass# DEF lapl_cyl ( V )
#
# Vect lapl. cyl. to cartes. mapping:
#
def lapl_cartes2cyl(V):
  U = [0 for i in xrange(3)]
  for i in range(0,3):
    ##U[i] = simplify( sympify( diffx_rho_phi(diffx_rho_phi(V[i]))+diffy_rho_phi(diffy_rho_phi(V[i]))+diff(diff(V[i],z),z) ) )
    U[i] = diffx_rho_phi(diffx_rho_phi(V[i]))+diffy_rho_phi(diffy_rho_phi(V[i]))+diff(diff(V[i],z),z)
  return U
pass# DEF lapl_cartes2cyl ( V )
#
# Vect. lapl. in gen. curv. coord. (V must have contravar. coord.):
#
def lapl_curv_ort_con( V ):
  P_h = prod( h )
  P_h1 = prod( h1 )
  U =                                                                   (
      [
       sum(
        diff( diff( V[tsum] * P_h, theta[tsum] ) * P_h1, theta[m] ) * h1[m]**2
        for tsum in xrange(3)
       )
       - sum(sum(sum(sum(
          LeviCivita( lsum+1, ksum+1, m+1 ) * LeviCivita( isum+1, jsum+1, ksum+1 )
          * diff
            (
             h[ksum]**2 * diff( V[jsum] * h[jsum]**2, theta[isum] ) * P_h1,
             theta[lsum]
            )
          for jsum in xrange(3) ) for isum in xrange(3) ) for ksum in xrange(3) ) for lsum in xrange(3) )
       * P_h1
       for m in xrange(3)
      ]                                                                 )
  return U
pass# DEF lapl_curv( V )
#
def lapl_curv_con( V_con, g_cov = gCurv_cov, g_con = gCurv_con, g_scal = gCurv, g_vars = theta ):
  U =                                             (
      [
       a-b for a,b in
       zip
       (
        grad_curv_con ( div_curv_from_con ( V_con, g_cov, g_con, g_vars ), g_cov, g_con, g_vars ),
        curl_curv_con ( dot_prod_l2r_tens_vect ( g_cov, curl_curv_con ( dot_prod_l2r_tens_vect ( g_cov, V_con ), g_cov, g_con, g_scal, g_vars ) ), g_cov, g_con, g_scal, g_vars )
       )
      ]                                           )
  return U
pass# DEF lapl_curv ( V )
#
# Matrix transponse:
#
def matr_trans( M ):
  MT = [[0 for i in xrange(3)] for i in xrange(3)]
  for i in range(0,3):
    for j in range(0,3):
      MT[i][j] = M[j][i]
    pass# FOR
  pass# FOR
  return MT
pass# DEF
#
# Substitution of bevel basis:
#
def bevel_basis_subst( auto_calculated ):
  auto_calcreeval = auto_calculated.subs( M11, 1/epsilon )
  auto_calcreeval = auto_calcreeval.subs( M22, 1/epsilon )
  auto_calcreeval = auto_calcreeval.subs( M12, 1/epsilon )
  auto_calcreeval = auto_calcreeval.subs( M33, lambda_11 )
  auto_calcreeval = auto_calcreeval.subs( M23, 0 )
  auto_calcreeval = auto_calcreeval.subs( M13, 0 )
  return auto_calcreeval
pass# DEF bevel_basis_subst
#
# Print matr. by row:
#
def print_matr_by_row_with_subst( M_ARG, f ):
  for i in range(0,3):
    print_row = 0
    for j in range(0,3):
      if ( "subs" in dir(M_ARG[i][j]) ):
        print_row = print_row + u_axis[j] * bevel_basis_subst( M_ARG[i][j] )
      else:
        print_row = print_row + u_axis[j] * M_ARG[i][j]
      pass# IF ( "subs" in dir(M_ARG[i][j]) )
    pass# FOR j in range(0,3)
    printTexExpressionInline( f, sympify( print_row ) )
  pass# FOR i in range(0,3)
pass# DEF print_matr_by_row_with_subst
#
#
#
#
#
#
#
getSimpleTimingData ( 'Precalculation timer' )
#
test = False
#print '\n\n\nDETERMINANT \n\n\n\n'
#pprint(simplify(powdenest(det3x3(M), force=True)))
#pprint( simplifySingleExpressioninMaxima(det3x3(M)))
#
if test:
  Mtest = conv_tens_bbb2eee_cyl                                           (
           curl_cyl_2ndord_tens___2
           (
            index_lower_1st_2ndord_tens_cyl( conv_tens_eee2bbb_cyl( M ) )
           )
                                                                          )
  #Mtest = curl_cyl_2ndord_tens___2                                       (
  #         index_lower_1st_2ndord_tens_cyl( conv_tens_eee2bbb_cyl( M ) )
  #                                                                       )
  #Mtest = index_lower_1st_2ndord_tens_cyl( conv_tens_eee2bbb_cyl( M ) )
  #Mtest = conv_tens_eee2bbb_cyl( M )
  #pprint                                              (
  # simplify
  # (
  #  cov_diff_2nd_ord_tens_single_comp( M, 1, 0, 0 ) -
  #  cov_diff_2nd_ord_tens_single_comp( M, 0, 0, 1 )
  # ),
  # wrap_line = True
  #                                                    )
  #pprint( Mtest, wrap_line=True )
  #print( latex( Mtest, mode = 'equation' ) )
else:
    # Skip main calculations if Fourier plane wave analysis is applied.
    #
    if ( not bFourPlaneWaveOnly ):
      #
      #
      #==============================_____MAIN CALCULATIONS_____===================================
      #
      # 1) k_0**2 * H[0..2] = curl( M[0..2][0..2] * (curl(H[0..2])) )
      #
      # 1.1) Right part
      #
      # 1.1.a) Cartesian coordinates
      #
      # M[0..2][0..2] * (curl(H[0..2])):
      if ( bFourierDomain ):
        M_curlH__ = dot_prod_l2r_tens_vect( M, curl_simple( Hfour[0], Hfour[1], Hfour[2] ) )
      else:
        M_curlH__ = dot_prod_l2r_tens_vect( M, curl_simple( Hcartes[0], Hcartes[1], Hcartes[2] ) )
      #
      # RHS:
      curl_RHS__ = curl_simple( M_curlH__[0], M_curlH__[1], M_curlH__[2] )
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
      ####curl_RHS = curl( f_MH )
      #
      #
      # 1.1.c) Cart. to cyl. coord.
      #
      # Dot prod. of 2nd order tens. and a curl of vect. field, but with spec. case of a curl:
      ##MH_ = dot_prod_l2r_tens_vect(M,curl_parallel_2_cyl(Hxyz))
      #
      # RHS:
      ##curl_RHS_ = curl_parallel_2_cyl(MH_)
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
      _Mcyl = convert_2nd_order_tens( M )
      #_Mcyl = M_dummy
      #
      # Convert M[0..2][0..2] from norm. covar. basis {e_rho; e_phi; e_z} to nonnorm. covar.      \
      # basis {b_rho; b_phi; b_z}:
      _Mbbb = conv_tens_eee2bbb_cyl( _Mcyl )
      #
      # Convert of curl(H[0..2]) from norm. covar. basis {e_rho; e_phi; e_z} to nonnorm. covar.   \
      # basis {b_rho; b_phi; b_z}:
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
      # Convert `M[0..2][0..2] * curl(H[0..2])' back to norm. covar. basis {e_rho; e_phi; e_z}    \
      # from nonnorm. covar. basis {b_rho; b_phi; b_z}:
      _McurlH_eee = conv_vect_bbb2eee_cyl( _M_curlH )
      #
      # RHS `curl( M[0..2][0..2] * curl(H[0..2]) )':
      _curl_RHS = curl( _McurlH_eee )
      #
      #
      # 1.1.e) Curv. coord. tens. transform.:
      #
      # Convert H[0..2] from contravar. to covar. basis:
      __Hcov = dot_prod_l2r_tens_vect ( gCurv_cov, __Hfour )
      #
      # curl(H):
      __curlH = curl_curv_con ( __Hcov )
      #
      # Convert curl ( H[0..2] ) from contravar. to covar. basis:
      __curlHcov = dot_prod_l2r_tens_vect ( gCurv_cov, __curlH )
      #
      # M[0..2][0..2] * curl(H[0..2]):
      __M_curlH = dot_prod_l2r_tens_vect ( __M, __curlHcov )
      #
      # Convert M[0..2][0..2] * curl ( H[0..2] ) from contravar. to covar. basis:
      __McurlH_cov = dot_prod_l2r_tens_vect ( gCurv_cov, __M_curlH )
      #
      # RHS `curl( M[0..2][0..2] * curl(H[0..2]) )':
      __curl_RHS = curl_curv_con ( __McurlH_cov )
      #
      #
      #
      # 1.2) Left part
      #
      # 1.2.a) Cartesian coordinates
      #
      # LHS:
      if ( bFourierDomain ):
        curl_LHS__ = prod_scal2vect( k_0**2, Hfour )
      else:
        curl_LHS__ = prod_scal2vect( k_0**2, Hcartes )
      #
      #
      ##### 1.2.b) Cyl. coord.
      #####
      ##### LHS:
      ####curl_LHS = prod_scal2vect( k_0**2, H )
      #####
      #####
      # 1.2.c) Cart. to cyl. coord.
      #
      # LHS:
      ##curl_LHS_ = prod_scal2vect(k_0**2,Hxyz)
      #
      #
      # 1.2.d) Cyl. coord. tens. transform.:
      #
      # LHS:
      _curl_LHS = prod_scal2vect( k_0**2, H )
      #
      #
      # 1.2.e) Curv. coord. tens. transform.:
      #
      # LHS:
      __curl_LHS = curl_LHS__
      #
      #
      #
      #
      # 2) ( N2[0..2][0..2] * nabla**2 + k_0**2 ) * H[0..2] =                 \
      #    (curl(M[0..2][0..2])) * (curl(H[0..2]))          +                 \
      #    ( N1[0..2][0..2] * nabla ) [X] (curl(H[0..2]))
      #
      # 2.1) Right part
      #
      # 2.1.a) Cartesian coordinates
      #
      # (curl(M[0..2][0..2])) * (curl(H[0..2])):
      if ( bFourierDomain ):
        curlM_curlH__ =                                              (
                        dot_prod_l2r_tens_vect
                        (
                         curl_parallel_l2r_tensor_2nd( M ),
                         curl_simple( Hfour[0], Hfour[1], Hfour[2] )
                        )
                                                                     )
      else:
        curlM_curlH__ =                                                    (
                        dot_prod_l2r_tens_vect
                        (
                         curl_parallel_l2r_tensor_2nd( M ),
                         curl_simple( Hcartes[0], Hcartes[1], Hcartes[2] )
                        )
                                                                           )
      #
      # ( N1[0..2][0..2] * nabla ) [X] (curl(H[0..2])):
      if ( bFourierDomain ):
        N1nabla_X_curlH__ = cross_prod_spec_case( N1, curl_simple( Hfour[0], Hfour[1], Hfour[2] ) )
      else:
        N1nabla_X_curlH__ =                                                    (
                            cross_prod_spec_case
                            (
                             N1,
                             curl_simple( Hcartes[0], Hcartes[1], Hcartes[2] )
                            )
                                                                               )
      #
      # RHS:
      lapl_RHS__ = [a+b for a,b in zip( N1nabla_X_curlH__, curlM_curlH__ )]
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
      ####lapl_RHS = [a+b-c for a,b,c in zip(curlM_curlH,divN1_x_curlH,div_N1curlH)]
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
      ##lapl_RHS_ = [a+b-c for a,b,c in zip(curlM_curlH_,divN1_x_curlH_,div_N1curlH_)]
      #
      #
      # 2.1.d) Cyl. coord. tens. transform.:
      ##########
      # M[0..2][0..2] in cyl. coord.:
      # _Mcyl = convert_2nd_order_tens(M)# DUPLICATE
      #
      # Convertion of M[0..2][0..2] from normalized covariant {e_rho; e_phi; e_z} to              \
      # nonnormalized covariant basis {b_rho; b_phi; b_z}:
      # _Mbbb = conv_tens_eee2bbb_cyl(_Mcyl)# DUPLICATE
      #
      # M[0..2][0..2] in dyadic representation `b_covar (X) b_contravar':
      _Mcovcon = index_lower_1st_2ndord_tens_cyl( _Mbbb )
      #
      # curl(H[0..2]):
      # _curlH = curl(H)# DUPLICATE
      #
      # Convert of curl(H[0..2]) from norm. covar. basis {e_rho; e_phi; e_z} to nonnorm. covar.   \
      # basis {b_rho; b_phi; b_z}:
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
      _N1cyl = convert_2nd_order_tens( N1 )
      #_N1cyl = N1_dummy
      #print "7\n\n"
      #pprint( _N1cyl, wrap_line=True )
      #
      # div(N1[0..2][0..2]):
      _divN1 = div_tens_cyl( _N1cyl )
      #
      # Convert `div(N1[0..2][0..2])' from norm. covar. basis {e_rho; e_phi; e_z} to nonnorm.     \
      # covar. basis {b_rho; b_phi; b_z}:
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
      # Convert `N1[0..2][0..2]' from norm. covar. basis {e_rho; e_phi; e_z} to nonnorm. covar.   \
      # basis {b_rho; b_phi; b_z}:
      _N1bbb = conv_tens_eee2bbb_cyl( _N1cyl )
      #print "6\n\n"
      #pprint( _N1bbb, wrap_line=True )
      #
      # Convert `N1[0..2][0..2]' from dyad. represent. of type `b_covar (X) b_covar' to dyad.     \
      # represent. of type `b_covar (X) b_contravar':
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
      _lapl_RHS = [a+b-c for a,b,c in zip( _curlMcurlH_eee, _divN1xcurlH_eee, _div_N1xcurlH )]
      #
      #
      # 2.1.e) Curv. coord. tens. transform.:
      #
      # (curl(M[0..2][0..2])) * (curl(H[0..2])):
      __curlM_curlH = dot_prod_l2r_tens_vect( __M, __curlHcov )
      #
      # div(N1[0..2][0..2]):
      __divN1 = div_tens_curv_con ( __N1 )
      #
      # div(N1[0..2][0..2]) [X] curl(H[0..2]):
      __divN1_x_curlH = cross_prod_curv_cov ( __divN1, __curlH )
      #
      # N1[0..2][0..2] [X] curl(H[0..2]):
      __N1_x_curlH = cross_prod_curv_concon_l2r_tens_vect( __N1, __curlH )
      #
      # div(N1[0..2][0..2] [X] curl(H[0..2])):
      __div_N1xcurlH = div_tens_curv_con ( __N1_x_curlH )
      #
      # Lapl. RHS:
      __lapl_RHS =                                                        (
                   [
                    a+b-c for a,b,c in
                    zip( __curlM_curlH, __divN1_x_curlH, __div_N1xcurlH )
                   ]                                                      )
      #
      #
      #
      # 2.2) Left part
      #
      # 2.2.a) Cartesian coordinates
      #
      # N2[0..2][0..2] * nabla**2 * H[0..2]:
      if ( bFourierDomain ):
        N2_nabla2H__ = dot_prod_l2r_tens_vect( N2, lapl( Hfour ) )
      else:
        N2_nabla2H__ = dot_prod_l2r_tens_vect( N2, lapl( Hcartes ) )
      #
      # k_0**2 * H[0..2]:
      if ( bFourierDomain ):
        k2_H__ = prod_scal2vect( k_0**2, Hfour )
      else:
        k2_H__ = prod_scal2vect( k_0**2, Hcartes )
      # LHS:
      lapl_LHS__ = [a+b for a,b in zip( N2_nabla2H__, k2_H__ )]
      lapl_LHS2RHS__ = N2_nabla2H__
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
      ##### k_0**2 * H[0..2]:
      ####k2_H = prod_scal2vect( k_0**2, H )
      ##### Sum:
      ####lapl_LHS = [a+b for a,b in zip( N2nabla2Hcyl, k2_H )]
      #####
      # 2.2.c) Cart. to cyl. coord.
      #
      # N2[0..2][0..2] * nabla**2 * H[0..2]:
      ##N2_nabla2H_ = dot_prod_l2r_tens_vect(N2,lapl_cartes2cyl(Hxyz))
      # k_0**2 * H[0..2] in cyl. coord.:
      ##k2_H_ = prod_scal2vect(k_0**2,Hxyz)
      # Sum:
      ##lapl_LHS_ = [a+b for a,b in zip( N2_nabla2H_, k2_H_ )]
      #
      #
      # 2.2.d) Cyl. coord. tens. transform.
      #
      # lapl(H[0..2]):
      _laplH = lapl_cyl( H )
      #
      # Convert lapl(H[0..2]) from norm. covar. basis {e_rho; e_phi; e_z} to nonnorm. covar. basis\
      # {b_rho; b_phi; b_z}:
      _laplH_bbb = conv_vect_eee2bbb_cyl( _laplH )
      #
      # Convert lapl(H[0..2]) from covar. basis to contravar. basis:
      _laplH_cov = index_lower_vect_cyl( _laplH_bbb )
      #
      _N2cyl = convert_2nd_order_tens( N2 )
      #_N2cyl = N2_dummy
      #
      # Convert N2 from norm. covar. basis {e_rho; e_phi; e_z} to nonnorm. covar. basis           \
      # {b_rho; b_phi; b_z}:
      _N2bbb = conv_tens_eee2bbb_cyl( _N2cyl )
      #
      # N2[0..2][0..2] * lapl(H[0..2]):
      _N2_laplH = dot_prod_l2r_tens_vect( _N2bbb, _laplH_cov )
      # Convert N2[0..2][0..2] * lapl(H[0..2]) back to norm. covar basis {e_rho; e_phi; e_z} from \
      # nonnorm. covar. basis {b_rho; b_phi; b_z}:
      _N2laplH_eee = conv_vect_bbb2eee_cyl( _N2_laplH )
      ##########
      # k_0**2 * H[0..2]:
      _k2_H = prod_scal2vect( k_0**2, H )
      ##########
      _lapl_LHS = [a+b for a,b in zip( _N2laplH_eee, _k2_H )]
      _lapl_LHS2RHS = _N2laplH_eee
      #
      #
      # 2.2.e) Curv. coord. tens. transform.:
      #
      # lapl(H[0..2]):
      __laplH = lapl_curv_con ( __Hfour )
      #
      # N2[0..2][0..2] * lapl(H[0..2]):
      __N2_laplH = dot_prod_l2r_tens_vect ( __N2, dot_prod_l2r_tens_vect(gCurv_cov, __laplH) )
      #
      # k_0**2 * H[0..2]:
      __k2_H = prod_scal2vect ( k_0**2, __Hfour )
      #
      # Lapl. LHS:
      __lapl_LHS = [a+b for a,b in zip( __N2_laplH, __k2_H )]
      __lapl_LHS2RHS = __N2_laplH
      #
      #
      #
      #
      # 3 ) delta = curl ( curl ( E[0..2] ) )        \
      #             - k**2 * M[0..2][0..2] * E[0..2]  
      #
      # 3.1) Right part
      #
      # 3.1.d) Cyl. coord. tens. transform.:
      ##########
      # M[0..2][0..2] in cyl. coord.:
      # _Mcyl = convert_2nd_order_tens(M)# DUPLICATE
      #
      # Convert M[0..2][0..2] from norm. covar. basis {e_rho; e_phi; e_z} \
      # to nonnorm. covar. basis {b_rho; b_phi; b_z}:                      
      # _Mbbb = conv_tens_eee2bbb_cyl( _Mcyl )# DUPLICATE
      #
      # Convert E[0..2] from covar. to contravar. basis:
      _Ecov = dot_prod_l2r_tens_vect ( g_covar, Efour )
      #
      # M[0..2][0..2] * E[0..2] in cyl. coord.:
      _MdotE = dot_prod_l2r_tens_vect ( _Mbbb, _Ecov )
      #
      # k**2 * M[0..2][0..2] * E[0..2] in cyl. coord.:
      _k_MdotE = prod_scal2vect ( k**2, _MdotE )
      #
      ##########
      # curl ( E[0..2] ):
      _curlE = curl_curv_con ( _Ecov, g_covar, g_contravar, g_scalar, rhophiz)
      #
      # curl ( E[0..2] ) from covar. to contravar. basis:
      _curlE_cov = dot_prod_l2r_tens_vect ( g_covar, _curlE )
      #
      # curl ( curl ( E[0..2] ) ):
      _curlcurlE = curl_curv_con ( _curlE_cov, g_covar, g_contravar, g_scalar, rhophiz )
      #
      ##########
      # LHS:
      _E_LHS = [a-b for a,b in zip ( _curlcurlE, _k_MdotE )]
      #
      #
      #
      #
      #
      #======================^^^^^MAIN CALCULATIONS^^^^^=======================
      #
      end_timestamp = timeit.default_timer()
      print 'Main calculations timer: '
      print end_timestamp - bgn_timestamp
      print '\n\n'
      bgn_timestamp = timeit.default_timer()
      #
      # The result of comparision between original expression and one with the laplacian:
      #
      #genform = numer( simplify( curl_LHS_[0] - curl_RHS_[0] - lapl_LHS_[0] + lapl_RHS_[0] ) )
      #
      genforms_LHS2RHS = [0 for i in xrange(6)]
      genforms_LHS2RHS_diff = [0 for i in xrange(3)]
      if ( bVarSeparate ):
        #
        # Using precise expr. with tens. curvilinear representation for var. separation.
        # Note that simplification has been moved to final calculations.
        #
        genform0 = sympify( _lapl_LHS[0] - _lapl_RHS[0] )# Rho
        genform1 = sympify( _lapl_LHS[1] - _lapl_RHS[1] )# Phi
        genform2 = sympify( _lapl_LHS[2] - _lapl_RHS[2] )# Z
        #
        for i in range(0,3):
          genforms_LHS2RHS[i] =  _lapl_RHS[i] - _lapl_LHS2RHS[i]
          genforms_LHS2RHS[i+3] = _curl_RHS[i]
        pass# FOR
        #
        genformA = sympify( _curl_LHS[0] - _curl_RHS[0] )# Rho
        genformB = sympify( _curl_LHS[1] - _curl_RHS[1] )# Phi
        genformC = sympify( _curl_LHS[2] - _curl_RHS[2] )# Z
      else:
        if ( bFourierDomain or bCartesian and not bCurvilinear ):
          #
          # Using cartes. coord. if Four. transform is specified.
          #
          genform0 = simplify( sympify( lapl_LHS__[0] - lapl_RHS__[0] ) )
          genform1 = simplify( sympify( lapl_LHS__[1] - lapl_RHS__[1] ) )
          genform2 = simplify( sympify( lapl_LHS__[2] - lapl_RHS__[2] ) )
          #
          for i in range(0,3):
            genforms_LHS2RHS[i] = lapl_RHS__[i] - lapl_LHS2RHS__[i]
            genforms_LHS2RHS[i+3] = curl_RHS__[i]
          pass# FOR
          #
          genformA = simplify( sympify( curl_LHS__[0] - curl_RHS__[0] ) )
          genformB = simplify( sympify( curl_LHS__[1] - curl_RHS__[1] ) )
          genformC = simplify( sympify( curl_LHS__[2] - curl_RHS__[2] ) )
        else:
          if ( bFourierDomain and bCurvilinear ):
            genform0 = simplify( sympify( __lapl_RHS[0] - __lapl_RHS[0] ) )
            genform1 = simplify( sympify( __lapl_LHS[1] - __lapl_RHS[1] ) )
            genform2 = simplify( sympify( __lapl_LHS[2] - __lapl_RHS[2] ) )
            #
            for i in range(0,3):
              genforms_LHS2RHS[i] = __lapl_RHS[i] - __lapl_LHS2RHS[i]
              genforms_LHS2RHS[i+3] = __curl_RHS[i]
            pass# FOR
            #
            genformA = simplify( sympify( __curl_LHS[0] - __curl_RHS[0] ) )
            genformB = simplify( sympify( __curl_LHS[1] - __curl_RHS[1] ) )
            genformC = simplify( sympify( __curl_LHS[2] - __curl_RHS[2] ) )
          else:
            genform0 = simplify( sympify( _lapl_LHS[0] - _lapl_RHS[0] ) )
            genform1 = simplify( sympify( _lapl_LHS[1] - _lapl_RHS[1] ) )
            genform2 = simplify( sympify( _lapl_LHS[2] - _lapl_RHS[2] ) )
            #
            for i in range(0,3):
              genforms_LHS2RHS[i] = _lapl_RHS[i] - _lapl_LHS2RHS[i]
              genforms_LHS2RHS[i+3] = _curl_RHS[i]
            pass# FOR
            #
            genformA = simplify( sympify( _curl_LHS[0] - _curl_RHS[0] ) )
            genformB = simplify( sympify( _curl_LHS[1] - _curl_RHS[1] ) )
            genformC = simplify( sympify( _curl_LHS[2] - _curl_RHS[2] ) )
          pass# IF ( bFourierDomain and bCurvilinear )
        pass# IF
      pass# IF
      for i in range(0,3):
        #
        # I should also try sum of two parts:
        #
        genforms_LHS2RHS_diff[i] = simplify( genforms_LHS2RHS[i] - genforms_LHS2RHS[i+3] )
      pass# FOR
      #
      end_timestamp = timeit.default_timer()
      print 'Vector simplification timer: '
      print end_timestamp - bgn_timestamp
      print '\n\n'
      bgn_timestamp = timeit.default_timer()
      #
      genforms = [genform0, genform1, genform2, genformA, genformB, genformC]
      Hrhophiz = [Hrho, Hphi, Hz]
      rhophiz_text = ["rho", "phi", "z"]
      Hrhophiz_text = ["Hrho", "Hphi", "Hz"]
      #
      #########################################################################
      #               Fourier method of variable separation                   #
      #########################################################################
      #
      if ( bVarSeparate ):
        f_Hrho_R = Function( 'Hrho_R' )
        f_Hrho_P = Function( 'Hrho_P' )
        f_Hrho_T = Function( 'Hrho_T' )
        f_Hrho_Z = Function( 'Hrho_Z' )
        f_Hphi_R = Function( 'Hphi_R' )
        f_Hphi_P = Function( 'Hphi_P' )
        f_Hphi_T = Function( 'Hphi_T' )
        f_Hphi_Z = Function( 'Hphi_Z' )
        f_Hz_R   = Function( 'Hz_R'   )
        f_Hz_P   = Function( 'Hz_P'   )
        f_Hz_T   = Function( 'Hz_T'   )
        f_Hz_Z   = Function( 'Hz_Z'   )
        #
        Hrho_R = f_Hrho_R( rho      )
        Hrho_P = f_Hrho_P( phi      )
        Hrho_T = f_Hrho_T( rho, phi )
        Hrho_Z = f_Hrho_Z(   z      )
        Hphi_R = f_Hphi_R( rho      )
        Hphi_P = f_Hphi_P( phi      )
        Hphi_T = f_Hphi_T( rho, phi )
        Hphi_Z = f_Hphi_Z(   z      )
        Hz_R   = f_Hz_R(   rho      )
        Hz_P   = f_Hz_P(   phi      )
        Hz_T   = f_Hz_T(   rho, phi )
        Hz_Z   = f_Hz_Z(     z      )
        #
        Hrhophiz_SV1 = [Hrho_T, Hrho_Z, Hphi_T, Hphi_Z, Hz_T, Hz_Z]
        Hrhophiz_SV2 = [Hrho_R, Hrho_P, Hphi_R, Hphi_P, Hz_R, Hz_P]
        Hrhophiz_T   = [Hrho_T,         Hphi_T,         Hz_T      ]
        #
        for i in range(0,3):
          for j in range(0,3):
            genforms[i] =                                     (
             genforms[i].subs
                         (
                          Hrhophiz[j],
                          Hrhophiz_SV1[j] * Hrhophiz_SV1[j+1]
                         )
                                                              )
          pass# FOR
        pass# FOR
        for i in range(0,3):
          for j in range(0,3):
            for k in range(0,3):
              genforms[i] =                                                         (
               genforms[i].subs
                           (
                            diff( Hrhophiz[j], rhophiz[k] ),
                            diff( Hrhophiz_SV1[j] * Hrhophiz_SV1[j+1], rhophiz[k] )
                           )
                                                                                    )
            pass# FOR
          pass# FOR
        pass# FOR
        for i in range(0,3):
          for j in range(0,3):
            for k in range(0,3):
              for l in range(0,3):
                genforms[i] =                                                           (
                 genforms[i].subs
                             (
                              diff
                              (
                               diff( Hrhophiz[j], rhophiz[k] ),
                               rhophiz[l]
                              ),
                              diff
                              (
                               diff( Hrhophiz_SV1[j] * Hrhophiz_SV1[j+1], rhophiz[k] ),
                               rhophiz[l]
                              )
                             )
                                                                                        )
              pass# FOR
            pass# FOR
          pass# FOR
        pass# FOR
        for i in range(0,3):
          genforms[i] = simplify( genforms[i] )
        pass# FOR
        #
        #for i in range(0,3):
        #  for j in range(0,3):
        #    genforms[i] = genforms[i].subs( Hrhophiz_T[j], Hrhophiz_SV2[j] * Hrhophiz_SV2[j+1] )
        #  genforms[i] = simplify( genforms[i] )
      pass# IF ( bVarSeparate )
      #
      #
      #
      end_timestamp = timeit.default_timer()
      print 'variable separation timer:'
      print end_timestamp - bgn_timestamp
      print '\n\n'
      bgn_timestamp = timeit.default_timer()
      #
      #########################################################################
      #^^^^^          Fourier method of variable separation              ^^^^^#
      #########################################################################
      #
      #
      # Substitute symbolic equations instead of partial derivatives.
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
      # Other substitutions:
      #
      if ( not bDoNotSubstitute ):
        #
        # Dummy matricies substitution:
        #
        #for i in range(0,3):
        #  for j in range(0,3):
        #    for k in range(0,3):
        #      genforms[i] = genforms[i].subs( M_dummy[j][k], _Mcyl0[j][k] )
        #      genforms[i] = genforms[i].subs( N1_dummy[j][k], _N1cyl0[j][k] )
        #      genforms[i] = genforms[i].subs( N2_dummy[j][k], _N2cyl0[j][k] )
        #  genforms[i] = simplify( genforms[i] )
        #
        # Eigenvalue substitution:
        #
        for i in range(0,3):
          genforms[i] = genforms[i].subs( Gamma, Gamma_s )
        pass# FOR
        searchTerm_Gamma = sympify( simplify( Gamma_s) )
      else:
        searchTerm_Gamma = sympify( simplify( Gamma ) )
        searchTerm_F = [0 for i in xrange ( 3 )]
        for i in range ( 0, 3 ) :
          searchTerm_F[i] = sympify ( simplify ( F[i] ) )
        pass# FOR i in range ( 0, 3 )
      pass# IF
      searchTerm_Gamma2 = sympify( simplify( searchTerm_Gamma**2 ) )
      #
      end_timestamp = timeit.default_timer()
      print 'substitution timer:'
      print end_timestamp - bgn_timestamp
      print '\n\n'
      bgn_timestamp = timeit.default_timer()
      #
      # Output:
      #
      # Collection of partial derivative symbolic aliases:
      #
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
      # Test output:
      #
      #pprint( genform0, wrap_line=True )
      #print( "\n\n" )
      #pprint( genform1, wrap_line=True )
      #print( "\n\n" )
      #pprint( genform2, wrap_line=True )
      #print( "\n\n" )
      #
      #genform_s = simplify( sympify( genform_collect_Hrhodphi2[S.One] ) )
      #
      if ( ( bDoNotCountDerivatives == False ) and ( bFourierDomain == False ) ):
        for l in range(0,3):
          print "\n"
          print "\n"
          print "\n"
          print "\n"
          print "\n"
          for i in range(0,3):
            for j in range(0,3):
              print "%s / %s :: " % (Hrhophiz_text[i], rhophiz_text[j])
              if ( bVarSeparate ):
                pprint( genforms[l].has( diff( Hrhophiz_SV1[i], rhophiz[j] ) ) )
                pprint( genforms[l+3].has( diff( Hrhophiz_SV1[i], rhophiz[j] ) ) )
              else:
                pprint( genforms[l].has( diff( Hrhophiz[i], rhophiz[j] ) ) )
                pprint( genforms[l+3].has( diff( Hrhophiz[i], rhophiz[j] ) ) )
              print "\n"
          for i in range(0,3):
            for j in range(0,3):
              for k in range(0,3):
                print "%s / %s %s :: " % (Hrhophiz_text[i], rhophiz_text[j], rhophiz_text[k])
                if ( bVarSeparate ):
                  pprint                                                                      (
                   genforms[l].has( diff( diff( Hrhophiz_SV1[i], rhophiz[j] ), rhophiz[k] ) )
                                                                                              )
                  pprint                                                                        (
                   genforms[l+3].has( diff( diff( Hrhophiz_SV1[i], rhophiz[j] ), rhophiz[k] ) )
                                                                                                )
                else:
                  pprint( genforms[l].has( diff( diff( Hrhophiz[i], rhophiz[j] ), rhophiz[k] ) ) )
                  pprint                                                                    (
                   genforms[l+3].has( diff( diff( Hrhophiz[i], rhophiz[j] ), rhophiz[k] ) )
                                                                                            )
                print "\n"
      pass# IF
      #
      if ( bFourierDomain and bRemovePlainWaveExponent ):
        x, y, z = [S.Zero for idx in xrange(3)]
        print "\nX, Y, Z set to zero.\n\n"
      pass# IF ( bFourierDomain and bRemovePlainWaveExponent )
      #
      if ( bFourierDomain == True ) :
        genform_collects = [0 for i in xrange(18)]
        genform_Gamma_collects = [0 for i in xrange(12)]
        for i in range(0,3):
          #
          # Sum of Gamma multipliers in diff. between expr. with lapl. and curl:
          #
          genform_Gamma_collect =                                                        (
                                  collect
                                  (
                                   expand( genforms_LHS2RHS_diff[i] ), searchTerm_Gamma,
                                   evaluate = False, exact = False
                                  )
                                                                                         )
          #
          # Same for expr. with lapl.:
          #
          eigenexp_lapl_Gamma_collect =                                           (
                                        collect
                                        (
                                         expand( genforms[i] ), searchTerm_Gamma,
                                         evaluate = False, exact = False
                                        )
                                                                                  )
          #
          # Same for expr. with curl:
          #
          eigenexp_curl_Gamma_collect =                                             (
                                        collect
                                        (
                                         expand( genforms[i+3] ), searchTerm_Gamma,
                                         evaluate = False, exact = False
                                        )
                                                                                    )
          #
          # First order terms:
          #
          if S.One in genform_Gamma_collect:
            for j in range(0,3):
              genform_F_collect =                                                 (
                                  collect
                                  (
                                   genform_Gamma_collect[S.One], searchTerm_F[j],
                                   evaluate = False, exact = False
                                  )
                                                                                  )
              if searchTerm_F[j] in genform_F_collect:
                auto_differenceEqTens[i][j] =                                       (
                                              simplify
                                              (
                                               removePWE(genform_F_collect[searchTerm_F[j]])
                                              )
                                                                                    )
                #pprint( auto_differenceEqTens[i][j] )
              else:
                auto_differenceEqTens[i][j] = 0
              pass# IF
              #
              if S.One in genform_F_collect:
                genform_F_collectWithin = genform_F_collect
                for kIdx in range(0,3):
                  if kIdx != j:
                    if S.One in genform_F_collectWithin:
                      genform_F_collectWithin =                                  (
                                                collect
                                                (
                                                 genform_F_collectWithin[S.One],
                                                 searchTerm_F[kIdx],
                                                 evaluate = False, exact = False
                                                )
                                                                                 )
                    pass# IF
                  pass# IF
                pass# FOR
                if S.One in genform_F_collectWithin:
                  print                                                                       (
                        '\n\nIn F[' + repr(i) + '] collection. There is S.One (zero-order)' +
                        ' F[x, y, z] term in zero-order term from first-order Gamma' +
                        ' collection. It contains: \n'
                                                                                              )
                  pprint( genform_F_collectWithin[S.One] )
                  print '\n\n'
                pass# IF
              pass# IF
            pass# FOR
          pass# IF
          #
          # Same for expr. with lapl.:
          #
          if S.One in eigenexp_lapl_Gamma_collect:
            for j in range(0,3):
              eigenexp_lapl_F_collect =                                                       (
                                        collect
                                        (
                                         eigenexp_lapl_Gamma_collect[S.One], searchTerm_F[j],
                                         evaluate = False, exact = False
                                        )
                                                                                              )
              if searchTerm_F[j] in eigenexp_lapl_F_collect:
                auto_EigenFTLaplPropagator[i][j] =                                         (
                                               simplify
                                               (
                                                removePWE(eigenexp_lapl_F_collect[searchTerm_F[j]])
                                               )
                                                                                           )
                #pprint( auto_differenceEqTens[i][j] )
              else:
                auto_EigenFTLaplPropagator[i][j] = 0
              pass# IF
            pass# FOR
          pass# IF
          #
          # Same for expr. with curl:
          #
          if S.One in eigenexp_curl_Gamma_collect:
            for j in range(0,3):
              eigenexp_curl_F_collect =                                                       (
                                        collect
                                        (
                                         eigenexp_curl_Gamma_collect[S.One], searchTerm_F[j],
                                         evaluate = False, exact = False
                                        )
                                                                                              )
              if searchTerm_F[j] in eigenexp_curl_F_collect:
                auto_EigenFTCurlPropagator[i][j] =                                         (
                                               simplify
                                               (
                                                removePWE(eigenexp_curl_F_collect[searchTerm_F[j]])
                                               )
                                                                                           )
                #pprint( auto_differenceEqTens[i][j] )
              else:
                auto_EigenFTCurlPropagator[i][j] = 0
              pass# IF
            pass# FOR
          pass# IF
          #
          # Same for expr. for electric field E:
          #
          for j in range(0,3):
            greensfunction_E_collect =                 \
             collect ( expand ( _E_LHS[i] ), searchTerm_F[j],     \
                       evaluate = False, exact = False \
                     )                                  
            if ( searchTerm_F[j] in greensfunction_E_collect ) :
              auto_greensfunction_E[i][j] = (
               simplify
               (
                removePWE
                (
                 greensfunction_E_collect[searchTerm_F[j]]
                )
               )                                           )
            else :
              auto_greensfunction_E[i][j] = 0
            pass# IF ( searchTerm_F[j] in greensfunction_E_collect )
          pass# FOR j in range(0,3)
          #
          # Second order terms:
          #
          if searchTerm_Gamma in genform_Gamma_collect:
            for j in range(0,3):
              genform_F_collect =                                                 (
                                  collect
                                  (
                                   genform_Gamma_collect[searchTerm_Gamma], searchTerm_F[j],
                                   evaluate = False, exact = False
                                  )
                                                                                  )
              if searchTerm_F[j] in genform_F_collect:
                auto_differenceEqGamma1Tens[i][j] =                                       (
                                                    simplify
                                                    (
                                                     removePWE(genform_F_collect[searchTerm_F[j]])
                                                    )
                                                                                          )
                #pprint( auto_differenceEqGamma1Tens[i][j] )
              else:
                auto_differenceEqGamma1Tens[i][j] = 0
              pass# IF
              if S.One in genform_F_collect:
                genform_F_collectWithin = genform_F_collect
                for kIdx in range(0,3):
                  if kIdx != j:
                    if S.One in genform_F_collectWithin:
                      genform_F_collectWithin =                                  (
                                                collect
                                                (
                                                 genform_F_collectWithin[S.One],
                                                 searchTerm_F[kIdx],
                                                 evaluate = False, exact = False
                                                )
                                                                                 )
                    pass# IF
                  pass# IF
                pass# FOR
                if S.One in genform_F_collectWithin:
                  print                                                                       (
                        '\n\nIn F[' + repr(i) + '] collection. There is S.One (zero-order)' +
                        ' F[x, y, z] term in first-order term from first-order Gamma' +
                        ' collection. It contains: \n'
                                                                                              )
                  pprint( genform_F_collectWithin[S.One] )
                  print '\n\n'
                pass# IF
              pass# IF
            pass# FOR
          pass# IF
          #
          # Same for expr. with lapl.:
          #
          if searchTerm_Gamma in eigenexp_lapl_Gamma_collect:
            for j in range(0,3):
              eigenexp_lapl_F_collect =                                                 (
                                        collect
                                        (
                                         eigenexp_lapl_Gamma_collect[searchTerm_Gamma],
                                         searchTerm_F[j],
                                         evaluate = False, exact = False
                                        )
                                                                                        )
              if searchTerm_F[j] in eigenexp_lapl_F_collect:
                auto_EigenFTLapltens1[i][j] =                                             (
                                          simplify
                                          (
                                           removePWE(eigenexp_lapl_F_collect[searchTerm_F[j]])
                                          )
                                                                                      )
                #pprint( auto_differenceEqTens[i][j] )
              else:
                auto_EigenFTLapltens1[i][j] = 0
              pass# IF
            pass# FOR
          pass# IF
          #
          # Same for expr. with curl:
          #
          if searchTerm_Gamma in eigenexp_curl_Gamma_collect:
            for j in range(0,3):
              eigenexp_curl_F_collect =                                                 (
                                        collect
                                        (
                                         eigenexp_curl_Gamma_collect[searchTerm_Gamma],
                                         searchTerm_F[j],
                                         evaluate = False, exact = False
                                        )
                                                                                        )
              if searchTerm_F[j] in eigenexp_curl_F_collect:
                auto_EigenFTCurltens1[i][j] =                                              (
                                               simplify
                                               (
                                                removePWE(eigenexp_curl_F_collect[searchTerm_F[j]])
                                               )
                                                                                           )
                #pprint( auto_differenceEqTens[i][j] )
              else:
                auto_EigenFTCurltens1[i][j] = 0
              pass# IF
            pass# FOR
          pass# IF
          #
          # Third order terms:
          #
          if searchTerm_Gamma**2 in genform_Gamma_collect:
            for j in range(0,3):
              genform_F_collect =                                                 (
                                  collect
                                  (
                                   genform_Gamma_collect[searchTerm_Gamma**2], searchTerm_F[j],
                                   evaluate = False, exact = False
                                  )
                                                                                  )
              if searchTerm_F[j] in genform_F_collect:
                auto_differenceEqGamma2Tens[i][j] =                                       (
                                                    simplify
                                                    (
                                                     removePWE(genform_F_collect[searchTerm_F[j]])
                                                    )
                                                                                          )
                #pprint( auto_differenceEqGamma2Tens[i][j] )
              else:
                auto_differenceEqGamma2Tens[i][j] = 0
              pass# IF
              #
              if S.One in genform_F_collect:
                genform_F_collectWithin = genform_F_collect
                for kIdx in range(0,3):
                  if kIdx != j:
                    if S.One in genform_F_collectWithin:
                      genform_F_collectWithin =                                  (
                                                collect
                                                (
                                                 genform_F_collectWithin[S.One],
                                                 searchTerm_F[kIdx],
                                                 evaluate = False, exact = False
                                                )
                                                                                 )
                    pass# IF
                  pass# IF
                pass# FOR
                if S.One in genform_F_collectWithin:
                  print                                                                       (
                        '\n\nIn F[' + repr(i) + '] collection. There is S.One (zero-order)' +
                        ' F[x, y, z] term in second-order term from first-order Gamma' +
                        ' collection. It contains: \n'
                                                                                              )
                  pprint( genform_F_collectWithin[S.One] )
                  print '\n\n'
                pass# IF
              pass# IF
            pass# FOR
          pass# IF
          #
          # Same for expr. with lapl.:
          #
          if searchTerm_Gamma**2 in eigenexp_lapl_Gamma_collect:
            for j in range(0,3):
              eigenexp_lapl_F_collect =                                                    (
                                        collect
                                        (
                                         eigenexp_lapl_Gamma_collect[searchTerm_Gamma**2],
                                         searchTerm_F[j],
                                         evaluate = False, exact = False
                                        )
                                                                                           )
              if searchTerm_F[j] in eigenexp_lapl_F_collect:
                auto_EigenFTLapltens2[i][j] =                                             (
                                          simplify
                                          (
                                           removePWE(eigenexp_lapl_F_collect[searchTerm_F[j]])
                                          )
                                                                                      )
                #pprint( auto_differenceEqTens[i][j] )
              else:
                auto_EigenFTLapltens2[i][j] = 0
              pass# IF
            pass# FOR
          pass# IF
          #
          # Same for expr. with curl:
          #
          if searchTerm_Gamma**2 in eigenexp_curl_Gamma_collect:
            for j in range(0,3):
              eigenexp_curl_F_collect =                                                    (
                                        collect
                                        (
                                         eigenexp_curl_Gamma_collect[searchTerm_Gamma**2],
                                         searchTerm_F[j],
                                         evaluate = False, exact = False
                                        )
                                                                                           )
              if searchTerm_F[j] in eigenexp_curl_F_collect:
                auto_EigenFTCurltens2[i][j] =                                              (
                                               simplify
                                               (
                                                removePWE(eigenexp_curl_F_collect[searchTerm_F[j]])
                                               )
                                                                                           )
                #pprint( auto_differenceEqTens[i][j] )
              else:
                auto_EigenFTCurltens2[i][j] = 0
              pass# IF
            pass# FOR
          pass# IF
          #
        pass# FOR i in range ( 0, 3 )
        #
        # PUT `chk_diff' HERE
        #
        chk_eigexp_lapl = [0 for i in xrange(3)]
        chk_eigexp_curl = [0 for i in xrange(3)]
        #
        chk_diff = [0 for i in xrange(3)]
        #
        if bCheckAutoCalculated:
          print 'Using auto-calculated values.'
          chk_eigexp_lapl0ord = dot_prod_l2r_tens_vect ( auto_EigenFTLaplPropagator, F )
          chk_eigexp_lapl1ord =                                                     (
                                prod_scal2vect
                                (
                                 Gamma,
                                 dot_prod_l2r_tens_vect( auto_EigenFTLapltens1, F )
                                )
                                                                                    )
          chk_eigexp_lapl2ord =                                                     (
                                prod_scal2vect
                                (
                                 Gamma**2,
                                 dot_prod_l2r_tens_vect( auto_EigenFTLapltens2, F )
                                )
                                                                                    )
          #
          chk_eigexp_curl0ord = dot_prod_l2r_tens_vect( auto_EigenFTCurlPropagator, F )
          chk_eigexp_curl1ord =                                                     (
                                prod_scal2vect
                                (
                                 Gamma,
                                 dot_prod_l2r_tens_vect( auto_EigenFTCurltens1, F )
                                )
                                                                                    )
          chk_eigexp_curl2ord =                                                     (
                                prod_scal2vect
                                (
                                 Gamma**2,
                                 dot_prod_l2r_tens_vect( auto_EigenFTCurltens2, F )
                                )
                                                                                    )
          #
          chk_diff0ord = dot_prod_l2r_tens_vect( auto_differenceEqTens, F )
          chk_diff1ord =                                                           (
                         prod_scal2vect
                         (
                          Gamma,
                          dot_prod_l2r_tens_vect( auto_differenceEqGamma1Tens, F )
                         )
                                                                                   )
          chk_diff2ord =                                                           (
                         prod_scal2vect
                         (
                          Gamma**2,
                          dot_prod_l2r_tens_vect( auto_differenceEqGamma2Tens, F )
                         )
                                                                                   )
        else:
          print 'Using predefined values.'
          chk_eigexp_lapl0ord = dot_prod_l2r_tens_vect( EigenFTLaplPropagator, F )
          chk_eigexp_lapl1ord =                                                (
                                prod_scal2vect
                                (
                                 Gamma,
                                 dot_prod_l2r_tens_vect( EigenFTLapltens1, F )
                                )
                                                                               )
          chk_eigexp_lapl2ord =                                                (
                                prod_scal2vect
                                (
                                 Gamma**2,
                                 dot_prod_l2r_tens_vect( EigenFTLapltens2, F )
                                )
                                                                               )
          #
          chk_diff0ord = dot_prod_l2r_tens_vect( differenceEqTens, F )
          chk_diff1ord =                                                      (
                         prod_scal2vect
                         (
                          Gamma,
                          dot_prod_l2r_tens_vect( differenceEqGamma1Tens, F )
                         )
                                                                              )
          chk_diff2ord =                                                      (
                         prod_scal2vect
                         (
                          Gamma**2,
                          dot_prod_l2r_tens_vect( differenceEqGamma2Tens, F )
                         )
                                                                              )
        pass# IF
        #
        for i in range(0,3):
          chk_eigexp_lapl[i] = (chk_eigexp_lapl0ord[i] + chk_eigexp_lapl1ord[i] +
                               chk_eigexp_lapl2ord[i])
          if ( bCheckAutoCalculated ) :
            chk_eigexp_curl[i] = (chk_eigexp_curl0ord[i] + chk_eigexp_curl1ord[i] +
                                 chk_eigexp_curl2ord[i])
          pass# IF ( bCheckAutoCalculated )
          chk_diff[i] = chk_diff0ord[i] + chk_diff1ord[i] + chk_diff2ord[i]
        pass# FOR
        #
        # WHAT IS THE GOAL OF THIS CALCULATIONS BELOW??????????
        #
        for i in range(0,6):
          # Collection of first-order terms:
          genform_Gamma_collects[i] = collect                 (
                                       expand( genforms[i] ),
                                       searchTerm_Gamma,
                                       evaluate = False,
                                       exact = False
                                                              )
          # Collection of second-order terms:
          genform_Gamma_collects[i + 6] = collect                 (
                                           expand( genforms[i] ),
                                           searchTerm_Gamma2,
                                           evaluate = False,
                                           exact = False
                                                                  )
          for j in range(0,3):
            genform_collects[3 * i + j] = collect                 (
                                           expand( genforms[i] ),
                                           F[j],
                                           evaluate = False,
                                           exact = True
                                                                  )
          pass# FOR
        pass# FOR
        #
        # WHAT IS THE GOAL OF THOSE CALCULATIONS ABOVE??????????
        #
        GF_E = Matrix ( auto_greensfunction_E )
        GF_E_inv = GF_E.inverse_ADJ()
        #GF_E_inv = GF_E
      pass# IF ( bFourierDomain == True )
      #
      end_timestamp = timeit.default_timer()
      print 'Search timer: '
      print end_timestamp - bgn_timestamp
      print '\n\n'
      bgn_timestamp = timeit.default_timer()
      #
      #########################################################################
      #                      PRINTING OUTPUT TO FILE                          #
      #########################################################################
      #
      f = open( 'test1.sympy.tex', 'w' )
      #
      f.write( '\\documentclass{article}\n' )
      f.write( '\\usepackage{amsmath}\n' )
      f.write( '\\usepackage{mathtools}\n' )
      f.write( '\\begin{document}\n' )
      f.write( '\\text{Results: }\\\n' )
      #
      if ( not bFourierDomain or bPrintAll ):
        text_varNamesSimple =                 (
                              [
                               'X[1] or RHO',
                               'X[2] or PHI',
                               'X[3] or Z'
                              ]
                                              )
        text_varNames =                                                 (
                        [
                         text_varNamesSimple[0] + ', e.g. for 1st raw',
                         text_varNamesSimple[1] + ', e.g. for 1st raw',
                         text_varNamesSimple[2] + ', e.g. for 1st raw'
                        ]
                                                                        )
        #
        f.write( '\\text{Raw output: }\n' )
        for i in range(0,3):
          if ( bForceExpandSimplified ):
            genforms[i] = expand( genforms[i] )
          pass# IF
          if ( bRemovePlainWaveExponent ):
            genforms[i] = removePWE(genforms[i])
            #genforms[i] = genforms[i] / plain_wave_exp
          pass# IF
          printTexExpressionInline( f, genforms[i] )
          if ( bFourierDomain ):
            printTexExpressionInline( f, expand( simplify( chk_eigexp_lapl[i] - genforms[i] ) ) )
            if ( bCheckAutoCalculated ) :
              printTexExpressionInline( f, expand( simplify( chk_eigexp_curl[i] - genforms[i] ) ) )
            pass# IF ( bCheckAutoCalculated )
          pass# IF
        pass# FOR
        #
        for i in range(0,3):
          if ( bForceExpandSimplified ):
            genforms[i+3] = expand( genforms[i+3] )
          pass# IF
          if ( bRemovePlainWaveExponent ):
            genforms[i+3] = removePWE(genforms[i+3])
            #genforms[i+3] = genforms[i+3] / plain_wave_exp
          pass# IF
          printTexExpressionInline( f, genforms[i+3] )
        pass# FOR
        #
        for i in range(0,3):
          if ( bForceExpandSimplified ):
            check_correct_difference = expand( genforms_LHS2RHS_diff[i] )
          else:
            check_correct_difference = genforms_LHS2RHS_diff[i]
          pass# IF
          if ( bRemovePlainWaveExponent ):
            check_correct_difference = removePWE(check_correct_difference)
            #check_correct_difference = check_correct_difference / plain_wave_exp
            diffs_comparision = simplify( chk_diff[i] - check_correct_difference )
          else:
            diffs_comparision = simplify( plain_wave_exp * chk_diff[i] - check_correct_difference )
          pass# IF
          printTexExpressionInline( f, check_correct_difference )
          printTexExpressionInline( f, diffs_comparision )
        pass# FOR
        #
        f.write( '\\text{Entire equations with Laplacian: }\n' )
        for i in range(0,3):
          #if ( bForceExpandSimplified ):
          #  genforms[i] = expand( genforms[i] )
          #pass# IF
          #if ( bRemovePlainWaveExponent ):
          #  genforms[i] = genforms[i] / plain_wave_exp
          #pass# IF
          f.write( '\\text{' + text_varNamesSimple[i] + ' component: }\n' )
          printTexExpressionInline( f, genforms[i] )
          if ( bFourierDomain ):
            f.write                                                                        (
               '\\text{Compar. of above eigen expr. with explicitly writen propagators ' +
               '(lapl): }\n'
                                                                                           )
            printTexExpressionInline                                (
             f,
             expand( simplify( chk_eigexp_lapl[i] - genforms[i] ) )
                                                                    )
            if ( bCheckAutoCalculated ) :
              f.write                                                                        (
                 '\\text{Compar. of above eigen expr. with explicitly writen propagators ' +
                 '(curl): }\n'
                                                                                             )
              printTexExpressionInline                                  (
               f,
               expand( simplify( chk_eigexp_curl[i] - genforms[i+3] ) )
                                                                        )
            pass# IF ( bCheckAutoCalculated )
          pass# IF
        pass# FOR
        #
        f.write( '\\text{Entire equations with Curl: }\n' )
        for i in range(0,3):
          #if ( bForceExpandSimplified ):
          #  genforms[i+3] = expand( genforms[i+3] )
          #pass# IF
          #if ( bRemovePlainWaveExponent ):
          #  genforms[i+3] = genforms[i+3] / plain_wave_exp
          #pass# IF
          f.write( '\\text{' + text_varNamesSimple[i] + ' component: }\n' )
          printTexExpressionInline( f, genforms[i+3] )
        pass# FOR
        #
        for i in range(0,3):
          f.write                                                    (
             '\\text{Diff. between eqs. with Lapl. and Curl. for ' +
             text_varNames[i] +
             ' if consider vectorial equation components: }\n'
                                                                     )
          if ( bForceExpandSimplified ):
            check_correct_difference = expand( genforms_LHS2RHS_diff[i] )
          else:
            check_correct_difference = genforms_LHS2RHS_diff[i]
          pass# IF
          if ( bRemovePlainWaveExponent ):
            check_correct_difference = removePWE(check_correct_difference)
            #check_correct_difference = check_correct_difference / plain_wave_exp
            diffs_comparision = simplify( chk_diff[i] - check_correct_difference )
          else:
            diffs_comparision = simplify( plain_wave_exp * chk_diff[i] - check_correct_difference )
          pass# IF
          printTexExpressionInline( f, check_correct_difference )
          f.write                                                                (
             '\\text{Compar. of above diff. with explicitly writen tens. for ' +
             text_varNames[i] +
             ' if consider vectorial equation componenets: }\n'
                                                                                 )
          printTexExpressionInline( f, diffs_comparision )
        pass# FOR
      pass#if
      #
      if ( bFourierDomain ):
        f.write( '\\text{Indexed equations are given below with lapl. ( indexes in ranges [0;2], [3;5], and [6;8] for F[1], F[2], and F[3] respectively ) and curl( indexes in ranges [9;11], [12;14], and [15;17] for F[1], F[2], and F[3] respectively ). }\\newline\n' )
        f.write( '\\text{Equations with Laplacian: }\\newline\n' )
        for i in range(0,6):
          for j in range(0,3):
            if ( 3 * i + j == 9 ):
              f.write( '\\text{Equations with Curl:}\\newline\n' )
            pass# IF
            genform_collect_temp = genform_collects[3 * i + j][F[j]]
            genform_collect_temp_exclude_Gamma =                                     (
                                             collect
                                             (
                                              genform_collects[3 * i + j][F[j]],
                                              searchTerm_Gamma,
                                              evaluate = False,
                                              exact = False
                                             )
                                                                                 )
            f.write( '\\text{Index of equation: '+repr(3 * i + j)+'}\\newline\n' )
            #printTexExpressionInline( f, 3 * i + j )
            f.write( '\\text{Component of vectorial equation: }\\newline\n' )
            printTexExpressionInline( f, xyz[i%3] )
            f.write( '\\text{Component of function F(O): }\\newline\n' )
            printTexExpressionInline( f, F[j] )
            printTexExpressionInline                          (
             f,
             expand ( simplify (
              removePWE(genform_collect_temp)
             ) )                                              )
            #printTexExpressionInline                                     (
            # f,
            # expand( simplify( genform_collect_temp / plain_wave_exp ) )
            #                                                             )
            f.write                                                       (
               '\\text{Terms withouth Gamma from ' + repr ( 3 * i + j ) +
               '\'th (previous, above) equation: }\\newline\n'
                                                                          )
            printTexExpressionInline                         (
             f,
             expand ( simplify (
              removePWE(genform_collect_temp_exclude_Gamma[S.One])
             ) )                                             )
            #printTexExpressionInline                                     (
            # f,
            # expand ( simplify (
            #  genform_collect_temp_exclude_Gamma[S.One] / plain_wave_exp
            # ) )                                                         )
          pass# FOR
        pass# FOR
        #
        f.write( '\\text{Components with Gamma ( Global Index with dot2 ), Gamma2 ( Global Index with dot3 ), and without it ( Global Index with dot1 ) from eq. with lapl. ( indexes in ranges [0;2] and [6;8] for first and second order respectively ) and curl( indexes in ranges [3;5] and [9;11] for first and second order respectively ):}\\newline\n' )
        for i in range(0,6):
          #
          # Separation of terms with Gamma and F for propagator calculation and moving terms with \
          # d/dz to the left hand side.
          #
          GammaClength = len ( genform_Gamma_collects[i] )
          f.write( '\\text{Length of Gamma['+repr(i)+']: '+repr(GammaClength)+'}\\newline\n' )
          #printTexExpressionInline( f, GammaClength )
          f.write( '\\text{Data in Gamma['+repr(i)+']: }\\newline\n' )
          # Print all zero-order terms from first-order collection:
          f.write( '\\text{Zero-order terms from first-order collection: }\\newline\n' )
          f.write( '\\text{Index: '+repr(i + 0.1)+'}\\newline\n' )
          #printTexExpressionInline( f, i + 0.1 )
          f.write( '    ' )
          printTexExpressionInline                (
           f,
           simplify
           (
            removePWE(genform_Gamma_collects[i][S.One])
           )                                      )
          #printTexExpressionInline                                       (
          # f,
          # simplify( genform_Gamma_collects[i][S.One] / plain_wave_exp )
          #                                                               )
          # Print all first-order terms from first-order collection:
          f.write( '\\text{First-order terms from first-order collection: }\\newline\n' )
          f.write( '\\text{Index: '+repr(i + 0.2)+'}\\newline\n' )
          #printTexExpressionInline( f, i + 0.2 )
          f.write( '    ' )
          printTexExpressionInline                                     (
           f,
           simplify
           (
            removePWE(genform_Gamma_collects[i][simplify(searchTerm_Gamma)])
           )                                                           )
          #printTexExpressionInline                                (
          # f,
          # simplify
          # (
          #  genform_Gamma_collects[i][simplify(searchTerm_Gamma)]
          #  / plain_wave_exp
          # )                                                      )
          if ( GammaClength > 2 ):
            # Print all second-order terms from first-order collection:
            f.write( '\\text{Second-order terms from first-order collection: }\\newline\n' )
            print                                                                                (
                  '\nPrinting all second-order terms from first-order collection i=' + repr(i) +
                  ' ...\n'
                                                                                                 )
            printTexExpressionInline                              (
             f,
             simplify
             (
              removePWE(genform_Gamma_collects[i][searchTerm_Gamma**2])
             )                                                    )
            #printTexExpressionInline                         (
            # f,
            # simplify
            # (
            #  genform_Gamma_collects[i][searchTerm_Gamma**2]
            #  / plain_wave_exp
            # )                                               )
          pass# IF
          if ( len ( genform_Gamma_collects[i + 6] ) > 2 ):
            # Print all `first-order' terms from second-order collection:
            f.write( '\\text{`First-order\' terms from second-order collection: }\\newline\n' )
            f.write( '\\text{Index: '+repr(i + 0.3)+'}\\newline\n' )
            #printTexExpressionInline( f, i + 0.3 )
            f.write( '    ' )
            print                                                                         (
                  '\nPrinting all `first-order\' terms from second-order collection i=' +
                  repr(i) + ' ...\n'
                                                                                          )
            printTexExpressionInline                                (
             f,
             simplify
             (
              removePWE(genform_Gamma_collects[i + 6][searchTerm_Gamma2])
             )                                                      )
            #printTexExpressionInline                           (
            # f,
            # simplify
            # (
            #  genform_Gamma_collects[i + 6][searchTerm_Gamma2]
            #  / plain_wave_exp
            # )                                                 )
          pass# IF
          # Collect all second-order terms within first-order collection:
          print '\nCollecting all second-order terms within first-order collection ...\n'
          tmp_CollectWithinCollection = collect                            (
                                         genform_Gamma_collects[i][S.One],
                                         searchTerm_Gamma2,
                                         exact = False,
                                         evaluate = False
                                        )
          # Print all collected `first-order' terms of second order from additonal second-order   \
          # collection:
          #f.write                                                                         (
          #   '\\text{`First-order\' terms of second order from additonal second-order ' +
          #   'collection: }\n'
          #                                                                                )
          #print                                                                         (
          #      '\nPrinting all collected `first-order\' terms of second order from ' +
          #      'additonal second-order collection ...\n'
          #                                                                              )
          #printTexExpressionInline                                    (
          # f,
          # simplify( tmp_CollectWithinCollection[searchTerm_Gamma2] )
          #                                                            )
        pass# FOR
      pass# IF
      u_axis = [u_i, u_j, u_k]
      f.write( '\\text{M}\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row = print_row + u_axis[j] * M[i][j]
        pass# FOR
        printTexExpressionInline( f, print_row )
      f.write( '\\text{N1}\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row = print_row + u_axis[j] * N1[i][j]
        printTexExpressionInline( f, print_row )
      f.write( '\\text{N2}\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row = print_row + u_axis[j] * N2[i][j]
        printTexExpressionInline( f, print_row )
      f.write( '\\text{EigenFTLaplPropagator}\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row = print_row + u_axis[j] * EigenFTLaplPropagator[i][j]
        printTexExpressionInline( f, print_row )
      f.write( '\\text{EigenFTLapltens1}\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row = print_row + u_axis[j] * EigenFTLapltens1[i][j]
        printTexExpressionInline( f, print_row )
      f.write( '\\text{EigenFTLapltens2}\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row = print_row + u_axis[j] * EigenFTLapltens2[i][j]
        printTexExpressionInline( f, print_row )
      f.write( '\\text{differenceEqTens}\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row = print_row + u_axis[j] * differenceEqTens[i][j]
        printTexExpressionInline( f, print_row )
      f.write( '\\text{auto differenceEqTens}\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row = print_row + u_axis[j] * auto_differenceEqTens[i][j]
        printTexExpressionInline( f, sympify( print_row ) )
      f.write( '\\text{Comparision of above two: }\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row = (print_row + u_axis[j] *
                      simplify(auto_differenceEqTens[i][j] - plain_wave_exp*differenceEqTens[i][j]))
        printTexExpressionInline( f, sympify( print_row ) )
      f.write( '\\text{differenceEqGamma1Tens}\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row = print_row + u_axis[j] * differenceEqGamma1Tens[i][j]
        printTexExpressionInline( f, print_row )
      f.write( '\\text{auto differenceEqGamma1Tens}\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row = print_row + u_axis[j] * auto_differenceEqGamma1Tens[i][j]
        printTexExpressionInline( f, sympify( print_row ) )
      f.write( '\\text{Comparision of above two: }\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row =                                                                            (
                      print_row + u_axis[j] *
                      simplify(auto_differenceEqGamma1Tens[i][j] - plain_wave_exp*differenceEqGamma1Tens[i][j])
                                                                                                 )
        printTexExpressionInline( f, sympify( print_row ) )
      f.write( '\\text{differenceEqGamma2Tens}\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row = print_row + u_axis[j] * differenceEqGamma2Tens[i][j]
        printTexExpressionInline( f, print_row )
      f.write( '\\text{auto differenceEqGamma2Tens}\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row = print_row + u_axis[j] * auto_differenceEqGamma2Tens[i][j]
        printTexExpressionInline( f, sympify( print_row ) )
      f.write( '\\text{Comparision of above two: }\n' )
      for i in range(0,3):
        print_row = 0
        for j in range(0,3):
          print_row = (print_row + u_axis[j] *
                      simplify(auto_differenceEqGamma2Tens[i][j] - plain_wave_exp*differenceEqGamma2Tens[i][j]))
        printTexExpressionInline( f, sympify( print_row ) )
      #
      # Reevaluating matricies A_1, A_2, B for the equation \
      #  ( gamma * A_1 + gamma**2 * A_2 ) * H = B * H        
      #
      #f.write( '\\text{Reevaluating matricies A_1, A_2, B for the equation ( gamma * A1 + gamma**2 * A2 ) * H = B * H}\\newline\n' )
      f.write( '\\text{Reevaluating matricies}\\newline\n' )
      f.write( '\\text{Auto-calculated difference between eq. with lapl. and curl in bevel principal axis basis of induced susceptibility tensor}\\newline\n' )
      f.write( '\\text{Tensor B}\n' )
      print_matr_by_row_with_subst( auto_differenceEqTens, f )
      f.write( '\\text{Tensor A1}\n' )
      print_matr_by_row_with_subst( auto_differenceEqGamma1Tens, f )
      f.write( '\\text{Tensor A2}\n' )
      print_matr_by_row_with_subst( auto_differenceEqGamma2Tens, f )
      f.write ( '\\text{equation for E components:}\\newline\n' )
      for i in range ( 0, 3 ) :
        printTexExpressionInline ( f, simplify ( _E_LHS[i] ) )
      pass# FOR i in range ( 0, 3 )
      f.write ( '\\text{Green\'s function linear operator:}\\newline\n' )
      for i in range ( 0, 3 ) :
        print_row = 0
        for j in range ( 0, 3 ) :
          print_row = ( print_row + u_axis[j] * simplify ( GF_E_inv[i, j] ) )
        pass# FOR j in range ( 0, 3 )
        printTexExpressionInline ( f, sympify ( print_row ) )
      pass# FOR i in range ( 0, 3 )
      f.write( '\\end{document}' )
      f.close()
      #
      # Skip all if calculating eigenvectors.
      #
      if ( bSolveEigenValue ):
        if eTensor4Eigenvalues == eTensors.DIFF:
          auto_calc_2ndord = auto_differenceEqGamma2Tens
          auto_calc_1stord = auto_differenceEqGamma1Tens
          auto_calc_0thord = auto_differenceEqTens
        pass# IF
        if eTensor4Eigenvalues == eTensors.LAPL:
          auto_calc_2ndord = auto_EigenFTLapltens2
          auto_calc_1stord = auto_EigenFTLapltens1
          auto_calc_0thord = auto_EigenFTLaplPropagator
        pass# IF
        if eTensor4Eigenvalues == eTensors.CURL:
          auto_calc_2ndord = auto_EigenFTCurltens2
          auto_calc_1stord = auto_EigenFTCurltens1
          auto_calc_0thord = auto_EigenFTCurlPropagator
        pass# IF
        if bCheckAutoCalculated:
          gamma_differenceEqGamma2Tens = prod_scal2matr( gamma, auto_calc_2ndord )
          A_2ndOrd = summ_tens( auto_calc_1stord, gamma_differenceEqGamma2Tens )
        else:
          gamma_differenceEqGamma2Tens = prod_scal2matr( gamma, differenceEqGamma2Tens )
          A_2ndOrd = summ_tens( differenceEqGamma1Tens, gamma_differenceEqGamma2Tens )
        pass# IF
        #invA = inverse_matr( A_2ndOrd )
        CA_matr = Cmatrix3x3( matr_trans( A_2ndOrd ) )
        DetA = simplify( powdenest( det3x3( A_2ndOrd ), force = True ) )
        getSimpleTimingData( 'Inverse matrix calculation timer' )
        bgn_timestamp_old = bgn_timestamp
        for i in range(0,3):
          for j in range(0,3):
            #invA[i][j] = simplify( powdenest( invA[i][j], force = True ), ratio = oo )
            CA_matr[i][j] = simplify( powdenest( CA_matr[i][j], force = True ) )
            end_timestamp = timeit.default_timer()
            print( 'Inverse matrix [' + repr(i) + '][' + repr(j) + '] simplification timer:' )
            print end_timestamp - bgn_timestamp
            print '\n\n'
            bgn_timestamp = timeit.default_timer()
          pass# FOR
        pass# FOR
        end_timestamp = timeit.default_timer()
        print 'Inverse matrix simplification timer:'
        print end_timestamp - bgn_timestamp_old
        print '\n\n'
        bgn_timestamp = timeit.default_timer()
        #Tens_2ndOrd = mult1_2nd_ord_tens( invA, differenceEqTens )
        if bCheckAutoCalculated:
          Tens_2ndOrd = mult1_2nd_ord_tens( CA_matr, auto_calc_0thord )
        else:
          Tens_2ndOrd = mult1_2nd_ord_tens( CA_matr, differenceEqTens )
        pass# IF
        I_gamma = prod_scal2matr( -gamma, Identity_M )
        DetA_x_IxGamma = prod_scal2matr( DetA, I_gamma )
        Tens_gamma = summ_tens( Tens_2ndOrd, I_gamma )
        getSimpleTimingData( 'Eigenexpression yield timer' )
        bgn_timestamp_old = bgn_timestamp
        for i in range(0,3):
          for j in range(0,3):
            #print '\n\nEigenexpression:\n'
            #pprint( Tens_gamma[i][j] )
            #print '\n\n'
            Tens_gamma[i][j] = simplify( powdenest( Tens_gamma[i][j], force=True ) )
            end_timestamp = timeit.default_timer()
            print( 'Eigenexpression [' + repr(i) + '][' + repr(j) + '] simplification timer:' )
            print end_timestamp - bgn_timestamp
            print '\n\n'
            bgn_timestamp = timeit.default_timer()
          pass# FOR
        pass# FOR
        end_timestamp = timeit.default_timer()
        print 'Eigenexpression simplification timer:'
        print end_timestamp - bgn_timestamp_old
        print '\n\n'
        bgn_timestamp = timeit.default_timer()
        #
        DET3X3_Tgamma = det3x3( Tens_gamma )
        #
        end_timestamp = timeit.default_timer()
        print 'DETERMINANT YIELDING TIMER:'
        print end_timestamp - bgn_timestamp
        print '\n\n'
        bgn_timestamp = timeit.default_timer()
        #
        if ( not bDoNotHandleDeterminant ) :
          #eigendet = simplify( powdenest( det3x3( Tens_gamma ), force = True ) )
          #eigendet = simplify( det3x3( Tens_gamma ) )
          #eigendet = powdenest( det3x3( Tens_gamma ), force = True )
          eigendet = expand( DET3X3_Tgamma )# HERE
          #eigendet = sympify( DET3X3_Tgamma )
        else :
          eigendet = DET3X3_Tgamma
        pass# IF ( not bDoNotHandleDeterminant )
        #
        end_timestamp = timeit.default_timer()
        print 'DETERMINANT HANDLING (expanding) TIMER:'
        print end_timestamp - bgn_timestamp
        print '\n\n'
        bgn_timestamp = timeit.default_timer()
        #
        #if ( not bDoNotHandleDeterminant ) :
        #  #eigendet = simplify( eigendet )# HERE
        #pass# IF ( not bDoNotHandleDeterminant )
        #
        end_timestamp = timeit.default_timer()
        print 'DETERMINANT HANDLING (simplification) TIMER:'
        print end_timestamp - bgn_timestamp
        print '\n\n'
        bgn_timestamp = timeit.default_timer()
        #
        gammaCollection = collect( eigendet, gamma, evaluate = False, exact = False )# HERE
        gammaPowers = gammaCollection.keys()# HERE
        #
        end_timestamp = timeit.default_timer()
        print 'Variable power separtion timer:'
        print end_timestamp - bgn_timestamp
        print '\n\n'
        bgn_timestamp = timeit.default_timer()
        #
        #a_6, a_5, a_4, a_3, a_2, a_1, a_0 = symbols                              (
        #                                     'a_6, a_5, a_4, a_3, a_2, a_1, a_0'
        #                                     ,
        #                                     real = True
        #                                                                         )
        #eq_test =(a_6 * gamma_TEST**6 + a_5 * gamma_TEST**5 +
        #          a_4 * gamma_TEST**4 + a_3 * gamma_TEST**3 +
        #          a_2 * gamma_TEST**2 + a_1 * gamma_TEST + a_0)
        #eigenvalues_test_6thOrd = solve( Eq( eq_test, 0 ), gamma_TEST )
        #eigenvalues = solve( Eq( eigendet, 0 ), gamma )
        #
        end_timestamp = timeit.default_timer()
        print 'Eigenvalues calculation timer:'
        print end_timestamp - bgn_timestamp
        print '\n\n'
        bgn_timestamp = timeit.default_timer()
        #
        f = open( 'test1.sympy.tex', 'w' )
        f.write( '\\documentclass{article}\n\\usepackage{amsmath}\n' )
        f.write( '\\usepackage{mathtools}\n\\begin{document}\n' )
        f.write( '\\text{Data. }\n' )
        #
        detA_simplified = simplify( powdenest( (DetA)**3, force = True ) )
        #
        end_timestamp = timeit.default_timer()
        print '3rd. pow. det. simplification timer:'
        print end_timestamp - bgn_timestamp
        print '\n\n'
        bgn_timestamp = timeit.default_timer()
        #
        for gamma_power in gammaPowers:
          if gamma_power in gammaCollection:
            f.write( '\\text{Gamma power: }\n' )
            print '\n'
            pprint( gamma_power )
            print '\n'
            printTexExpressionInline( f, gamma_power )
            f.write( '\\text{Gamma constant multiplier: }\n' )
            gammaConstMult = simplify( gammaCollection[gamma_power] )
            end_timestamp = timeit.default_timer()
            print 'Gamma const. mult. simplification timer:'
            print end_timestamp - bgn_timestamp
            print '\n\n'
            bgn_timestamp = timeit.default_timer()
            printTexExpressionInline( f, gammaConstMult )
            #f.write( '\\text{Gamma constant divided over det3x3: }\n' )
            #gammaConst_div_DetA = simplify( gammaConstMult / detA_simplified )
            #end_timestamp = timeit.default_timer()
            #print 'Gamma const. mult. divided over det. simplification timer:'
            #print end_timestamp - bgn_timestamp
            #print '\n\n'
            #bgn_timestamp = timeit.default_timer()
            #printTexExpressionInline( f, gammaConst_div_DetA )
            #f.write( '\\text{Gamma constant multiplier factorized: }\n' )
            #printTexExpressionInline( f, factor( gammaConstMult ) )
            #end_timestamp = timeit.default_timer()
            #print 'Gamma const. mult. factorization timer:'
            #print end_timestamp - bgn_timestamp
            #print '\n\n'
            #bgn_timestamp = timeit.default_timer()
          pass# IF
        pass# FOR
        #
        #for evalue in eigenvalues:
        #  printTexExpressionInline( f, evalue )
        #pass# FOR
        #
        f.write( '\\end{document}' )
        f.close()
        end_timestamp = timeit.default_timer()
        print 'Print timer:'
        print end_timestamp - bgn_timestamp
        print '\n\n'
        bgn_timestamp = timeit.default_timer()
        #
        #print 'Eigevalues:'
        #pprint( eigenvalues )
        #pprint( factor( eigendet, gaussian = True ) )
        #print 'Matrix determinant:'
        #pprint( eigendet )
      pass# IF ( bSolveEigenValue )
      #
      #########################################################################
      #^^^^^                 PRINTING OUTPUT TO FILE                     ^^^^^#
      #########################################################################
      #
      end_timestamp = timeit.default_timer()
      print 'print timer:'
      print end_timestamp - bgn_timestamp
      print '\n\n'
      bgn_timestamp = timeit.default_timer()
    else:
      Hdummy[2] =                                                            (
                  -( I/( k_0*muT[2][2] ) * ( diff( Ey, x ) - diff( Ex, y ) ) -
                  1/muT[2][2] * summ2( muT, Hdummy, 2 ) )
                                                                             )
      Edummy[2] =                                                              (
                  I/( k_0*epsilonT[2][2] ) * ( diff( Hy, x ) - diff( Hx, y ) ) -
                  1/epsilonT[2][2] * summ2( epsilonT, Edummy, 2 )
                                                                               )
      dEx_dz = diff( Edummy[2], x ) + I*k_0*summ3( muT, Hdummy, 1 )
      dEy_dz = diff( Edummy[2], y ) - I*k_0*summ3( muT, Hdummy, 0 )
      dHx_dz = diff( Hdummy[2], x ) - I*k_0*summ3( epsilonT, Edummy, 1 )
      dHy_dz = diff( Hdummy[2], y ) + I*k_0*summ3( epsilonT, Edummy, 0 )
      #
      #########################################################################
      #                      PRINTING OUTPUT TO FILE                          #
      #########################################################################
      #
      f = open( 'test.sympy.tex', 'w' )
      f.write( '\\documentclass{article}\n\\usepackage{amsmath}\n' )
      f.write( '\\usepackage{mathtools}\n\\begin{document}\n' )
      f.write                                                         (
         '\\text{Fourier plane wave decomposition:}\n' +
         '\\begin{equation}\n\\begin{gathered}\n\\begin{multlined}\n'
                                                                      )
      f.write( latex( expand(simplify(simplify( dEx_dz )-simplify( dHx_dz ))), mode='inline' ) )
      f.write                                                         (
         '\n\\end{multlined}\n\\end{gathered}\n\\end{equation}\n' +
         '\\begin{equation}\n\\begin{gathered}\n\\begin{multlined}\n'
                                                                      )
      f.write( latex( expand(simplify(simplify( dEy_dz )-simplify( dHy_dz ))), mode='inline' ) )
      #f.write                                                         (
      #   '\n\\end{multlined}\n\\end{gathered}\n\\end{equation}\n' +
      #   '\\begin{equation}\n\\begin{gathered}\n\\begin{multlined}\n'
      #                                                                )
      #f.write( latex( simplify( dHx_dz ), mode='inline' ) )
      #f.write                                                         (
      #   '\n\\end{multlined}\n\\end{gathered}\n\\end{equation}\n' +
      #   '\\begin{equation}\n\\begin{gathered}\n\\begin{multlined}\n'
      #                                                                )
      #f.write( latex( simplify( dHy_dz ), mode='inline' ) )
      f.write( '\n\\end{multlined}\n\\end{gathered}\n\\end{equation}\n\\end{document}' )
      f.close()
      #
      #########################################################################
      #^^^^^                 PRINTING OUTPUT TO FILE                     ^^^^^#
      #########################################################################
      #
      end_timestamp = timeit.default_timer()
      print 'simplify & print timer:'
      print end_timestamp - bgn_timestamp
      print '\n\n'
    pass# IF
pass# IF
print 'v1.4.build04'
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
