#!/usr/bin/octave -qf
#
pkg load symbolic
###############################################################################
#                                   SYMBOLS                                   #
###############################################################################
#
syms z x r;
simplify (x + z + 1 / (x - 2*z) + r / (x^2));
########syms w k w_0 t t_0 chi_1221 chi_1122 chi_1212 r_p a_p t_p
########syms tau infSymb x y z
#########
#######################################################################################
#########^^^^^                              SYMBOLS                              ^^^^^#
#######################################################################################
#########
#########
#######################################################################################
#########                                  CONSTANTS                                  #
#######################################################################################
#########
########eTransProfiles.TEM = 0
########eTransProfiles.TE = 1
########eTransProfiles.TM = 2
#########
#######################################################################################
#########^^^^^                             CONSTANTS                             ^^^^^#
#######################################################################################
#########
#########
#######################################################################################
#########                                   SETTINGS                                  #
#######################################################################################
#########
########bSwitchExec2LibOnly   = False
########bQuasiSteadyState     =  True
########bLaguerreSymbolic     = False
########bSubstituteValues     =  True
########bCheckEquals          = False
########bDoNotEvaluateDiffs   =  True
########bUsePump2             = False
########bPump1QWP             =  True# TODO # HERE
########bPump2QWP             =  True# TODO
########bProbeQWP             =  True# TODO
########bIsNumericConstants   =  True# Is w_0, k, chi_1221, chi_1122, chi_1212 an     \
#########                            # arbitrary numeric constants?
########bSimplifyEachEquation = False#
#########
########iRotatePump1Pump2 = 0# 90 180 270 # TODO # HERE
########iRotatePump1Probe = 0# 90 180 270 # TODO
#########
########fPump2ZDisplacement = 0.0# TODO # HERE
########fProbeZDisplacement = 0.0# TODO
#########
########eSwitch_probe_TE2TM = eTransProfiles.TEM
########eSwitch_pump_TE2TM  = eTransProfiles.TEM
#########
#######################################################################################
#########^^^^^                              SETTINGS                             ^^^^^#
#######################################################################################
#########
#########
#######################################################################################
#########                                  CONSTANTS                                  #
#######################################################################################
#########
########wavelength = 630e-9
#########
########if ( bIsNumericConstants )
########  k = 2 * pythmath.pi / wavelength
########  if ( bSimplifyEachEquation )
########    k = simplify ( k )
########  endif# IF ( bSimplifyEachEquation )
########  w_0 = 1e-3
########  chi_1chloronaphtallene = ( 1.632 * 3.04e+14 * 488e-9 ) / ( 4 * pythmath.pi )
########  chi_1chloronaphtallene_diff = chi_1chloronaphtallene - chi_1chloronaphtallene * 0.05
########  nrand_vals = stdnormal_rnd ( 3 )
########  chi_1221 = chi_1chloronaphtallene_diff + chi_1chloronaphtallene * 0.1 * nrand_vals(0)
########  chi_1122 = chi_1chloronaphtallene_diff + chi_1chloronaphtallene * 0.1 * nrand_vals(1)
########  chi_1212 = chi_1chloronaphtallene_diff + chi_1chloronaphtallene * 0.1 * nrand_vals(2)
########  # HERE!!!02.01.15-22:10!!! Assign correct values for 1-cloronaphtallene or  \
########  # anethole                                                                   
########  # DONE !!!02.10.15-23:10!!!
########endif# IF ( bIsNumericConstants )#
#########
#######################################################################################
#########^^^^^                             CONSTANTS                             ^^^^^#
#######################################################################################
#########
#########x = rho * cos ( phi )
#########y = rho * sin ( phi )
#########
########z_r = k * w_0 .^ 2 / 2
########if ( bSimplifyEachEquation )
########  z_r = simplify ( z_r )
########endif# IF ( bSimplifyEachEquation )
##
#rho = sqrt ( x .^ 2 + y .^ 2 )
#if ( bSimplifyEachEquation )
#  rho = simplify ( rho )
#pass# IF ( bSimplifyEachEquation )
##
#def r_vector( z_displ )
#  r_vector_res = sqrt( rho**2 + ( z + z_displ )**2 )
#  if ( bSimplifyEachEquation )
#    r_vector_res = simplify( r_vector_res )
#  pass# IF ( bSimplifyEachEquation )
#  #
#  return r_vector_res
#pass# DEF r_vector ( z_displ )
##
#f = 1 / ( k * w_0 )
#if ( bSimplifyEachEquation )
#  f = simplify ( f )
#pass# IF ( bSimplifyEachEquation )
##
#phi = atan( y / x )
#if ( bSimplifyEachEquation ) :
#  phi = simplify ( phi )
#pass# IF ( bSimplifyEachEquation )
##
#L_ra = Function( 'L_ra' )
#L = L_ra( r_p, a_p, t_p )
#
#
#
#
#
#
#
#
#
#
#
#example octave script
#syms s
#arg_list = argv ();
#num = str2int(arg_list{1});
#printf ("Name of Octave script: ", program_name ());
#simplify ((Sin(s))^2 + (Cos(s))^2)
#tic();
#for i=1:num
#   a(i) = i;
# endfor
# elapsed = toc();
# printf("Elapsed time: %.4f seconds", elapsed);
