#!/usr/bin/env python
#
import timeit
import os
from sympy import latex
from sympy import sympify
#
#
###############################################################################
#                                     DATA                                    #
###############################################################################
#
maximacodebgn =                                                       (
            "batchload ( \"../maxima_srcs/common/utils.mac\" );\n" +
            "/*simp : true;*/\n" +
            "/*ratsimpexpons : true;*/\n" +
            "/*logsimp : true;*/\n" +
            "/*factorflag : true;*/\n" +
            "/*simpsum : false;*/\n" +
            "/*facexpand : false;*/\n" +
            "the_main_routine () := block\n" +
            "(\n" +
            "  redirect_stream : openw ( \"redirect2sympy.txt\" ),\n"
                                                                      )
maximacodeend =                                                               (
            "  close ( redirect_stream ),\n" +
            "  redirect_done_stream : openw ( \"redirect2sympy.done.txt\" ),\n" +
            "  printf ( redirect_done_stream, \"DONE\" ),\n" +
            "  close ( redirect_done_stream ),\n" +
            "  return ( done )\n" +
            ");\n" +
            "the_main_routine()$"
                                                                              )
#
###############################################################################
#^^^^^                                DATA                               ^^^^^#
###############################################################################
#
def enum(**enums):
  return type('Enum', (), enums)
#
total_timestamp = timeit.default_timer()
bgn_timestamp = timeit.default_timer()
end_timestamp = timeit.default_timer()
#
#
###############################################################################
#                                  FUNCTIONS                                  #
###############################################################################
#
def simplifySingleExpressioninMaxima ( expr ) :
  f = open ( "redirect2maxima.max", "w" )
  f.write ( maximacodebgn )
  #f.write                                                          (
  #   "printf ( redirect_stream, \"~a\", SAGElikeFullSimplify ( 0" +
  #   " ) ),\n"
  #                                                                 )
  f.write                                                          (
     "  printf ( redirect_stream, \"~a\", SAGElikeLiteSimplify ( " +
     ((repr(expr)).replace ( "**", "^" )).replace( "I", "%i" ) +
     " ) ),\n"
                                                                   )
  f.write ( maximacodeend )
  f.close ()
  os.system ( "maxima -q --very-quiet -b ../maxima_srcs/common/wrapper.max" )
  idx_wait_counter = 1
  while ( not os.path.isfile("redirect2sympy.done.txt") ) :
    ++idx_wait_counter
    if ( idx_wait_counter % 1000 == 0 ) :
      print ( "\nWaiting for Maxima simplified output ...\n" )
    pass# IF ( idx_wait_counter % 1000 == 0 )
  pass# WHILE ( not os.path.isfile("redirect2sympy.done.txt") )
  r = open ( "redirect2sympy.txt", "r" )
  expr2sympify = r.readline ()
  expr2sympify = expr2sympify.replace ( "^", "**" )
  expr2sympify = expr2sympify.replace ( "%e", "(exp(1))" )
  expr2sympify = expr2sympify.replace ( "%i", "(I)" )
  simplified_res = sympify ( expr2sympify )
  r.close ()
  os.remove( "redirect2sympy.done.txt" )
  return simplified_res
  #return 0
pass# DEF simplifySingleExpressioninMaxima ( expr )
#
def printMaximaRedirectMessage () :
  print                                                                                   (
        '\n\nRedirecting internal representation of simplified variable to maxima batch ' +
        'file\n\n'
                                                                                          )
pass# DEF printMaximaRedirectMessage (  )
#
###############################################################################
#^^^^^                             FUNCTIONS                             ^^^^^#
###############################################################################
#
#
def getSimpleTimingData( timingMessageText, bGetTotal = False ) :
  global end_timestamp
  global bgn_timestamp
  end_timestamp = timeit.default_timer()
  print timingMessageText + ':  '
  if ( bGetTotal ) :
    print end_timestamp - total_timestamp
  else :
    print end_timestamp - bgn_timestamp
  bgn_timestamp = timeit.default_timer()
pass# DEF getSimpleTimingData( timingMessageText )
#
def beginSimpleTiming():
  getSimpleTimingData( 'Simple test timer' )
pass# DEF
#
def getTotalTime():
  getSimpleTimingData( 'Total timer', bGetTotal = True )
pass# DEF
#
def printTexExpressionInline( texfile, expr ):
  texfile.write( '\\begin{equation}\n\\begin{split}\n' )
  texfile.write( latex( expr, mode = 'inline' ).replace ( '$', '' ) )
  texfile.write( '\n' )
  texfile.write( '\\end{split}\n\\end{equation}\n' )
pass# DEF
#
