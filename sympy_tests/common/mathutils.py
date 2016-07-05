#!/usr/bin/env python
#
Ident_M = [[1, 0, 0], [0, 1, 0], [0, 0, 1] ]
#
# Dot prod.:
#
def dot_product ( U, V, bForceSimplify = False ) :
  ##return simplify( sympify( U[0] * V[0] + U[1] * V[1] + U[2] * V[2] ) )
  dot_product_res = U[0] * V[0] + U[1] * V[1] + U[2] * V[2]
  if ( bForceSimplify ) :
    return simplify ( dot_product_res )
  else :
    return dot_product_res
pass# DEF dot_product ( U, V )
#
# Dot prod. from l2r of vect. with 2nd. ord. tens.:
#
def dot_prod_l2r_tens_vect ( M, V, bForceSimplify = False ) :
  R = [0 for i in xrange(3)]
  U = [0 for i in xrange(3)]
  for i in range ( 0, 3 ) :
    for j in range ( 0, 3 ) :
      R[j] = M[i][j]
    pass# FOR j in range ( 0, 3 )
    U[i] = dot_product( R, V, bForceSimplify )
  pass# FOR i in range ( 0, 3 )
  #
  return U
pass# DEF dot_prod_l2r_tens_vect ( M, V )
#
