#include <stdlib.h>
#include "quadqags2d_internals.hpp"


VariablesAndFunction::VariablesAndFunction () : \
  y (0.0),                                      \
  lower_boundary (0.0),                         \
  upper_boundary (0.0),                         \
  params (NULL),                                \
  function (NULL)                                
{}

VariablesAndFunction::VariablesAndFunction      \
(                                               \
  const VariablesAndFunction &vafReference      \
) :                                             \
  y (vafReference.y),                           \
  lower_boundary (vafReference.lower_boundary), \
  upper_boundary (vafReference.upper_boundary), \
  params (vafReference.params),                 \
  function (vafReference.function)               
{}
