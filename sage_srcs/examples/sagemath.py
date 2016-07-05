r"""
Conversion between Sage Symbolic Ring and Maxima via ECL library interface
"""
###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2009 Nils Bruin <nbruin@sfu.ca>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################


from sage.libs.ecl import *
import sage.symbolic.expression
from sage.symbolic.ring import is_SymbolicVariable
from sage.symbolic.ring import SR

ecl_eval("(require 'maxima)")
ecl_eval("(in-package :maxima)")

car=EclObject("car")
cdr=EclObject("cdr")
caar=EclObject("caar")
cadadr=EclObject("cadadr")
meval=EclObject("meval")
NIL=EclObject("NIL")
ratdisrep=EclObject("ratdisrep")

#first we define the dictionary in text form because that is easier to edit

op_sage_to_max = {
    sage.symbolic.expression.operator.abs : "MABS",
    sage.symbolic.expression.operator.add : "MPLUS",
    sage.symbolic.expression.operator.div : "MQUOTIENT",
    sage.symbolic.expression.operator.eq : "MEQUAL",
    sage.symbolic.expression.operator.ge : "MGEQP",
    sage.symbolic.expression.operator.gt : "MGREATERP",
    sage.symbolic.expression.operator.le : "MLEQP",
    sage.symbolic.expression.operator.lt : "MLESSP",
    sage.symbolic.expression.operator.mul : "MTIMES",
    sage.symbolic.expression.operator.ne : "MNOTEQUAL",
    sage.symbolic.expression.operator.neg : "MMINUS",
    sage.symbolic.expression.operator.pow : "MEXPT",
    sage.symbolic.expression.operator.or_ : "MOR",
    sage.symbolic.expression.operator.and_ : "MAND",
    sage.functions.trig.acos : "%ACOS",
    sage.functions.trig.acot : "%ACOT",
    sage.functions.trig.acsc : "%ACSC",
    sage.functions.trig.asec : "%ASEC",
    sage.functions.trig.asin : "%ASIN",
    sage.functions.trig.atan : "%ATAN",
    sage.functions.trig.cos : "%COS",
    sage.functions.trig.cot : "%COT",
    sage.functions.trig.csc : "%CSC",
    sage.functions.trig.sec : "%SEC",
    sage.functions.trig.sin : "%SIN",
    sage.functions.trig.tan : "%TAN",
    sage.functions.log.exp : "%EXP",
    sage.functions.log.ln : "%LOG",
    sage.functions.log.log : "%LOG",
    sage.functions.other.factorial : "MFACTORIAL",
    sage.functions.other.erf : "%ERF",
    sage.calculus.calculus._limit : "$LIMIT",
    sage.calculus.calculus.dummy_diff : "$DIFF",
    sage.calculus.calculus._integrate : "$INTEGRATE"
}
#we compile the dictionary
op_sage_to_max = dict([(k,EclObject(op_sage_to_max[k])) for k in op_sage_to_max])

#and also construct the dictionary in the other direction
op_max_to_sage = dict([(op_sage_to_max[k],k) for k in op_sage_to_max])

def add_vararg(*args):
    S=0
    for a in args:
        S=S+a
    return S

def mul_vararg(*args):
    P=1
    for a in args:
        P=P*a
    return P

def sage_rat(x,y):
    return x/y

op_max_to_sage[EclObject("MPLUS")]=add_vararg
op_max_to_sage[EclObject("MTIMES")]=mul_vararg
op_max_to_sage[EclObject("RAT")]=sage_rat

sym_sage_to_max={}
sym_max_to_sage={}

def symbol_factory(packager,prefix):
    r"""
    Produce a symbol generator
    
    EXAMPLES::

        sage: from sagemax import *
        ;;; Loading #P"/usr/local/sage/4.1.2/local/lib/ecl/maxima.fas"
        sage: f=symbol_factory(SR,"myname")
        sage: [f.next() for i in [1..10]]  

        [myname1,
         myname2,
         myname3,
         myname4,
         myname5,
         myname6,
         myname7,
         myname8,
         myname9,
         myname10]
    
    """
    i=1
    while True:
        yield packager(prefix+str(i))
        i +=1

def mrat_to_sage(expr):
    r"""
    Convert a maxima MRAT expression to Sage SR
    
    Maxima has an optimised representation for multivariate rational expressions.
    The easiest way to translate those to SR is by first asking maxima to give
    the generic representation of the object. That is what RATDISREP does in
    maxima.
    """
    return max_to_sage(meval(EclObject([[ratdisrep],expr])))

special_max_to_sage={
    EclObject("MRAT") : mrat_to_sage
}

maxop=symbol_factory(EclObject,"SAGE-OP-")
maxsym=symbol_factory(EclObject,"SAGE-SYM-")
sageop=symbol_factory(sage.calculus.calculus.function,"maxima_op_")
sagesym=symbol_factory(sage.calculus.calculus.var,"maxima_sym_")

def max_read(s):
    r"""
    Maxima's reader
    
    Convert a string into a maxima object via maxima's reader. Only use for
    development purposes, since the macro used in this routine suffers from
    performance loss with heavy usage. Does give good insight into internal
    representation of maxima objects, though.
    
    EXAMPLES::
    
        sage: from sagemax import *
        ;;; Loading #P"/usr/local/sage/4.1.2/local/lib/ecl/maxima.fas"
        sage: max_read("integral(sin(x),x)")
        <ECL: (($INTEGRAL) ((%SIN) $X) $X)>

    """
    return cadadr(EclObject("#$%s$"%s))

def pyobject_to_max(obj):
    r"""
    Translate a python object into a maxima object
    
    Mainly just a wrapper around EclObject, but needs to be in place because
    some objects might be translated better into maxima than just into lisp
    (think vectors and matrices).
    """
    return EclObject(obj)

def sage_to_max(expr):
    r"""
    Convert a symbolic ring expression to maxima
    
    EXAMPLES::

        sage: from sagemax import *
        ;;; Loading #P"/usr/local/sage/4.1.2/local/lib/ecl/maxima.fas"
        sage: from sage.calculus.calculus import _integrate
        sage: I=_integrate(cos(x),x)
        sage: sage_to_max(I)
        <ECL: (($INTEGRATE) ((%COS) SAGE-SYM-1) SAGE-SYM-1)>
        sage: meval(sage_to_max(I))
        <ECL: ((%SIN SIMP) SAGE-SYM-1)>
        sage: max_to_sage(meval(sage_to_max(I)))
        sin(x)
        
    This process has defined a mapping from the sage SR element x to maxima::
    
        sage: sym_sage_to_max[x]
        <ECL: SAGE-SYM-1>

    And this mapping exists in the other direction too (this is why EclObjects
    should be quickly hashable)::
    
        sage: sym_max_to_sage[sym_sage_to_max[x]]
        x

    Expressions that do not have a special meaning get translated relatively
    robustly. For instance, formal derivatives have not been implemented yet::
    
        sage: f=SFunction('f')
        sage: L=sage_to_max(derivative(f(x),x))
        sage: L
        <ECL: ((SAGE-OP-2) SAGE-SYM-1)>

    As you can see, the nature of the derivative is not translated at all.
    However, the reverse mapping unpacks that again::
    
        sage: max_to_sage(L)
        D[0](f)(x)
    """
    
    global op_sage_to_max, op_max_to_sage
    global sym_sage_to_max, sym_max_to_sage
    op = expr.operator()
    if op:
        if not (op in op_sage_to_max):
            op_max=maxop.next()
            op_sage_to_max[op]=op_max
            op_max_to_sage[op_max]=op
        return EclObject(([op_sage_to_max[op]], [sage_to_max(o) for o in expr.operands()]))
    elif is_SymbolicVariable(expr):
        if not expr in sym_sage_to_max:
            sym_max=maxsym.next()
            sym_sage_to_max[expr]=sym_max
            sym_max_to_sage[sym_max]=expr
        return sym_sage_to_max[expr]
    else:
        return pyobject_to_max(expr.pyobject())

def max_to_sage(expr):
    r"""
    Convert a maxima expression to sage symbolic ring
    
    EXAMPLES::

        sage: from sagemax import *
        ;;; Loading #P"/usr/local/sage/4.1.2/local/lib/ecl/maxima.fas"
        sage: from sage.calculus.calculus import _integrate
        sage: I=_integrate(cos(x),x)
        sage: sage_to_max(I)
        <ECL: (($INTEGRATE) ((%COS) SAGE-SYM-1) SAGE-SYM-1)>
        sage: meval(sage_to_max(I))
        <ECL: ((%SIN SIMP) SAGE-SYM-1)>
        sage: max_to_sage(meval(sage_to_max(I)))
        sin(x)
    
        
    """
    global op_sage_to_max, op_max_to_sage
    global sym_sage_to_max, sym_max_to_sage
    if expr.consp():
        op_max=caar(expr)
        if op_max in special_max_to_sage:
            return special_max_to_sage[op_max](expr)
        if not(op_max in op_max_to_sage):
            op=sageop.next()
            op_max_to_sage[op_max]=op
            op_sage_to_max[op]=op_max
        op=op_max_to_sage[op_max]
        max_args=cdr(expr)
        args=[]
        while not(max_args.nullp()):
            args.append(max_to_sage(car(max_args)))
            max_args=cdr(max_args)
        return op(*args)
    elif expr.symbolp():
        if not(expr in sym_max_to_sage):
            sym=sagesym.next()
            sym_max_to_sage[expr]=sym
            sym_sage_to_max[sym]=expr
        sym=sym_max_to_sage[expr]
        return sym
    else:
        return expr.python()
