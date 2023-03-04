"""
Symbolic Integration
"""


#*****************************************************************************
#       Copyright (C) 2013-2023 Svetlin Tassev <stassev@alum.mit.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

########################################################################
from sage.all import *
#print (erf(1.))
#from sage.symbolic.relation import solve
#from sage.calculus.calculus import limit
#from sage.calculus.functional import simplify
#from sage.misc.functional import integrate
#from sage.symbolic.assumptions import assume
from multiprocessing import Pool,TimeoutError
from time import time
#from sage.symbolic.expression import Expression
#from sage.calculus.functional import diff
#from sage.matrix.constructor import matrix
#from sage.misc.functional import det
#from sage.functions.log import Function_exp,Function_log1
#exp = Function_exp()
#log = Function_log1()
#from sage.functions.trig import arctan
#from sage.functions.hyperbolic import tanh
#from sage.symbolic.constants import pi
#from sage.calculus.var import var
#from sage.rings.integer import Integer
#from sage.functions.generalized import sign
#from sage.rings.infinity import Infinity

def sigmoid_logistic(x):
    """
    A sigmoid function.
    """
    return 1/(1+exp(-x))

def sigmoid_atan(x):
    """
    A sigmoid function.
    """
    return (arctan(x)/pi+Integer(1)/Integer(2))

def sigmoid_tanh(x):
    """
    A sigmoid function.
    """
    return (tanh(x)+1)/2

def sigmoid_gd(x):
    """
    A sigmoid function.
    """
    return 2*arctan(tanh(x/2))/pi+Integer(1)/Integer(2)

cdef _jacobian_and_substitution(ranges,new_vars,simplify,inf_var,sigmoid=[sigmoid_logistic],use_limits=True):
    """
    Take integration intervals and the new variable 
    names that live on the unit hypercube and compute the Jacobian 
    of the transformation between the integration variables and the 
    hypercube variables, as well as the transformation itself 
    between them. Passing a simplify value of True, will try to 
    simplify the results. See ``symbolic_multidim_integral``.
    """
    

    X = []
    v = []
    tosub = []
    transform = []
    
    for i in range(len(ranges)):
            if len(sigmoid)>1:
                T=sigmoid[i]
            else:
                T=sigmoid[0]
            x = ranges[i][0]
            y = new_vars[i]
            v.append(x)
            a = ranges[i][1]# a and b can explicitly depend on the integration 
                            # variables. we eliminate them below.
            b = ranges[i][2]
            if (use_limits):
                if str(abs(a))=='+Infinity':
                    a=inf_var*sign(a)
                if str(abs(b))=='+Infinity':
                    b=inf_var*sign(b)
            try:
                if str(abs(a))!='+Infinity' and str(abs(b))!='+Infinity':
                    try:
                        X.append((x-a)/(b-a))
                    except(ZeroDivisionError):
                        return 0,0  # NOTE: The function returns zero for zero 
                        # intervals (i.e. when in fact 1/jac=0). So, in this case 
                        # the integrator must return 0. When in fact jac=0, we raise 
                        # an exception below.
                    tosub.append((x,y*(b-a)+a))
                    transform.append(y*(b-a)+a)
                elif str(abs(a))!='+Infinity':
                    if b>0:
                        X.append(1-T(-x)/T(-a))
                        Tinv=x.subs(solve(1-T(-x)/T(-a)==y,x,solution_dict=True)[0])
                        tosub.append((x, Tinv))
                        transform.append(Tinv)
                    else:
                        X.append(1-T(x)/T(a))
                        Tinv=x.subs(solve(1-T(x)/T(a)==y,x,solution_dict=True)[0])
                        tosub.append((x, Tinv))
                        transform.append(Tinv)
                elif str(abs(b))!='+Infinity':
                    if a>0:
                        X.append(T(-x)/T(-b))
                        Tinv=x.subs(solve(T(-x)/T(-b)==y,x,solution_dict=True)[0])
                        tosub.append((x, Tinv))
                        transform.append(Tinv)
                    else:
                        X.append(T(x)/T(b))
                        Tinv=x.subs(solve(T(x)/T(b)==y,x,solution_dict=True)[0])
                        tosub.append((x, Tinv))
                        transform.append(Tinv)
                elif str(a)=='-Infinity' and (str(b)=='+Infinity' or str(b)=='Infinity'):
                    X.append(T(x))
                    Tinv=x.subs(solve(T(x)==y,x,solution_dict=True)[0])
                    tosub.append((x, Tinv))
                    transform.append(Tinv)
                elif str(b)=='-Infinity' and (str(a)=='+Infinity' or str(a)=='Infinity'):
                    X.append(T(-x))
                    Tinv=x.subs(solve(T(-x)==y,x,solution_dict=True)[0])
                    tosub.append((x, Tinv))
                    transform.append(Tinv)
                else:
                    raise Exception("Error: Not sure how to handle the integration ranges.") 
            except:
                print "Error: Not sure how to handle the integration ranges."
                raise # Exception("Error: Not sure how to handle the integration ranges.") 
            
    jac = det(matrix([[diff(xx,vv) for xx in X] for vv in v]))
    
    if (jac == 0): #So, when in fact jac=0, we raise an exception
        raise Exception('Error: Zero jacobian! Check your integration range.')
    
    tosub = dict(tosub)
    
    # This will eliminate variables for all integration ranges which make any sense. 
    # An example which does    make sense: \int^1_0 dx \int^1_x dy f(x,y)
    # An example which doesn't make sense: \int^y_0 dx \int^1_x dy f(x,y)
    for i in range(len(ranges)):
        tmp = []
        for t in transform:
            t1 = t.subs(tosub)
            if (simplify):
                t1 = t1.simplify_full()
            tmp.append(t1)
        transform = tmp
        
        jac = jac.subs(tosub)
        if (simplify):
            jac = jac.simplify_full()

            
    
    for i in transform:
        for vv in v:
            if i.has(vv):
                raise  Exception('Error: Cannot eliminate variables! Check your integration range.')
    
    if (jac == 0): #So, when in fact jac=0, we raise an exception
        raise Exception('Error: Zero jacobian! Check your integration range.')
        
    return jac,dict(zip(v,transform))




def symbolic_multidim_integral(func, *ranges, 
                         verbose=1, 
                         dimension_limit=0, time_limit=5,
                         simplify_func=False,sigmoid=[sigmoid_logistic],
                         dummy_var_prefix='_X_',
                         use_limits=False,
                         algorithm='Default'):
    r'''
    Calculate the symbolic multidimensional integral of a function. 
    
    If the integral cannot be obtained within the constraints imposed by 
    ``dimension_limit`` and ``time_limit``, then
    a partial result is provided if available. The partial result 
    is the result of doing as many of the integrals in the 
    multidimensional integral as possible. Any remaining integrals are over 
    dummy variables labelled ``_X_i`` by default.
    
    The symbolic integration is done by first changing variables so that 
    the new integration ranges span the unit hypercube. That transformation is 
    linear for proper integrals. For improper integrals, the transformation is
    non-linear, and is done using a sigmoid function. So, if you suspect the 
    integral is not converging because of that transformation, you may 
    want to try converting your integral to a proper integrals, and then
    take the limit to infinity of the endpoinds of the integration ranges.
        
    INPUT:
    
    -   ``func`` -- the integrand.
    -   ``ranges`` -- the integration intervals in each variable as 
        tuples, e.g. (x,0,10),(y,x^2,200). Each interval is of the form 
        (variable name, lower bound, upper bound). It can include 
        expressions or numbers. The functional interdependence 
        of the intervals is sorted out internally, and so for the 
        example above one can equivalently pass (y,x^2,200),(x,0,10). 
        The integration variable names must match the function argument 
        names. 
    -   ``verbose`` -- integer between 0 and 3 (default:1). Sets the code
        verbosity level. All critical warnings are displayed when 
        ``verbose`` >=1. The calculation proceeds silently for ``verbose`` =0.
    -   ``time_limit`` -- number (default:5; but for fast 
        machines setting this to more than 1 may not result in 
        improvements) of seconds to attempt the integral after any
        simplifications are finished (see ``simplify`` below). Note that 
        this time limit is only approximate.
    -   ``dimension_limit`` -- an integer (default:0). The 
        lowest dimension (number of integration variables) to which the 
        symbolic integrator should attempt to reduce the integral. 
    -   ``simplify_func`` -- boolean (default: False). If 
        set to True, the code tries to symbolically simplify the 
        integrand after the transformation of variables to the unit
        hypercube. This could slow down the calculation significantly.
    -   ``sigmoid`` -- a list of sigmoid functions 
        (default: ``[sigmoid_logistic]``). 
        The sigmoid
        functions must be monotonic and have limits of 0 (for argument going
        to -Infinity) and 1 (for argument going to +Infinity). They are used
        to map improper integrals to integrals on the unit hypercube. Not 
        used for finite integration ranges. If more than one function is 
        supplied, then each function will be applied
        separately as needed to the corresponding integration dummy variable in 
        the ``ranges`` list.  Thus, both lists must have the same length.
        The module supplies a few sigmoid functions: ``sigmoid_gd``, 
        ``sigmoid_logistic``, ``sidmoid_atan``, ``sigmoid_tanh``.
        If one function is supplied, then that is applied to all dummy
        variables. Not used if ``use_limits`` is set to ``True``.
    -   ``dummy_var_prefix`` -- a string (default: ``_X_``). It is used
        when the calculation cannot perform all integrals. It is the 
        prefix added to the dummy integration variable names for
        the left-over integrals. See ``ivars`` below.
    -   ``use_limits`` -- a bool (default: ``False``). Whether to use
        a limit in handling improper integrals. If set to ``False``, then
        the sigmoidal map of the integration ranges to the unit hypercube
        is used. See ``sigmoid`` above.
    -   ``algorithm`` -- to be passed on to ``integrate``. One of ``Default``,
        ``maxima``, ``sympy``, ``giac``.
    
    OUTPUT:
        
    If all integrals are successfully performed, the function returns:
        
        ``integration result``
        
        where:
        
        ``integration result`` is an expression or a number.
    
    If the only some of the integrals can be performed, then the 
    result is a tuple containing: 
        
        ``[integration result, ivars]``

        where:

        ``integration result`` is the result of the partial integration. 
        The dummy variables of the remaining integrals to be performed are
        listed in ``ivars``.
        
        ``ivars`` contains the variables over which one still has to 
        integrate the result. The range of integration for those variables
        is the unit hypercube, so it is [0,1] for each. Each dummy variable 
        name is of the form ``dummy_var_prefix+str(i)``, with ``i`` some
        integer.
    
        
    
    EXAMPLES:
    
    Let us integrate `x\times y+z` over the interval `(0,2)` in all variables::
    
        sage: from multidim_integration.symbolic_definite_integration import symbolic_multidim_integral
        sage: y,z = var('y z')
        sage: symbolic_multidim_integral(x*y+z,(x,0,2),(y,0,2),(z,0,2))
        16
    
    The order in which the integration ranges are supplied does not matter.
    The integrator sorts the functional dependence of the dummy variables
    internally::
        
        sage: from multidim_integration.symbolic_definite_integration import symbolic_multidim_integral
        sage: y=var('y')
        sage: symbolic_multidim_integral(x,(x,0,y),(y,0,1))
        1/6
        sage: symbolic_multidim_integral(x,(y,0,1),(x,0,y))
        1/6
    
    Here is an example, where we supply a function::
    
        sage: from multidim_integration.symbolic_definite_integration import symbolic_multidim_integral
        sage: y,z = var('y z')
        sage: f(x,y)=x*sin(y)
        sage: symbolic_multidim_integral(f,(y,0,200*x),(x,0,1))
        -1/40000*cos(200) - 1/200*sin(200) + 20001/40000
        
        
    If one of the integration intervals is of measure zero, then the result is always zero::
    
        sage: symbolic_multidim_integral(sqrt(-y),(x,0,1),(y,2*x,2*x)) 
        0
    
    Let us find the value of `\pi`::
        
        sage: y,z = var('y z')
        sage: symbolic_multidim_integral(1,(x,-1,1),(y,-sqrt(1-x^2),sqrt(1-x^2))) 
        pi
        sage: 3/4*symbolic_multidim_integral(1,(x,-1,1),(y,-sqrt(1-x^2),sqrt(1-x^2)),(z,-sqrt(1-x^2-y^2),sqrt(1-x^2-y^2)))
        pi
    
    Here is an improper 1-D integral::
    
        sage: symbolic_multidim_integral(exp(-x^2),(x,0,Infinity))
        1/2*sqrt(pi)

    And here is an improper integral using a different sigmoid function to map
    the infinite range to the unit interval::
    
        sage: from multidim_integration.symbolic_definite_integration import sigmoid_atan,sigmoid_tanh,sigmoid_logistic
        sage: symbolic_multidim_integral(exp(-x^2-x),(x,0,Infinity),sigmoid=[sigmoid_atan],simplify_func=True,use_limits=False)
        -1/2*sqrt(pi)*(erf(1/2) - 1)*e^(1/4)
        sage: symbolic_multidim_integral(exp(-y^2-x^2),(x,0,Infinity),(y,0,Infinity),sigmoid=[sigmoid_atan],simplify_func=True,use_limits=False)
        1/4*pi
        sage: symbolic_multidim_integral(exp(-y^2-x^2),(x,0,Infinity),(y,0,Infinity),sigmoid=[sigmoid_logistic],simplify_func=True,use_limits=False)
        1/4*pi
        sage: symbolic_multidim_integral(exp(-y^2-x^2),(x,0,Infinity),(y,0,Infinity),sigmoid=[sigmoid_tanh],simplify_func=True,use_limits=False)
        1/4*pi

    And here we use two different sigmoidal functions to help the integrator::
    
        sage: from multidim_integration.symbolic_definite_integration import symbolic_multidim_integral
        sage: from multidim_integration.symbolic_definite_integration import sigmoid_atan
        sage: from multidim_integration.symbolic_definite_integration import sigmoid_logistic
        sage: y = var('y')
        sage: symbolic_multidim_integral(exp(-x^2-y),(y,x,Infinity),(x,0,Infinity),sigmoid=[sigmoid_logistic,sigmoid_atan],use_limits=False)
        -1/2*sqrt(pi)*erf(1/2)*e^(1/4) + 1/2*sqrt(pi)*e^(1/4)
        
    The same result can be achieved by using limits and avoid performing
    the improper integral::
        
        sage: limit(symbolic_multidim_integral(exp(-x^2-y),(y,x,t),(x,0,t)),t=Infinity)
        -1/2*sqrt(pi)*erf(1/2)*e^(1/4) + 1/2*sqrt(pi)*e^(1/4)
        
    Or alternatively we can use the internal limit function which creates a dummy ``_X__INF`` which needs to be sent to infinity::
    
        sage: symbolic_multidim_integral(exp(-x^2-y),(y,x,Infinity),(x,0,Infinity),use_limits=True)
        -1/2*sqrt(pi)*erf(1/2)*e^(1/4) + 1/2*sqrt(pi)*erf(_X__INF + 1/2)*e^(1/4) - 1/2*sqrt(pi)*erf(_X__INF)*e^(-_X__INF)
        sage: limit(-1/2*sqrt(pi)*erf(1/2)*e^(1/4) + 1/2*sqrt(pi)*erf(_X__INF + 1/2)*e^(1/4) - 1/2*sqrt(pi)*erf(_X__INF)*e^(-_X__INF),_X__INF=Infinity)
        -1/2*sqrt(pi)*erf(1/2)*e^(1/4) + 1/2*sqrt(pi)*e^(1/4)

    Let us do an improper 3-dimensional integral:: 
    
        sage: a,b,c,t=var('a b c t')
        sage: assume(a>0,c>0,b>0)
        sage: symbolic_multidim_integral(exp(-x^2/2/a^2-y^2/2/b^2-z^2/2/c^2)/(2*pi)^(3/2)/(a*b*c),(x,-t,t),(y,-t,t),(z,-t,t))
        erf(1/2*sqrt(2)*t/a)*erf(1/2*sqrt(2)*t/b)*erf(1/2*sqrt(2)*t/c)
        sage: limit(erf(1/2*sqrt(2)*t/a)*erf(1/2*sqrt(2)*t/b)*erf(1/2*sqrt(2)*t/c),t=Infinity)
        1
        sage: # Using explicitly infinite limits in the integral does not produce a closed-form expression.
        sage: # The reason is that we map infinite ranges to the unit hypercube through a non-linear map, which may produce
        sage: # intractable integrands.
        
    This is an example of an integral which cannot be fully performed.
    The result has one remaining integration to be done in the dummy
    variable ``_X_2``.
    
        sage: symbolic_multidim_integral(x*sin(y*z),(x,0,1),(y,0,2000*x),(z,0,10))
        [1/4000000*(2000000*_X_2^2 - 2000*_X_2*sin(2000*_X_2) - cos(2000*_X_2))/_X_2^3 + 1/4000000/_X_2^3, [_X_2]]
        
    When the integral is only partially successful, we can then integrate
    the remaining integrals numerically::
    
        sage: symbolic_multidim_integral(x*sin(x^x*y),(y,0,x),(x,0,1))
        [-(_X_1^(-_X_1 - 1)*cos(_X_1^(_X_1 + 1)) - _X_1^(-_X_1 - 1))*_X_1^2, [_X_1]]
        sage: # The remaining 1-D integral can be done numerically:
        sage: numerical_integral(-(_X_1^(-_X_1 - 1)*cos(_X_1^(_X_1 + 1)) - _X_1^(-_X_1 - 1))*_X_1^2,0,1)
        (0.10224975449206862, 1.1352003169936412e-15)
        
    Sometimes the default internal symbolic integrator may require
    additional assumptions that will be spit out before the output as 
    warnings if the verbosity level is set higher than 1. Note warnings are spit
    out by each integration thread. Since only one of those threads will provide
    a final answer, some warnings are not applicable::
    
        sage: from multidim_integration.symbolic_definite_integration import symbolic_multidim_integral
        sage: y,z = var('y z')
        sage: symbolic_multidim_integral(x*y/z^2,(x,0,1),(y,x,2*x),(z,2*y,5.43*x^2),verbose=3) # LONG OUTPUT
    
    To make progress, one can do the multiple integral, one integral at a time,
    feeding assumptions at each step by using ``assume()``. Or one can try
    a different integration algorithm, which may be more successful::
    
        sage: symbolic_multidim_integral(x*y/z^2,(x,0,1),(y,x,2*x),(z,2*y,5.43*x^2),algorithm='sympy')
        31/1086
        sage: symbolic_multidim_integral(x*y/z^2,(x,0,1),(y,x,2*x),(z,2*y,5.43*x^2),algorithm='giac')
        31/1086
    
    
    TESTS:
    
    Make sure the integral over a constant integrand works::
    
        sage: from sage.libs.cuba.cuba import *
        sage: y,z = var('y z')
        sage: symbolic_multidim_integral(1,(x,0,2),(y,0,2),(z,0,2)) 
        8
        sage: # This also reduces to a constant after the transformation of 
        sage: # variables to the unit cube. To see that explicitly, run 
        sage: # with verbose=3
        sage: symbolic_multidim_integral(1,(x,0,2),(y,0,2),(z,x,2+x))
        8

    
    Make sure the 1-dim integral also works::
    
        sage: symbolic_multidim_integral(x,(x,0,2))
        2

    

        
    ALGORITHM:
    
    Here is a sketch of the algorithm.
    
    1. The code transforms the variables of integration to new 
    variables over the unit hypercube. It computes the Jacobian of 
    the transformation as well.
    
    2. If the integrand is constant over the unit hypercube, the 
    constant is returned as it is the result of the integration.
    
    3. Copies of the integral are then dispatched to multiple threads. 
    Each thread attempts to perform the multiple integrals in different
    order. Whichever calculation succeeds first, returns the result. If 
    not all integrals can be done, then the code returns the result with 
    fewest number of remaining integrals to be done. The dummy variables of 
    integration in the remaining integrals are labeled with ``_X_i``. 
    The prefix ``_X_`` can be changed by supplying a 
    custom ``dummy_var_prefix``.
    
    
    .. NOTE::
    
        -   When modifying the code, make sure that your 
            modifications do not result in speed regressions.
        
    
    AUTHORS:
    
    - Svetlin Tassev (2013-2023)
        
   '''
    
    ranges = list(ranges)
    
    r_vars = [v[0] for v in ranges]
    num_vars = len(r_vars)
    if (len(sigmoid)!=1 and len(sigmoid)!=num_vars):
        raise Exception("sigmoid should be a list with 1 or "+str(num_vars)+" elements.")
    
    # Make sure this is a float and not some symbolic number
    # Needed for when we check for constant integrand; or for passing
    # to the compiled function, when the jacobian is not a constant.

    
    new_vars=[]
    for i in range(len(ranges)):
        new_vars.append(var(dummy_var_prefix+str(i)))

    inf_var=var(dummy_var_prefix+'_INF')
    
    # Need to map integration range onto the unit hypercube for Cuba
    jac,tosub=_jacobian_and_substitution(ranges,new_vars,simplify_func,inf_var,sigmoid=sigmoid,use_limits=use_limits)
    if jac==0: # Note that _jacobian_and_substitution returning zero is a special 
               # case of zero intervals. So, the integral is zero.
        if (verbose>=2):
            print "One or more of the integration intervals is zero. Thus, the integral is zero."
        return 0
    jac=1/jac
    try:
        f_novars=func(*func.arguments()) # in case a function without arguments was passed
    except:
        f_novars=func

    if hasattr(f_novars,'subs'):
        fjac=f_novars.subs(tosub)*jac
    else:
        fjac=f_novars*jac
    
    if (simplify_func):
        try:
            fjac=fjac.full_simplify()
        except:
            fjac=simplify(fjac)
    #print jac
    #print f_novars
    #print func
    
    #handle constant case
    try:
        fjac.n() # check for an error.
        if (verbose>=2):
            print "The integrand is constant over the unit hypercube."
            print "The integral is "+str(fjac)
        return fjac
    except (AttributeError,TypeError):
        pass
    
    if (verbose>=3):
        print "The function to be integrated over the unit hypercube spanned by _X_i is:"
        print "-------------------"
        tmp=str(fjac)
        #tmp=repl("(.*)\|-->","",tmp)
        print tmp
        print "-------------------"
    
    
    if num_vars > dimension_limit:
        if verbose >= 2:
            print "Starting symbolic integration."
        #print new_vars
        #print fjac
        X = _multidim_analytical_integration_unit_cube(fjac,new_vars, time_limit, 
                                                dimlimit=dimension_limit,simplify_func=simplify_func,algorithm=algorithm,verbose=verbose)
        if (simplify_func):
            try:
                X[1]=X[1].full_simplify()
            except:   
                X[1]=simplify(X[1])
        if X[0] == 0:
            if (verbose >= 2):
                print "Symbolic calculation was successful!"
                print "The integrand is constant over the unit hypercube."
                print "The integral is "+str(X[1])
            if (use_limits):
                #print X[1]
                #print inf_var
                #print limit(X[1],inf_var=Infinity)
                X[1]=limit(X[1],inf_var=Infinity,taylor=True)
                return X[1]
            else:
                return X[1]
        else:
            num_vars=X[0]
            if (verbose >= 2):
                print "Symbolic calculation reduced dimensionality of integral to: "+str(X[0])
            if (verbose >= 3):
                print "The resulting function over the unit hypercube is:"
                print "----------------------"
                print str(X[1])
                print "----------------------" 
    else:
        X=[num_vars,fjac]
    
    lv = list(X[1].variables())
    str_lv=[str(v) for  v in lv]
    lv_out=[]
    str_vars=[str(r) for r in r_vars]
    for i in range(len(str_lv)):
        if (dummy_var_prefix in str_lv[i]) and not('_INF' in str_lv[i]):
            lv_out.append(lv[i])
        if (str_lv[i] in str_vars):
            lv_out.append(lv[i])
    return [X[1],lv_out]


from sympy import nsimplify

def _integ(toint):
    """
    Do a 1-d symbolic integration using integrate(). 
    Take a tuple containing ``(integrand,integration variable, 
    list of free variables,algorithm,verbosity)`` .The result is either an Expression or 
    a None. If the result contains an integrate() (e.g. Maxima did 
    not succeed) it returns a None, as it is useless for numerics, 
    and this allows for automatically terminating further 
    integration attempts for this particilar variable order by its 
    caller ``_multidim_analytical_integration_unit_cube`` . See that 
    function for more. This is tested by all the examples in the 
    rest of this file which use symbolic integration. Here is one 
    test for the sake of it::
    
        sage: from multidim_integration.symbolic_definite_integration import _integ
        sage: y,z = var('y z')
        sage: _integ((x^2*y*sin(z),x,[y,z]))
        1/3*y*sin(z)
    
    
    """
    for v in toint[2]: # The analytical integration is done over the unit hypercube, 
                       # so make that explicit to Sage:
        assume(v > 0)
        assume(v < 1)
    
    algorithm=toint[3]
    
    l = len(toint[0].variables())
    
    try:
        if algorithm=='Default':
            r = integrate(nsimplify(toint[0]),toint[1],Integer(0),Integer(1))
        else:
            r = integrate(nsimplify(toint[0]),toint[1],Integer(0),Integer(1),algorithm=algorithm)
    except Exception as e:
        if toint[4]>1:
            print e
        return None
        
    if not isinstance(r,Expression):
        return r
    
    if l > len(r.variables()): # This check ensures that r does not contain 
                               # integrate() expressions. 
                               # Those are useless for numerics, and moreover, 
                               # if passed through the 
                               # pools, result sporadically in segfaults!!!
            return r
    else:
        return None



########################################################################
    
    

def  _multidim_analytical_integration_unit_cube(func,ivars, timelimit, dimlimit=0,simplify_func=False,algorithm='Default',verbose=0):
    """ 
    Return the analytical multidimensional integral of a 
    function over the unit hypercube. 
    
    May return only a partially integrated function, if the integral 
    cannot be integrated in all integration variables, either in 
    principle (by Maxima), or within the ``timelimit`` and 
    ``dimlimit`` .
    
    This function should be considered as a preprocessor for 
    ``multidim_integration``.
    
    INPUT:
    
    -   ``func`` -- the integrand, which must be an Expression. The 
        integral is assumed to be over all ``func.variables()`` of 
        ``func``, and therefore the result ultimately should be numeric. 
        However, the analytical integral is not always possible to be 
        reduced to a number. In such cases, the result is a function 
        which is to be integrated over the unit hypercube. In many 
        cases, the resulting function depends on fewer variables than 
        ``func``, and therefore the dimensionality of the integral is 
        reduced.
    
    -   ``ivars`` -- a list containing the variables to be integrated over
        the unit hypercube.
    
    -   ``timelimit`` -- a number. It gives the number of seconds 
        integration should be attempted before a result is returned. 
   
    -   ``dimlimit``  -- an integer (default:0). It gives the 
        minimum number of integration variables to which the integrator 
        should try to reduce the integrand. See ``multidim_integration``.
   
    -   ``simplify_func`` -- a bool (default:False). Whether to try
        to simplify intermediate results.
   
    OUTPUT: 
    
    A tuple containing:
    
    ``(ndim,result)``
    
    where ``result`` is either a number (for which ``ndim`` =0), in 
    which case the integral has been fully performed; or an 
    Expression, which is to be integrated over the unit hypercube in 
    all ``result.variables()`` . The latter result is obtained when 
    only some of the integrals are possible to do analytically by 
    Sage.
    
    
    .. WARNING::
    
        -   The second warning in the Warning section of 
            ``multidim_integration`` applies to this function as well.
        -   This function should mostly be used as a low-level 
            version of ``multidim_integration``. Use that function 
            instead as most debugging checks are done there. Also, note 
            that the examples in ``multidim_integration`` test this 
            function as well.
   
    ALGORITHM:
      
    Take a function and try to integrate it symbolically on the 
    unit hypercube within timelimit by running many forks in 
    parallel. Each fork has a different integration variable order. 
    We do that because whether an integration succeeds or not and how fast 
    depends on that order. Compare the time it takes to do:
    
    integrate(integrate(2000*y^2*sin(2000*x*y),(y,0,1)),(x,0,1)) 
    
    to the time it takes to perform:
    
    integrate(integrate(2000*y^2*sin(2000*x*y),(x,0,1)),(y,0,1)) 
    
    
    EXAMPLES::
        
        sage: from multidim_integration.symbolic_definite_integration import _multidim_analytical_integration_unit_cube
        sage: y,z,w=var('y z w')
        sage: f=2000*y^2*sin(2000*x*y)
        sage: # May fail for very slow computers within set timelimit:
        sage: _multidim_analytical_integration_unit_cube(f,[x,y],10) # long time (1sec). 
        (0, -1/4000000*cos(2000) - 1/2000*sin(2000) + 2000001/4000000)
        sage: # Integrating just in y
        sage: _multidim_analytical_integration_unit_cube(f,[y],10)
        (0, -1/2000000*((2000000*x^2 - 1)*cos(2000*x) - 2000*x*sin(2000*x))/x^3 - 1/2000000/x^3)
        sage: # Integrating just in x
        sage: _multidim_analytical_integration_unit_cube(f,[x],10)
        (0, -y^2*(cos(2000*y)/y - 1/y))
        sage: # For the next example, one should get one of the answers below depending on
        sage: # which proceeds faster -- the integral in x or in y. The result below is 
        sage: # twice as fast as the one above as we do not reduce the integral to zero 
        sage: # integration variables. 
        sage: _multidim_analytical_integration_unit_cube(f,[x,y],10,dimlimit=1) # random output long time
        (1, -y^2*(cos(2000*y)/y - 1/y))
        (1, -1/2000000*((2000000*x^2 - 1)*cos(2000*x) - 2000*x*sin(2000*x))/x^3 - 1/2000000/x^3)
        sage: f=20000*y^2*sin(20000*x*y*z)
        sage: # Output below may have x instead of z, depending on which fork 
        sage: # finishes first.
        sage: _multidim_analytical_integration_unit_cube(f,[x,y,z],10) # random output long time
        (1, 1/400000000*(200000000*z^2 - 20000*z*sin(20000*z) - cos(20000*z))/z^3 + 1/400000000/z^3)
        sage: # The calculation below is 100x faster than the above example, as the  
        sage: # integration terminates immediately when the integrand is reduced to 
        sage: # one integration variable. 
        sage: _multidim_analytical_integration_unit_cube(f,[x,y,z],10,dimlimit=1) # random output, long time
        (1, 1/400000000*(200000000*x^2 - 20000*x*sin(20000*x) - cos(20000*x))/x^3 + 1/400000000/x^3)
        sage: f=20000*y^2*sin(20000*x^x*y^y*z^z)
        sage: # The result below cannot be integrated in any of its variables:
        sage: _multidim_analytical_integration_unit_cube(f,[x,y,z],10) # long time
        (3, 20000*y^2*sin(20000*x^x*y^y*z^z))
        sage: f=20000*y^2*sin(20000*x*y^y*z^z) 
        sage: # If we set dimlimit=2 below, the computation is 300 times faster as further
        sage: # reduction of the integral (below that dimlimit) cannot be done.
        sage: # May fail for very slow computers within the timelimit:
        sage: _multidim_analytical_integration_unit_cube(f,[x,y,z],10,dimlimit=2) # 
        (2, -y^2*(cos(20000*y^y*z^z)/(y^y*z^z) - 1/(y^y*z^z)))
        sage: f=2000*y^2*sin(2000*x*y*sin(z))
        sage: _multidim_analytical_integration_unit_cube(f,[x,y,z],10)  # long time . random output
        (1, -1/2000000*(2000000*(2*cos(2*z) - 1)*cos(4*z) - 1000000*cos(4*z)^2 - 4000000*cos(2*z)^2 + 1000*(cos(4*z) - 2*cos(2*z) + 1)*cos(3*z + 2000*sin(z)) - (cos(4*z) - 2*cos(2*z) + 1)*cos(2*z + 2000*sin(z)) - 1000*(cos(4*z) - 2*cos(2*z) + 1)*cos(z + 2000*sin(z)) + 1000*(cos(4*z) - 2*cos(2*z) + 1)*cos(-z + 2000*sin(z)) - (cos(4*z) - 2*cos(2*z) + 1)*cos(-2*z + 2000*sin(z)) - 1000*(cos(4*z) - 2*cos(2*z) + 1)*cos(-3*z + 2000*sin(z)) - 1000000*sin(4*z)^2 + 4000000*sin(4*z)*sin(2*z) - 4000000*sin(2*z)^2 + 1000*(sin(4*z) - 2*sin(2*z))*sin(3*z + 2000*sin(z)) - (sin(4*z) - 2*sin(2*z))*sin(2*z + 2000*sin(z)) - 1000*(sin(4*z) - 2*sin(2*z))*sin(z + 2000*sin(z)) - 1000*(sin(4*z) - 2*sin(2*z))*sin(-z + 2000*sin(z)) + (sin(4*z) - 2*sin(2*z))*sin(-2*z + 2000*sin(z)) + 1000*(sin(4*z) - 2*sin(2*z))*sin(-3*z + 2000*sin(z)) + 4000000*cos(2*z) - 1000000)/(cos(4*z)^2*sin(z) + 4*cos(2*z)^2*sin(z) + sin(4*z)^2*sin(z) - 4*sin(4*z)*sin(2*z)*sin(z) + 4*sin(2*z)^2*sin(z) - 2*(2*cos(2*z)*sin(z) - sin(z))*cos(4*z) - 4*cos(2*z)*sin(z) + sin(z)) - 1/1000000*(cos(4*z)*cos(2*z) - 2*cos(2*z)^2 + sin(4*z)*sin(2*z) - 2*sin(2*z)^2 + cos(2*z))/(cos(4*z)^2*sin(z) + 4*cos(2*z)^2*sin(z) + sin(4*z)^2*sin(z) - 4*sin(4*z)*sin(2*z)*sin(z) + 4*sin(2*z)^2*sin(z) - 2*(2*cos(2*z)*sin(z) - sin(z))*cos(4*z) - 4*cos(2*z)*sin(z) + sin(z)))
        

    .. NOTE::
    
        -   When modifying the code, make sure that the function can 
            return numerical answers (if any) as they become available, 
            not waiting for all forks to finish. See, for example, the 
            first example above. Also, make sure all forks are killed in 
            that case before the function returns. And test extensively 
            to make sure that no segfaults appear randomly especially by 
            issuing one and the same command in rapid succession.
        
        -   When modifying the code, make sure that your 
            modifications do not result in speed regressions.
    
    
    
    AUTHORS:
    
    - Svetlin Tassev (2013-2023)
    
        
    """
    
    
    
    # _name__ == '__main__':

     
    result=func
    n_vars = len(ivars)
    
    if n_vars <= dimlimit:
       return [n_vars, result]
    
    it = []
    pool = []
    ns = []
    v_arr = []

    start = time()

    try:
        pool.append( Pool(processes=n_vars) )
    
        it.append( pool[0].imap_unordered(_integ, zip([func]*n_vars,ivars,[ivars]*n_vars,[algorithm]*n_vars,[verbose]*n_vars)) )
        ns.append( n_vars )
        v_arr.append( n_vars )

        nthreads = 1
        next_thread = 0
        nmax = n_vars
        
        while (nthreads > 0) and ((time() - start) < timelimit):
            if ns[next_thread] > 0:
                try:
                    func = it[next_thread].next(timeout = 0.01)
                    
                    
                    ns[next_thread] -= 1
                    if func != None:
                        
                        if not isinstance(func, Expression):
                            func.n() # will trigger exception if it fails
                            result = func
                            nthreads = 0
                            n_vars = 0
                            break
                        
                        if (simplify_func):
                            try:
                                func=func.full_simplify()
                            except:
                                func=simplify(func)
                            
                        ivars1 = set(func.variables()).intersection(set(ivars)) 
                        n_vars1 = len(ivars1)



                        if n_vars1 <= dimlimit:
                            result = func
                            n_vars = n_vars1
                            nthreads = 0
                            break
                        
                        nthreads += 1
                        if nmax > n_vars1:
                            nmax = n_vars1
                            n_vars = n_vars1
                            result = func
                           
                        
                        v_arr.append( n_vars1 )
                        pool.append( Pool(processes=n_vars1) )
                        it.append( pool[nthreads-1].imap_unordered(_integ, zip([func]*n_vars1,ivars1,[ivars1]*n_vars1,[algorithm]*n_vars1,[verbose]*n_vars1)) )
                        ns.append( n_vars1 )
                        
                    next_thread = (next_thread+1) % nthreads
                except (TimeoutError):
                    next_thread = (next_thread+1) % nthreads
            else:
                pool.pop(next_thread)
                it.pop(next_thread)
                ns.pop(next_thread)
                v_arr.pop(next_thread)
                nthreads -= 1 
                if nthreads==0:
                    break
                next_thread = next_thread % nthreads
    finally:
        for i in pool:
            if not isinstance(i,int):
                 i.terminate()
                 i.join()
    return [n_vars, result]

