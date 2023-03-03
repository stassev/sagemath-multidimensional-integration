# Multidimensional Integration in SageMath
This is a python package that provides both symbolic and numerical integration in [SageMath](https://www.sagemath.org/) over non-trivial integration domains.

To compile the code, run:

python setup.py build_ext --inplace 

Then you can import the functions in Sage and follow the examples in the .pyx files. The manuals for the two main functions are provided below.

**Note:** This code is a modification of the original code published here [https://github.com/sagemath/sage/issues/15492](https://github.com/sagemath/sage/issues/15492). The code is refactored a lot from that original, and no longer relies on the CUBA library but instead uses the internal SageMath numerical integration methods.

**Author**: Svetlin Tassev

**License**: GPLv3 or higher

## `symbolic_multidim_integral`

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
        (variable name, lower bound, upper bound). It can only include 
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
    
        sage: symbolic_multidim_integral(exp(-x^2),(x,0,Infinity),use_limits=False)
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
        
        
## `numerical_multidim_integral`

    Calculate the numerical multidimensional integral of a function along
    with an error estimate. Non-trivial integration domains are supported
    with ranges that depend on some of the dummy integration variables.
    
    As an option, symbolic integration can 
    be attempted before the numerical integration, which can be 
    especially useful when the numerical integration is slow to 
    converge to a requested precision. The result is always 
    converted to a floating point number. 
    
    INPUT:
    
    -   ``func`` -- the integrand. It is the same input type as that for 
        ``monte_carlo_integral`` and/or ``numerical_integral`` as those
        are called internally.
        
    -   ``symbolic`` -- bool (default: True). Whether to attempt symbolic
        integration. Even if symbolic integration cannot integrate over all
        integration variables, it may be able to perform some of them, thus
        reducing the dimensionality of the numerical integral. Highly recommended.
    
    -   ``verbose``, ``dimension_limit``, ``time_limit``, ``sigmoid``,
        ``simplify_func``: parameters to be passed on to ``symbolic_multidim_integral``. 
        See that function for help. 
        
    -   ``calls``, ``algorithmN``: parameters to be passed on to
        ``monte_carlo_integral``.
        
    -   ``algorithm1``, ``eps_abs``, ``eps_rel``, ``rule``, 
        ``calls`` (=``max_points``): parameters 
        to be passed on to ``numerical_integral``. 
        
    OUTPUT:

    A tuple whose first component is the answer and whose second
    component is an error estimate. If the error estimate is zero, then
    most probably the symbolic calculation was successful in fully 
    evaluating the integral.

    EXAMPLES:
    
        sage: from multidim_integration.numerical_integration import numerical_multidim_integral
        sage: y,z = var('y z')
        sage: numerical_multidim_integral(x*y+z,(x,0,2),(y,0,2),(z,0,2))
        (16.0000000000000, 0.000000000000000)
        
        sage: numerical_multidim_integral(1,(x,-1,1),(y,-sqrt(1-x^2),sqrt(1-x^2)),(z,-sqrt(1-x^2-y^2),sqrt(1-x^2-y^2)),xl_embed=[-1,-1,-1],xu_embed=[1,1,1],symbolic=False)
        (4.192101965800699, 0.007228186070794007)
        sage: numerical_multidim_integral(1,(x,-1,1),(y,-sqrt(1-x^2),sqrt(1-x^2)),(z,-sqrt(1-x^2-y^2),sqrt(1-x^2-y^2)),symbolic=True)
        (4.18879020478639, 0.000000000000000)
        
        sage: numerical_multidim_integral(exp(-x^2-y^2),(x,-Infinity,Infinity),(y,-Infinity,Infinity),symbolic=True)
        (3.1415926534957728, 5.064472029536488e-07)
        sage: # Using symbolic=False above results in tiny numbers, as the integrand is practically never sampled near the origin.
        
    Now let us integrate `xy/z^2` for `x` in the interval `(0,1)`, for `y`
    in ther interval `(x,2x)` and for `z` in the interval `(2y,5.43x^2)`::
    
        sage: numerical_multidim_integral(x*y/z^2,(x,0,1),(y,x,2*x),(z,2*y,5.43*x^2),calls=1e5,time_limit=0.1)
        (0.028538489254568394, 4.6051349513320865e-06)
        sage: f(x,y,z)=x*y/z^2
        sage: numerical_multidim_integral(f,(x,0,1),(y,x,2*x),(z,2*y,5.43*x^2),calls=1e6,time_limit=0.1)
        (0.028544687794419644, 1.013962345476828e-06)
        
    Integrating a python function is also possible::
        
        sage: y,z = var('y z')
        sage: def g(x,y):
        ....:     if (x^2+y^2<1.0):
        ....:         return 1.0
        ....:     else:
        ....:         return 0.0
        ....: 
        sage: numerical_multidim_integral(g,(x,-1,1),(y,-1,1))
        (3.14306762453169, 0.0014065408252111505)
