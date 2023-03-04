
#*****************************************************************************
#       Copyright (C) 2013-2023 Svetlin Tassev <stassev@alum.mit.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.all import *
from multidim_integration.symbolic_definite_integration import *

def numerical_multidim_integral(func, *ranges,xl_embed=[],xu_embed=[],
                         symbolic=True,
                         verbose=1, 
                         dimension_limit=0, time_limit=5,sigmoid=[sigmoid_logistic],
                         simplify_func=False,algorithmS='Default',
                         calls=1e4,
                         algorithmN='vegas',
                         algorithm1='qags',eps_abs=1.e-6,eps_rel=1.e-6,rule=6):
    r'''
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
        
    -   ``ranges`` -- the integration intervals in each variable as 
        tuples, e.g. (x,0,10),(y,x^2,200). Each interval is of the form 
        (variable name, lower bound, upper bound). It can include 
        expressions or numbers. The functional interdependence 
        of the intervals is sorted out internally, and so for the 
        example above one can equivalently pass (y,x^2,200),(x,0,10). 
        The integration variable names must match the function argument 
        names. 
        
    -   ``x[l|u]_embed`` -- lists of [lower|upper] embedding limits. In 
		cases when the integrand is not transformed internally to the unit
		hypercube, one needs to specify an outer hypercube in which the
		intervals in ``ranges`` reside. So, for example, for a ranges 
		(x,0,1),(y,0,x), one may need to to supply xl_embed=[0,0], 
		xu_embed=[1,1], which are limits which contain the integration range.
		The function will notify the user when x*_embed are needed.
        
    -   ``symbolic`` -- bool (default: True). Whether to attempt symbolic
        integration. Even if symbolic integration cannot integrate over all
        integration variables, it may be able to perform some of them, thus
        reducing the dimensionality of the numerical integral. Highly recommended.
    
    -   ``verbose``, ``dimension_limit``, ``time_limit``, ``sigmoid``,
        ``simplify_func``, ``algorithmS``: parameters to be passed on to ``symbolic_multidim_integral``. 
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
	
	Having a lower bound above the upper bound is also treated correctly::
	
		sage: numerical_multidim_integral(x,(x,1,0))
		(-0.500000000000000, 0.0)
		sage: numerical_multidim_integral(x,(x,0,1))
		(0.500000000000000, 0.0)
		sage: numerical_multidim_integral(x,(x,1,0),xl_embed=[0],xu_embed=[1],symbolic=False)
		(-0.49999990229069513, 2.58729524973818e-07)
		sage: numerical_multidim_integral(x,(x,0,1),xl_embed=[0],xu_embed=[1],symbolic=False)
		(0.49999990229069513, 2.58729524973818e-07)

    
    AUTHORS:
    
    - Svetlin Tassev (2013-2023)
        
        
    '''
    if (symbolic):
        try: 
            f=symbolic_multidim_integral(func,*ranges,dummy_var_prefix='_X_',use_limits=False,verbose=verbose, 
                            dimension_limit=dimension_limit, time_limit=time_limit,sigmoid=sigmoid,
                            simplify_func=simplify_func,algorithm=algorithmS)
            if f==None:
                raise
            if type(f)!=list:
                return (f.n(),0.0)
            else:
                if len(f[1])>1:
                    return monte_carlo_integral(f[0],[0]*len(f[1]),[1]*len(f[1]),calls=calls,algorithm=algorithmN)
                else:
                    return numerical_integral(f[0],0,1,algorithm=algorithm1,eps_abs=eps_abs,eps_rel=eps_rel,rule=rule,max_points=calls)
        except:
            None
    
    xl=[v[1] for v in ranges]
    xu=[v[2] for v in ranges]
    
    # Check if ranges are constant. Do that by attempting the integral.
    
    try:
        if len(ranges)==1:
            numerical_integral(f[0],xl[0],xu[0],algorithm=algorithm1,eps_abs=eps_abs,eps_rel=eps_rel,rule=rule,max_points=calls) # This should take care of all remaining 1D integrals.
        else:
            return monte_carlo_integral(func,xl,xu,calls=calls,algorithm=algorithmN)
    except:
        None

    # Construct integrand embedded in x[l|u]_embed that is 0 outside the
    # integration ranges. Use heaviside(x) to do that.
    
    var_list=[v[0] for v in ranges]
    H=1
    for i in range(len(var_list)):
        H*=(heaviside(var_list[i]-xl[i])*heaviside(-var_list[i]+xu[i])    # If xl<x<xu
			- heaviside(-var_list[i]+xl[i])*heaviside(var_list[i]-xu[i])) # If xl>x>xu, hence the minus sign
    
    if len(xl_embed)==0 or len(xu_embed)==0:
        raise Exception("xl_embed and xu_embed must be set.")
    
    try:
        return monte_carlo_integral(H*func(*var_list),xl_embed,xu_embed,calls=calls,algorithm=algorithmN)
    except:
        return monte_carlo_integral(H*func,xl_embed,xu_embed,calls=calls,algorithm=algorithmN)
    


