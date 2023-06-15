# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 10:15:34 2023

@author: txjam
Numerical Methods for finding Roots
"""
import numpy as np
import matplotlib.pyplot as plt
#%%
def bisrf(F, xl,xu, eps = 0.01, maxit = 100):
    '''
    Root finder by bisection
    
    Inputs:
    F = function evaluated
    xl = lower bracket
    xu = upper bracket
    eps = interval size desired
    iter = iterations
    
    
    '''
    err = xu - xl
    iter = 0
    detect = False
    if maxit > iter:
        while err > eps:
                xm = (xl+xu)/2
                '''
                Test if [xl,xm] is good bracket
                If yes new bracket is [xl,xm] if not then [xm,xu]
                '''
                if F(xl) * F(xm) <= 0 :
                    xu = xm
                elif F(xm) * F(xu) < 0:
                    xl = xm
                else:
                    print('Error: No root between these two points')
                    detect = True
                    break
                err = xu - xl
                iter = iter + 1
        if detect == False:        
            print(f'The root approx. is {(xu+xl)/2} and it took {iter} iterations.')
    return(xl,xu, iter)
#%%
#Example function
def f(x):
    return x**3 - 5*x**2 + x + 2.5
def df(x):
    return 3*x**2 - 10*x + 1
xc = np.linspace(-2.5,7) #Creating vector for x-coordinates
plt.plot(xc, f(xc))
plt.plot(xc, np.zeros(xc.shape)) #plotting x = 0
bisrf(f,0,2)
#%%
def falseposrf(F, xl,xu, eps = 0.01, maxiter = 100):
    '''
    Root finder by False position
    
    Inputs:
    F = function evaluated
    xl = lower bracket
    xu = upper bracket
    eps = interval size desired
    iter = iterations
    
    
    '''
    err = xu - xl
    iter = 0
    detect = False
    if maxiter > iter:
        while err > eps:
                xr = xu - (F(xu)*(xl-xu))/(F(xl)-F(xu))
                '''
                Test if [xl,xm] is good bracket
                If yes new bracket is [xl,xm] if not then [xm,xu]
                '''
                if F(xl) * F(xr) <= 0 :
                    xu = xr
                elif F(xr) * F(xu) < 0:
                    xl = xr
                else:
                    print('Error: No root between these two points')
                    detect = True
                    break
                err = xu - xl
                iter = iter + 1
        if detect == False:        
            print(f'The root approx. is {(xu+xl)/2} and it took {iter} iterations.')
    return(xl,xu, iter)
#%%
falseposrf(f, 0, 2)
xc1 = np.linspace(0.9,0.92) #Creating vector for x-coordinates
plt.plot(xc1, f(xc1))
plt.plot(xc1, np.zeros(xc1.shape)) #plotting x = 0
#Much slower than bisection for this function
#%%
def fixprf(F,eps = 0.05, x0 = 0., maxit = 100.):
    '''
    Fixed Point iteration method for finding roots
    Turn f(x) = 0 into g(x) = x
    Note that g'(x) must be <1 or this method will diverge
    F = function
    eps = error threshold
    x0 = starting value
    it = current iteration number
    maxit = maximum # of iterations
    '''
    #initialize
    x1 = F(x0)
    err = np.abs(x1 - x0)
    it = 0
    while err > eps and it >= maxit:
        x1 = F(x0)
        err = np.abs(x1 - x0)
        x0 = x1
        it += 1
    print(f'Root estimate = {x1}. This took {it} iterations.')
#%%
#Given function x**4 - x - 10 use:
def F(x):
    return x**4 - x - 10
def g(x):
    return (x + 10)**(1/4)
fixprf(g)
xc2 = np.linspace(1.7,1.9)
plt.plot(xc2, F(xc2))
plt.plot(xc2, np.zeros(xc2.shape)) #plotting x = 0
#%%
def newtrf(F,dF, eps = 0.05, x0 = 1., maxit = 100.):
    '''
    Newton method for finding roots
    F = function
    dF = derivative of function
    eps = error threshold
    x0 = starting value
    maxit = maximum # of iterations
    '''
    #initialize
    it = 0
    x1 = x0 - F(x0) / dF(x0)
    err = np.abs(x1 - x0)

    while err > eps and it >= maxit:
        x1 = x0 - F(x0)/dF(x0)
        err = np.abs(x1 - x0)
        x0 = x1
        it += 1
    print(f'Root estimate = {x1}. This took {it} iterations.')
    return x1,it
#%%
newtrf(f, df)
#Consistent with previous results
#%%
def secrf(F, eps = 0.05, x0 = 1, x1 = 2, maxit = 20.):
    '''
    Secant Method for root approx.
    Does require two starting values

    F = function
    eps = error threshold
    x0 = starting value 1
    x1 = starting value 2
    x2 = approx. root value
    it = iteration #
    maxit = maximum # of iterations
    '''
    #initialize
    it = 0
    err = 2*eps

    while err > eps and it < maxit:
        '''
        print(f' x0 in {it}iteration = {x0}')
        print(f' x1 in {it}iteration = {x1}')
        Used for testing
        '''

        x2 = x1 - (F(x1)*(x1-x0) / (F(x1)-F(x0) ))
        #print(f' x2 in {it}iteration = {x2}')

        err = np.abs(x2 - x1)
        #print(f' err in {it}iteration = {err}')

        x0 = x1
        x1 = x2
        it += 1
    print(f'Root estimate = {x1}. This took {it} iterations.')
    return x2,it
#%%
secrf(f)
