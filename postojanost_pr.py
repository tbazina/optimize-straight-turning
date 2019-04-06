# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 19:05:57 2015

@author: tomislav
"""

import numpy as np
from scipy.optimize import minimize, differential_evolution

def postojanost(x, v, f, ap):
    return x[0]/(v**x[1])/(f**x[2])/(ap**x[3])

def fun(x, t, v, f, ap):
    sum_n = t*(v**x[1])*(f**x[2])*(ap**x[3])/x[0]
    sum_m = (np.sum(sum_n, 1) - 1.)**2.
    return np.sum(sum_m)
    
def fun_jac(x, t, v, f, ap):
    sum_n0 = -t*(v**x[1])*(f**x[2])*(ap**x[3])/(x[0]**2.)
    sum_m0 = 2.*np.sum(sum_n0, 1)
    dfdx0 = np.sum(sum_m0)
    
    sum_n1 = t*(v**x[1])*np.log(v)*(f**x[2])*(ap**x[3])/x[0]
    sum_m1 = 2.*np.sum(sum_n1, 1)
    dfdx1 = np.sum(sum_m1)

    sum_n2 = t*(v**x[1])*(f**x[2])*np.log(f)*(ap**x[3])/x[0]
    sum_m2 = 2.*np.sum(sum_n2, 1)
    dfdx2 = np.sum(sum_m2)

    sum_n3 = t*(v**x[1])*(f**x[2])*(ap**x[3])*np.log(ap)/x[0]
    sum_m3 = 2.*np.sum(sum_n3, 1)
    dfdx3 = np.sum(sum_m3)
    
    return np.array([dfdx0, dfdx1, dfdx2, dfdx3])
    
def postojanost_min(x, v, f, ap):
    return x[0]/(v**x[1])/(f**x[2])/(ap**x[3])
    
def postojanost_min_jac(x, v, f, ap):
    dfdx0 = 1./(v**x[1])/(f**x[2])/(ap**x[3])
    
    dfdx1 = -x[0]*np.log(v)/(v**x[1])/(f**x[2])/(ap**x[3])
    
    dfdx2 = -x[0]*np.log(f)/(v**x[1])/(f**x[2])/(ap**x[3])
    
    dfdx3 = -x[0]*np.log(ap)/(v**x[1])/(f**x[2])/(ap**x[3])
    
    return np.array([dfdx0, dfdx1, dfdx2, dfdx3])
    
def postojanost_max(x, v, f, ap):
    return 500.-x[0]/(v**x[1])/(f**x[2])/(ap**x[3])
    
def postojanost_max_jac(x, v, f, ap):
    dfdx0 = -1./(v**x[1])/(f**x[2])/(ap**x[3])
    
    dfdx1 = x[0]*np.log(v)/(v**x[1])/(f**x[2])/(ap**x[3])
    
    dfdx2 = x[0]*np.log(f)/(v**x[1])/(f**x[2])/(ap**x[3])
    
    dfdx3 = x[0]*np.log(ap)/(v**x[1])/(f**x[2])/(ap**x[3])
    
    return np.array([dfdx0, dfdx1, dfdx2, dfdx3])
    
t = np.array([[43.6, 3.8], [18.38, 9.8], [18.33, 15.19], [14.87, 4.99],
              [20.22, 14.82], [24.55, 13.61], [9.76, 9.99], [10.29, 6.38]])

ap = np.array([[0.6, 0.9], [0.6, 0.9], [0.7, 0.8], [0.6, 0.9], [0.6, 0.7],
               [0.7, 0.8], [0.7, 0.6], [0.6, 0.9]])

f = np.array([[0.14, 0.14], [0.2, 0.16], [0.2, 0.14], [0.16, 0.14], [0.16, 0.2],
              [0.14, 0.16], [0.2, 0.16], [0.2, 0.14]])

v = np.array([[180, 225], [190, 200], [185, 195], [210, 220], [200, 180], 
              [190, 185], [210, 220], [225, 195]], dtype=np.float64)

bnds = ((270.**3, 270**6), (1e-6, 100.), (1e-6, 100.), (1e-6, 100.))

cons = ({'type': 'ineq', 'fun': lambda x: x[1]-x[2],
         'jac': lambda x: np.array([0., 1., -1., 0.])
         },
        {'type': 'ineq', 'fun': lambda x: x[1]-x[3],
         'jac': lambda x: np.array([0., 1., 0., -1.])
         })
for i in np.arange(np.size(t)):
    cons += ({'type': 'ineq', 'fun': postojanost_min, 'jac': postojanost_min_jac,
             'args': (np.ndarray.flatten(v)[i], np.ndarray.flatten(f)[i], 
                      np.ndarray.flatten(ap)[i])},
            {'type': 'ineq', 'fun': postojanost_max, 'jac': postojanost_max_jac,
             'args': (np.ndarray.flatten(v)[i], np.ndarray.flatten(f)[i], 
                      np.ndarray.flatten(ap)[i])})

res = minimize(fun, (270**5., 0.1, 0.3, 1.), args=(t, v, f, ap), 
               method='SLSQP', bounds=bnds, constraints=cons,
               options={'disp': True, 'maxiter': 1000, 'ftol': 1e-15})

x = np.array([res.x[0]**(1./res.x[1]), 1./res.x[1], res.x[2]/res.x[1],
              res.x[3]/res.x[1]])

resga = differential_evolution(fun, bnds, args=(t, v, f, ap), strategy='rand1bin',
                               tol=1e-12, disp=True)

xga = np.array([resga.x[0]**(1./resga.x[1]), 1./resga.x[1], resga.x[2]/resga.x[1],
              resga.x[3]/resga.x[1]])

import matplotlib.pyplot as plt

mi = 'minimum'
ma = 'maksimum'
sr = 'srednje'
arr = 'array'
              
vplt = {mi: np.min(v), ma: np.max(v), sr: np.mean(v)}
vplt[arr] = np.linspace(vplt[mi], vplt[ma], 1000)

fplt = {mi: np.min(f), ma: np.max(f), sr: np.mean(f)}
fplt[arr] = np.linspace(fplt[mi], fplt[ma], 1000)

applt = {mi: np.min(ap), ma: np.max(ap), sr: np.mean(ap)}
applt[arr] = np.linspace(applt[mi], applt[ma], 1000)

plt.figure()
plt.plot(vplt[arr], postojanost(res.x, vplt[arr], fplt[sr], applt[sr]))

plt.figure()
plt.plot(fplt[arr], postojanost(res.x, vplt[sr], fplt[arr], applt[sr]))

plt.figure()
plt.plot(applt[arr], postojanost(res.x, vplt[sr], fplt[sr], applt[arr]))

plt.figure()
plt.plot(vplt[arr], postojanost(res.x, vplt[arr], fplt[arr], applt[arr]))
plt.plot(vplt[sr], postojanost(res.x, vplt[sr], fplt[sr], applt[sr]), 'k+')
plt.plot(vplt[arr], postojanost(np.array([273.**5., 5., 0.15/0.2, 1.]),
         vplt[arr], fplt[arr], applt[arr]))

plt.plot(np.ndarray.flatten(v), np.ndarray.flatten(postojanost(res.x, v, f, ap)),
         'bo')
         
plt.plot(np.ndarray.flatten(v), 
         np.ndarray.flatten(postojanost(np.array([273.**5., 5., 0.15/0.2, 1.]),
                                        v, f, ap)), 'gs')

#plt.show()