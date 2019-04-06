# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 07:51:17 2015

@author: tomislav
"""

import numpy as np
from scipy.optimize import differential_evolution, minimize, basinhopping

def postojanost(x, v, f, ap):
    return x[0]/(v**x[1])/(f**x[2])/(ap**x[3])

def vr_pop_v_n(Dv, D1, v, f, ip):
    return (Dv**2+D1**2)*np.pi/4e3/v/f*ip
    
def brzina_rezanja(promjer, brzina_vrtnje):
    """
    
    """
    return promjer*np.pi*brzina_vrtnje/1e3

def brzina_vrtnje(promjer, brzina_rezanja):
    """
    
    """
    return brzina_rezanja*1e3/promjer/np.pi
    
def strojno_vrijeme_n(duljina_prolaza, brzina_vrtnje, posmak, broj_prolaza):
    """
    
    """
    return duljina_prolaza/brzina_vrtnje/posmak*broj_prolaza
    
def promjer_rez(v, n):
    return v*1e3/np.pi/n
    
def strojno_vrijeme_pop_n(Dv, n, f, ip):
    return Dv/2./f/n*ip
    
def strojno_vrijeme_prs_n(Dv, Du, f, n, ip):
    return (Dv-Du)/2./f/n*ip
    
def vsr_n(Dv, n):
    return 0.7*Dv*np.pi*n/1e3
    
def Dsr(Dv, Du):
    return np.sqrt(Dv**2/2.+Du**2/2.)
    
class Strojno_vrijeme(object):
    def __init__(self, tip, const, *args):
        self.tip = tip
        self.const = const
        
        if (tip == 'str') and (const == 'v'):
            k, v, f, ap, D, z1, z2 = args
            self.t = k*np.pi*D*np.abs(z2-z1)/1e3/v/f
            self.v = v
            self.f = f
            self.ap = ap
            self.D = D
            self.L = np.abs(z2-z1)
            
        if (tip == 'str') and (const == 'n'):
            k, n, f, ap, D, z1, z2 = args
            self.t = k*np.abs(z2-z1)/n/f
            self.v = brzina_rezanja(D, n)
            self.f = f
            self.ap = ap
            self.D = D
            self.L = np.abs(z2-z1)
            
        if (tip == 'fce') and (const == 'v'):
            k, v, f, ap, x1, x2 = args
            self.t = k*np.pi*np.abs((x2**2.-x1**2.)/4.)/1e3/v/f
            self.v = v
            self.f = f
            self.ap = ap
            self.L = np.abs(x2-x1)/2.
            self.D = Dsr(x1, x2)
            
        if (tip == 'fce') and (const == 'n'):
            k, n, f, ap, x1, x2 = args
            self.t = k*np.abs(x2-x1)/2./n/f
            self.v = brzina_rezanja(Dsr(x1, x2), n)
            self.f = f
            self.ap = ap
            self.L = np.abs(x2-x1)/2.
            self.D = Dsr(x1, x2)
            
        if (tip == 'tpr') and (const == 'v'):
            k, v, f, ap, x1, x2, z1, z2 = args
            fi = np.arctan(np.abs(x2-x1)/2./np.abs(z2-z1))
            self.fi = fi
            self.t = k*np.pi*np.abs((x2**2.-x1**2.)/4./np.sin(fi))/1e3/v/f
            self.v = v
            self.f = f
            self.ap = ap
            self.L = np.abs(x2-x1)/2./np.sin(fi)
            self.D = Dsr(x1, x2)
            
        if (tip == 'tpr') and (const == 'n'):
            k, n, f, ap, x1, x2, z1, z2 = args
            fi = np.arctan(np.abs(x2-x1)/2./np.abs(z2-z1))
            self.fi = fi
            self.t = k*np.abs((x2-x1)/2./np.sin(fi))/n/f
            self.v = brzina_rezanja(Dsr(x1, x2), n)
            self.f = f
            self.ap = ap
            self.L = np.abs(x2-x1)/2./np.sin(fi)
            self.D = Dsr(x1, x2)
            
        if (tip == 'rad') and (const == 'v'):
            k, v, f, ap, x1, x2, z1, z2, ra, xc, zc = args
            if z1 == zc:
                fi1 = np.arctan(np.inf)
            else:
                fi1 = np.arctan(np.abs(x1-xc)/2./np.abs(z1-zc))
            if z2 == zc:
                fi2 = np.arctan(np.inf)
            else:
                fi2 = np.arctan(np.abs(x2-xc)/2./np.abs(z2-zc))
            self.fi1 = fi1
            self.fi2 = fi2
            self.t = (k*np.pi*ra*np.abs(xc*(fi2-fi1)-ra*(np.cos(fi2)-np.cos(fi1)))
                      /500./v/f)
            self.v = v
            self.f = f
            self.ap = ap
            self.D = Dsr(x1, x2)
            self.ra = ra
                      
        if (tip == 'rad') and (const == 'n'):
            k, n, f, ap, x1, x2, z1, z2, ra, xc, zc = args
            if z1 == zc:
                fi1 = np.arctan(np.inf)
            else:
                fi1 = np.arctan(np.abs(x1-xc)/2./np.abs(z1-zc))
            if z2 == zc:
                fi2 = np.arctan(np.inf)
            else:
                fi2 = np.arctan(np.abs(x2-xc)/2./np.abs(z2-zc))
            self.fi1 = fi1
            self.fi2 = fi2
            self.t = k*ra*np.abs(fi2-fi1)/n/f
            self.v = brzina_rezanja(Dsr(x1, x2), n)
            self.f = f
            self.ap = ap
            self.D = Dsr(x1, x2)
            self.ra = ra
            

k = 67.
T = {1: ('str', 'v', k, 100., 0.2, 0.815, 18., 0., -28.),
     2: ('str', 'v', k, 100., 0.2, 1., 16., 0., -28.),
     3: ('str', 'v', k, 100., 0.2, 1., 14., 0., -28.),
     4: ('str', 'v', k, 100., 0.2, 0.4, 13.2, -24.5, -28.),
     5: ('str', 'n', k, 2500., 0.1, 1., 12., 0., -11.),
     6: ('tpr', 'n', k, 2500., 0.1, 0.9175, 9.33, 11., 0., -11.),
     7: ('tpr', 'n', k, 2500., 0.1, 0.25, 8.33, 11., 0., -11.),
     8: ('str', 'v', k, 100., 0.1, 0.52, 12.96, -11.52, -24.5),
     9: ('str', 'n', k, 2500., 0.1, 0.5, 12.2, -24.5, -28.),
     10: ('fce', 'v', k, 100., 0.1, 0.3, 19.63, promjer_rez(100., 2500.)),
     11: ('fce', 'n', k, 2500., 0.1, 0.3, promjer_rez(100., 2500.), 0.),
     12: ('tpr', 'v', k, 100., 0.1, 0.25, 16.8, 19.63, -28., -28.52)}


for i in T.iterkeys():
    T[i] = Strojno_vrijeme(*T[i])

t, v, f, ap = np.array([]), np.array([]), np.array([]), np.array([])

for i in T.itervalues():
    t = np.append(t, i.t)
    v = np.append(v, i.v)
    f = np.append(f, i.f)
    ap = np.append(ap, i.ap)
    
def sum_kv(x, t, v, f, ap):
    sum_n = t*(v**x[1])*(f**x[2])*(ap**x[3])/x[0]
    sum_m = (np.sum(sum_n) - 1.)**2.
    return np.sum(sum_m)
    

cons = ({'type': 'ineq', 'fun': lambda x: x[1]-x[2],
         'jac': lambda x: np.array([0., 1., -1., 0.])},
        {'type': 'ineq', 'fun': lambda x: x[1]-x[3],
         'jac': lambda x: np.array([0., 1., 0., -1.])},
        {'type': 'ineq', 'fun': lambda x: x[2]-x[3],
         'jac': lambda x: np.array([0., 0., 1., -1.])},
        {'type': 'ineq', 'fun': lambda x: x[0]-x[1],
         'jac': lambda x: np.array([1., -1., 0., 0.])},
        {'type': 'ineq', 'fun': lambda x: x[0]**(1./x[1])-100.},
        {'type': 'ineq', 'fun': lambda x: 500.-x[0]**(1./x[1])}
         )

class MyBounds(object):
    def __init__(self, bnds = [(1e8, 5e16), (0.1, 10.), (0.1, 10.), (0.1, 10.)]):
        xmax, xmin = np.array([]), np.array([])
        for i in bnds:
            xmax = np.append(xmax, i[1])
            xmin = np.append(xmin, i[0])
        self.bnds = bnds
        self.xmax = xmax
        self.xmin = xmin
    def __call__(self, **kwargs):
        x = kwargs['x_new']
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return (tmax and tmin)

class MyStep(object):
    def __init__(self, stepsize=5e-1):
        self.stepsize = stepsize
    def __call__(self, x):
        s = self.stepsize
        x[0] += np.random.uniform(-1e11*s, 1e11*s)
        x[1:] += np.random.uniform(-s, s, x[1:].shape)
        return x
        

bnds = MyBounds().bnds
mybnds = MyBounds()
mystep = MyStep()

#res = minimize(sum_kv, (1e11, 1., 0.5, 1.), args=(t, v, f, ap), 
#               method='SLSQP', bounds=bnds, constraints=cons,
#               options={'disp': True, 'maxiter': 1000, 'ftol': 1e-15})
result = np.array([])

def callback_f(x, convergence):
    print x
    if (postojanost(x, v.min(), f.min(), ap.min()) <= 300.) and (
        postojanost(x, v.max(), f.max(), ap.max()) >= 1.):
            print x
        

#res = differential_evolution(sum_kv, bnds, args=(t, v, f, ap), strategy='best1bin',
#                             maxiter=1000, tol=1e-6, callback=callback_f,
#                             disp=True, polish=True)

res = basinhopping(sum_kv, (1e14, 0.5, 0.5, 0.5), niter=200, accept_test=mybnds,
                   take_step=mystep, 
                   minimizer_kwargs={'method': 'SLSQP', 'args': (t, v, f, ap),
                                     'constraints': cons,
                                     'options': {'maxiter': 1000.
                                                 }}, disp=True)

x = np.array([res.x[0]**(1./res.x[1]), 1./res.x[1], res.x[2]/res.x[1],
              res.x[3]/res.x[1]])

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
plt.plot(vplt[arr], postojanost(res.x, vplt[arr], fplt[arr], applt[arr]))
plt.plot(vplt[sr], postojanost(res.x, vplt[sr], fplt[sr], applt[sr]), 'k+')

plt.plot(np.ndarray.flatten(v), np.ndarray.flatten(postojanost(res.x, v, f, ap)),
         'bo')
         
plt.show()
