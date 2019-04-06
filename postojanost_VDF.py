# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:28:47 2015

@author: tomislav
"""

import numpy as np
from scipy.optimize import basinhopping, differential_evolution
from func import *

k = np.array([138., 110.])
VDF = ([('fce', 'v', k[0]*3, 160., 0.2, 0.4, 34.4, promjer_rez(160., 2400.)),
        ('fce', 'n', k[0]*3, 2400., 0.2, 0.4, promjer_rez(160., 2400.), 0.),
        ('str', 'v', k[0], 160., 0.25, 1.3, 31.8, 0., -15.95),
        ('str', 'v', k[0], 160., 0.25, 1.6, 28.6, 0., -5.),
        ('str', 'v', k[0], 160., 0.25, 0.15, 28.3, 0., -5.),
        ('tpr', 'v', k[0], 160., 0.25, 0.325, 28.3, 30.8, -5., -5.4),
        ('str', 'v', k[0], 160., 0.25, 0.5, 30.8, -5.4, -15.95),
        ('tpr', 'v', k[0], 160., 0.2, 0.313, 30.8, 29.55, -10., -10.2),
        ('str', 'v', k[0], 160., 0.2, 0.625, 29.55, -10.2, -15.95),
        ('tpr', 'v', k[0], 160., 0.2, 0.313, 29.55, 28.3, -10.2, -10.4),
        ('str', 'v', k[0], 160., 0.2, 0.625, 28.3, -10.4, -15.95),
        ('str', 'v', k[0], 160., 0.2, 0.4, 27.5, 0., -5.),
        ('tpr', 'v', k[0], 160., 0.2, 0.4, 27.5, 30., -5., -5.4),
        ('str', 'v', k[0], 160., 0.2, 0.4, 30., -5.4, -10.),
        ('tpr', 'v', k[0], 160., 0.2, 0.4, 30., 27.5, -10., -10.4),
        ('str', 'v', k[0], 160., 0.2, 0.4, 27.5, -10.4, -16.),
        ('fce', 'v', k[0], 160., 0.2, 0.05, 27.5, 34.4),
        ('fce', 'v', k[0]*3, 150., 0.15, 0.3, 34.4, promjer_rez(150., 2400.)),
        ('fce', 'n', k[0]*3, 2400., 0.15, 0.3, promjer_rez(150., 2400.), 0.)],
       [('fce', 'v', k[1]*3, 160., 0.2, 0.4, 34.4, promjer_rez(160., 2400.)),
        ('fce', 'n', k[1]*3, 2400., 0.2, 0.4, promjer_rez(160., 2400.), 0.),
        ('str', 'v', k[1], 160., 0.27, 1.3, 31.8, 0., -15.95),
        ('str', 'v', k[1], 160., 0.27, 1.6, 28.6, 0., -5.),
        ('str', 'v', k[1], 160., 0.27, 0.15, 28.3, 0., -5.),
        ('tpr', 'v', k[1], 160., 0.27, 0.325, 28.3, 30.8, -5., -5.4),
        ('str', 'v', k[1], 160., 0.27, 0.5, 30.8, -5.4, -15.95),
        ('tpr', 'v', k[1], 160., 0.22, 0.313, 30.8, 29.55, -10., -10.2),
        ('str', 'v', k[1], 160., 0.22, 0.625, 29.55, -10.2, -15.95),
        ('tpr', 'v', k[1], 160., 0.22, 0.313, 29.55, 28.3, -10.2, -10.4),
        ('str', 'v', k[1], 160., 0.22, 0.625, 28.3, -10.4, -15.95),
        ('str', 'v', k[1], 160., 0.22, 0.4, 27.5, 0., -5.),
        ('tpr', 'v', k[1], 160., 0.22, 0.4, 27.5, 30., -5., -5.4),
        ('str', 'v', k[1], 160., 0.22, 0.4, 30., -5.4, -10.),
        ('tpr', 'v', k[1], 160., 0.22, 0.4, 30., 27.5, -10., -10.4),
        ('str', 'v', k[1], 160., 0.22, 0.4, 27.5, -10.4, -16.),
        ('fce', 'v', k[1], 160., 0.22, 0.05, 27.5, 34.4),
        ('fce', 'v', k[1]*3, 150., 0.22, 0.3, 34.4, promjer_rez(150., 2400.)),
        ('fce', 'n', k[1]*3, 2400., 0.22, 0.3, promjer_rez(150., 2400.), 0.)])


t = StrojnoVrijeme(VDF)
T = PostojanostAlata()


def sum_kv(x=np.array([], dtype=np.longfloat)):
    return (((t.t/T(x, t.vc, t.f, t.ap)).sum(1)-1.)**2.).sum(0)


def sum_kv_jac(x):
    dfdx0 = (2.*((t.t/T(x, t.vc, t.f, t.ap)).sum(1)-1.) *
             (t.t/T(x, t.vc, t.f, t.ap)*(-x[1]/x[0])).sum(1)).sum(0)
    dfdx1 = (2.*((t.t/T(x, t.vc, t.f, t.ap)).sum(1)-1.) *
             (t.t/T(x, t.vc, t.f, t.ap)*np.log(t.vc/x[0])).sum(1)).sum(0)
    dfdx2 = (2.*((t.t/T(x, t.vc, t.f, t.ap)).sum(1)-1.) *
             (t.t/T(x, t.vc, t.f, t.ap)*np.log(t.f)).sum(1)).sum(0)
    dfdx3 = (2.*((t.t/T(x, t.vc, t.f, t.ap)).sum(1)-1.) *
             (t.t/T(x, t.vc, t.f, t.ap)*np.log(t.ap)).sum(1)).sum(0)
    return np.array([dfdx0, dfdx1, dfdx2, dfdx3])


bnds = [(1., 600.), (1., 3.6), (0.2, 3.), (0.2, 3.)]
# plotf = PlotSumKv(sum_kv, bnds)
# init = (plotf.Cm, plotf.cv, plotf.cf, plotf.ca)


call = CallBackDE(10, 4)
"""
res = basinhopping(sum_kv, init, niter=10000, accept_test=mybnds,
                   take_step=mystep, T=0.01, interval=50, callback=minima,
                   minimizer_kwargs={'method': 'SLSQP',
                                     'bounds': bnds,
                                     'jac': sum_kv_jac,
                                     'options': {'maxiter': 1000,
                                                 # 'disp': True
                                                 }}, disp=True)
"""

res = differential_evolution(sum_kv, bnds, strategy='best1bin', disp=True,
                             popsize=40, maxiter=30000, tol=1e-6,
                             mutation=(0.5, 1.), recombination=0.6,
                             polish=True)

# SaveToFile(res, 'VDF.pkl')
plot = PlotT(res, t.vc, t.f, t.ap)

plt.show()
