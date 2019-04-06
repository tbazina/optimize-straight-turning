# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 20:42:28 2015

@author: tomislav
"""

import numpy as np
from scipy.optimize import basinhopping, differential_evolution
from func import *

k = 67.
HiCell = ([('str', 'v', k, 100., 0.2, 0.815, 18., 0., -28.),
           ('str', 'v', k, 100., 0.2, 1., 16., 0., -28.),
           ('str', 'v', k, 100., 0.2, 1., 14., 0., -28.),
           ('str', 'v', k, 100., 0.2, 0.4, 13.2, -24.5, -28.),
           ('str', 'n', k, 2500., 0.1, 1., 12., 0., -11.),
           ('tpr', 'n', k, 2500., 0.1, 0.9175, 9.33, 11., 0., -11.),
           ('tpr', 'n', k, 2500., 0.1, 0.25, 8.33, 11., 0., -11.),
           ('str', 'v', k, 100., 0.1, 0.52, 12.96, -11.52, -24.5),
           ('str', 'n', k, 2500., 0.1, 0.5, 12.2, -24.5, -28.),
           ('fce', 'v', k, 100., 0.1, 0.3, 19.63, promjer_rez(100., 2500.)),
           ('fce', 'n', k, 2500., 0.1, 0.3, promjer_rez(100., 2500.), 0.),
           ('tpr', 'v', k, 100., 0.1, 0.25, 16.8, 19.63, -28., -28.52)],)

t = StrojnoVrijeme(HiCell)
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

bnds = [(180., 230.), (1.6, 3.4), (0.2, 1.5), (0.2, 1.5)]
# bnds = [(168., 175.), (3.9, 4.1), (0.46, 0.58), (0.35, 0.6)]
# plotf = PlotSumKv(sum_kv, bnds)
# init = (plotf.Cm, plotf.cv, plotf.cf, plotf.ca)

"""
class MyStep(object):
    def __init__(self, stepsize=0.6e-1):
        self.stpsize = stepsize

    def __call__(self, x):
        s = self.stpsize
        x[0] += np.random.uniform(-10.*s, 10.*s)
        x[1] += np.random.uniform(-s, s)
        x[2:] += np.random.uniform(-s, s, x[2:].shape)
        return x

mybnds = GraniceBasHop(bnds)
mystep = MyStep()

rezultati = []

for i in xrange(1):

    res = basinhopping(sum_kv, init, niter=30, accept_test=mybnds,
                       take_step=mystep, T=1e-18, interval=50,
                       minimizer_kwargs={'method': 'Newton-CG',
                                         'jac': sum_kv_jac,
                                         'options': {'maxiter': 1000,
                                                     'xtol': 1e-8
                                                     }}, disp=True)
    mybnds = GraniceBasHop(bnds)
    mystep = MyStep()
    rezultati.append(res)

for i in rezultati:
    print i.x, i.fun
"""

res = differential_evolution(sum_kv, bnds, strategy='best1bin', disp=True,
                             popsize=40, maxiter=30000, tol=1e-6,
                             mutation=(0.5, 1.), recombination=0.6,
                             polish=True)

plot = PlotT(res, t.vc, t.f, t.ap)
plt.show()
