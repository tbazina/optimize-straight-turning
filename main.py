# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 12:04:13 2015

@author: tomislav
"""

from kivy.app import App
from kivy.uix.screenmanager import Screen, ScreenManager
from kivy.lang import Builder
from kivy.properties import StringProperty, AliasProperty, NumericProperty

import numpy as np
from scipy.optimize import basinhopping, minimize
import matplotlib as mpl
mpl.use(u'GTK3Agg')

import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from func import GraniceBasHop

Builder.load_file('optim.kv')


class ConsScreen(Screen):
    pass


class PathScreen(Screen):
    def StockCall(self):
        if (self.ids.d0.text != '') and (self.ids.l0.text != ''):
            Stck(self.ids.d0.text, self.ids.l0.text)

    def ResetContour(self):
        Cnt.cnt = {'Straight': [], 'Face': [], 'Taper': [], 'Circular': []}
        Stck(Stck.d0, Stck.l0)
        for i in ('x1', 'x2', 'z1', 'z2', 'xc', 'zc', 'ra'):
            self.ids[i].text = ''

    def AddContour(self):
        if self.ids.typ.text != 'Choose':
            coord = [self.ids.typ.text]
            for i in ('x1', 'x2', 'z1', 'z2', 'xc', 'zc', 'ra'):
                if self.ids[i].disabled is False:
                    if self.ids[i].text == '':
                        return
                    coord.append(self.ids[i].text)
            Cnt(coord)

    def NextCheck(self):
        for key, value in Cnt.cnt.iteritems():
            if key == 'Circular':
                for i in value:
                    if i[1] >= Stck.d0:
                        self.ids.nextpath.disabled = False
            elif key == 'Straight':
                for i in value:
                    if i[-1] >= Stck.l0:
                        self.ids.nextpath.disabled = False
            elif key == 'Taper':
                for i in value:
                    if i[1] >= Stck.d0:
                        self.ids.nextpath.disabled = False
            elif key == 'Face':
                for i in value:
                    if i[1] >= Stck.d0:
                        self.ids.nextpath.disabled = False

    def SpinnerText(self, typ):
        if typ == 'Straight':
            self.ids.x1.disabled = True
            self.ids.x2.disabled = False
            self.ids.z1.disabled = False
            self.ids.z2.disabled = False
            self.ids.xc.disabled = True
            self.ids.zc.disabled = True
            self.ids.ra.disabled = True
        if typ == 'Face':
            self.ids.x1.disabled = False
            self.ids.x2.disabled = False
            self.ids.z1.disabled = True
            self.ids.z2.disabled = False
            self.ids.xc.disabled = True
            self.ids.zc.disabled = True
            self.ids.ra.disabled = True
        if typ == 'Taper':
            self.ids.x1.disabled = False
            self.ids.x2.disabled = False
            self.ids.z1.disabled = False
            self.ids.z2.disabled = False
            self.ids.xc.disabled = True
            self.ids.zc.disabled = True
            self.ids.ra.disabled = True
        if typ == 'Circular':
            self.ids.x1.disabled = False
            self.ids.x2.disabled = False
            self.ids.z1.disabled = False
            self.ids.z2.disabled = False
            self.ids.xc.disabled = False
            self.ids.zc.disabled = False
            self.ids.ra.disabled = False


class ParamScreen(Screen):
    pass


class OptimScreen(Screen):
    def switching_function(self, *args):
        if self.ids.radiotime.active:
            sm.current = 'mptime_start'
        elif self.ids.radiocost.active:
            sm.current = 'mpcost_start'

    def on_radio_active(self):
        if self.ids.radiotime.active:
            self.ids.nextoptim.bind(on_press=self.switching_function)
            self.ids.nextoptim.disabled = False
        elif self.ids.radiocost.active:
            self.ids.nextoptim.bind(on_press=self.switching_function)
            self.ids.nextoptim.disabled = False
        else:
            self.ids.nextoptim.disabled = True


class MpTimeStart(Screen):
    pass


class MpCostStart(Screen):
    def set_cost(self):
        for key, val in self.ids.items():
            if (val.text != '') and (key != 'nextcost'):
                Cst(key, float(val.text))

        keys = ['ko', 'kt']
        if all(k in Cst.c for k in keys):
            self.ids.nextcost.disabled = False


class MpTimeInput(Screen):
    def switching_function(self, *args):
        if sm.get_screen('optim').ids.radiotime.active:
            sm.current = 'mptime_start'
        if sm.get_screen('optim').ids.radiocost.active:
            sm.current = 'mpcost_start'

    def next_switching_function(self, *args):
        if sm.get_screen('optim').ids.radiotime.active:
            sm.current = 'mptool_life'
        if sm.get_screen('optim').ids.radiocost.active:
            sm.current = 'mptool_life'

    def SetTime(self):
        for key, val in self.ids.items():
            if (val.text != '') and (key != 'nexttime'):
                Time(key, val.text)

    def NextCheck(self):
        if ('tl' in Time.t) and ('ts' in Time.t) and ('tc' in Time.t):
            self.ids.nexttime.disabled = False


class MpToolLife(Screen):
    def switching_function(self, *args):
        if sm.get_screen('optim').ids.radiotime.active:
            sm.current = 'mptime_input'
        if sm.get_screen('optim').ids.radiocost.active:
            sm.current = 'mptime_input'

    def next_switching_function(self, *args):
        if sm.get_screen('optim').ids.radiotime.active:
            sm.current = 'mpbounds_screen'
        if sm.get_screen('optim').ids.radiocost.active:
            sm.current = 'mpbounds_screen'

    def SetConstants(self):
        for key, val in self.ids.items():
            if (val.text != '') and (key != 'nexttool'):
                Tool(key, val.text)

    def NextCheck(self):
        if (('c' in Tool.cons) and ('kv' in Tool.cons) and
           ('kf' in Tool.cons) and ('ka' in Tool.cons)):
            self.ids.nexttool.disabled = False


class MpBoundsScreen(Screen):
    def switching_function(self, *args):
        if sm.get_screen('optim').ids.radiotime.active:
            sm.current = 'mptool_life'
        if sm.get_screen('optim').ids.radiocost.active:
            sm.current = 'mptool_life'

    def next_switching_function(self, *args):
        if sm.get_screen('optim').ids.radiotime.active:
            sm.current = 'mpsurfrel_screen'
        if sm.get_screen('optim').ids.radiocost.active:
            sm.current = 'mpsurfrel_screen'

    def SetBounds(self):
        keys = ['vr', 'fr', 'ar', 'epsr', 'm', 'vf', 'ff', 'af', 'epsf', 'tl']
        for i in keys:
            minvalue = self.ids[i+'min'].text
            maxvalue = self.ids[i+'max'].text
            if (minvalue != '') and (maxvalue != ''):
                Bounds(i, float(minvalue), float(maxvalue))

        if all(k in Bounds.bnds for k in keys):
            self.ids.nextbounds.disabled = False


class MpSurfRelScreen(Screen):
    def switching_function(self, *args):
        if sm.get_screen('optim').ids.radiotime.active:
            sm.current = 'mpbounds_screen'
        if sm.get_screen('optim').ids.radiocost.active:
            sm.current = 'mpbounds_screen'

    def next_switching_function(self, *args):
        if sm.get_screen('optim').ids.radiotime.active:
            sm.current = 'mpforce_screen'
        if sm.get_screen('optim').ids.radiocost.active:
            sm.current = 'mpforce_screen'

    def set_parameters(self, *args):
        for key, val in self.ids.items():
            if (val.text != '') and (key != 'nextsurf'):
                Param(key, float(val.text))

        keys = ['k1', 'k2', 'k3', 'reps', 'Ra', 'l1', 'l2']
        if all(k in Param.p for k in keys):
            self.ids.nextsurf.disabled = False


class MpForceScreen(Screen):
    def switching_function(self, *args):
        if sm.get_screen('optim').ids.radiotime.active:
            sm.current = 'mpsurfrel_screen'
        if sm.get_screen('optim').ids.radiocost.active:
            sm.current = 'mpsurfrel_screen'

    def next_switching_function(self, *args):
        if sm.get_screen('optim').ids.radiotime.active:
            sm.current = 'mpoptimization_screen'
        if sm.get_screen('optim').ids.radiocost.active:
            sm.current = 'mpoptimization_screen'

    def set_parameters(self, *args):
        for key, val in self.ids.items():
            if (val.text != '') and (key != 'nextforce'):
                Frc(key, float(val.text))

        keys = ['Kc', 'mc', 'kappa', 'Fmax', 'Ps', 'etas']
        if all(k in Frc.p for k in keys):
            self.ids.nextforce.disabled = False


class MpOptimizationScreen(Screen):
    def switching_function(self, *args):
        if sm.get_screen('optim').ids.radiotime.active:
            sm.current = 'mpforce_screen'
        if sm.get_screen('optim').ids.radiocost.active:
            sm.current = 'mpforce_screen'

    def next_switching_function(self, *args):
        if sm.get_screen('optim').ids.radiotime.active:
            sm.current = 'mprun_screen'
        if sm.get_screen('optim').ids.radiocost.active:
            sm.current = 'mprun_screen'

    def set_parameters(self, *args):
        for key, val in self.ids.items():
            if (val.text != '') and (key != 'nextoptimization'):
                if key in ['T', 'stepsize', 'ftol']:
                    Param(key, float(val.text))
                else:
                    Param(key, int(val.text))

        keys = ['niter', 'T', 'stepsize', 'interval', 'ftol', 'maxiter']
        if all(k in Param.p for k in keys):
            self.ids.nextoptimization.disabled = False


class MpRunScreen(Screen):
    output = StringProperty()
    result = StringProperty()
    button_text = StringProperty('Start Optimization')

    def switching_function(self, *args):
        if sm.get_screen('optim').ids.radiotime.active:
            sm.current = 'mpoptimization_screen'
        if sm.get_screen('optim').ids.radiocost.active:
            sm.current = 'mpoptimization_screen'

    def run_optimization(self, *args):
        RunOpt()

sm = ScreenManager()
sm.add_widget(PathScreen(name='path'))
sm.add_widget(OptimScreen(name='optim'))
sm.add_widget(MpTimeStart(name='mptime_start'))
sm.add_widget(MpCostStart(name='mpcost_start'))
sm.add_widget(MpTimeInput(name='mptime_input'))
sm.add_widget(MpToolLife(name='mptool_life'))
sm.add_widget(MpBoundsScreen(name='mpbounds_screen'))
sm.add_widget(MpSurfRelScreen(name='mpsurfrel_screen'))
sm.add_widget(MpForceScreen(name='mpforce_screen'))
sm.add_widget(MpOptimizationScreen(name='mpoptimization_screen'))
sm.add_widget(MpRunScreen(name='mprun_screen'))
sm.add_widget(ParamScreen(name='param'))
sm.add_widget(ConsScreen(name='cons'))


class Stock(object):
    def __init__(self):
        self.d0 = 1.
        self.l0 = 1.

    def __call__(self, d0, l0):
        self.d0 = float(d0)
        self.l0 = float(l0)
        fig = plt.figure(1, figsize=(6, 4))
        plt.clf()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Z')
        ax.set_ylabel('D')
        ax.grid(True, which='both')
        ax.plot([0, 0, self.l0], [0, self.d0, self.d0], c='g', lw=4,
                label='Stock')
        ax.set_xlim(self.l0, -5)
        ax.set_ylim(0, self.d0*1.2)
        ax.legend()
        self.ax = ax
        plt.savefig('contour.jpg', dpi=180, format='jpg')

Stck = Stock()


class Contour(object):
    def __init__(self):
        self.cnt = {'Straight': [], 'Face': [], 'Taper': [], 'Circular': []}

    def __call__(self, coord):
        x = np.array(coord[1:], dtype=np.float)
        self.cnt[coord[0]].append(x)
        if coord[0] == 'Straight':
            Stck.ax.plot(x[1:], np.array([x[0], x[0]]), c='b', lw=4)
        elif coord[0] == 'Face':
            Stck.ax.plot(np.array([x[-1], x[-1]]), x[:-1], c='b', lw=4)
        elif coord[0] == 'Taper':
            Stck.ax.plot(x[-2:], x[:2], c='b', lw=4)
        elif coord[0] == 'Circular':
            theta1 = np.degrees(np.arctan(np.divide((x[0]-x[-3])/2,
                                                    np.abs(x[2]-x[-2]))))
            theta2 = np.degrees(np.arctan(np.divide((x[1]-x[-3])/2,
                                                    np.abs(x[3]-x[-2]))))
            circ = Arc(xy=x[5:3:-1], width=x[-1]*2, height=x[-1]*4,
                       theta1=theta1+90, theta2=theta2+90, color='b', lw=4)
            Stck.ax.add_patch(circ)
        plt.savefig('contour.jpg', dpi=180, format='jpg')

Cnt = Contour()


class TimeInput(object):
    def __init__(self):
        self.t = {}

    def __call__(self, name, value, *args):
        self.t[name] = float(value)

    def tsf(self, x, *args):
        m = RunOpt.args[0]
        return Time.t['ts'] * (m + 1.)

    def trf(self, x, *args):
        m, L, d0 = RunOpt.args[:3]
        return np.pi*L*(m*d0-2.*x[2]*(m-1.))/1000./x[0]/x[1]

    def tff(self, x, *args):
        L, df = RunOpt.args[1:4:2]
        return np.pi*L*(df+2.*x[5])/1000./x[3]/x[4]

Time = TimeInput()


class CostInput(object):
    def __init__(self):
        self.c = {}

    def __call__(self, name, value, *args):
        self.c[name] = value

Cst = CostInput()


class ToolLife(object):
    def __init__(self):
        self.cons = {}

    def __call__(self, name, value, *args):
        self.cons[name] = float(value)

    def Tlf(self, x, *args):
        return ((Time.trf(x, args) + Time.tff(x, args)) / (Time.trf(x, args) /
                self.Trf(x, args) + Time.tff(x, args) / self.Tff(x, args)))

    def Trf(self, x, *args):
        C, kv, kf, ka = RunOpt.args[4:8]
        return C**kv * x[0]**(-kv) * x[1]**(-kf) * x[2]**(-ka)

    def Tff(self, x, *args):
        C, kv, kf, ka = RunOpt.args[4:8]
        return C**kv * x[3]**(-kv) * x[4]**(-kf) * x[5]**(-ka)

    def Tlminf(self, x, *args):
        Tlmin = RunOpt.args[8]
        return self.Tlf(x, args) - Tlmin

    def Tlmaxf(self, x, *args):
        Tlmax = RunOpt.args[9]
        return Tlmax - self.Tlf(x, args)

Tool = ToolLife()


class ParameterBounds(object):
    def __init__(self):
        self.bnds = {}

    def __call__(self, name, minvalue, maxvalue, *args):
        self.bnds[name] = (minvalue, maxvalue)

Bounds = ParameterBounds()


class Parameters(object):
    def __init__(self):
        self.p = {}

    def __call__(self, name, value, *args):
        self.p[name] = value

Param = Parameters()


class Force(object):
    def __init__(self):
        self.p = {}

    def __call__(self, name, value, *args):
        self.p[name] = float(value)

    def Frf(self, x, *args):
        kappa, Kc, mc = RunOpt.args[21:24]
        return x[1]**(1.-mc)*x[2]*Kc/np.sin(np.radians(kappa))**mc

    def Fff(self, x, *args):
        kappa, Kc, mc = RunOpt.args[21:24]
        return x[4]**(1.-mc)*x[5]*Kc/np.sin(np.radians(kappa))**mc

Frc = Force()


class WorkConditions(object):
    def epsrminf(self, x, *args):
        epsrmin = RunOpt.args[10]
        return x[2] / x[1] - epsrmin

    def epsrmaxf(self, x, *args):
        epsrmax = RunOpt.args[11]
        return epsrmax - x[2] / x[1]

    def epsfminf(self, x, *args):
        epsfmin = RunOpt.args[12]
        return x[5] / x[4] - epsfmin

    def epsfmaxf(self, x, *args):
        epsfmax = RunOpt.args[13]
        return epsfmax - x[5] / x[4]

    def vrelf(self, x, *args):
        k1 = RunOpt.args[14]
        return x[3] - k1 * x[0]

    def frelf(self, x, *args):
        k2 = RunOpt.args[15]
        return x[1] - k2 * x[4]

    def arelf(self, x, *args):
        k3 = RunOpt.args[16]
        return x[2] - k3 * x[5]

    def surfacef(self, x, *args):
        reps, Ra = RunOpt.args[17:19]
        return np.sqrt(32.*reps*Ra) - x[4]

    def geometryf(self, x, *args):
        m = RunOpt.args[0]
        d0, df = RunOpt.args[2:4]
        return df + 2. * m * x[2] + 2. * x[5] - d0

    def forcerf(self, x, *args):
        Fmax = RunOpt.args[24]
        return Fmax - Frc.Frf(x)

    def forceff(self, x, *args):
        Fmax = RunOpt.args[24]
        return Fmax - Frc.Fff(x)

    def powerrf(self, x, *args):
        Ps, etas = RunOpt.args[19:21]
        return Ps*etas - Frc.Frf(x)*x[0]/60e3

    def powerff(self, x, *args):
        Ps, etas = RunOpt.args[19:21]
        return Ps*etas - Frc.Fff(x)*x[3]/60e3

WrkCond = WorkConditions()


class RunOptimization(object):
    def __init__(self):
        self.res = {}
        self.mystep = self.MyStep()
        self.x = self.dummyclass()

    class MyStep(object):
        def __init__(self, stepsize=1e-1):
            self.stepsize = stepsize

        def __call__(self, x):
            s = self.stepsize
            x[::3] += np.random.uniform(-500.*s, 500.*s, x[::3].shape)
            x[1::3] += np.random.uniform(-s, s, x[1::3].shape)
            x[2::3] += np.random.uniform(-15.*s, 15.*s, x[2::3].shape)
            return x

    class CallBackFunction(object):
        def __init__(self):
            self.step = 1
            self.text = ''

        def __call__(self, x, f, accept):
            txt = '{}: x={}, f={}, accept={}\n'.format(
                self.step, x, f, accept)
            self.text += txt
            self.step += 1

    class dummyclass(object):
        def __init__(self):
            self.fun = np.float('inf')
            self.success = False

    def __call__(self):
        self.x = self.dummyclass()
        if sm.get_screen('optim').ids.radiotime.active:
            Time.t['tc*'] = Time.t['tc']
        if sm.get_screen('optim').ids.radiocost.active:
            Time.t['tc*'] = Time.t['tc'] + Cst.c['kt'] / Cst.c['ko']
        self.niter = Param.p['niter']
        self.stepsize = Param.p['stepsize']
        self.T = Param.p['T']
        self.interval = Param.p['interval']
        self.ftol = Param.p['ftol']
        self.maxiter = Param.p['maxiter']
        self.bnds = [Bounds.bnds['vr'], Bounds.bnds['fr'], Bounds.bnds['ar'],
                     Bounds.bnds['vf'], Bounds.bnds['ff'], Bounds.bnds['af']]
        self.init = np.mean(self.bnds, axis=1)
        self.mybnds = GraniceBasHop(self.bnds)
        d0 = Stck.d0
        df = Cnt.cnt['Straight'][0][0]
        ar = Bounds.bnds['ar']
        af = Bounds.bnds['af']
        afmax = np.min((af[1], ar[1]/Param.p['k3']))
        self.mlim = np.arange(np.max((np.int((d0-df-2.*afmax)/2./ar[1]+0.9999),
                                      np.int(Bounds.bnds['m'][0]))),
                              np.min((np.int((d0-df-2.*af[0])/2./ar[0]+0.9999),
                                      np.int(Bounds.bnds['m'][1]+1))))
        self.callback = self.CallBackFunction()
        j = 0
        for i in self.mlim:
            self.args = np.array(
                (i, np.diff(Cnt.cnt['Straight'][0][1:])+Param.p['l1'] +
                 Param.p['l2'], Stck.d0, Cnt.cnt['Straight'][0][0],
                 Tool.cons['c'], Tool.cons['kv'], Tool.cons['kf'],
                 Tool.cons['ka'], Bounds.bnds['tl'][0], Bounds.bnds['tl'][1],
                 Bounds.bnds['epsr'][0], Bounds.bnds['epsr'][1],
                 Bounds.bnds['epsf'][0], Bounds.bnds['epsf'][1],
                 Param.p['k1'], Param.p['k2'], Param.p['k3'], Param.p['reps'],
                 Param.p['Ra'], Frc.p['Ps'], Frc.p['etas'], Frc.p['kappa'],
                 Frc.p['Kc'], Frc.p['mc'], Frc.p['Fmax']
                 ), dtype=np.float64)

            self.cons = ({'type': 'ineq', 'fun': Tool.Tlminf},
                         {'type': 'ineq', 'fun': Tool.Tlmaxf},
                         {'type': 'ineq', 'fun': WrkCond.epsrminf},
                         {'type': 'ineq', 'fun': WrkCond.epsrmaxf},
                         {'type': 'ineq', 'fun': WrkCond.epsfminf},
                         {'type': 'ineq', 'fun': WrkCond.epsfmaxf},
                         {'type': 'ineq', 'fun': WrkCond.vrelf},
                         {'type': 'ineq', 'fun': WrkCond.frelf},
                         {'type': 'ineq', 'fun': WrkCond.arelf},
                         {'type': 'ineq', 'fun': WrkCond.surfacef},
                         {'type': 'eq', 'fun': WrkCond.geometryf},
                         {'type': 'ineq', 'fun': WrkCond.forcerf},
                         {'type': 'ineq', 'fun': WrkCond.forceff},
                         {'type': 'ineq', 'fun': WrkCond.powerrf},
                         {'type': 'ineq', 'fun': WrkCond.powerff}
                         )

            self.mystep = self.MyStep(self.stepsize)
            self.mybnds = GraniceBasHop(self.bnds)

            self.res[j] = basinhopping(
                            ObjFun.f, self.init, niter=self.niter, T=self.T,
                            take_step=self.mystep, disp=True,
                            accept_test=self.mybnds, interval=self.interval,
                            callback=self.callback,
                            minimizer_kwargs={'method': 'SLSQP',
                                              'bounds': self.bnds,
                                              'constraints': self.cons,
                                              'options': {'ftol': self.ftol,
                                                          'maxiter':
                                                              self.maxiter}
                                              })
            if (self.res[j].lowest_optimization_result.success and
               (self.res[j].fun < self.x.fun)):
                self.x = self.res[j]
                self.x.x = np.append(self.x.x, i)
            j += 1
        sqp = self.x.lowest_optimization_result
        resprint = ('{}\n{}\nFunction value: {}\nFunction evaluations: {}\n'
                    'Jacobian evaluations: {}\nIterations: {}\n\nRoughing:\nm '
                    '= {} passes\nvr = {} m/min\nfr = {} mm/rev\nar = {} mm\n'
                    '\nFinishing:\nvf = {} m/min\nff = {} mm/rev\naf = {} mm'
                    ).format(
                    sqp.message, self.x.message[0], sqp.fun, self.x.nfev,
                    self.x.njev, self.x.nit, self.x.x[6], sqp.x[0], sqp.x[1],
                    sqp.x[2], sqp.x[3], sqp.x[4], sqp.x[5])
        sm.get_screen('mprun_screen').output = self.callback.text
        sm.get_screen('mprun_screen').result = resprint
        sm.get_screen('mprun_screen').button_text = ('Optimization finished '
                                                     '(Rerun)')

RunOpt = RunOptimization()


class ObjectiveFunction(object):
    def f(self, x, *args):
        return (Time.t['tl'] + Time.tsf(x, args) + Time.trf(x, args) +
                Time.tff(x, args) + Time.t['tc*'] * (Time.trf(x, args) +
                Time.tff(x, args)) / Tool.Tlf(x, args))

ObjFun = ObjectiveFunction()


class OptimizeApp(App):
    def build(self):
        return sm


if __name__ == '__main__':
    OptimizeApp().run()
