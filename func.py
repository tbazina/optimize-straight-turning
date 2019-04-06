# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 20:04:03 2015

@author: tomislav
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import cPickle as pickle

pltfnt = {'family': 'sans-serif',
          'weight': 'normal',
          'size': 6}
fntdict = {'fontsize': 15,
           'fontweight': 'bold'}
plt.rc('font', **pltfnt)


def brzina_rezanja(promjer, brzina_vrtnje):
    """
    brzina_rezanja [m/min]
    \tpromjer [mm]\n
    \tbrzina_vrtnje [okr/min]
    """
    return promjer*np.pi*brzina_vrtnje/1e3


def brzina_vrtnje(promjer, brzina_rezanja):
    """
    brzina_vrtnje [okr/min]
    \tpromjer [mm]\n
    \tbrzina_rezanja [m/min]
    """
    return brzina_rezanja*1e3/promjer/np.pi


def volumen_strugotine(presjek_strugotine, brzina_rezanja):
    """
    volumen_strugotine [mm3/min]
    \tpresjek_strugotine [mm2]\n
    \tbrzina_rezanja [m/min]
    """
    return 1e3*presjek_strugotine*brzina_rezanja


def masa_strugotine(volumen_strugotine, gustoca_materijala):
    """
    masa_strugotine [kg/min]
    \tvolumen_strugotine [mm3/min]\n
    \tgustoca_materijala [kg/m3]
    """
    return 1e-9*volumen_strugotine*gustoca_materijala


def glavna_sila_rezanja(fs, A, Kgamma, Kv, Ka, Ki):
    """
    glavna_sila_rezanja [N]
    \tfs - specificna sila rezanja [N/mm2]\n
    \tA - presjek strugotine [mm2]
    """
    return fs*A*Kgamma*Kv*Ka*Ki


def presjek_strugotine(dubina_rezanja, posmak):
    """
    presjek_strugotine [mm2]
    \tdubina_rezanja [mm]\n
    \tposmak [mm/okr]
    """
    return dubina_rezanja*posmak


def debljina_strugotine(posmak, prisloni_kut):
    """
    debljina_strugotine [mm]
    \tposmak [mm/okr]\n
    \tprisloni_kut [°]
    """
    return posmak*np.sin(np.radians(prisloni_kut))


def spec_sila_rezanja(fs1x1, debljina_strugotine, eksponent_Kienzlea):
    """
    spec_sila_rezanja [N/mm2]
    \tfs1x1 [N/mm2]\n
    \tdebljina_strugotine [mm]\n
    \teksponent_Kienzlea [/]
    """
    return fs1x1/(debljina_strugotine)**eksponent_Kienzlea


def snaga_tokarenje(glavna_sila, brzina_rezanja):
    """
    snaga_tokarenje [kW]
    \tglavna_sila [N]\n
    \tbrzina_rezanja [m/min]
    """
    return glavna_sila*brzina_rezanja/60e3


def snaga_stroja(snaga_rezanja, stupanj_iskoristivosti):
    """
    snaga_stroja [kW]
    \tsnaga_rezanja [kW]\n
    \tstupanj_iskoristivosti [/]
    """
    return snaga_rezanja/stupanj_iskoristivosti


def Ra(posmak, radijus_alata):
    """
    teorijska_hrapavost [mm]
    \tposmak [mm/okr]\n
    \tradijus_alata [mm]
    """
    return posmak**2./32./radijus_alata


def promjer_rez(brzina_rezanja, brzina_vrtnje):
    """
    promjer_rez [mm]
    \tbrzina_rezanja [m/min]\n
    \tbrzina_vrtnje [okr/min]
    """
    return brzina_rezanja*1e3/np.pi/brzina_vrtnje


def Dsr(Dv, Du):
    return np.sqrt(Dv**2/2.+Du**2/2.)


class PostojanostAlata(object):
    """
    Postojanost [min]
    potrebno prilagoditi za optimiranje
    """
    def __call__(self, *args):
        if args.__len__() == 4:
            Cm, cv, cf, ca = args[0]
            vc, f, ap = args[1:]
            Ct = Cm**cv
            kv = -cv
            kf = -cf
            ka = -ca
            # self.Cm = Cm
            # self.cv = cv
            # self.cf = cf
            # self.ca = ca
        elif args.__len__() == 7:
            vc, f, ap, Ct, kv, kf, ka = args
        else:
            print 'ERROR'
        # self.vc = vc
        # self.f = f
        # self.ap = ap
        # self.Ct = Ct
        # self.kv = kv
        # self.kf = kf
        # self.ka = ka
        self.T = Ct*(vc**kv)*(f**kf)*(ap**ka)
        return self.T


class PostojanostAlataGrad(object):
    """
    Postojanost gradijent
    """
    def __call__(self, *args):
        if args.__len__() == 4:
            Cm, cv, cf, ca = args[0]
            vc, f, ap = args[1:]
            Ct = Cm**cv
            kv = -cv
            kf = -cf
            ka = -ca
            # self.Cm = Cm
            # self.cv = cv
            # self.cf = cf
            # self.ca = ca
        elif args.__len__() == 7:
            vc, f, ap, Ct, kv, kf, ka = args
        else:
            print 'ERROR'
        # self.vc = vc
        # self.f = f
        # self.ap = ap
        # self.Ct = Ct
        # self.kv = kv
        # self.kf = kf
        # self.ka = ka
        self.Tv = Ct*kv*(vc**(kv-1.))*(f**kf)*(ap**ka)
        self.Tf = Ct*kf*(vc**kv)*(f**(kf-1.))*(ap**ka)
        self.Tap = Ct*ka*(vc**kv)*(f**kf)*(ap**(ka-1.))
        return self.Tv, self.Tf, self.Tap


class StrojnoVrijeme(object):
    """
    StrojnoVrijeme [min]
    """
    def __init__(self, thist):
        alati = thist.__len__()
        rezimi = np.array([], dtype=np.longfloat)
        for i in thist.__iter__():
            rezimi = np.append(rezimi, i.__len__())

        tuk = np.zeros([alati, rezimi.max()], dtype=np.longfloat)
        vuk = np.zeros([alati, rezimi.max()], dtype=np.longfloat)
        fuk = np.zeros([alati, rezimi.max()], dtype=np.longfloat)
        apuk = np.zeros([alati, rezimi.max()], dtype=np.longfloat)

        a = 0
        for i in thist.__iter__():
            b = 0
            for j in i.__iter__():
                tip, const = j[0:2]

                if (tip == 'str') and (const == 'v'):
                    k, v, f, ap, D, z1, z2 = j[2:]
                    tuk[a][b] = k*np.pi*D*np.abs(z2-z1)/1e3/v/f
                    vuk[a][b] = v

                if (tip == 'str') and (const == 'n'):
                    k, n, f, ap, D, z1, z2 = j[2:]
                    tuk[a][b] = k*np.abs(z2-z1)/n/f
                    v = brzina_rezanja(D, n)
                    vuk[a][b] = v

                if (tip == 'fce') and (const == 'v'):
                    k, v, f, ap, x1, x2 = j[2:]
                    tuk[a][b] = k*np.pi*np.abs((x2**2.-x1**2.)/4.)/1e3/v/f
                    vuk[a][b] = v

                if (tip == 'fce') and (const == 'n'):
                    k, n, f, ap, x1, x2 = j[2:]
                    tuk[a][b] = k*np.abs(x2-x1)/2./n/f
                    v = brzina_rezanja(Dsr(x1, x2), n)
                    vuk[a][b] = v

                if (tip == 'tpr') and (const == 'v'):
                    k, v, f, ap, x1, x2, z1, z2 = j[2:]
                    fi = np.arctan(np.abs(x2-x1)/2./np.abs(z2-z1))
                    tuk[a][b] = k*np.pi*np.abs((x2**2.-x1**2.)/4./np.sin(fi))/1e3/v/f
                    vuk[a][b] = v

                if (tip == 'tpr') and (const == 'n'):
                    k, n, f, ap, x1, x2, z1, z2 = j[2:]
                    fi = np.arctan(np.abs(x2-x1)/2./np.abs(z2-z1))
                    tuk[a][b] = k*np.abs((x2-x1)/2./np.sin(fi))/n/f
                    v = brzina_rezanja(Dsr(x1, x2), n)
                    vuk[a][b] = v

                if (tip == 'rad') and (const == 'v'):
                    k, v, f, ap, x1, x2, z1, z2, ra, xc, zc = j[2:]
                    if z1 == zc:
                        fi1 = np.arctan(np.inf)
                    else:
                        fi1 = np.arctan(np.abs(x1-xc)/2./np.abs(z1-zc))
                    if z2 == zc:
                        fi2 = np.arctan(np.inf)
                    else:
                        fi2 = np.arctan(np.abs(x2-xc)/2./np.abs(z2-zc))
                    tuk[a][b] = (k*np.pi*ra*np.abs(xc*(fi2-fi1) -
                                 ra*(np.cos(fi2)-np.cos(fi1)))/500./v/f)
                    vuk[a][b] = v

                if (tip == 'rad') and (const == 'n'):
                    k, n, f, ap, x1, x2, z1, z2, ra, xc, zc = j[2:]
                    if z1 == zc:
                        fi1 = np.arctan(np.inf)
                    else:
                        fi1 = np.arctan(np.abs(x1-xc)/2./np.abs(z1-zc))
                    if z2 == zc:
                        fi2 = np.arctan(np.inf)
                    else:
                        fi2 = np.arctan(np.abs(x2-xc)/2./np.abs(z2-zc))
                    tuk[a][b] = k*ra*np.abs(fi2-fi1)/n/f
                    v = brzina_rezanja(Dsr(x1, x2), n)
                    vuk[a][b] = v

                fuk[a][b] = f
                apuk[a][b] = ap
                b += 1

            a += 1

        self.alat = alati
        self.rezim = rezimi
        self.t = np.ma.masked_equal(tuk, 0.)
        self.vc = np.ma.masked_equal(vuk, 0.)
        self.f = np.ma.masked_equal(fuk, 0.)
        self.ap = np.ma.masked_equal(apuk, 0.)

    def __call__(self, tip, const, *args):
        if (tip == 'str') and (const == 'v'):
            v, f, D, z1, z2 = args
            return np.pi*D*np.abs(z2-z1)/1e3/v/f

        if (tip == 'str') and (const == 'n'):
            n, f, D, z1, z2 = args
            return np.abs(z2-z1)/n/f

        if (tip == 'fce') and (const == 'v'):
            v, f, x1, x2 = args
            return np.pi*np.abs((x2**2.-x1**2.)/4.)/1e3/v/f

        if (tip == 'fce') and (const == 'n'):
            n, f, x1, x2 = args
            return np.abs(x2-x1)/2./n/f

        if (tip == 'tpr') and (const == 'v'):
            v, f, x1, x2, z1, z2 = args
            fi = np.arctan(np.abs(x2-x1)/2./np.abs(z2-z1))
            return np.pi*np.abs((x2**2.-x1**2.)/4./np.sin(fi))/1e3/v/f

        if (tip == 'tpr') and (const == 'n'):
            n, f, x1, x2, z1, z2 = args
            fi = np.arctan(np.abs(x2-x1)/2./np.abs(z2-z1))
            return np.abs((x2-x1)/2./np.sin(fi))/n/f

        if (tip == 'rad') and (const == 'v'):
            v, f, x1, x2, z1, z2, ra, xc, zc = args
            if z1 == zc:
                fi1 = np.arctan(np.inf)
            else:
                fi1 = np.arctan(np.abs(x1-xc)/2./np.abs(z1-zc))
            if z2 == zc:
                fi2 = np.arctan(np.inf)
            else:
                fi2 = np.arctan(np.abs(x2-xc)/2./np.abs(z2-zc))
            return (np.pi*ra*np.abs(xc*(fi2-fi1)-ra*(np.cos(fi2)-np.cos(fi1)))
                    / 500./v/f)

        if (tip == 'rad') and (const == 'n'):
            n, f, x1, x2, z1, z2, ra, xc, zc = args
            if z1 == zc:
                fi1 = np.arctan(np.inf)
            else:
                fi1 = np.arctan(np.abs(x1-xc)/2./np.abs(z1-zc))
            if z2 == zc:
                fi2 = np.arctan(np.inf)
            else:
                fi2 = np.arctan(np.abs(x2-xc)/2./np.abs(z2-zc))
            return ra*np.abs(fi2-fi1)/n/f


class GraniceBasHop(object):
    """
    Granice pri Basinhopping metodi trazenja minimuma
    """
    def __init__(self, bnds):
        xmax, xmin = np.array([]), np.array([])
        for i in bnds:
            xmax = np.append(xmax, i[1])
            xmin = np.append(xmin, i[0])
        self.xmax = xmax
        self.xmin = xmin

    def __call__(self, **kwargs):
        x = kwargs['x_new']
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return (tmax and tmin)


class CalBaBasHop(object):
    """
    Spemanje pronadenih minimuma
    """
    def __init__(self, n, x_size):
        mi = np.empty((n, x_size+1))
        mi[:, -1].fill(np.inf)
        self.minima = mi

    def __call__(self, x, f, accept):
        if (bool(np.any(f < self.minima[:, -1])) and bool(accept)):
            self.minima[self.minima[:, -1].argmax()] = np.append(x, f)


class CallBackDE(object):
    """
    Spemanje pronadenih minimuma
    """
    def __init__(self, n, x_size):
        mi = np.empty((n, x_size+1))
        mi[:, -1].fill(0.)
        self.minima = mi

    def __call__(self, xk, convergence):
        x, val = xk, convergence
        if bool(np.any(val > self.minima[:, -1])):
            self.minima[self.minima[:, -1].argmin()] = np.append(x, val)


class PlotT(object):
    """
    Plot postojanosti alata
    """
    def __init__(self, res, vc, f, ap):
        mi = 'minimum'
        ma = 'maksimum'
        sr = 'srednje'
        arr = 'array'

        nplt = 500

        vplt = {mi: vc.min(), ma: vc.max(), sr: vc.mean()}
        vplt[arr] = np.linspace(vplt[mi], vplt[ma], nplt)

        fplt = {mi: f.min(), ma: f.max(), sr: f.mean()}
        fplt[arr] = np.linspace(fplt[mi], fplt[ma], nplt)

        applt = {mi: ap.min(), ma: ap.max(), sr: ap.mean()}
        applt[arr] = np.linspace(applt[mi], applt[ma], nplt)

        t = PostojanostAlata()

        self.d3 = plt.figure(figsize=(22, 11), dpi=300)
        ax = self.d3.add_subplot(111, projection='3d')
        ax.set_title('Postojanost alata pri razlicitim dubinama rezanja',
                     fntdict)
        ax.set_xlabel(r'Brzina rezanja $[\frac{m}{min}]$', fntdict)
        ax.set_ylabel(r'Posmak $[\frac{mm}{okr}]$', fntdict)
        ax.set_zlabel(r'Postojanost alata $[min]$', fntdict)
        V, F = np.meshgrid(vplt[arr], fplt[arr])
        Ap = np.empty(V.shape)
        Postmax = np.array([])
        for i in applt[mi], applt[sr], applt[ma]:
            Ap[:] = i
            Post = t(res.x, V, F, Ap)
            Postmax = np.append(Postmax, Post.max())
            N = Post/Postmax.max()
            ax.plot_surface(V, F, Post, linewidth=0, facecolors=cm.jet_r(N),
                            antialiased=False, shade=False)


class PlotSumKv(object):
    """
    3D plot funkcije
    """
    def __init__(self, sum_kv, bnds):

        Cm, cv, cf, ca = {}, {}, {}, {}
        mi = 'minimum'
        ma = 'maksimum'
        sr = 'srednje'
        arr = 'array'

        (Cm[mi], Cm[ma]), (cv[mi], cv[ma]), (cf[mi], cf[ma]), (ca[mi], ca[ma]) = bnds
        nplt = 50

        for i in Cm, cv, cf, ca:
            i[arr] = np.linspace(i[mi], i[ma], nplt)
            i[sr] = np.mean(np.array([i[mi], i[ma]]))

        self.Cm = Cm[sr]
        self.cv = cv[sr]
        self.cf = cf[sr]
        self.ca = ca[sr]

        fntdict = {'fontsize': 11}

        fig = plt.figure(1, figsize=(12, 8), dpi=99)
        sp = fig.add_subplot(111)
        sp.spines['top'].set_color('none')
        sp.spines['bottom'].set_color('none')
        sp.spines['left'].set_color('none')
        sp.spines['right'].set_color('none')
        sp.tick_params(labelcolor='w', top='off', bottom='off', left='off',
                       right='off')
        sp.set_title(r'Ovisnost sume kvadrata o $C_m$, $c_v$, $c_f$ i $c_a$',
                     fntdict)
        ax = fig.add_subplot(2, 3, 1, projection='3d')
        ax.set_xlabel(r'$C_m$', fntdict)
        ax.set_ylabel(r'$c_v$', fntdict)
        ax.set_zlabel(r'$f(C_m,$ $c_v)$', fntdict)
        X, Y = np.meshgrid(Cm[arr], cv[arr], sparse=True, copy=False)
        Sum = np.empty((X+Y).shape)
        it = np.nditer([X, Y], flags=['multi_index'])
        while not it.finished:
            x = np.array([it[0], it[1], cf[sr], ca[sr]])
            Sum[it.multi_index] = sum_kv(x)
            it.iternext()
        ax.plot_surface(X, Y, Sum, linewidth=0, cmap=cm.coolwarm,
                        antialiased=False, shade=False)

        ax = fig.add_subplot(2, 3, 2, projection='3d')
        ax.set_xlabel(r'$C_m$', fntdict)
        ax.set_ylabel(r'$c_f$', fntdict)
        ax.set_zlabel(r'$f(C_m,$ $c_f)$', fntdict)
        X, Y = np.meshgrid(Cm[arr], cf[arr], sparse=True, copy=False)
        Sum = np.empty((X+Y).shape)
        it = np.nditer([X, Y], flags=['multi_index'])
        while not it.finished:
            x = np.array([it[0], cv[sr], it[1], ca[sr]])
            Sum[it.multi_index] = sum_kv(x)
            it.iternext()
        ax.plot_surface(X, Y, Sum, linewidth=0, cmap=cm.coolwarm,
                        antialiased=False, shade=False)

        ax = fig.add_subplot(2, 3, 3, projection='3d')
        ax.set_xlabel(r'$C_m$', fntdict)
        ax.set_ylabel(r'$c_a$', fntdict)
        ax.set_zlabel(r'$f(C_m,$ $c_a)$', fntdict)
        X, Y = np.meshgrid(Cm[arr], ca[arr], sparse=True, copy=False)
        Sum = np.empty((X+Y).shape)
        it = np.nditer([X, Y], flags=['multi_index'])
        while not it.finished:
            x = np.array([it[0], cv[sr], cf[sr], it[1]])
            Sum[it.multi_index] = sum_kv(x)
            it.iternext()
        ax.plot_surface(X, Y, Sum, linewidth=0, cmap=cm.coolwarm,
                        antialiased=False, shade=False)

        ax = fig.add_subplot(2, 3, 4, projection='3d')
        ax.set_xlabel(r'$c_v$', fntdict)
        ax.set_ylabel(r'$c_f$', fntdict)
        ax.set_zlabel(r'$f(c_v,$ $c_f)$', fntdict)
        X, Y = np.meshgrid(cv[arr], cf[arr], sparse=True, copy=False)
        Sum = np.empty((X+Y).shape)
        it = np.nditer([X, Y], flags=['multi_index'])
        while not it.finished:
            x = np.array([Cm[sr], it[0], it[1], ca[sr]])
            Sum[it.multi_index] = sum_kv(x)
            it.iternext()
        ax.plot_surface(X, Y, Sum, linewidth=0, cmap=cm.coolwarm,
                        antialiased=False, shade=False)

        ax = fig.add_subplot(2, 3, 5, projection='3d')
        ax.set_xlabel(r'$c_v$', fntdict)
        ax.set_ylabel(r'$c_a$', fntdict)
        ax.set_zlabel(r'$f(c_v,$ $c_a)$', fntdict)
        X, Y = np.meshgrid(cv[arr], ca[arr], sparse=True, copy=False)
        Sum = np.empty((X+Y).shape)
        it = np.nditer([X, Y], flags=['multi_index'])
        while not it.finished:
            x = np.array([Cm[sr], it[0], cf[sr], it[1]])
            Sum[it.multi_index] = sum_kv(x)
            it.iternext()
        ax.plot_surface(X, Y, Sum, linewidth=0, cmap=cm.coolwarm,
                        antialiased=False, shade=False)

        ax = fig.add_subplot(2, 3, 6, projection='3d')
        ax.set_xlabel(r'$c_f$', fntdict)
        ax.set_ylabel(r'$c_a$', fntdict)
        ax.set_zlabel(r'$f(c_f,$ $c_a)$', fntdict)
        X, Y = np.meshgrid(cf[arr], ca[arr], sparse=True, copy=False)
        Sum = np.empty((X+Y).shape)
        it = np.nditer([X, Y], flags=['multi_index'])
        while not it.finished:
            x = np.array([Cm[sr], cv[sr], it[0], it[1]])
            Sum[it.multi_index] = sum_kv(x)
            it.iternext()
        ax.plot_surface(X, Y, Sum, linewidth=0, cmap=cm.coolwarm,
                        antialiased=False, shade=False)

        self.fig = fig.tight_layout()


def SaveToFile(data, file_name):
    """
    Spremanje podataka
    """
    output = open(file_name, 'wb')
    pickle.dump(data, output)
    output.close


def ReadFromFile(file_name):
    """
    Ucitavanje podataka
    """
    pkl_file = open(file_name, 'rb')
    data = pickle.load(pkl_file)
    pkl_file.close()
    return data
