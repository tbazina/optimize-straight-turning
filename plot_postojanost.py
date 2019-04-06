# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 15:42:57 2016

@author: tomislav
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from func import *

majorLocator = MultipleLocator(5)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(2.5)
t = PostojanostAlata()

vg_a = np.linspace(90., 120., 100)
vg = 104.45
vf_a = np.linspace(150., 180., 100)
vf = 172.53
fg_a = np.linspace(0.25, 0.35, 100)
fg = 0.35
ff_a = np.linspace(0.15, 0.25, 100)
ff = 0.233
ag_a = np.linspace(2., 3.5, 100)
ag = 3.284
af_a = np.linspace(0.8, 2., 100)
af = 0.933

res = ReadFromFile('VDF.pkl')

fig = plt.figure(figsize=(8, 11), dpi=300)
fig.suptitle('Postojanost alata u ovisnosti o parametrima obrade', fontsize=14)

sp1 = fig.add_subplot(221)
sp1.plot(vg_a, t(res.x, vg_a, fg, ag), label='Gruba obrada', lw=3)
sp1.plot(vf_a, t(res.x, vf_a, ff, af), label='Fina obrada', lw=3)
sp1.set_ylabel(u'Postojanost alata $T$ [min]')
sp1.set_xlabel(u'Brzina rezanja $v$ [m/min]')
sp1.set_title(u'$T_R(v_R, f_R = 0.35, a_R = 3.284)$, '
              '$T_F(v_F, f_F = 0.233, a_F = 0.933)$')
sp1.legend()
sp1.grid(which=u'both')
sp1.yaxis.set_major_locator(MultipleLocator(2.5))
sp1.yaxis.set_minor_locator(MultipleLocator(1.25))
sp1.xaxis.set_minor_locator(MultipleLocator(5))

sp2 = fig.add_subplot(222)
sp2.plot(fg_a, t(res.x, vg, fg_a, ag), label='Gruba obrada', lw=3)
sp2.plot(ff_a, t(res.x, vf, ff_a, af), label='Fina obrada', lw=3)
sp2.set_ylabel(u'Postojanost alata $T$ [min]')
sp2.set_xlabel(u'Posmak $f$ [mm/okr]')
sp2.set_title(u'$T_R(v_R = 104.45, f_R, a_R = 3.284)$, '
              '$T_F(v_F = 172.53, f_F, a_F = 0.933)$')
sp2.legend()
sp2.grid(which=u'both')
sp2.yaxis.set_major_locator(MultipleLocator(5))
sp2.yaxis.set_minor_locator(MultipleLocator(2.5))
sp2.xaxis.set_minor_locator(MultipleLocator(0.025))


sp3 = fig.add_subplot(223)
sp3.plot(ag_a, t(res.x, vg, fg, ag_a), label='Gruba obrada', lw=3)
sp3.plot(af_a, t(res.x, vf, ff, af_a), label='Fina obrada', lw=3)
sp3.set_ylabel(u'Postojanost alata $T$ [min]')
sp3.set_xlabel(u'Dubina rezanja $a$ [mm]')
sp3.set_title(u'$T_R(v_R = 104.45, f_R = 0.35, a_R)$, '
              '$T_F(v_F = 172.53, f_F = 0.233, a_F)$')
sp3.legend()
sp3.grid(which=u'both')
sp3.yaxis.set_major_locator(MultipleLocator(2.5))
sp3.yaxis.set_minor_locator(MultipleLocator(1.25))
sp3.xaxis.set_minor_locator(MultipleLocator(0.25))

plt.savefig('post_parametri.png', dpi=300, bbox_inches='tight',
            transparent=False)
