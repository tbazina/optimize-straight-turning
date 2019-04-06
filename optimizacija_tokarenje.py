# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 15:13:48 2015

@author: tomislav
"""


import numpy as np
from scipy.optimize import minimize

"""
******************************************************************************
                                    Tokarenje
******************************************************************************
"""

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
    
def duljina_prolaza(duljina_zahvata, duljina_ulaza, duljina_izlaza):
    """
    
    """
    return duljina_zahvata+duljina_ulaza+duljina_izlaza
    
def volumen_strugotine(presjek_strugotine, brzina_rezanja):
    """
    
    """
    return 1e3*presjek_strugotine*brzina_rezanja
    
def masa_strugotine(volumen_strugotine, gustoca_materijala):
    """
    
    """
    return 1e-9*volumen_strugotine*gustoca_materijala
    
def glavna_sila_rezanja(fs, A, Kgamma, Kv, Ka, Ki):
    """
    
    """
    return fs*A*Kgamma*Kv*Ka*Ki
    
def presjek_strugotine(dubina_rezanja, posmak):
    """
    
    """
    return dubina_rezanja*posmak
    
def dubina_rezanja(pocetni_promjer, zavrsni_promjer, dodatak_fino, broj_prolaza):
    """
    
    """
    return ((pocetni_promjer-zavrsni_promjer)/2.-dodatak_fino)/broj_prolaza
    
def debljina_strugotine(posmak, prisloni_kut):
    """
    
    """
    return posmak*np.sin(np.radians(prisloni_kut))
    
def spec_sila_rezanja(fs1x1, debljina_strugotine, eksponent_Kienzlea):
    """
    
    """
    return fs1x1/(debljina_strugotine)**eksponent_Kienzlea
    
def snaga_tokarenje(glavna_sila, brzina_rezanja):
    """
    
    """
    return glavna_sila*brzina_rezanja/60e3
    
def snaga_stroja(snaga_rezanja, stupanj_iskoristivosti):
    """
    
    """
    return snaga_rezanja/stupanj_iskoristivosti
    
def teorijska_hrapavost(posmak, radijus_alata):
    """
    
    """
    return posmak**2./8./radijus_alata
    
def posmak_Ra(prosjecno_odstupanje_profila, radijus_alata):
    """
    
    """
    return np.sqrt(32.*prosjecno_odstupanje_profila*radijus_alata)
    
def broj_prolaza(pocetni_promjer, zavrsni_promjer, dodatak_fino, dubina_rezanja):
    """
    
    """
    return ((pocetni_promjer-zavrsni_promjer)/2.-dodatak_fino)/dubina_rezanja
    
def postojanost_alata(vc, f, ap, Ct, kv, kf, ka):
    """
    
    """
    return Ct*vc**kv*f**kf*ap**ka
    
def strojno_vrijeme(d, l, v, f, ap, ip):
    """
    
    """
    return np.pi*l*(ip*d-2.*ap*(ip-1.))/v/f/1000.

"""
******************************************************************************
                     Optimiranje parametara obrade
******************************************************************************
Proizvodac plocice: TaeguTec
Tip: DNMG 150504 MT
Grade: TT5100 - CVD coated

Materijal obratka: Č.3990 - čelik za automate
Stroj: HiCell
Promjer šipke: 30 mm
Završni promjer: 22h10

x[0] = brzina rezanja pri gruboj obradi
x[1] = posmak pri gruboj obradi
x[2] = brzina rezanja pri finoj obradi
x[3] = posmak pri finoj obradi
"""

"""
--------Podaci
"""

mi = 'minimum'
ma = 'maksimum'
gr = 'grubo'
fi = 'fino'
s = 'stroj'
na = 'namjestanje'
pro = 'promjena'
pr = 'prolaz'
ul = 'ulaz'
iz = 'izlaz'
za = 'zahvat'
pri = 'prisloni'

d = {0: 30., fi: 21.958}
l = {ul: 1., iz: 1., za: 40.}
l[pr] = duljina_prolaza(l[za], l[ul], l[iz])
v = {gr: {mi: 70., ma: 350.}, fi: {mi: 70., ma: 350.}}
f = {gr: {mi: 0.15, ma: 0.4}, fi: {mi: 0.01, ma: 0.4}}
ap = {gr: {mi: 0.8, ma: 4.}, fi: 0.225}
Ra = 0.0032
r = {'alat': 0.4}
Ct = 273.**5
k = {'v': -5., 'f': -0.75, 'a': -1., 's': 1.3}
P = {s: 10.}
eta = 0.7
n = {s: {ma: 2500.}}
t = {na: 1./60., pro: 1.5}
fs = {'1x1': 1680.}
z1 = 0.28
kp = {pri: 93.}
F = {ma: 7000.}
ksi = {mi: 4., ma: 16.}

"""
--------Geometrijsko ograničenje
"""

ip = {mi: np.ceil(broj_prolaza(d[0], d[fi], ap[fi], ap[gr][ma])),
      ma: np.floor(broj_prolaza(d[0], d[fi], ap[fi], ap[gr][mi]))}
      
ip['arr'] = np.arange(ip[mi], ip[ma]+1., 1.)
ap['arr'] = dubina_rezanja(d[0], d[fi], ap[fi], ip['arr'])

"""
--------Jedinično vrijeme izrade
"""

def vrijeme_izrade(x, ts, tc, ds, df, ip, ap, apf, l, Ct, kv, kf, ka):
    """
    
    """
    return ((ip+1.)*ts+np.pi*l/1000./x[0]/x[1]*(ip*ds-2.*ap*(ip-1.))*
            (1.+tc/Ct/x[0]**kv/x[1]**kf/ap**ka)+np.pi*l*(df+2.*apf)/1000./x[2]/
            x[3]*(1.+tc/Ct/x[2]**kv/x[3]**kf/apf**ka))
            
def vrijeme_izrade_jac(x, ts, tc, ds, df, ip, ap, apf, l, Ct, kv, kf, ka):
    """
    
    """
    dfdx0 = (-np.pi*l/1000./x[0]**2./x[1]*(ip*ds-2.*ap*(ip-1.))*(1.+
            tc*(kv+1.)/Ct/x[0]**kv/x[1]**kf/ap**ka))
    
    dfdx1 = (-np.pi*l/1000./x[0]/x[1]**2.*(ip*ds-2.*ap*(ip-1.))*(1.+
            tc*(kf+1.)/Ct/x[0]**kv/x[1]**kf/ap**ka))
    
    dfdx2 = (-np.pi*l*(df+2.*apf)/1000./x[2]**2./x[3]*(1.+
            tc*(kv+1.)/Ct/x[2]**kv/x[3]**kf/apf**ka))
    
    dfdx3 = (-np.pi*l*(df+2.*apf)/1000./x[2]/x[3]**2.*(1.+
            tc*(kf+1.)/Ct/x[2]**kv/x[3]**kf/apf**ka))
    
    return np.array([dfdx0, dfdx1, dfdx2, dfdx3])
            
"""
--------Ograničenje brzine vrtnje alatnog stroja
"""

if v[fi][ma] > brzina_rezanja(d[fi], n[s][ma]):
    v[fi][ma] = brzina_rezanja(d[fi], n[s][ma])
    
if v[gr][ma] > brzina_rezanja(d[fi]+2.*ap[fi], n[s][ma]):
    v[gr][ma] = brzina_rezanja(d[fi]+2.*ap[fi], n[s][ma])
    
"""
--------Ograničenja zbog hrapavosti površine
"""

if f[fi][ma] > posmak_Ra(Ra, r['alat']):
    f[fi][ma] = posmak_Ra(Ra, r['alat'])
    
"""
--------Ograničenje zbog vitkosti strugotine
"""

for i in [mi, ma]:
    if f[fi][i] > (ap[fi]/ksi[mi]):
        f[fi][i] = (ap[fi]/ksi[mi])
    elif f[fi][i] < (ap[fi]/ksi[ma]):
        f[fi][i] = (ap[fi]/ksi[ma])
        
def vitkost_min(x, ap, ksi_min):
    """
    
    """
    return np.array([ap/ksi_min-x[1]])
    
def vitkost_min_jac(x, ap, ksi_min):
    """
    
    """
    return np.array([0., -1., 0., 0.])
    
def vitkost_max(x, ap, ksi_max):
    """
    
    """
    return np.array([x[1]-ap/ksi_max])
    
def vitkost_max_jac(x, ap, ksi_max):
    """
    
    """
    return np.array([0., 1., 0., 0.])
    
"""
--------Ograničenje postojanosti alata
Minimalna postojanost alata mora biti veća od vremena potrebnog za obradu
jednog komada
"""

def postojanost_grubo(x, ds, df, ip, ap, apf, l, Ct, kv, kf, ka):
    """
    
    """
    return np.array([Ct*x[0]**kv*x[1]**kf*ap**ka-np.pi*l/1000.*(1./x[0]/x[1]*
            (ip*ds-2.*ap*(ip-1.))+(df+2.*apf)/x[2]/x[3])])
            
def postojanost_grubo_jac(x, ds, df, ip, ap, apf, l, Ct, kv, kf, ka):
    """
    
    """
    dfdx0 = (Ct*kv*x[0]**(kv-1.)*x[1]**kf*ap**ka+np.pi*l/1000./x[0]**2./x[1]*
            (ip*ds-2.*ap*(ip-1.)))
            
    dfdx1 = (Ct*kf*x[0]**kv*x[1]**(kf-1.)*ap**ka+np.pi*l/1000./x[0]/x[1]**2.*
            (ip*ds-2.*ap*(ip-1.)))
            
    dfdx2 = np.pi*l*(df+2.*apf)/1000./x[2]**2./x[3]
    
    dfdx3 = np.pi*l*(df+2.*apf)/1000./x[2]/x[3]**2.
    
    return np.array([dfdx0, dfdx1, dfdx2, dfdx3])
    
def postojanost_fino(x, ds, df, ip, ap, apf, l, Ct, kv, kf, ka):
    """
    
    """
    return np.array([Ct*x[2]**kv*x[3]**kf*apf**ka-np.pi*l/1000.*(1./x[0]/x[1]*
            (ip*ds-2.*ap*(ip-1.))+(df+2.*apf)/x[2]/x[3])])
            
def postojanost_fino_jac(x, ds, df, ip, ap, apf, l, Ct, kv, kf, ka):
    """
    
    """
    dfdx0 = np.pi*l/1000./x[0]**2./x[1]*(ip*ds-2.*ap*(ip-1.))
    
    dfdx1 = np.pi*l/1000./x[0]/x[1]**2.*(ip*ds-2.*ap*(ip-1.))
    
    dfdx2 = (Ct*kv*x[2]**(kv-1.)*x[3]**kf*apf**ka+
            np.pi*l*(df+2.*apf)/1000./x[2]**2./x[3])
    
    dfdx3 = (Ct*kf*x[2]**kv*x[3]**(kf-1.)*apf**ka+
            np.pi*l*(df+2.*apf)/1000./x[2]/x[3]**2.)
    
    return np.array([dfdx0, dfdx1, dfdx2, dfdx3])

"""
--------Ograničenje zbog sile na alatu
Sila na alatu ne smije bit veća od maksimalne dozvoljene koju propisuje
proizvođač alata
"""

def sila_alat_grubo(x, ap, Fmax, kp, ks, fs1x1, z1):
    """
    
    """
    return np.array([Fmax*np.sin(np.radians(kp))**z1/ks/fs1x1-x[1]**(1.-z1)*ap])
    
def sila_alat_grubo_jac(x, ap, Fmax, kp, ks, fs1x1, z1):
    """
    
    """
    return np.array([0., (z1-1.)*x[1]**(-z1)*ap, 0., 0.])
    
"""
--------Ograničenje zbog instalirane snage stroja
"""

def snaga_grubo(x, ap, Ps, eta, kp, fs1x1, z1):
    """
    
    """
    return np.array([60e3*Ps*eta*np.sin(np.radians(kp))**z1/fs1x1-
                    x[0]*x[1]**(1.-z1)*ap])

def snaga_grubo_jac(x, ap, Ps, eta, kp, fs1x1, z1):
    """
    
    """
    return np.array([-x[1]**(1.-z1)*ap, (z1-1.)*x[0]*x[1]**(-z1)*ap, 0., 0.])

"""
--------Ukupna lista granica
"""

bnds = ((v[gr][mi], v[gr][ma]), (f[gr][mi], f[gr][ma]), (v[fi][mi], v[fi][ma]),
        (f[fi][mi], f[fi][ma]))
    
"""
--------Ukupna lista ograničenja
"""

cons = {}
j = 0
for i in ip['arr']:
    cons[j] = ({'type': 'ineq', 'fun': vitkost_min, 'jac': vitkost_min_jac,
                'args': (ap['arr'][j], ksi[mi])},
               {'type': 'ineq', 'fun': vitkost_max, 'jac': vitkost_max_jac,
                'args': (ap['arr'][j], ksi[ma])},
               {'type': 'ineq', 'fun': postojanost_grubo,
                'jac': postojanost_grubo_jac,
                'args': (d[0], d[fi], i, ap['arr'][j], ap[fi], l[pr], Ct, 
                         k['v'], k['f'], k['a'])},
               {'type': 'ineq', 'fun': postojanost_fino,
                'jac': postojanost_fino_jac,
                'args': (d[0], d[fi], i, ap['arr'][j], ap[fi], l[pr], Ct,
                         k['v'], k['f'], k['a'])},
               {'type': 'ineq', 'fun': sila_alat_grubo,
                'jac': sila_alat_grubo_jac,
                'args': (ap['arr'][j], F[ma], kp[pri], k['s'], fs['1x1'], z1)},
               {'type': 'ineq', 'fun': snaga_grubo, 'jac': snaga_grubo_jac,
                'args': (ap['arr'][j], P[s], eta, kp[pri], fs['1x1'], z1)}
              )
    j+=1

"""
--------Početne vrijednosti
"""

init = (v[gr][mi], f[gr][mi], v[fi][mi], f[fi][mi])

"""
--------Traženje minimuma funkcije
"""

res = {}
class dummyclass():
    fun = np.float('inf')
    success = False
x = dummyclass()

j = 0
for i in ip['arr']:
    res[j] = minimize(vrijeme_izrade, init, args=(t[na], t[pro], d[0], d[fi], i,
             ap['arr'][j], ap[fi], l[pr], Ct, k['v'], k['f'], k['a']),
             method='SLSQP', jac=vrijeme_izrade_jac, bounds=bnds,
             constraints=cons[j], options={'disp': True, 'maxiter': 100})
    if (res[j].success == True) and (res[j].fun < x.fun):
        x = res[j]
        x.x = np.append(x.x, np.array([i, ap['arr'][j]]))
    j+=1

"""
--------Ispis rezultata
"""

print '\n\n\t\tSuccess:  {}'.format(x.success)


T = {gr: postojanost_alata(x.x[0], x.x[1], x.x[5], Ct, k['v'], k['f'], k['a']),
     fi: postojanost_alata(x.x[2], x.x[3], ap[fi], Ct, k['v'], k['f'], k['a'])}

t[s] = {gr: strojno_vrijeme(d[0], l[pr], x.x[0], x.x[1], x.x[5], x.x[4]),
        fi: strojno_vrijeme(d[fi]+2.*ap[fi], l[pr], x.x[2], x.x[3], ap[fi], 1.)}

n[gr] = brzina_vrtnje(d[fi]+2.*ap[fi], x.x[0])
n[fi] = brzina_vrtnje(d[fi], x.x[2])

A = {gr: presjek_strugotine(x.x[5], x.x[1]),
     fi: presjek_strugotine(ap[fi], x.x[3])}

V = {gr: volumen_strugotine(A[gr], x.x[0]),
     fi: volumen_strugotine(A[fi], x.x[2])}
     
F[gr] = glavna_sila_rezanja(spec_sila_rezanja(fs['1x1'],
        debljina_strugotine(x.x[1], kp[pri]), z1), A[gr], 1., 1., 1., 1.)
        
F[fi] = glavna_sila_rezanja(spec_sila_rezanja(fs['1x1'],
        debljina_strugotine(x.x[3], kp[pri]), z1), A[fi], 1., 1., 1., 1.)

P[gr] = snaga_tokarenje(F[gr], x.x[0])
P[fi] = snaga_tokarenje(F[fi], x.x[2])

P[mi] = {gr: snaga_stroja(P[gr], eta), fi: snaga_stroja(P[fi], eta)}


print '\nGruba obrada:\n'
print '\tBroj prolaza: {:.1f}'.format(x.x[4])
print '\tDubina rezanja u svakom prolazu: {:.3f} mm'.format(x.x[5])
print '\tBrzina rezanja: {:.0f} m/min'.format(x.x[0])
print '\tBrzina vrtnje na najmanjem promjeru: {:.0f} okr/min'.format(n[gr])
print '\tPosmak: {:.3f} mm/okr'.format(x.x[1])
print '\tPostojanost alata: {:.3f} min'.format(T[gr])
print '\tOdnos ap/f: {:.2f}'.format(x.x[5]/x.x[1])
print '\tGlavno strojno vrijeme: {:.1f} s'.format(t[s][gr]*60.)
print '\tPovršina presjeka strugotine: {:.2f} mm2'.format(A[gr])
print '\tVolumen skinute strugotine: {:.1f} mm3/min'.format(V[gr])
print '\tGlavna sila rezanja: {:.1f} N'.format(F[gr])
print '\tSnaga potrebna za tokarenje: {:.2f} kW'.format(P[gr])
print '\tMinimalna potrebna snaga stroja: {:.2f} kW'.format(P[mi][gr])


print '\nFina obrada:\n'
print '\tDodatak za finu obradu: {:.3f} mm'.format(ap[fi])
print '\tBrzina rezanja: {:.0f} m/min'.format(x.x[2])
print '\tBrzina vrtnje: {:.0f} okr/min'.format(n[fi])
print '\tPosmak: {:.3f} mm/okr'.format(x.x[3])
print '\tPostojanost alata: {:.3f} min'.format(T[fi])
print '\tOdnos ap/f: {:.2f}'.format(ap[fi]/x.x[3])
print '\tGlavno strojno vrijeme: {:.1f} s'.format(t[s][fi]*60.)
print '\tPovršina presjeka strugotine: {:.2f} mm2'.format(A[fi])
print '\tVolumen skinute strugotine: {:.1f} mm3/min'.format(V[fi])
print '\tGlavna sila rezanja: {:.1f} N'.format(F[fi])
print '\tSnaga potrebna za tokarenje: {:.2f} kW'.format(P[fi])
print '\tMinimalna potrebna snaga stroja: {:.2f} kW'.format(P[mi][fi])






