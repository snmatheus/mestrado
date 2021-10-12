#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 00:20:04 2020

@author: administrador
"""

import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy.special import gamma
import matplotlib.pyplot as plt
import SAIDA as OUT
from datetime import datetime
from constants import path_dados

start = datetime.now()

lats = OUT.lats
lons = OUT.lons

def ln(x):
    from math import log
    return log(x)

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def weib(v,k,c):
    return (k/c) * (v/c)**(k-1) * np.exp(-(v/c)**k)

def DPE(v):
    return 0.5*1.225*v**3

fll = open("./ZEE_pontos/SSE1000_TEC.csv")
ijs, latlons = [], []
for latlon in fll.readlines()[1:]:
    lns = float(latlon.replace("\n","").split(",")[1])
    lts = float(latlon.replace("\n","").split(",")[0])
    i = np.where(lons == lns)
    j = np.where(lats == lts)
    ijs.append((j[0][0], i[0][0]))
    latlons.append((lts, lns))
    # 
arq = "../DADOS/NPY/30ANOS/ROCI_medio_geral.npy"
tmeans_r = np.load(arq)
means_r = []
for ii, jj in ijs:
    means_r.append(tmeans_r[ii, jj])
medias_rad = pd.DataFrame({"latlons": latlons, "dp": means_r})


mlatlon, medias = [], {}
devpads = {}
for ano in range (1990, 2020):
    ano = str(ano)
    medias_ws100 = pd.read_csv(PATH_DADOS+"MD_DPS/medias_ws100_"+ano)
    for linha in medias_ws100.iterrows():
        cont = linha[0]
        if cont not in medias:
            medias[cont] = []
        mlatlon.append(linha[1].latlons)
        medias[cont].append(linha[1].dp)
    
    devpads_ws100 = pd.read_csv(PATH_DADOS+"MD_DPS/devpads_ws100_"+ano)
    for linha in devpads_ws100.iterrows():
        cont = linha[0]
        if cont not in devpads:
            devpads[cont] = []
        devpads[cont].append(linha[1].dp)

mmed, mstd = [], []
for i in range(0, len(devpads)):
    mmed.append(np.mean(medias[i]))
    mstd.append(np.mean(devpads[i]))

devpad = np.array(mstd)
med = np.array(mmed)
#%%
dados = pd.read_csv('../DADOS/power15M.csv', sep=";")

k = (devpad/med)**-1.086
c = (med/gamma(1+(1/k)))
#f = [weib(vm, k, c) for vm in np.arange(0, 28, .1)]

PEs = [(0.5*1.225*v**3) * weib(v,k,c) * (8760/1000) for v in np.arange(0, 28)]
WAEP = np.sum(PEs)

SPD = np.load("../DADOS/NPY/30ANOS/ROCI_medio_geral.npy")
As = 1000
Ttsh = 12*365
Lr = 0.9
SAEP = SPD * As * Ttsh * Lr/1000

wpd_rating = 15
Pwn = wpd_rating * 8760/1000 #MWh

spd_rating = 15
Psn = spd_rating * 8760/1000 #MWh

AE = Pwn+Psn

CF = (WAEP+SAEP)/(AE)

print('k', k)
print('c', c)
print('TAEP: ', TAEP)
print('AE: ', AE)
print('CF: ', CF)
#%%

OUT.plota_mapa2




print ('- finished! Tempo:', datetime.now() - start)

# plt.rcParams["figure.figsize"]=[8, 5]

# plt.rcParams["figure.figsize"]=[8, 4]
# plt.title("Power Curve", size=18)
# plt.ylabel('Turbine Power (MW)', size=14)
# plt.xlabel('wind speed (m/s)', size=14)

# plt.ylim(dados['Power (MW)'].min(), 16)
# plt.xlim([0, 27])
# plt.yticks(np.arange(1,17,2))
# plt.xticks(np.arange(0,28,1))

# plt.plot(dados['Wind (m/s)'], dados['Power (MW)'], label='IEC Class IB', color='blue')
# plt.legend(loc='best')
# plt.savefig("../SAIDAS/POWER_CURVE.png", bbox_inches='tight')
# plt.show()



# plt.title("Weibull", size=18)
# plt.xlabel('wind speed (m/s)', size=14)
# plt.ylabel('Probability', size=14)

# xvalues = np.arange(0, 28, .1)
# xlim = [0, 27]
# plt.xlim(xlim)
# plt.xticks(np.arange(0, 28))
# ylim = [0, np.max(f)]
# plt.ylim(ylim)
# yvalues = np.arange(0, np.max(f)+0.1, .05)
# yticks = [str(round(pd*100))+"%" for pd in yvalues]
# plt.yticks(yvalues, yticks)

# plt.plot(np.arange(0, 28, .1), f, linewidth=1.5, color="red" , label="k = "+str(round(k,1))+", c = "+str(round(c,1)))
# plt.legend(fontsize=16)
# plt.savefig("../SAIDAS/WEIBULL.png", bbox_inches='tight')
# plt.show()