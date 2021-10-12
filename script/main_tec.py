import os
import math
import numpy as np
import pandas as pd
from output import plota_bat
from datetime import datetime
from scipy.special import gamma
from constants import path_pos, path_points
from utils import load_variable, extrapolate, climatology


def filter_by_points(df_base, df):
    points_base = df_base[['lat', 'lon']].values.tolist()
    points = df[['lat', 'lon']].values.tolist()
    loutput = []
    for latlon in points:
        if latlon in points_base:
            loutput.append(latlon)
    df_output = pd.DataFrame(loutput, columns=['lon','lat'])
    return df_output


def filter_by_points2(lats, lons, data, df):
    mdata = data * np.nan
    points = df[['lat', 'lon']].values.tolist()
    for i, lon in enumerate(lons):
        for j, lat in enumerate(lats):
            if [lat, lon] in points:
                mdata[j, i] = data[j, i]
                points.remove([lat, lon])
    return mdata  


def make_mask(lats, lons, data, df):
    mdata = data * np.nan
    points = df[['lat', 'lon']].values.tolist()
    for i, lon in enumerate(lons):
        for j, lat in enumerate(lats):
            if [lat, lon] in points:
                mdata[j, i] = 1
                points.remove([lat, lon])
    return mdata  


def weibull(v, k, c):
    w = np.ones(c.shape) * np.nan
    for i in range(w.shape[0]):
        for j in range(w.shape[1]):
            k_ = k[i, j]
            c_ = c[i, j]
            w[i, j] = (k_/c_) * (v/c_)**(k_-1) * np.exp(-(v/c_)**k_)
            # if not np.isnan(w[i, j]):
            #     print(k_, c_, w[i, j])
    return w

k_ = 0.3966
c_ = 1.0177
v = 0

def wpd(v):
    return 0.5*1.225*v**3


def printa(number_points_20m, number_points_2050m, number_points_50100m,
           number_points_1001000m, number_points_zee, turbines_per_points,
           panel_area_per_points):
    print('Número de pontos de grade:')
    print('\t0 - 20 m:', number_points_20m)
    print('\t20 - 50 m :', number_points_2050m)
    print('\t50 - 100 m:', number_points_50100m)
    print('\t100 - 1000 m:', number_points_1001000m)
    print('\tZEE:', number_points_zee)
    print('Número de turbinas:', )
    print('\t0 - 20 m:', number_points_20m * turbines_per_points)
    print('\t20 - 50 m :', number_points_2050m * turbines_per_points)
    print('\t50 - 100 m:', number_points_50100m * turbines_per_points)
    print('\t100 - 1000 m:', number_points_1001000m * turbines_per_points)
    print('\tZEE:', number_points_zee * turbines_per_points)
    print('\nÁrea coberta por painéis:', )
    print('\t0 - 20 m:', number_points_20m * panel_area_per_points)
    print('\t20 - 50 m :', number_points_2050m * panel_area_per_points)
    print('\t50 - 100 m:', number_points_50100m * panel_area_per_points)
    print('\t100 - 1000 m:', number_points_1001000m * panel_area_per_points)
    print('\tZEE:', number_points_zee * panel_area_per_points)



#%%
times = pd.date_range(start='1990', end='2020', freq='Y')
lats = load_variable('save.lats')
lons = load_variable('save.lons')
spd_means = list(load_variable('save.spd_means').values())

nne_mean = np.load('../dados/K_C/NNE_MED_WS100.npy')
sse_mean = np.load('../dados/K_C/SSE_STD_WS100.npy')
nne_std = np.load('../dados/K_C/NNE_STD_WS100.npy')
sse_std = np.load('../dados/K_C/SSE_MED_WS100.npy')   
data_turbine = pd.read_csv('../dados/power15M.csv', names=['v','power'])
dic_turbine = dict(zip(data_turbine.v, data_turbine.power))

df_zee_tec = pd.read_csv(path_points + 'ZEE1000_TEC.csv')
df_zee_1000 = pd.read_csv(path_points + 'ZEE1000.csv')

#%%
# Diâmetro do rotor:
rotor_diameter = 240
# Distância entre torres dada pelo valor de 7x o diâmetro do rotor:
turbines_distance = 7 * rotor_diameter

'''
 A área disponível para os painéis solares é equivalente a área circular
 ao redor de cada turbina com metade do raio de varredura da turbina(metade 
 de 7D. Sendo assim, é calculada a partir de um raio equivalente à 1/4 da
 distância entre as turbinas.
'''
sweep_radius = turbines_distance/4
panel_area = math.pi * (sweep_radius**2)

# Quantidade de painéis de 2 metros de área ao redor de 1 turbina:
panels_per_turbines = int(panel_area/2)

# Quantidade de turbinas de 7D de distância ao redor de cada ponto de grade:
turbines_per_points = 16 * 16

# Quantidade de painéis por pontos de grade:
panels_per_points = panels_per_turbines * turbines_per_points

# Quantidade de pontos de grade dentro da ZEE, considerando critérios técnicos:
grid_points_in_zee = len(df_zee_tec.index)

#%%
mask = np.ones((len(lats), len(lons))) * np.nan
masked_tec = make_mask(lats, lons, mask , df_zee_tec)

clim_wind_mean = np.nansum([nne_mean, sse_mean], axis=0) * masked_tec
clim_wind_std = np.nansum([nne_std, sse_std], axis=0) * masked_tec

clim_spd_means = climatology(spd_means) * masked_tec

TNH = 8760/1000


#%% --- Eólica
k = (clim_wind_std/clim_wind_mean)**-1.086
c = (clim_wind_mean/gamma(1+(1/k)))

# Wind Annual Energy Production
AEPw = masked_tec * np.nan
for v in range(1, 25):
    EP = dic_turbine[v] * weibull(v, k, c) * TNH
    AEPw = np.nansum([AEPw, EP], axis=0)
    # print(AEPw[106][58])

# Wind Available Energy
Pnw = float(data_turbine.power.mode()) # Potência nominal
AEw = masked_tec * Pnw * TNH

# Wind Capacity Factor
CFw = (AEPw/AEw) * 100


#%% --- Solar
# Silicon Amorphous cells
nSA = 0.104 # eficiencia
PnSA = 38.3 # fator de capacidade
# CdTe thin film cells
nCTF = 0.172 # eficiencia
PnCTF = 63 # fator de capacidade
# Silicon Crystalline cells
nSC = 0.268 # eficiencia
PnSC = 98.7 # fator de capacidade

npv = nSC # eficiência escolhida: Silicon Crystalline cells
Apv = 2 # Um painel tem 2 metros de área

# Solar Annual Energy Production:
Ttsh = TNH/2
AEPs = clim_spd_means * Apv * npv * Ttsh

# Solar Available Energy:
Pns = PnSC # potência nominal escolhida: Silicon Crystalline cells
AEs = masked_tec * Pns * Apv * TNH

# Solar Capacity Factor
CFs = (AEPs/AEs) * 100


#%% --- Wind + Solar
# Annual Energy Production
AEP = AEPw + AEPs
# Capacity Factor
CF = ((AEPw + AEPs)/(AEPw + AEPs)) * 100

#%%
print('Wind')
print(AEPw[106][58], AEw[106][58], CFw[106][58])
print('Solar')
print(AEPs[106][58], AEs[106][58], CFs[106][58])
print('Wind - Solar')
print(AEP[106][58], CF[106][58])
#%%
plota_bat(lats, lons, AEPw, 'AEPw', 'Tecnico/', 'lime')
plota_bat(lats, lons, AEPs, 'AEPs', 'Tecnico/', 'lime')
plota_bat(lats, lons, AEP, 'AEP', 'Tecnico/', 'lime')

plota_bat(lats, lons, CFw, 'CFw', 'Tecnico/', 'lime')
plota_bat(lats, lons, CFs, 'CFs', 'Tecnico/', 'lime')
plota_bat(lats, lons, CF, 'CF', 'Tecnico/', 'lime')