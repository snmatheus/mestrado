import os
import math
import numpy as np
import pandas as pd
import netCDF4 as nc
from output import plota_bat
from datetime import datetime
from scipy.special import gamma
from constants import path_pos, path_points, path_dados
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


def make_mask(lats, lons, df):
    mask = np.ones((len(lats), len(lons))) * np.nan
    points = df[['lat', 'lon']].values.tolist()
    for i, lon in enumerate(lons):
        for j, lat in enumerate(lats):
            if [lat, lon] in points:
                mask[j, i] = 1
                points.remove([lat, lon])
    return mask


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
wpd_means = list(load_variable('save.wpd_means').values())
clim_spd_means = climatology(spd_means)
clim_wpd_means = climatology(wpd_means)

clim_wnd100_means = load_variable('clim_wnd100_means.pickle')
clim_wnd100_std = load_variable('clim_wnd100_std.pickle')

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
masked_tec = make_mask(lats, lons, df_zee_tec)

clim_spd_means = clim_spd_means * masked_tec

clim_wpd_means = clim_wpd_means * masked_tec
clim_wnd100_means = clim_wnd100_means * masked_tec
clim_wnd100_std = clim_wnd100_std * masked_tec

#%% --- Eólica
TNH = 8760

k = (clim_wnd100_std/clim_wnd100_means)**-1.086
c = (clim_wnd100_means/gamma(1+(1/k)))

# Wind Annual Energy Production
AEPw = masked_tec * np.nan
for v in range(1, 25):
    EP = dic_turbine[v] * weibull(v, k, c) * TNH/1000
    AEPw = np.nansum([AEPw, EP], axis=0)
    # print(AEPw[106][58])

# Wind Available Energy
Pnw = float(data_turbine.power.mode()) # Potência nominal
AEw = masked_tec * Pnw * TNH/1000

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
Apv = 1 # Um painel tem 2 metros de área

# Solar Annual Energy Production:
Ttsh = TNH/2
AEPs = clim_spd_means * Apv * npv * Ttsh/1000

# Solar Available Energy:
Pns = PnSC # potência nominal escolhida: Silicon Crystalline cells
AEs = masked_tec * Pns * Apv * TNH/1000

# Solar Capacity Factor
CFs = (AEPs/AEs) * 100

#%% --- Wind + Solar
# Annual Energy Production
AEP = AEPw + AEPs
# Capacity Factor
CF = ((AEPw + AEPs)/(AEw + AEs)) * 100

#%%


#%%
plota_bat(lats, lons, AEPw, 'AEPw', 'Tecnico/', 'lime', 'NNE')
plota_bat(lats, lons, AEPw, 'AEPw', 'Tecnico/', 'lime', 'SSE')
plota_bat(lats, lons, CFw, 'CFw', 'Tecnico/', 'lime', 'NNE')
plota_bat(lats, lons, CFw, 'CFw', 'Tecnico/', 'lime', 'SSE')

plota_bat(lats, lons, AEPs, 'AEPs', 'Tecnico/', 'lime', 'NNE')
plota_bat(lats, lons, AEPs, 'AEPs', 'Tecnico/', 'lime', 'SSE')
plota_bat(lats, lons, CFs, 'CFs', 'Tecnico/', 'lime', 'NNE')
plota_bat(lats, lons, CFs, 'CFs', 'Tecnico/', 'lime', 'SSE')

plota_bat(lats, lons, AEP, 'AEP', 'Tecnico/', 'lime', 'NNE')
plota_bat(lats, lons, AEP, 'AEP', 'Tecnico/', 'lime', 'SSE')
plota_bat(lats, lons, CF, 'CF', 'Tecnico/', 'lime', 'NNE')
plota_bat(lats, lons, CF, 'CF', 'Tecnico/', 'lime', 'SSE')
