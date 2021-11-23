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


# % % --- Constantes
TNH = 8760
Ttsh = 4380 # TNH/2
WMW = 1000000**-1
MWTW = 1000000**-1

lats = load_variable('save.lats')
lons = load_variable('save.lons')
df_zee_tec = pd.read_csv(path_points + 'ZEE1000_TEC.csv')
grid_points_in_zee = len(df_zee_tec.index) # Quantidade de pontos de grade dentro da ZEE (critérios técnicos)

# % % --- Turbinas
rotor_diameter = 240 # m, diâmetro do rotor
turbines_distance = 7 * rotor_diameter # Distância entre torres dada pelo valor de 7x o diâmetro do rotor
At = math.pi * (rotor_diameter/2)**2 # Turbine sweep area 
Cp = 0.98 # Generator efficiency
data_turbine = pd.read_csv('../dados/power15M.csv', names=['v','power'])
dic_turbine = dict(zip(data_turbine.v, data_turbine.power))
Pnw = float(data_turbine.power.mode()) # 15 MW = Potência nominal de cada turbina

# % % --- Painéis - Silicon Crystalline cells
npv = 0.268 # = eficiência
Pns = 0.00022110215053763443 # MW [98.7/446400] = Potência nominal de cada painel
As = 2

# % % Interação com Turbinas e Painés
'''
 A área disponível para os painéis solares é equivalente a área circular
 ao redor de cada turbina com metade do raio de varredura da turbina(metade 
 de 7D. Sendo assim, é calculada a partir de um raio equivalente à 1/4 da
 distância entre as turbinas.
'''
sweep_radius = turbines_distance/4
panel_area = math.pi * (sweep_radius**2)

panels_per_turbines = int(panel_area/2) # Quantidade de painéis de 2 metros de área ao redor de 1 turbina
turbines_per_points = 16 * 16 # Quantidade de turbinas de 7D de distância ao redor de cada ponto de grade
panels_per_points = panels_per_turbines * turbines_per_points # Quantidade de painéis por pontos de grade

# % % Funções
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


def min(array):
    return np.nanmin(array)


def max(array):
    return np.nanmax(array)


def filter_by_points(df_base, df):
    points_base = df_base[['lat', 'lon']].values.tolist()
    points = df[['lat', 'lon']].values.tolist()
    loutput = []
    for latlon in points:
        if latlon in points_base:
            loutput.append(latlon)
    df_output = pd.DataFrame(loutput, columns=['lat','lon'])
    return df_output


def make_mask(df):
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
    return w

#%% Leitura de variáveis climatológicas (vento, radiação, densidades de potência eólica e solar)
wpd_means = list(load_variable('save.wpd_means').values())
spd_means = list(load_variable('save.spd_means').values())
clim_wpd_means = climatology(wpd_means)
clim_spd_means = climatology(spd_means)
clim_wnd100_means = load_variable('clim_wnd100_means.pickle')
clim_wnd100_std = load_variable('clim_wnd100_std.pickle')

#%% Definindo e aplicando a máscara de critérios técnicos aos campos
masked_tec = make_mask(df_zee_tec)
clim_spd_means = clim_spd_means * masked_tec
clim_wpd_means = clim_wpd_means * masked_tec
clim_wnd100_means = clim_wnd100_means * masked_tec
clim_wnd100_std = clim_wnd100_std * masked_tec

#%% --- Potêncial técnico eólico
k = (clim_wnd100_std/clim_wnd100_means)**-1.086
c = (clim_wnd100_means/(gamma(1+(1/k))))

print('mean - Mínimo: {}, Máximo; {}'.format(min(clim_wnd100_means), max(clim_wnd100_means)))
print('std - Mínimo: {}, Máximo; {}'.format(min(clim_wnd100_std), max(clim_wnd100_std)))
print('k - Mínimo: {}, Máximo; {}'.format(min(k), max(k)))
print('c - Mínimo: {}, Máximo; {}'.format(min(c), max(c)))

EP = [dic_turbine[v] * weibull(v, k, c) for v in range(0, 25)] # MW, Produced Energy
AEPw = masked_tec * turbines_per_points * np.nansum(EP, axis=0) * TNH # MW * h = MWh, Wind Annual Energy Production
AEw = masked_tec * turbines_per_points * (Pnw * TNH) # MW * h = MWh, Wind Available Energy
CFw = (AEPw/AEw) * 100 # %, Wind Capacity Factor

#%% --- Potêncial técnico solar
AEPs = masked_tec * panels_per_points * (clim_spd_means * As * Ttsh * npv) * WMW  # W/m² * m² * h * 1e-06 = MWh, Solar Annual Energy Production
AEs = masked_tec * panels_per_points * (Pns * Ttsh) # MW * h = MWh, Solar Available Energy
CFs = (AEPs/AEs) * 100 # %, Solar Capacity Factor

#%% --- Potêncial técnico eólico e solar
AEP = AEPw + AEPs # MWh, Annual Energy Production
CF = (CFw + CFs)/2 # %, Capacity Factor

#%% Convertendo para GWh
AEPw = AEPw * MWTW
AEPs = AEPs * MWTW
AEP = AEP * MWTW

#%% Criação dos mapas
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


#%%
df_20 = pd.read_csv(path_points + '20_points.csv')
df_2050 = pd.read_csv(path_points + '2050_points.csv')
df_50100 = pd.read_csv(path_points + '50100_points.csv')
df_1001000 = pd.read_csv(path_points + '1001000_points.csv')

fdf_20 = filter_by_points(df_zee_tec, df_20)
fdf_2050 = filter_by_points(df_zee_tec, df_2050)
fdf_50100 = filter_by_points(df_zee_tec, df_50100)
fdf_1001000 = filter_by_points(df_zee_tec, df_1001000)

number_points_20m = len(fdf_20.index)
number_points_2050m = len(fdf_2050.index)
number_points_50100m = len(fdf_50100.index)
number_points_1001000m = len(fdf_1001000.index)
number_points_zee = len(df_zee_tec.index)

panel_area_per_points = panels_per_points * As

printa(number_points_20m, number_points_2050m, number_points_50100m,
           number_points_1001000m, number_points_zee, turbines_per_points,
           panel_area_per_points)
    

for bat, df in zip(['20', '2050', '50100', '1001000', 'zee'],
                   [df_20, df_2050, df_50100, df_1001000, df_zee_tec]):
    masked_tec = make_mask(df)
    # plota_bat(lats, lons, AEPw * masked_tec, 'AEPw_' + bat, 'Tecnico/', 'lime', 'NNE')
    # plota_bat(lats, lons, AEPw * masked_tec, 'AEPw_' + bat, 'Tecnico/', 'lime', 'SSE')
    # plota_bat(lats, lons, AEPs * masked_tec, 'AEPs_' + bat, 'Tecnico/', 'lime', 'NNE')
    # plota_bat(lats, lons, AEPs * masked_tec, 'AEPs_' + bat, 'Tecnico/', 'lime', 'SSE')    
    # plota_bat(lats, lons, AEP * masked_tec, 'AEP_' + bat, 'Tecnico/', 'lime', 'NNE')
    # plota_bat(lats, lons, AEP * masked_tec, 'AEP_' + bat, 'Tecnico/', 'lime', 'SSE')

    sum_AEPw = np.nansum(AEPw * masked_tec)
    sum_AEPs = np.nansum(AEPs * masked_tec)
    sum_AEP = np.nansum(AEP * masked_tec)

    mean_CFw = np.nanmean(CFw * masked_tec)
    mean_CFs = np.nanmean(CFs * masked_tec)
    mean_CF = np.nanmean(CF * masked_tec)

    print('\nBatimetria:', bat)
    print('\tAEPw', round(sum_AEPw), 'TWh')
    print('\tAEPs', round(sum_AEPs), 'TWh')
    print('\tAEP', round(sum_AEP), 'TWh\n')

    print('\tmean_CFw', round(mean_CFw), '%')
    print('\tmean_CFs', round(mean_CFs), '%')
    print('\tmean_CF', round(mean_CF), '%\n')



# AEPw -> abaixo de 50.000
# CF aceitável > 45% 
# kg/m³ * m³/s³ => (kg*m³)/(m³*s³) => (kg*m*m²)/(s*m³*s²) => (kg*m/s²) * m²/(s*m³) => N * m²/(s*m³) =>
#=> N * m/(s*m²) => (N*m/s) * (1/m²) => (J/s) * (1/m²) => W * (1/m²) => W/m², Wind Power Density
