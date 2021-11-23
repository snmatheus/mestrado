import os
import math
import numpy as np
import pandas as pd
from scipy.special import gamma
from output import plota_bat_m, plota_ts_tec
from constants import path_pos, path_points, path_dados
from utils import load_variable, extrapolate, climatology_m, find_nearest


# % % --- Constantes
TNH = 744
Ttsh = TNH/2
WMW = 1000000**-1
MWGW = 1000**-1

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
spd_means = list(load_variable('save.m_spd_means').values())
wpd_means = list(load_variable('save.m_wpd_means').values())
m_wpd_means = climatology_m(wpd_means)
m_spd_means = climatology_m(spd_means)
m_wnd100_means = np.array(list(load_variable('m_wnd100_means.pickle').values()))
m_wnd100_std = np.array(list(load_variable('m_wnd100_std.pickle').values()))

#%% Definindo e aplicando a máscara de critérios técnicos aos campos
masked_tec = make_mask(df_zee_tec)
m_spd_means = m_spd_means * masked_tec
m_wpd_means = m_wpd_means * masked_tec
m_wnd100_means = m_wnd100_means * masked_tec
m_wnd100_std = m_wnd100_std * masked_tec

#%% --- Potêncial técnico eólico
k = (m_wnd100_std/m_wnd100_means)**-1.086
c = (m_wnd100_means/(gamma(1+(1/k))))

print('mean - Mínimo: {}, Máximo; {}'.format(min(m_wnd100_means), max(m_wnd100_means)))
print('std - Mínimo: {}, Máximo; {}'.format(min(m_wnd100_std), max(m_wnd100_std)))
print('k - Mínimo: {}, Máximo; {}'.format(min(k), max(k)))
print('c - Mínimo: {}, Máximo; {}'.format(min(c), max(c)))

EP = [dic_turbine[v] * weibull(v, k, c) for v in range(0, 25)] # MW, Produced Energy
AEPw = masked_tec * turbines_per_points * np.nansum(EP, axis=0) * TNH # MW * h = MWh, Wind Annual Energy Production
AEw = masked_tec * turbines_per_points * (Pnw * TNH) # MW * h = MWh, Wind Available Energy
CFw = (AEPw/AEw) * 100 # %, Wind Capacity Factor

#%% --- Potêncial técnico solar
AEPs = masked_tec * panels_per_points * (m_spd_means * As * Ttsh * npv) * WMW  # W/m² * m² * h * 1e-06 = MWh, Solar Annual Energy Production
AEs = masked_tec * panels_per_points * (Pns * Ttsh) # MW * h = MWh, Solar Available Energy
CFs = (AEPs/AEs) * 100 # %, Solar Capacity Factor

#%% --- Potêncial técnico eólico e solar
AEP = AEPw + AEPs # MWh, Annual Energy Production
CF = (CFw + CFs)/2 # %, Capacity Factor

#%% Convertendo para GWh
AEPw = AEPw * MWGW
AEPs = AEPs * MWGW
AEP = AEP * MWGW

#%%

latlons_regions = {
    'North': [-51.25, -46, -0.75, 5.25],
    'Northeast': [-46, -33.25, -17.75, 0.75],
    'Southeast': [-44.5, -36, -25.5, -17.75],
    'South': [-53, -45, -34.5, -25.5],
}

def filter_region(data, region):
    i1 = find_nearest(lons, latlons_regions[region][0])
    i2 = find_nearest(lons, latlons_regions[region][1])
    j1 = find_nearest(lats, latlons_regions[region][3])
    j2 = find_nearest(lats, latlons_regions[region][2])
    return data[:, j1:j2, i1:i2]

generations = pd.read_csv('../dados/Gerações.csv')

for region in list(latlons_regions.keys()):
    mean_AEP = [np.nanmean(filter_region(AEP, region)[i]) for i in range(0, 12)]
    plota_ts_tec(generations, mean_AEP, 'AEP', 'Tecnico/m_ts_', region)

#%% Criação dos mapas
plota_bat_m(lats, lons, AEP, latlons_regions, 'AEP', 'Tecnico/m_')
plota_bat_m(lats, lons, CF, latlons_regions, 'CF', 'Tecnico/m_')

