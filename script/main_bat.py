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


def weibull(v,k,c):
    return (k/c) * (v/c)**(k-1) * np.exp(-(v/c)**k)


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
data_turbine = pd.read_csv('../dados/power15M.csv',
                           sep=';', names=['v','power'])

df_20 = pd.read_csv(path_points + '20_points.csv')
df_2050 = pd.read_csv(path_points + '2050_points.csv')
df_50100 = pd.read_csv(path_points + '50100_points.csv')
df_1001000 = pd.read_csv(path_points + '1001000_points.csv')
df_zee_tec = pd.read_csv(path_points + 'ZEE1000_TEC.csv')
df_zee_1000 = pd.read_csv(path_points + 'ZEE1000.csv')

#%%
clim_spd_means = climatology(spd_means)

rotor_diameter = 240/1000 # km
sweep_diameter = turbines_distance = 7 * rotor_diameter
sweep_radius = sweep_diameter/2
panel_area = math.pi * (sweep_radius/2)**2

turbines_per_points = 16 * 16
panel_area_per_points = round(panel_area * turbines_per_points)

fdf_20 = filter_by_points(df_zee_1000, df_20)
fdf_2050 = filter_by_points(df_zee_1000, df_2050)
fdf_50100 = filter_by_points(df_zee_1000, df_50100)
fdf_1001000 = filter_by_points(df_zee_1000, df_1001000)

number_points_20m = len(fdf_20.index)
number_points_2050m = len(fdf_2050.index)
number_points_50100m = len(fdf_50100.index)
number_points_1001000m = len(fdf_1001000.index)
number_points_zee = len(df_zee_tec.index)

wind_mean = np.nansum([nne_mean, sse_mean])
wind_std = np.nansum([nne_std, sse_std])

#%%
k = (wind_std/wind_mean)**-1.086
c = (wind_mean/gamma(1+(1/k)))
 
nominal_power = float(data_turbine.power.mode())

TNH = 8760/1000

PEs = [row.power * weibull(row.v, k, c) * TNH
       for index, row in data_turbine.iterrows()]
AEPw = np.sum(PEs)
AEw = nominal_power * TNH
CFw = AEPw/AEw

As = 1000
Ttsh = TNH/2
Lr = 0.9
AEPs = spd_means * As * Lr * Ttsh * TNH
AEs = spd_means * As * TNH
CFs = AEPs/AEs

CF = CFw + CFs

#%%
mclim_wpd_means = filter_by_points2(lats, lons, clim_wpd_means, df_zee_tec)
plota_bat(lats, lons, mclim_wpd_means,
          'Clim_Average_WPD_ZEE_TEC', '', 'lime')

CF = filter_by_points2(lats, lons, CF, df_zee_1000)
plota_bat(lats, lons, mclim_wpd_means,
          'Clim_Average_WPD_ZEE_1000', '', 'lime')

print('- finished! Tempo:', datetime.now() - start)
