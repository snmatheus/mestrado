from output import plota_mapa, plota_mapa_anual, plota_area, plota_ts
import os
import math
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime
from constants import path_wnd, path_rad, wpdmin, spdmin, path_pos
from utils import load_variable, save_variable, wpd, extrapolate, climatology


start = datetime.now()


lats = []
lons = []

u_means = {}
v_means = {}
wnd10_means = {}
wnd100_means = {}
wpd_means = {}
spd_means = {}
solar_availabilitys = {}
wind_availabilitys = {}

wcss = {}
scws = {}
sinergys = {}

if os.path.isfile(path_pos + 'save.lats'):
    lats = load_variable('save.lats')
    lons = load_variable('save.lons')
    u_means = load_variable('save.u_means')
    v_means = load_variable('save.v_means')
    wnd10_means = load_variable('save.wnd10_means')
    wnd100_means = load_variable('save.wnd100_means')
    wpd_means = load_variable('save.wpd_means')
    spd_means = load_variable('save.spd_means')
    solar_availabilitys = load_variable('save.solar_availabilitys')
    wind_availabilitys = load_variable('save.wind_availabilitys')
    wcss = load_variable('save.wcss')
    scws = load_variable('save.scws')
    sinergys = load_variable('save.sinergys')
    times = pd.date_range(start='1990', end='2020', freq='Y')

else:
    w_file_list = [file_name for file_name in sorted(os.listdir(path_wnd))
                   if file_name.endswith('.nc')]
    r_file_list = [file_name for file_name in sorted(os.listdir(path_rad))
                   if file_name.endswith('.nc')]

    for wnd_file, rad_file in zip(w_file_list, r_file_list):
        year = wnd_file[7:11]
        print(year)
        wnd_ds = xr.open_dataset(path_wnd + wnd_file)
        rad_ds = xr.open_dataset(path_rad + rad_file)

        times = wnd_ds['time'].data
        tnh = len(times)

        if len(lats) == 0:
            lats = wnd_ds['latitude'].data
            lons = wnd_ds['longitude'].data

        r_data = rad_ds['ssrd'][:, :, :].data
        spd_data = r_data/3600

        u_data = wnd_ds['u10'][:, :, :].data
        v_data = wnd_ds['v10'][:, :, :].data

        wnd10_data = (u_data**2 + v_data**2)**0.5
        wnd100_data = extrapolate(wnd10_data)
        wpd_data = wpd(wnd100_data)

        rad_ds.close()
        wnd_ds.close()

        u_means[year] = np.mean(u_data, axis=0)
        v_means[year] = np.mean(v_data, axis=0)
        wnd10_means[year] = np.mean(wnd10_data, axis=0)
        wnd100_means[year] = np.mean(wnd100_data, axis=0)

        spd_means[year] = np.mean(spd_data, axis=0)
        wpd_means[year] = np.mean(wpd_data, axis=0)

        nh_spd_available = spd_data.copy()
        nh_spd_available[spd_data < spdmin] = 0
        nh_spd_available[spd_data >= spdmin] = 1
        nh_spd_available = np.sum(nh_spd_available, axis=0)
        solar_availability = (nh_spd_available/tnh) * 100
        solar_availabilitys[year] = solar_availability

        nh_wpd_available = wpd_data.copy()
        nh_wpd_available[wpd_data < wpdmin] = 0
        nh_wpd_available[wpd_data >= wpdmin] = 1
        nh_wpd_available = np.sum(nh_wpd_available, axis=0)
        wind_availability = (nh_wpd_available/tnh) * 100
        wind_availabilitys[year] = wind_availability

        nh_wpd_ge_spd_le = np.zeros(wpd_data.shape)
        nh_wpd_ge_spd_le[(wpd_data >= wpdmin) & (spd_data < spdmin)] = 1
        nh_wpd_ge_spd_le = np.sum(nh_wpd_ge_spd_le, axis=0)
        wcss[year] = nh_wpd_ge_spd_le/tnh * 100

        nh_wpd_le_spd_ge = np.zeros(wpd_data.shape)
        nh_wpd_le_spd_ge[(wpd_data < wpdmin) & (spd_data >= spdmin)] = 1
        nh_wpd_le_spd_ge = np.sum(nh_wpd_le_spd_ge, axis=0)
        scws[year] = nh_wpd_le_spd_ge/tnh * 100

        nh_wpd_ge_spd_ge = np.zeros(wpd_data.shape)
        nh_wpd_ge_spd_ge[(wpd_data >= wpdmin) & (spd_data >= spdmin)] = 1
        nh_wpd_ge_spd_ge = np.sum(nh_wpd_ge_spd_ge, axis=0)
        sinergys[year] = nh_wpd_ge_spd_ge/tnh * 100

    save_variable(lats, 'save.lats')
    save_variable(lons, 'save.lons')
    save_variable(u_means, 'save.u_means')
    save_variable(v_means, 'save.v_means')
    save_variable(wnd10_means, 'save.wnd10_means')
    save_variable(wnd100_means, 'save.wnd100_means')
    save_variable(wpd_means, 'save.wpd_means')
    save_variable(spd_means, 'save.spd_means')
    save_variable(solar_availabilitys, 'save.solar_availabilitys')
    save_variable(wind_availabilitys, 'save.wind_availabilitys')
    save_variable(wcss, 'save.wcss')
    save_variable(scws, 'save.scws')
    save_variable(sinergys, 'save.sinergys')

clim_u_means = climatology(u_means)
clim_v_means = climatology(v_means)
clim_wnd10_means = climatology(wnd10_means)
clim_wnd100_means = climatology(wnd100_means)
clim_wpd_means = climatology(wpd_means)
clim_spd_means = climatology(spd_means)
clim_wind_availabilitys = climatology(wind_availabilitys)
clim_solar_availabilitys = climatology(solar_availabilitys)
clim_wcss = climatology(wcss)
clim_scws = climatology(scws)
clim_sinergys = climatology(sinergys)

# %%
plota_area(lats, lons, 'Time_Series_points', 'ts/')
wpd_points, spd_points = plota_ts(lats, lons, times, wpd_means, spd_means,
                                  '%Y', 'Ano', 'TS_Anual', 'ts/')
np.save('TS.Anual.npy', [times, wpd_points, spd_points])

# %%

plota_mapa_anual(lats, lons, wnd10_means,
                 'Annual_Average_Wind_Speed_at_10_m', 'Anual/', 'lime',
                 u_means, v_means)

plota_mapa_anual(lats, lons, wnd100_means,
                 'Annual_Average_Wind_Speed_at_100_m', 'Anual/', 'lime',
                 u_means, v_means)

plota_mapa_anual(lats, lons, wpd_means, 'Annual_Average_WPD', 'Anual/', 'lime')
plota_mapa_anual(lats, lons, spd_means, 'Annual_Average_SPD', 'Anual/', 'lime')

plota_mapa_anual(lats, lons, wind_availabilitys,
                 'Annual_WPD_Availability', 'Anual/', 'lime')
plota_mapa_anual(lats, lons, solar_availabilitys,
                 'Annual_SPD_Availability', 'Anual/', 'lime')

plota_mapa_anual(lats, lons, wcss, 'Annual_Wind_Complements_Solar', 'Anual/',
                 'lime')
plota_mapa_anual(lats, lons, scws, 'Annual_Solar_Complements_Wind', 'Anual/',
                 'lime')
plota_mapa_anual(lats, lons, sinergys, 'Annual_Wind_and_Solar_Synergy',
                 'Anual/', 'lime')

# %%
plota_mapa(lats, lons, clim_wnd10_means,
           'Clim_Average_Wind_Speed_at_10_m', '30anos/', 'lime',
           clim_u_means, clim_v_means)

plota_mapa(lats, lons, clim_wnd100_means,
           'Clim_Average_Wind_Speed_at_100_m', '30anos/', 'lime',
           clim_u_means, clim_v_means)

plota_mapa(lats, lons, clim_wpd_means, 'Clim_Average_WPD', '30anos/', 'lime')
plota_mapa(lats, lons, clim_spd_means, 'Clim_Average_SPD', '30anos/', 'lime')

plota_mapa(lats, lons, clim_wind_availabilitys,
           'Clim_WPD_Availability', '30anos/', 'lime')

plota_mapa(lats, lons, clim_solar_availabilitys,
           'Clim_SPD_Availability', '30anos/', 'lime')

plota_mapa(lats, lons, clim_wcss, 'Clim_Wind_Complements_Solar', '30anos/',
           'lime')

plota_mapa(lats, lons, clim_scws, 'Clim_Solar_Complements_Wind', '30anos/',
           'lime')

plota_mapa(lats, lons, clim_sinergys, 'Clim_Wind_and_Solar_Synergy', '30anos/',
           'lime')

print('- finished! Tempo:', datetime.now() - start)
