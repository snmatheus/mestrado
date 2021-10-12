from output import plota_mapa_hourly, plota_area, plota_ts
import os
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime
from constants import path_wnd, path_rad, wpdmin, spdmin, path_pos, \
    selected_hours
from utils import load_variable, save_variable, wpd, extrapolate


start = datetime.now()

lats = load_variable('save.lats')
lons = load_variable('save.lons')

if os.path.isfile(path_pos + 'save.h_wnd10_means'):
    h_u_means = load_variable('save.h_u_means')
    h_v_means = load_variable('save.h_v_means')

    h_wnd10_means = load_variable('save.h_wnd10_means')
    h_wnd100_means = load_variable('save.h_wnd100_means')

    h_wpd_means = load_variable('save.h_wpd_means')
    h_spd_means = load_variable('save.h_spd_means')

    h_wpd_availabilitys = load_variable('save.h_wpd_availabilitys')
    h_spd_availabilitys = load_variable('save.h_spd_availabilitys')

    h_wpd_ge_spd_les = load_variable('save.h_wpd_ge_spd_les')
    h_wpd_le_spd_ges = load_variable('save.h_wpd_le_spd_ges')
    h_wpd_ge_spd_ges = load_variable('save.h_wpd_ge_spd_ges')
    times = pd.date_range(start='1990-01-01T00:00', end='1990-01-01T23:00', freq='H')

else:
    h_u_means = {}
    h_v_means = {}
    h_wnd10_means = {}
    h_wnd100_means = {}
    h_wpd_means = {}
    h_spd_means = {}
    h_wpd_availabilitys = {}
    h_spd_availabilitys = {}
    h_wpd_ge_spd_les = {}
    h_wpd_le_spd_ges = {}
    h_wpd_ge_spd_ges = {}

    shape = (1, len(lats), len(lons))
    tnh = 0

    w_file_list = [file_name for file_name in sorted(os.listdir(path_wnd))
                   if file_name.endswith('.nc')]
    r_file_list = [file_name for file_name in sorted(os.listdir(path_rad))
                   if file_name.endswith('.nc')]

    for wnd_file, rad_file in zip(w_file_list, r_file_list):

        year = wnd_file[7:11]

        wnd_ds = xr.open_dataset(path_wnd + wnd_file)
        rad_ds = xr.open_dataset(path_rad + rad_file)

        times64 = wnd_ds['time'].data
        times = np.array([pd.Timestamp(time).to_pydatetime()
                         for time in times64])

        times2 = [[time, time.hour] for time in times]
        df_time = pd.DataFrame(times2, columns=['datetime', 'hour'])

        for H in selected_hours:
            indexes = df_time[df_time.hour == H]

            time = times[indexes.index]

            r_data = rad_ds['ssrd'][indexes.index].data
            spd_data = r_data/3600

            u_data = wnd_ds['u10'][indexes.index].data
            v_data = wnd_ds['v10'][indexes.index].data
            wnd10_data = (u_data**2 + v_data**2)**0.5

            tnh = len(time)

            print('{:%Y-%m-%dT%Hhr} {:%Y-%m-%dT%Hhr}'.format(time[0], time[-1]))

            spd_data = r_data/3600
            wnd10_data = (u_data**2 + v_data**2)**0.5

            wnd100_data = extrapolate(wnd10_data)
            wpd_data = wpd(wnd100_data)

            nh_wpd_available = wpd_data.copy()
            nh_wpd_available[wpd_data < wpdmin] = 0
            nh_wpd_available[wpd_data >= wpdmin] = 1

            nh_spd_available = spd_data.copy()
            nh_spd_available[spd_data < spdmin] = 0
            nh_spd_available[spd_data >= spdmin] = 1

            nh_wpd_ge_spd_le = np.zeros(wpd_data.shape)
            nh_wpd_ge_spd_le[(wpd_data >= wpdmin) &
                             (spd_data < spdmin)] = 1

            nh_wpd_le_spd_ge = np.zeros(wpd_data.shape)
            nh_wpd_le_spd_ge[(wpd_data < wpdmin) &
                             (spd_data >= spdmin)] = 1

            nh_wpd_ge_spd_ge = np.zeros(wpd_data.shape)
            nh_wpd_ge_spd_ge[(wpd_data >= wpdmin) &
                             (spd_data >= spdmin)] = 1

            rad_ds.close()
            wnd_ds.close()

            h_u = np.mean(u_data, axis=0)
            h_v = np.mean(v_data, axis=0)

            h_wnd10 = np.mean(wnd10_data, axis=0)
            h_spd = np.mean(spd_data, axis=0)

            h_wpd_available = np.sum(nh_wpd_available, axis=0)/tnh * 100
            h_spd_available = np.sum(nh_spd_available, axis=0)/tnh * 100

            h_wpd_ge_spd_le = np.sum(nh_wpd_ge_spd_le, axis=0)/tnh * 100
            h_wpd_le_spd_ge = np.sum(nh_wpd_le_spd_ge, axis=0)/tnh * 100
            h_wpd_ge_spd_ge = np.sum(nh_wpd_ge_spd_ge, axis=0)/tnh * 100

            if H not in h_wnd10_means:
                h_u_means[H] = h_u.reshape(shape)
                h_v_means[H] = h_v.reshape(shape)

                h_wnd10_means[H] = h_wnd10.reshape(shape)
                h_spd_means[H] = h_spd.reshape(shape)

                h_wpd_availabilitys[H] = h_wpd_available.reshape(shape)
                h_spd_availabilitys[H] = h_spd_available.reshape(shape)

                h_wpd_ge_spd_les[H] = h_wpd_ge_spd_le.reshape(shape)
                h_wpd_le_spd_ges[H] = h_wpd_le_spd_ge.reshape(shape)
                h_wpd_ge_spd_ges[H] = h_wpd_ge_spd_ge.reshape(shape)

            else:
                h_u_means[H] = np.vstack((h_u_means[H], h_u.reshape(shape)))
                h_v_means[H] = np.vstack((h_v_means[H], h_v.reshape(shape)))

                h_wnd10_means[H] = np.vstack((h_wnd10_means[H],
                                              h_wnd10.reshape(shape)))

                h_spd_means[H] = np.vstack(
                    (h_spd_means[H], h_spd.reshape(shape)))

                h_wpd_availabilitys[H] = np.vstack((h_wpd_availabilitys[H],
                                                    h_wpd_available
                                                    .reshape(shape)))

                h_spd_availabilitys[H] = np.vstack((h_spd_availabilitys[H],
                                                    h_spd_available
                                                    .reshape(shape)))

                h_wpd_ge_spd_les[H] = np.vstack((h_wpd_ge_spd_les[H],
                                                 h_wpd_ge_spd_le
                                                .reshape(shape)))

                h_wpd_le_spd_ges[H] = np.vstack((h_wpd_le_spd_ges[H],
                                                 h_wpd_le_spd_ge
                                                .reshape(shape)))

                h_wpd_ge_spd_ges[H] = np.vstack((h_wpd_ge_spd_ges[H],
                                                 h_wpd_ge_spd_ge
                                                .reshape(shape)))

    for H in selected_hours:
        h_wnd100_means[H] = extrapolate(h_wnd10_means[H])
        h_wpd_means[H] = wpd(h_wnd100_means[H])

    save_variable(h_u_means, 'save.h_u_means')
    save_variable(h_v_means, 'save.h_v_means')
    save_variable(h_wnd10_means, 'save.h_wnd10_means')
    save_variable(h_wnd100_means, 'save.h_wnd100_means')
    save_variable(h_wpd_means, 'save.h_wpd_means')
    save_variable(h_spd_means, 'save.h_spd_means')
    save_variable(h_wpd_availabilitys, 'save.h_wpd_availabilitys')
    save_variable(h_spd_availabilitys, 'save.h_spd_availabilitys')
    save_variable(h_wpd_ge_spd_les, 'save.h_wpd_ge_spd_les')
    save_variable(h_wpd_le_spd_ges, 'save.h_wpd_le_spd_ges')
    save_variable(h_wpd_ge_spd_ges, 'save.h_wpd_ge_spd_ges')

hourly_u_means = {}
hourly_v_means = {}

hourly_wnd10_means = {}
hourly_wnd100_means = {}

hourly_wpd_means = {}
hourly_spd_means = {}

hourly_wpd_availabilitys = {}
hourly_spd_availabilitys = {}

hourly_wpd_ge_spd_les = {}
hourly_wpd_le_spd_ges = {}
hourly_wpd_ge_spd_ges = {}

for H in selected_hours:
    hourly_u_means[H] = np.mean(h_u_means[H], axis=0)
    hourly_v_means[H] = np.mean(h_v_means[H], axis=0)

    hourly_wnd10_means[H] = np.mean(h_wnd10_means[H], axis=0)
    hourly_wnd100_means[H] = np.mean(h_wnd100_means[H], axis=0)

    hourly_wpd_means[H] = np.mean(h_wpd_means[H], axis=0)
    hourly_spd_means[H] = np.mean(h_spd_means[H], axis=0)

    hourly_wpd_availabilitys[H] = np.mean(h_wpd_availabilitys[H], axis=0)
    hourly_spd_availabilitys[H] = np.mean(h_spd_availabilitys[H], axis=0)

    hourly_wpd_ge_spd_les[H] = np.mean(h_wpd_ge_spd_les[H], axis=0)
    hourly_wpd_le_spd_ges[H] = np.mean(h_wpd_le_spd_ges[H], axis=0)
    hourly_wpd_ge_spd_ges[H] = np.mean(h_wpd_ge_spd_ges[H], axis=0)

# %%
plota_area(lats, lons, 'Time_Series_points', 'ts/')
wpd_points, spd_points = plota_ts(lats, lons, times, hourly_wpd_means,
                                  hourly_spd_means, '%H', 'Hora (UTC)',
                                  'TS_Horario', 'ts/')
np.save('TS.Horario.npy', [times, wpd_points, spd_points])

# %%

plota_mapa_hourly(lats, lons, hourly_wnd10_means,
                  'Hourly_Average_Wind_Speed_at_10_m', 'Hourly/', 'lime',
                  hourly_u_means, hourly_v_means)

plota_mapa_hourly(lats, lons, hourly_wnd100_means,
                  'Hourly_Average_Wind_Speed_at_100_m', 'Hourly/', 'lime',
                  hourly_u_means, hourly_v_means)

plota_mapa_hourly(lats, lons, hourly_wpd_means,
                  'Hourly_Average_WPD', 'Hourly/', 'lime')

plota_mapa_hourly(lats, lons, hourly_spd_means,
                  'Hourly_Average_SPD', 'Hourly/', 'lime')

plota_mapa_hourly(lats, lons, hourly_wpd_availabilitys,
                  'Hourly_WPD_Availability', 'Hourly/', 'lime')

plota_mapa_hourly(lats, lons, hourly_spd_availabilitys,
                  'Hourly_SPD_Availability', 'Hourly/', 'lime')

plota_mapa_hourly(lats, lons, hourly_wpd_ge_spd_les,
                  'Hourly_Wind_Complements_Solar', 'Hourly/', 'lime')

plota_mapa_hourly(lats, lons, hourly_wpd_le_spd_ges,
                  'Hourly_Solar_Complements_Wind', 'Hourly/', 'lime')

plota_mapa_hourly(lats, lons, hourly_wpd_ge_spd_ges,
                  'Hourly_Wind_and_Solar_Synergy', 'Hourly/', 'lime')

print('- finished! Tempo:', datetime.now() - start)
