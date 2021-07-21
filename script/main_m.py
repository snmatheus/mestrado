from output import plota_mapa_monthly, plota_area, plota_ts
import os
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime
from constants import path_wnd, path_rad, wpdmin, spdmin, path_pos, idxes, \
                      selected_months
from utils import load_variable, save_variable, wpd, extrapolate


start = datetime.now()

lats = load_variable('save.lats')
lons = load_variable('save.lons')

if os.path.isfile(path_pos + 'save.m_wnd10_means'):
    m_u_means = load_variable('save.m_u_means')
    m_v_means = load_variable('save.m_v_means')

    m_wnd10_means = load_variable('save.m_wnd10_means')
    m_wnd100_means = load_variable('save.m_wnd100_means')

    m_wpd_means = load_variable('save.m_wpd_means')
    m_spd_means = load_variable('save.m_spd_means')

    m_wpd_availabilitys = load_variable('save.m_wpd_availabilitys')
    m_spd_availabilitys = load_variable('save.m_spd_availabilitys')

    m_wpd_ge_spd_les = load_variable('save.m_wpd_ge_spd_les')
    m_wpd_le_spd_ges = load_variable('save.m_wpd_le_spd_ges')
    m_wpd_ge_spd_ges = load_variable('save.m_wpd_ge_spd_ges')
    times = pd.date_range(start='1990-01', end='1991-01', freq='M')

else:
    m_u_means = {}
    m_v_means = {}
    m_wnd10_means = {}
    m_wnd100_means = {}
    m_wpd_means = {}
    m_spd_means = {}
    m_wpd_availabilitys = {}
    m_spd_availabilitys = {}
    m_wpd_ge_spd_les = {}
    m_wpd_le_spd_ges = {}
    m_wpd_ge_spd_ges = {}

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

        for m in selected_months:

            i1, i2 = idxes[m-1]

            if len(times) == 8784:
                if m >= 2:
                    i2 += 24
                if m > 2:
                    i1 += 24

            time = times[i1: i2]
            r_data = rad_ds['ssrd'][i1: i2].data
            u_data = wnd_ds['u10'][i1: i2].data
            v_data = wnd_ds['v10'][i1: i2].data

            tnh = len(time)

            print(time[0].strftime('%Y-%m-%dT%H:00'),
                  time[-1].strftime('%Y-%m-%dT%H:00'))

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

            m_u = np.mean(u_data, axis=0)
            m_v = np.mean(v_data, axis=0)

            m_wnd10 = np.mean(wnd10_data, axis=0)
            m_spd = np.mean(spd_data, axis=0)

            m_wpd_available = np.sum(nh_wpd_available, axis=0)/tnh * 100
            m_spd_available = np.sum(nh_spd_available, axis=0)/tnh * 100

            m_wpd_ge_spd_le = np.sum(nh_wpd_ge_spd_le, axis=0)/tnh * 100
            m_wpd_le_spd_ge = np.sum(nh_wpd_le_spd_ge, axis=0)/tnh * 100
            m_wpd_ge_spd_ge = np.sum(nh_wpd_ge_spd_ge, axis=0)/tnh * 100

            if m not in m_wnd10_means:
                m_u_means[m] = m_u.reshape(shape)
                m_v_means[m] = m_v.reshape(shape)

                m_wnd10_means[m] = m_wnd10.reshape(shape)
                m_spd_means[m] = m_spd.reshape(shape)

                m_wpd_availabilitys[m] = m_wpd_available.reshape(shape)
                m_spd_availabilitys[m] = m_spd_available.reshape(shape)

                m_wpd_ge_spd_les[m] = m_wpd_ge_spd_le.reshape(shape)
                m_wpd_le_spd_ges[m] = m_wpd_le_spd_ge.reshape(shape)
                m_wpd_ge_spd_ges[m] = m_wpd_ge_spd_ge.reshape(shape)

            else:
                m_u_means[m] = np.vstack((m_u_means[m], m_u.reshape(shape)))
                m_v_means[m] = np.vstack((m_v_means[m], m_v.reshape(shape)))

                m_wnd10_means[m] = np.vstack((m_wnd10_means[m],
                                              m_wnd10.reshape(shape)))

                m_spd_means[m] = np.vstack(
                    (m_spd_means[m], m_spd.reshape(shape)))

                m_wpd_availabilitys[m] = np.vstack((m_wpd_availabilitys[m],
                                                    m_wpd_available
                                                    .reshape(shape)))

                m_spd_availabilitys[m] = np.vstack((m_spd_availabilitys[m],
                                                    m_spd_available
                                                    .reshape(shape)))

                m_wpd_ge_spd_les[m] = np.vstack((m_wpd_ge_spd_les[m],
                                                 m_wpd_ge_spd_le
                                                .reshape(shape)))

                m_wpd_le_spd_ges[m] = np.vstack((m_wpd_le_spd_ges[m],
                                                 m_wpd_le_spd_ge
                                                .reshape(shape)))

                m_wpd_ge_spd_ges[m] = np.vstack((m_wpd_ge_spd_ges[m],
                                                 m_wpd_ge_spd_ge
                                                .reshape(shape)))
    # % %
    for m in selected_months:
        m_wnd100_means[m] = extrapolate(m_wnd10_means[m])
        m_wpd_means[m] = wpd(m_wnd100_means[m])

    # % %
    save_variable(m_u_means, 'save.m_u_means')
    save_variable(m_v_means, 'save.m_v_means')
    save_variable(m_wnd10_means, 'save.m_wnd10_means')
    save_variable(m_wnd100_means, 'save.m_wnd100_means')
    save_variable(m_wpd_means, 'save.m_wpd_means')
    save_variable(m_spd_means, 'save.m_spd_means')
    save_variable(m_wpd_availabilitys, 'save.m_wpd_availabilitys')
    save_variable(m_spd_availabilitys, 'save.m_spd_availabilitys')
    save_variable(m_wpd_ge_spd_les, 'save.m_wpd_ge_spd_les')
    save_variable(m_wpd_le_spd_ges, 'save.m_wpd_le_spd_ges')
    save_variable(m_wpd_ge_spd_ges, 'save.m_wpd_ge_spd_ges')

# % %
monthly_u_means = {}
monthly_v_means = {}

monthly_wnd10_means = {}
monthly_wnd100_means = {}

monthly_wpd_means = {}
monthly_spd_means = {}

monthly_wpd_availabilitys = {}
monthly_spd_availabilitys = {}

monthly_wpd_ge_spd_les = {}
monthly_wpd_le_spd_ges = {}
monthly_wpd_ge_spd_ges = {}

for m in selected_months:
    monthly_u_means[m] = np.mean(m_u_means[m], axis=0)
    monthly_v_means[m] = np.mean(m_v_means[m], axis=0)

    monthly_wnd10_means[m] = np.mean(m_wnd10_means[m], axis=0)
    monthly_wnd100_means[m] = np.mean(m_wnd100_means[m], axis=0)

    monthly_wpd_means[m] = np.mean(m_wpd_means[m], axis=0)
    monthly_spd_means[m] = np.mean(m_spd_means[m], axis=0)

    monthly_wpd_availabilitys[m] = np.mean(m_wpd_availabilitys[m], axis=0)
    monthly_spd_availabilitys[m] = np.mean(m_spd_availabilitys[m], axis=0)

    monthly_wpd_ge_spd_les[m] = np.mean(m_wpd_ge_spd_les[m], axis=0)
    monthly_wpd_le_spd_ges[m] = np.mean(m_wpd_le_spd_ges[m], axis=0)
    monthly_wpd_ge_spd_ges[m] = np.mean(m_wpd_ge_spd_ges[m], axis=0)

# %%
plota_area(lats, lons, 'Time_Series_points', 'ts/')
wpd_points, spd_points = plota_ts(lats, lons, times, monthly_wpd_means,
                                  monthly_spd_means, '%b', 'Mês', 'TS_Mensal',
                                  'ts/')
np.save('TS.Mensal.npy', [times, wpd_points, spd_points])

# %%

plota_mapa_monthly(lats, lons, monthly_wnd10_means,
                   'Monthly_Average_Wind_Speed_at_10_m', 'Monthly/', 'lime',
                   monthly_u_means, monthly_v_means)

plota_mapa_monthly(lats, lons, monthly_wnd100_means,
                   'Monthly_Average_Wind_Speed_at_100_m', 'Monthly/', 'lime',
                   monthly_u_means, monthly_v_means)

plota_mapa_monthly(lats, lons, monthly_wpd_means,
                   'Monthly_Average_WPD', 'Monthly/', 'lime')

plota_mapa_monthly(lats, lons, monthly_spd_means,
                   'Monthly_Average_SPD', 'Monthly/', 'lime')

plota_mapa_monthly(lats, lons, monthly_wpd_availabilitys,
                   'Monthly_WPD_Availability', 'Monthly/', 'lime')

plota_mapa_monthly(lats, lons, monthly_spd_availabilitys,
                   'Monthly_SPD_Availability', 'Monthly/', 'lime')

plota_mapa_monthly(lats, lons, monthly_wpd_ge_spd_les,
                   'Monthly_Wind_Complements_Solar', 'Monthly/', 'lime')

plota_mapa_monthly(lats, lons, monthly_wpd_le_spd_ges,
                   'Monthly_Solar_Complements_Wind', 'Monthly/', 'lime')

plota_mapa_monthly(lats, lons, monthly_wpd_ge_spd_ges,
                   'Monthly_Wind_and_Solar_Synergy', 'Monthly/', 'lime')


print('- finished! Tempo:', datetime.now() - start)
