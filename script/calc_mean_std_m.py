#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import xarray as xr
import pandas as pd
from constants import path_wnd, selected_months, idxes
from utils import save_variable, load_variable, climatology, extrapolate

lats = load_variable('save.lats')
lons = load_variable('save.lons')

try:
    m_wnd100_mean = load_variable('m_wnd100_means.pickle')
    m_wnd100_std = load_variable('m_wnd100_std.pickle')

except:
    m_wnd100_mean = {}
    m_wnd100_std = {}
    w_file_list = [file_name for file_name in sorted(os.listdir(path_wnd)) if file_name.endswith('.nc')]

    for m in selected_months:
        i1, i2 = idxes[m-1]
        m_total = None

        for wnd_file in w_file_list:
            wnd_ds = xr.open_dataset(path_wnd + wnd_file)
            times64 = wnd_ds['time'].data
            times = np.array([pd.Timestamp(time).to_pydatetime() for time in times64])

            if len(times) == 8784: # ano bissexto
                if m >= 2:
                    i2 += 24
                if m > 2:
                    i1 += 24

            time = times[i1: i2]
            print(time[0].strftime('%Y-%m-%dT%H:00'), time[-1].strftime('%Y-%m-%dT%H:00'))

            u_data = wnd_ds['u10'][i1: i2].data
            v_data = wnd_ds['v10'][i1: i2].data

            m_wnd10 = (u_data**2 + v_data**2)**0.5
            m_wnd100 = extrapolate(m_wnd10)

            if m_total is None:
                m_total = m_wnd100
            else:
                m_total = np.vstack((m_total, m_wnd100))

            wnd_ds.close()

        m_wnd100_mean[m] = np.nanmean(m_total, axis=0)
        m_wnd100_std[m] = np.nanstd(m_total, axis=0, ddof=1)

save_variable(m_wnd100_mean, 'm_wnd100_means.pickle')
save_variable(m_wnd100_std, 'm_wnd100_std.pickle')
