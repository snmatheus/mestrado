#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import xarray as xr
from constants import path_wnd
from utils import save_variable, load_variable, climatology, extrapolate

lats = load_variable('save.lats')
lons = load_variable('save.lons')

try:
    wnd10_means = load_variable('wnd10_means.pickle')
    wnd100_means = load_variable('wnd100_means.pickle')
    wnd10_std = load_variable('wnd10_std.pickle')
    wnd100_std = load_variable('wnd100_std.pickle')

    clim_wnd10_means = load_variable('clim_wnd10_means.pickle')
    clim_wnd100_means = load_variable('clim_wnd100_means.pickle')
    clim_wnd10_std = load_variable('clim_wnd10_std.pickle')
    clim_wnd100_std = load_variable('clim_wnd100_std.pickle')

except:
    wnd10_means = {}
    wnd100_means = {}
    wnd10_std = {}
    wnd100_std = {}

    w_file_list = [file_name for file_name in sorted(os.listdir(path_wnd))
                   if file_name.endswith('.nc')]

    for wnd_file in w_file_list:
        year = wnd_file[7:11]
        print(year)
        wnd_ds = xr.open_dataset(path_wnd + wnd_file)

        u_data = wnd_ds['u10'][:, :, :].data
        v_data = wnd_ds['v10'][:, :, :].data

        wnd10_data = (u_data**2 + v_data**2)**0.5
        wnd100_data = extrapolate(wnd10_data)

        wnd10_means[year] = np.nanmean(wnd10_data, axis=0)
        wnd100_means[year] = np.nanmean(wnd100_data, axis=0)

        wnd10_std[year] = np.nanstd(wnd10_data, axis=0, ddof=1)
        wnd100_std[year] = np.nanstd(wnd100_data, axis=0, ddof=1)

        save_variable(wnd10_means, 'wnd10_means.pickle')
        save_variable(wnd100_means, 'wnd100_means.pickle')
        save_variable(wnd10_std, 'wnd10_std.pickle')
        save_variable(wnd100_std, 'wnd100_std.pickle')

        wnd_ds.close()


clim_wnd10_means = climatology(list(wnd10_means.values()))
clim_wnd100_means = climatology(list(wnd100_means.values()))
clim_wnd10_std = climatology(list(wnd10_std.values()))
clim_wnd100_std = climatology(list(wnd100_std.values()))

save_variable(clim_wnd10_means, 'clim_wnd10_means.pickle')
save_variable(clim_wnd100_means, 'clim_wnd100_means.pickle')
save_variable(clim_wnd10_std, 'clim_wnd10_std.pickle')
save_variable(clim_wnd100_std, 'clim_wnd100_std.pickle')
