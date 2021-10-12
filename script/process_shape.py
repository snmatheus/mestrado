#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pickle
import numpy as np
import pandas as pd
import shapefile as shp
import cartopy.crs as ccrs
import cartopy.feature as cft
from utils import find_nearest
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from shapely.geometry import Point
import matplotlib.ticker as mticker
from cartopy.io.shapereader import Reader
from shapely.geometry.polygon import Polygon
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from constants import pos_siglas, pos_siglas2, path_output, meses, path_shp, \
                      latlon_points, wpdmin, spdmin, path_points

proj = ccrs.PlateCarree()
RES = 0.25
LATi, LATf = -36, 8
LONi, LONf = -53, -23

bounds = [LONi, LONf, LATi, LATf]

path_0_20 = path_shp + 'Batimetria/0-20.shp'
path_20_50 = path_shp + 'Batimetria/20-50.shp'
path_50_100 = path_shp + 'Batimetria/50-100.shp'
path_100_1000 = path_shp + 'Batimetria/100-1000.shp'
path_zee = path_shp + 'ZEE/EEZ_land_v1.shp'
path_brasil = path_shp + 'Brasil_degradado/Brasil_degradado.shp'

bat_0_20 = Reader(path_0_20).geometries()
bat_20_50 = Reader(path_20_50).geometries()
bat_50_100 = Reader(path_50_100).geometries()
bat_100_1000 = Reader(path_100_1000).geometries()
shp_zee = Reader(path_zee).geometries()
shp_brasil = Reader(path_brasil).geometries()

# SHP_18KM = Reader('shps/18km/distancia-18km.shp').geometries()
# SHP_DIF = Reader('/discolocal/marcolino/mestrado/SCRIPT/shps/Oceano_corte_ZEE/Oceano_ZEE.shp').geometries()
# SHP_1000 = Reader('shps/Batimetria_1000m/Batimetria_1000m.shp').geometries()


def moldura(ax):
    gl = ax.gridlines(crs=proj, draw_labels=True, linewidth=1, color='gray',
                      alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(-20, -51, -5))
    gl.top_labels = False
    gl.xlines = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylocator = mticker.FixedLocator(np.arange(5, -36, -5))
    gl.right_labels = False
    gl.ylines = False


def plota_mapa_tecnico_regiao(points, filename):
    plt.rcParams['figure.figsize'] = [15, 15]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    ax.set_extent(bounds)
    ax.set_facecolor('lightblue')

    for _, row in points.iterrows():
        ax.plot(row.lon, row.lat, '.', color='navy', markersize=1, zorder=6,
                transform=proj)

    ax.add_geometries(shp_zee, proj, zorder=1, facecolor='snow',
                      edgecolor='black', linewidth=0.2)

    ax.add_geometries(shp_brasil, proj, zorder=1.5, facecolor='silver',
                      edgecolor='black', linewidth=0.2)

    ax.add_geometries(bat_100_1000, proj, zorder=2, facecolor='ivory',
                      edgecolor='black', linewidth=0.2)

    ax.add_geometries(bat_50_100, proj, zorder=3, facecolor='azure',
                      edgecolor='black', linewidth=0.2)

    ax.add_geometries(bat_20_50, proj, zorder=4, facecolor='azure',
                      edgecolor='black', linewidth=0.2)

    ax.add_geometries(bat_0_20, proj, zorder=5, facecolor='azure',
                      edgecolor='black', linewidth=0.2)


    # ax.add_geometries(SHP_UC, proj, facecolor='#32CD32',
    #                   edgecolor='#32CD32', zorder=3, linewidth=0.5)
    # ax.add_geometries(SHP_18KM, proj, facecolor='#4c1130',
    #                   edgecolor='black', linewidth=0.5, zorder=3)

    # cmap = ListedColormap(['#e06969','#e06969'])
    # ax.contourf(lons, lats, ws7, cmap=cmap, transform=proj, zorder=3)

    # cmap = ListedColormap(['red','red'])
    # ax.contourf(lons, lats, ROCI113, cmap=cmap, transform=proj, zorder=3)

    moldura(ax)
    plt.savefig(path_output + filename + '.png', dpi=300)
    plt.show()


def read_shp2df(path):
    sf = shp.Reader(path)
    fields = [x[0] for x in sf.fields][1:]
    records = [list(i) for i in sf.records()]
    shps = [s.points for s in sf.shapes()]
    df = pd.DataFrame(columns=fields, data=records)
    df = df.assign(coords=shps)
    return df


def zee_filter(points):
    df_zee = read_shp2df(path_zee).sum()    
    zee_polygon = Polygon(df_zee.coords)
    filtered_points = []
    for lonlat in points:
        point = Point(lonlat)
        if point.within(zee_polygon):
            filtered_points.append(lonlat)
    return filtered_points


def brasil_filter(points):
    df_brasil = read_shp2df(path_brasil)
    filtered_points = []
    for lonlat in points:
        point = Point(lonlat)
        add = True
        for i, row in df_brasil.iterrows():
            state_polygon = Polygon(row.coords)
            if point.within(state_polygon):
               add = False
               break
        if add:
            filtered_points.append(lonlat)
    return filtered_points


def bat_filter(path, points, aux):
    print('Filtrando:', aux)
    df = read_shp2df(path).sum()    
    polygon = Polygon(df.coords)
    filtered_points = []
    for _, row in points.iterrows():
        point = Point((row.lon, row.lat))
        if point.within(polygon):
            filtered_points.append((row.lon, row.lat))
    return pd.DataFrame(filtered_points, columns=['lon','lat'])


def save_dict(filename, data):
    file = open("./" + filename + ".pkl", "wb")
    pickle.dump(data, file)


def read_dict(filename):
    file = open("./" + filename + ".pkl", "rb")
    data = pickle.load(file)
    return data


#%%
lats = np.arange(LATi, LATf+RES, RES)
lons = np.arange(LONi, LONf+RES, RES)

lonlats = []
for lon in lons:
    for lat in lats:
        lonlats.append((lon,lat))




