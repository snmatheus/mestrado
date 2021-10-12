#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np

import cartopy.crs as ccrs
import cartopy.feature as cft
from cartopy.io.shapereader import Reader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.ticker as mticker

from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import ListedColormap

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from constants import pos_siglas, pos_siglas2, path_output, meses, \
    latlon_points, wpdmin, spdmin, path_shp
from utils import find_nearest

proj = ccrs.PlateCarree()

path_0_20 = path_shp + 'Batimetria/0-20.shp'
path_20_50 = path_shp + 'Batimetria/20-50.shp'
path_50_100 = path_shp + 'Batimetria/50-100.shp'
path_100_1000 = path_shp + 'Batimetria/100-1000.shp'
path_zee = path_shp + 'ZEE/EEZ_land_v1.shp'
path_brasil = path_shp + 'Brasil_degradado/Brasil_degradado.shp'

LATi, LONi = -30.38012, -30.38012
LATf, LONf = -30.625, -30.625
bounds = [(LONi, LONf, LATi, LATf)]


def ponto(ax, lat, lon, texto, color, mk='x'):
    ax.plot(lon, lat, mk, markerfacecolor=color,
            markeredgecolor=color, markersize=5)
    ax.text(lon-0.8, lat+0.5, texto, color=color)


def legenda(ax, lat, lon):
    ax.text(lon, lat,
            'Exclusive Economic Zone (EEZ)\n' +
            'Exclusion of areas with:\n' +
            '    Bathymetry < 1 km\n' +
            '    Annual average wind speed < 7 m/s\n' +
            '    Annual average SSRD < 113 W/m²\n' +
            '    Distante from the coastline > 18 km\n' +
            '    Environmental protection\n',

            zorder=6, fontsize=10, ha='left', va='center', family='fantasy',
            bbox=dict(boxstyle='square', ec='black', fc='white'))

    # ax.plot(lon+.3, lat+3.98, '_', markerfacecolor='#f818ea',
    #         markeredgecolor='k', markersize=8, zorder=7)
    # ax.plot(lon+.3, lat+3, 's', markerfacecolor='#f818ea',
    #         markeredgecolor='k', markersize=8, zorder=7)
    # ax.plot(lon+.3, lat+2, 's', markerfacecolor='red', markeredgecolor='k',
    #         markersize=8, zorder=7)
    # ax.plot(lon+.3, lat-0.71, 's', markerfacecolor='#e06969',
    #         markeredgecolor='k', markersize=8, zorder=7)
    # ax.plot(lon+.3, lat-1.59, 's', markerfacecolor='#741b47',
    #         markeredgecolor='k', markersize=8, zorder=7)
    # ax.plot(lon+.3, lat-2.48, 's', markerfacecolor='#6fa8dc',
    #         markeredgecolor='k', markersize=8, zorder=7)


def siglas(ax, r):
    if r == 'NNE' or r == 'BR':
        for k, v in pos_siglas.items():
            ax.text(v[1], v[0], k, color='#4F4F4F', zorder=14)

    if r == 'SSE' or r == 'BR':
        for k, v in pos_siglas2.items():
            ax.text(v[1], v[0], k, color='#4F4F4F', zorder=14)

    if r == 'NNE':
        ax.text(-36.2, 5, 'Atlantic Ocean', color='#4F4F4F',
                fontstyle='italic', fontsize=20, zorder=14)

    elif r == 'SSE' or r == 'BR':
        ax.text(-37.2, -36.5, 'Atlantic Ocean', color='#4F4F4F',
                fontstyle='italic', fontsize=20, zorder=14)


def continente(ax, batimetria=False, zee_color='grey'):
    bat_0_20 = Reader(path_0_20).geometries()
    bat_20_50 = Reader(path_20_50).geometries()
    bat_50_100 = Reader(path_50_100).geometries()
    bat_100_1000 = Reader(path_100_1000).geometries()
    shp_zee = Reader(path_zee).geometries()
    shp_brasil = Reader(path_brasil).geometries()
    # ax.add_geometries(shp_zee, proj, linewidth=1,
    #                   facecolor='none', edgecolor=zee_color, zorder=4.9)
    ax.add_geometries(shp_zee, proj, zorder=13, facecolor='none',
                      edgecolor='black', linewidth=0.2)

    ax.add_geometries(shp_brasil, proj, zorder=13, facecolor='white',
                      edgecolor='gray', linewidth=0.2)
    if batimetria:

        ax.add_geometries(bat_100_1000, proj, zorder=13, facecolor='none',
                          edgecolor='magenta', linewidth=0.2)
    
        ax.add_geometries(bat_50_100, proj, zorder=13, facecolor='none',
                          edgecolor='green', linewidth=0.2)
    
        ax.add_geometries(bat_20_50, proj, zorder=13, facecolor='none',
                          edgecolor='blue', linewidth=0.2)
    
        ax.add_geometries(bat_0_20, proj, zorder=13, facecolor='none',
                          edgecolor='red', linewidth=0.2)

    # states = cft.NaturalEarthFeature(category='cultural',
    #                                  name='admin_1_states_provinces',
    #                                  scale='10m',
    #                                  facecolor='silver')
    # ax.add_feature(states, edgecolor='grey', linewidth=.5, zorder=4.9)

    # ax.coastlines(resolution='10m', color='grey', linewidth=.25)


def moldura(ax):
    gl = ax.gridlines(crs=proj, draw_labels=True, linewidth=1,
                      color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(-20, -51, -5))
    gl.top_labels = False
    gl.xlines = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylocator = mticker.FixedLocator(np.arange(5, -36, -5))
    gl.right_labels = False
    gl.ylines = False


def moldura2(ax, hori, vert):
    gl = ax.gridlines(crs=proj, draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(-50, -22, 8))
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylocator = mticker.FixedLocator(np.arange(5, -36, -5))
    gl.xlines = False
    gl.ylines = False
    gl.xlabel_style = {'size': 14, 'color': 'black'}
    gl.ylabel_style = {'size': 14, 'color': 'black'}

    if vert and hori:
        gl.top_labels = False
        gl.right_labels = False
    elif vert and not hori:
        gl.top_labels = False
        gl.bottom_labels = False
        gl.right_labels = False
    elif not vert and hori:
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = False
    elif not vert and not hori:
        gl.top_labels = False
        gl.bottom_labels = False
        gl.right_labels = False
        gl.left_labels = False


def plota_area(lats, lons, fname, dir_):
    plt.rcParams['figure.figsize'] = [10, 10]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    bounds = [(lons.min(), lons.max(), lats.min(), lats.max())]
    ax.set_extent(*bounds, crs=proj)

    pt = 1

    for lat, lon in latlon_points.items():
        plt.text(lon, lat+.5, '#{:02d}'.format(pt), c='blue',
                 transform=proj, fontsize=10)

        plt.plot(lon, lat, 's', c='blue', markersize=4,
                 transform=proj)
        pt += 1

    continente(ax, 'red')
    siglas(ax, 'BR')
    moldura(ax)

    plt.savefig(path_output + dir_ + fname + '.png', bbox_inches='tight')
    plt.show()


def plota_ts(lats, lons, times, wpd, spd, format, xlabel, fname, dir_):    
    pt = 1
    for lat, lon in latlon_points.items():
        ilat = find_nearest(lats, lat)
        ilon = find_nearest(lons, lon)

        wpd_points, spd_points = [], []

        for key in wpd.keys():
            wpd_points.append(wpd[key][ilat, ilon])
            spd_points.append(spd[key][ilat, ilon])

        xs = np.arange(len(times))

        fig, ax = plt.subplots(figsize=(18, 4))

        lns1 = ax.bar(xs-0.15, wpd_points, width=.25, color='#0000FF80',
                      align='center', label='WPD')

        lns2 = ax.bar(xs+0.15, spd_points, width=.25, color='#ff000080',
                      align='center', label='SPD')

        lns3 = ax.plot(xs, [wpdmin] * len(times), '--', color='#0000FF',
                       label='WPDmin (%s)' % wpdmin)

        lns4 = ax.plot(xs, [spdmin] * len(times), '--', color='#ff0000',
                       label='SPDmin (%s)' % spdmin)

        ax.set_xlabel(xlabel, fontsize=14)
        ax.set_ylabel('W/m²', fontsize=14)

        plt.xticks(xs, [dt.strftime(format) for dt in times])
        plt.xlim([xs[0], xs[-1]])

        levels1 = np.arange(0, 1001, 100)

        ax.set_ylim([levels1[0], levels1[-1]])

        ax.legend((lns1[0], lns2[0], lns3[0], lns4[0]),
                  ('WPD', 'SPD', 'WPDmin = %s' % wpdmin,
                   'SPDmin = %s' % spdmin), loc=1)

        ax.grid()

        props = dict(facecolor='blue')
        ax.text(0.007, .97, 'Ponto #{:02d}'.format(pt), fontsize=16, color='w',
                bbox=props, transform=ax.transAxes, verticalalignment='top')

        plt.savefig('{}{}{}_p{}.png'.format(path_output, dir_, fname, pt),
                    bbox_inches='tight', dpi=300)
        plt.show()

        pt += 1

    return wpd_points, spd_points


def plota_ts2(lats, lons, times, wpd, spd, format, xlabel, fname, dir_):
    pt = 1
    for lat, lon in latlon_points.items():
        ilat = find_nearest(lats, lat)
        ilon = find_nearest(lons, lon)

        wpd_points, spd_points = [], []

        for key in wpd.keys():
            wpd_points.append(wpd[key][ilat, ilon])
            spd_points.append(spd[key][ilat, ilon])

        fig, ax1 = plt.subplots(figsize=(18, 4))

        ax2 = ax1
        # ax2 = ax1.twinx()
        ax1.fill_between(times, wpd_points, color='#0000FF55', linewidth=0.5,
                         edgecolor='#0000FF')

        ax2.fill_between(times, spd_points, color='#ff000055', linewidth=0.5,
                         edgecolor='#ff0000')

        lns1 = ax1.plot(times, wpd_points, '-', color='#0000FF', linewidth=0.5,
                        label='WPD')
        lns2 = ax2.plot(times, spd_points, '-', color='#ff0000', linewidth=0.5,
                        label='SPD')

        lns3 = ax1.plot(times, [wpdmin] * len(times), '--', color='#0000FF',
                        label='WPDmin (%s)' % wpdmin)

        lns4 = ax2.plot(times, [spdmin] * len(times), '--', color='#ff0000',
                        label='SPDmin (%s)' % spdmin)

        ax1.set_xlabel(xlabel, fontsize=14)
        ax1.set_ylabel('W/m²', fontsize=14)

        plt.xticks(times, [dt.strftime(format) for dt in times])
        plt.xlim([times[0], times[-1]])

        levels1 = np.arange(0, 1001, 100)

        ax1.set_ylim([levels1[0], levels1[-1]])

        lns = lns1 + lns2 + lns3 + lns4
        labs = [ln.get_label() for ln in lns]
        ax1.legend(lns, labs, loc=0)

        ax1.grid()

        plt.title('Ponto # {:02d}'.format(pt), fontsize=16)
        plt.savefig('{}{}{}_p{}.png'.format(path_output, dir_, fname, pt),
                    bbox_inches='tight', dpi=300)
        plt.show()

        pt += 1


def plota_mapa(lats, lons, data, fname, dir_, zee_color='red', u=[], v=[]):

    plt.rcParams['figure.figsize'] = [10, 10]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    bounds = [(lons.min(), lons.max(), lats.min(), lats.max())]
    ax.set_extent(*bounds, crs=proj)

    continente(ax, True, 'black')
    siglas(ax, 'BR')
    moldura(ax)

    if 'Wind_Speed' in fname:
        levels = np.arange(1, 12.1, 1)
        cmap, label = 'afmhot_r', 'm/s'

    elif 'Capacity' in fname or 'Availability' in fname or \
         'Complements' in fname or 'Synergy' in fname:
        levels = np.arange(0, 100.1, 10)
        cmap, label = 'inferno', '%'

    elif 'WPD' in fname:
        levels = np.arange(0, 1101, 100)
        cmap, label = 'afmhot_r', 'W/m²'

    elif 'SPD' in fname:
        levels = np.arange(0, 351, 50)
        cmap, label = 'afmhot_r', 'W/m²'

    elif 'Scale' in fname:
        levels = np.arange(0, 15, 1)
        cmap, label = 'hot_r', 'm/s'

    elif 'Shape' in fname:
        levels = np.arange(0, 22, 2)
        cmap, label = 'hot_r', ''

    elif 'Total Annual Energy Production' in fname:
        levels = np.arange(0, 1101, 50)
        cmap, label = 'hot_r', 'GWh'

    cf = ax.contourf(lons, lats, data, levels=levels, cmap=cmap,
                     transform=proj, zorder=3)

    if len(u) != 0:
        ax.streamplot(lons, lats, u, v, density=[0.9, 2], zorder=3)

    axins = inset_axes(ax, width='100%', height='100%', loc='lower left',
                       borderpad=0, bbox_to_anchor=(0.03, 0.35, 0.06, 0.4),
                       bbox_transform=ax.transAxes,
                       axes_kwargs={'title': label})

    fig.colorbar(cf, cax=axins)

    plt.savefig(path_output + dir_ + fname + '.png', bbox_inches='tight', dpi=300)
    plt.show()


def plota_mapa_anual(lats, lons, data, fname, dir_, zee_color, u=[], v=[]):

    plt.rcParams['figure.figsize'] = [30, 30]

    def configura(ax, data, u, v, ano, hori, vert):
        continente(ax, zee_color)

        ax.set_title(ano, size=18)

        cf = ax.contourf(lons, lats, data[ano], levels=levels, cmap=cmap,
                         transform=proj, zorder=2)

        if len(u) != 0:
            ax.set_extent([lons.min(), lons.max(), lats.min(), lats.max()])
            ax.streamplot(lons, lats, u[ano], v[ano], density=[1, 2], zorder=2)

        moldura2(ax, hori, vert)

        return cf

    if 'Wind_Speed' in fname:
        levels = np.arange(1, 12.1, 1)
        cmap, label = 'afmhot_r', 'm/s'

    elif 'Capacity' in fname or 'Availability' in fname or \
         'Complements' in fname or 'Synergy' in fname:
        levels = np.arange(0, 100.1, 10)
        cmap, label = 'inferno', '%'

    elif 'WPD' in fname:
        levels = np.arange(0, 1101, 100)
        cmap, label = 'afmhot_r', 'W/m²'

    elif 'SPD' in fname:
        levels = np.arange(0, 351, 50)
        cmap, label = 'afmhot_r', 'W/m²'

    fig = plt.figure()
    gs = fig.add_gridspec(6, 5, hspace=0.11, wspace=-0.82, left=0.0, right=0.8)
    ax1 = fig.add_subplot(gs[0, 0], projection=proj)
    ax2 = fig.add_subplot(gs[0, 1], projection=proj)
    ax3 = fig.add_subplot(gs[0, 2], projection=proj)
    ax4 = fig.add_subplot(gs[0, 3], projection=proj)
    ax5 = fig.add_subplot(gs[0, 4], projection=proj)
    ax6 = fig.add_subplot(gs[1, 0], projection=proj)
    ax7 = fig.add_subplot(gs[1, 1], projection=proj)
    ax8 = fig.add_subplot(gs[1, 2], projection=proj)
    ax9 = fig.add_subplot(gs[1, 3], projection=proj)
    ax10 = fig.add_subplot(gs[1, 4], projection=proj)
    ax11 = fig.add_subplot(gs[2, 0], projection=proj)
    ax12 = fig.add_subplot(gs[2, 1], projection=proj)
    ax13 = fig.add_subplot(gs[2, 2], projection=proj)
    ax14 = fig.add_subplot(gs[2, 3], projection=proj)
    ax15 = fig.add_subplot(gs[2, 4], projection=proj)
    ax16 = fig.add_subplot(gs[3, 0], projection=proj)
    ax17 = fig.add_subplot(gs[3, 1], projection=proj)
    ax18 = fig.add_subplot(gs[3, 2], projection=proj)
    ax19 = fig.add_subplot(gs[3, 3], projection=proj)
    ax20 = fig.add_subplot(gs[3, 4], projection=proj)
    ax21 = fig.add_subplot(gs[4, 0], projection=proj)
    ax22 = fig.add_subplot(gs[4, 1], projection=proj)
    ax23 = fig.add_subplot(gs[4, 2], projection=proj)
    ax24 = fig.add_subplot(gs[4, 3], projection=proj)
    ax25 = fig.add_subplot(gs[4, 4], projection=proj)
    ax26 = fig.add_subplot(gs[5, 0], projection=proj)
    ax27 = fig.add_subplot(gs[5, 1], projection=proj)
    ax28 = fig.add_subplot(gs[5, 2], projection=proj)
    ax29 = fig.add_subplot(gs[5, 3], projection=proj)
    ax30 = fig.add_subplot(gs[5, 4], projection=proj)

    configura(ax1, data, u, v, '1990', False, True)
    configura(ax2, data, u, v, '1991', False, False)
    configura(ax3, data, u, v, '1992', False, False)
    configura(ax4, data, u, v, '1993', False, False)
    configura(ax5, data, u, v, '1994', False, False)
    configura(ax6, data, u, v, '1995', False, True)
    configura(ax7, data, u, v, '1996', False, False)
    configura(ax8, data, u, v, '1997', False, False)
    configura(ax9, data, u, v, '1998', False, False)
    configura(ax10, data, u, v, '1999', False, False)
    configura(ax11, data, u, v, '2000', False, True)
    configura(ax12, data, u, v, '2001', False, False)
    configura(ax13, data, u, v, '2002', False, False)
    configura(ax14, data, u, v, '2003', False, False)
    configura(ax15, data, u, v, '2004', False, False)
    configura(ax16, data, u, v, '2005', False, True)
    configura(ax17, data, u, v, '2006', False, False)
    configura(ax18, data, u, v, '2007', False, False)
    configura(ax19, data, u, v, '2008', False, False)
    configura(ax20, data, u, v, '2009', False, False)
    configura(ax21, data, u, v, '2010', False, True)
    configura(ax22, data, u, v, '2011', False, False)
    configura(ax23, data, u, v, '2012', False, False)
    configura(ax24, data, u, v, '2013', False, False)
    configura(ax25, data, u, v, '2014', False, False)
    configura(ax26, data, u, v, '2015', True, True)
    configura(ax27, data, u, v, '2016', True, False)
    configura(ax28, data, u, v, '2017', True, False)
    configura(ax29, data, u, v, '2018', True, False)
    cf = configura(ax30, data, u, v, '2019', True, False)

    axins = inset_axes(ax30, width='100%', height='100%',
                       loc='lower left', borderpad=0,
                       bbox_to_anchor=(-4.2, -.3, 5.2, 0.11),
                       bbox_transform=ax30.transAxes,
                       axes_kwargs={'title': label})
    fig.colorbar(cf, cax=axins, orientation='horizontal')

    plt.savefig(path_output + dir_ + fname + '.png', bbox_inches='tight')
    plt.show


def plota_mapa_monthly(lats, lons, data, fname, dir_, zee_color='red',
                       u=[], v=[]):

    plt.rcParams['figure.figsize'] = [20, 20]

    def configura(ax, data, month, hori, vert):

        continente(ax, zee_color)

        ax.set_title(meses[month], size=18)

        cf = ax.contourf(lons, lats, data[month], levels=levels, cmap=cmap,
                         transform=proj, zorder=2)

        if len(u) != 0:
            ax.set_extent([lons.min(), lons.max(), lats.min(), lats.max()])
            ax.streamplot(lons, lats, u[month], v[month], density=[0.5, 1])

        moldura2(ax, hori, vert)

        return cf

    if 'Wind_Speed' in fname:
        levels = np.arange(1, 12.1, 1)
        cmap, label = 'afmhot_r', 'm/s'

    elif 'Capacity' in fname or 'Availability' in fname or \
         'Complements' in fname or 'Synergy' in fname:
        levels = np.arange(0, 100.1, 10)
        cmap, label = 'inferno', '%'

    elif 'WPD' in fname:
        levels = np.arange(0, 1101, 100)
        cmap, label = 'afmhot_r', 'W/m²'

    elif 'SPD' in fname:
        levels = np.arange(0, 351, 50)
        cmap, label = 'afmhot_r', 'W/m²'

    fig = plt.figure()
    gs = fig.add_gridspec(4, 3, hspace=0.11, wspace=-0.75, left=0.0, right=0.8)
    ax1 = fig.add_subplot(gs[0, 0], projection=proj)
    ax2 = fig.add_subplot(gs[0, 1], projection=proj)
    ax3 = fig.add_subplot(gs[0, 2], projection=proj)
    ax4 = fig.add_subplot(gs[1, 0], projection=proj)
    ax5 = fig.add_subplot(gs[1, 1], projection=proj)
    ax6 = fig.add_subplot(gs[1, 2], projection=proj)
    ax7 = fig.add_subplot(gs[2, 0], projection=proj)
    ax8 = fig.add_subplot(gs[2, 1], projection=proj)
    ax9 = fig.add_subplot(gs[2, 2], projection=proj)
    ax10 = fig.add_subplot(gs[3, 0], projection=proj)
    ax11 = fig.add_subplot(gs[3, 1], projection=proj)
    ax12 = fig.add_subplot(gs[3, 2], projection=proj)

    configura(ax1, data, 1, False, True)
    configura(ax2, data, 2, False, False)
    configura(ax3, data, 3, False, False)
    configura(ax4, data, 4, False, True)
    configura(ax5, data, 5, False, False)
    configura(ax6, data, 6, False, False)
    configura(ax7, data, 7, False, True)
    configura(ax8, data, 8, False, False)
    configura(ax9, data, 9, False, False)
    configura(ax10, data, 10, True, True)
    configura(ax11, data, 11, True, False)
    cf = configura(ax12, data, 12, True, False)

    axins = inset_axes(ax12, width='100%', height='100%',
                       loc='lower left', borderpad=0,
                       bbox_to_anchor=(-2.3, -0.3, 3.2, 0.1),
                       bbox_transform=ax12.transAxes,
                       axes_kwargs={'title': label})

    fig.colorbar(cf, cax=axins, orientation='horizontal')

    plt.savefig(path_output + dir_ + fname + '.png', bbox_inches='tight')
    plt.show()


def plota_mapa_hourly(lats, lons, data, fname, dir_, zee_color, u=[], v=[]):

    plt.rcParams['figure.figsize'] = [20, 20]

    def configura(ax, data, hora, hori, vert):

        continente(ax, zee_color)

        ax.set_extent([lons.min(), lons.max(), lats.min(), lats.max()])
        ax.set_title('{} UTC'.format(hora), size=18)
        cf = ax.contourf(lons, lats, data[hora], levels=levels, cmap=cmap,
                         transform=proj, zorder=2)

        if len(u) != 0:
            ax.set_extent([lons.min(), lons.max(), lats.min(), lats.max()])
            ax.streamplot(lons, lats, u[hora], v[hora], density=[0.5, 1])

        moldura2(ax, hori, vert)

        return cf

    if 'Wind_Speed' in fname:
        levels = np.arange(1, 12.1, 1)
        cmap, label = 'afmhot_r', 'm/s'

    elif 'Capacity' in fname or 'Availability' in fname or \
         'Complements' in fname or 'Synergy' in fname:
        levels = np.arange(0, 100.1, 10)
        cmap, label = 'inferno', '%'

    elif 'WPD' in fname or 'SPD' in fname:
        levels = np.arange(0, 1101, 100)
        cmap, label = 'afmhot_r', 'W/m²'

    fig = plt.figure()
    # gs = fig.add_gridspec(4, 6, hspace=0.15, wspace=0.05)
    gs = fig.add_gridspec(4, 6, hspace=0.11, wspace=-0.15, left=0.0, right=0.8)
    ax1 = fig.add_subplot(gs[0, 0], projection=proj)
    ax2 = fig.add_subplot(gs[0, 1], projection=proj)
    ax3 = fig.add_subplot(gs[0, 2], projection=proj)
    ax4 = fig.add_subplot(gs[0, 3], projection=proj)
    ax5 = fig.add_subplot(gs[0, 4], projection=proj)
    ax6 = fig.add_subplot(gs[0, 5], projection=proj)
    ax7 = fig.add_subplot(gs[1, 0], projection=proj)
    ax8 = fig.add_subplot(gs[1, 1], projection=proj)
    ax9 = fig.add_subplot(gs[1, 2], projection=proj)
    ax10 = fig.add_subplot(gs[1, 3], projection=proj)
    ax11 = fig.add_subplot(gs[1, 4], projection=proj)
    ax12 = fig.add_subplot(gs[1, 5], projection=proj)
    ax13 = fig.add_subplot(gs[2, 0], projection=proj)
    ax14 = fig.add_subplot(gs[2, 1], projection=proj)
    ax15 = fig.add_subplot(gs[2, 2], projection=proj)
    ax16 = fig.add_subplot(gs[2, 3], projection=proj)
    ax17 = fig.add_subplot(gs[2, 4], projection=proj)
    ax18 = fig.add_subplot(gs[2, 5], projection=proj)
    ax19 = fig.add_subplot(gs[3, 0], projection=proj)
    ax20 = fig.add_subplot(gs[3, 1], projection=proj)
    ax21 = fig.add_subplot(gs[3, 2], projection=proj)
    ax22 = fig.add_subplot(gs[3, 3], projection=proj)
    ax23 = fig.add_subplot(gs[3, 4], projection=proj)
    ax24 = fig.add_subplot(gs[3, 5], projection=proj)

    configura(ax1, data, 0, False, True)
    configura(ax2, data, 1, False, False)
    configura(ax3, data, 2, False, False)
    configura(ax4, data, 3, False, False)
    configura(ax5, data, 4, False, False)
    configura(ax6, data, 5, False, False)
    configura(ax7, data, 6, False, True)
    configura(ax8, data, 7, False, False)
    configura(ax9, data, 8, False, False)
    configura(ax10, data, 9, False, False)
    configura(ax11, data, 10, False, False)
    configura(ax12, data, 11, False, False)
    configura(ax13, data, 12, False, True)
    configura(ax14, data, 13, False, False)
    configura(ax15, data, 14, False, False)
    configura(ax16, data, 15, False, False)
    configura(ax17, data, 16, False, False)
    configura(ax18, data, 17, False, False)
    configura(ax19, data, 18, True, True)
    configura(ax20, data, 19, True, False)
    configura(ax21, data, 20, True, False)
    configura(ax22, data, 21, True, False)
    configura(ax23, data, 22, True, False)
    cf = configura(ax24, data, 23, True, False)

    axins = inset_axes(ax24, width='100%', height='100%',
                       loc='lower left', borderpad=0,
                       bbox_to_anchor=(-5.6, -0.5, 6.5, 0.2),
                       bbox_transform=ax24.transAxes,
                       axes_kwargs={'title': label})
    fig.colorbar(cf, cax=axins, orientation='horizontal')

    plt.savefig(path_output + dir_ + fname + '.png', bbox_inches='tight')


def plota_mapa_tecnico_regiao(pontos_dentro, bar_dentro, bounds=None):
    plt.rcParams['figure.figsize'] = [15, 7]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    if bounds:
        ax.set_extent(*bounds, crs=proj)

    ax.text(-36.0, -34.4, 'Atlantic Ocean', style='italic',
            weight='bold', size=12, c='white', zorder=4)

    for lon, lat in pontos_dentro[:]:
        ax.plot(lon, lat, '.', color='red', markersize=1,
                transform=proj, zorder=5)

    gl = ax.gridlines(crs=proj, draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(-20, -51, -5))
    gl.top_labels = False
    gl.xlines = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylocator = mticker.FixedLocator(np.arange(5, -36, -5))
    gl.right_labels = False
    gl.ylines = False

    plt.savefig('../SAIDAS/Mapa_Tecnico.png', bbox_inches='tight')


def plota_power_curve(dados, f, k, c):
    plt.rcParams['figure.figsize'] = [8, 5]

    plt.title('Weibull', size=18)
    plt.xlabel('wind speed (m/s)', size=14)
    plt.ylabel('Probability', size=14)

    xvalues = np.arange(0, 28, .1)
    xlim = [0, 27]
    plt.xlim(xlim)
    plt.xticks(np.arange(0, 28))
    ylim = [0, np.max(f)]
    plt.ylim(ylim)
    yvalues = np.arange(0, np.max(f)+0.1, .05)
    yticks = [str(round(pd*100))+'%' for pd in yvalues]
    plt.yticks(yvalues, yticks)

    plt.plot(np.arange(0, 28, .1), f, linewidth=1.5, color='red',
             label='k = '+str(round(k, 1))+', c = '+str(round(c, 1)))
    plt.legend(fontsize=16)
    plt.savefig('../SAIDAS/WEIBULL.png', bbox_inches='tight')
    plt.show()

    plt.rcParams['figure.figsize'] = [8, 4]
    plt.title('Power Curve', size=18)
    plt.ylabel('Turbine Power (MW)', size=14)
    plt.xlabel('wind speed (m/s)', size=14)

    plt.ylim(dados['Power (MW)'].min(), 16)
    plt.xlim([0, 27])
    plt.yticks(np.arange(1, 17, 2))
    plt.xticks(np.arange(0, 28, 1))

    plt.plot(dados['Wind (m/s)'], dados['Power (MW)'],
             label='IEC Class IB', color='blue')
    plt.legend(loc='best')
    plt.savefig('../SAIDAS/POWER_CURVE.png', bbox_inches='tight')
    plt.show()


def plota_grid():
    plt.rcParams['figure.figsize'] = [10, 10]

    LATi, LONi = -30, -30
    LATf, LONf = -30.75, -30.75
    bounds = [(LONi, LONf, LATi, LATf)]
    area_turbina = (7*240)/(111.12*1000)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    ax.set_extent(*bounds, crs=proj)

    for ln in np.arange(LONi-0.25, LONf, -0.25):
        for lt in np.arange(LATi-0.25, LATf, -0.25):
            pt = plt.plot(ln, lt, 'o', color='red', markeredgecolor='black',
                          markerfacecolor='red', markersize=10)

    for ln in np.arange(-30-0.125-area_turbina, -30.5+0.125, -area_turbina):
        for lt in np.arange(-30-0.125-area_turbina, -30.5+0.125, -area_turbina):
            t1 = plt.plot(ln, lt, marker='1',
                          markeredgecolor='blue', markersize=6)

    for ln in np.arange(-30.25-0.125-area_turbina, -30.75+0.125, -area_turbina):
        for lt in np.arange(-30.25-0.125-area_turbina, -30.75+0.125, -area_turbina):
            t2 = plt.plot(ln, lt, marker='1',
                          markeredgecolor='purple', markersize=6)

    pt = mlines.Line2D([], [], color='red', marker='o', markersize=10,
                       markeredgecolor='black', label='Grid points')
    t1 = mlines.Line2D([], [], color='blue', marker='1',
                       markersize=15, label='Example Turbines 1')
    t2 = mlines.Line2D([], [], color='purple', marker='1',
                       markersize=15, label='Example Turbines 2')

    plt.legend(handles=[pt, t1, t2])

    gl = ax.gridlines(crs=proj, draw_labels=True,
                      linewidth=1, color='gray', alpha=1, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(LONi, LONf, -0.25))
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylocator = mticker.FixedLocator(np.arange(LATi, LATf, -0.25))

    plt.savefig('../SAIDAS/grade.png', bbox_inches='tight')


def plota_gridzoom():
    plt.rcParams['figure.figsize'] = [10, 10]

    area_turbina = (7*240)/(111.12*1000)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    ax.set_extent(*bounds, crs=proj)

    # ax.add_geometries(SHP_DIF, proj, facecolor='#bee8ff', edgecolor='none', linewidth=0.5)

    lts, lns = [], []

    for ln in np.arange(-30.25-0.125-area_turbina, -30.75+0.125, -area_turbina):
        for lt in np.arange(-30.25-0.125-area_turbina, -30.75+0.125, -area_turbina):
            t2 = plt.plot(ln, lt, marker='1',
                          markeredgecolor='purple', markersize=6)
            # ax.add_patch(mpatches.Circle(xy=[ln, lt],
            #             linewidth=0.5, radius=area_turbina*0.5, color='blue',
            #             fill=False, transform=proj, zorder=3))
            lns.append(ln)
            lts.append(lt)
    df = pd.DataFrame({'lon': lns, 'lat': lts})
    df.to_csv('../SAIDAS/grade_zoom.csv')
    # pt = mlines.Line2D([], [], color='red', marker='o', markersize=10, markeredgecolor='black', label='Grid points')
    # t1 = mlines.Line2D([], [], color='blue', marker='1', markersize=15, label='Example Turbines 1')
    # t2 = mlines.Line2D([], [], color='purple', marker='1', markersize=15, label='Example Turbines 2')
    # plt.legend(handles=[pt, t1, t2])
    plt.savefig('../SAIDAS/grade_zoom.png', bbox_inches='tight')


def plota_bat(lats, lons, data, fname, dir_, zee_color='red'):

    plt.rcParams['figure.figsize'] = [10, 10]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    cmap = ListedColormap(['#bee8ff'])
    levels2 = [-1000, -998]
    norm = BoundaryNorm(levels2, 1)

    cf = ax.contourf(lons, lats, -999*(1**data), levels=levels2, cmap=cmap,
                     norm=norm, transform=proj, zorder=13)

    bounds = [(lons.min(), lons.max(), lats.min(), lats.max())]
    ax.set_extent(*bounds, crs=proj)

    continente(ax, True, 'black')
    siglas(ax, 'BR')
    moldura(ax)
    ax.set(facecolor = "orange")

    if 'Wind_Speed' in fname:
        vmin, vmax = 0, 12.1
        cmap, label = 'afmhot_r', 'm/s'

    elif 'Availability' in fname or 'Complements' in fname \
         or 'Synergy' in fname or 'CF' in fname:
        vmin, vmax = 0, 100
        cmap, label = 'inferno', '%'

    elif 'WPD' in fname:
        vmin, vmax = 0, 1101
        cmap, label = 'afmhot_r', 'W/m²'

    elif 'SPD' in fname:
        vmin, vmax = 0, 351
        cmap, label = 'afmhot_r', 'W/m²'

    elif 'AEP' in fname:
        vmin, vmax = 0, 1101
        cmap, label = 'hot_r', 'MWh'

    vmin, vmax = data[data>0].min(), data[data>0].max()
    levels = MaxNLocator(nbins=10).tick_values(vmin, vmax)
    cmap = plt.get_cmap(cmap)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    data[data<=0] = np.nan
    cf = ax.pcolormesh(lons, lats, data, cmap=cmap, norm=norm, transform=proj,
                       zorder=13)

    axins = inset_axes(ax, width='100%', height='100%', loc='lower left',
                       borderpad=0, bbox_to_anchor=(0.03, 0.35, 0.06, 0.4),
                       bbox_transform=ax.transAxes,
                       axes_kwargs={'title': label})

    fig.colorbar(cf, cax=axins)

    plt.savefig(path_output + dir_ + fname + '.png', bbox_inches='tight', dpi=300)
    plt.show()
