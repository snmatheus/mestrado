#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 22:40:43 2020

@author: administrador
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

plt.rcParams["figure.figsize"]=[10,10]

meses = ["Jan","Fev","Mar","Abr","Mai","Jun","Jul","Ago","Set","Out","Nov","Dez",]
mesoes = ["JAN","FEV","MAR","ABR","MAI","JUN","JUL","AGO","SET","OUT","NOV","DEZ",]
lats = np.load("../DADOS/NPY/lats.npy")
lons = np.load("../DADOS/NPY/lons.npy")

def plota_sat(mes, vel, u, v, lats, lons):
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    from cartopy.io.shapereader import Reader
    import matplotlib.ticker as mticker

    plt.rcParams["figure.figsize"]=[10,10]    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    LATi, LONi=0, -49
    LATf, LONf=-18, -30

    LATi, LONi=30, -42
    LATf, LONf=-50, -11

    bounds = [(LONi, LONf, LATi, LATf)]
    ax.set_extent(*bounds, crs=ccrs.PlateCarree())
 
    cf = ax.contourf(lons, lats, vel, levels=np.arange(0, np.nanmax(vel)+0.5,0.5), cmap="viridis", transform=ccrs.PlateCarree(), zorder=2)
    ax.set_title('Vento a 10 m (m/s) e Dir. Vento (setas)',loc='left')
    cb = fig.colorbar(cf, orientation='horizontal', extend='max', aspect=65, shrink=0.5, pad=0.05, extendrect='True')
    cb.set_label('Velocidade')
    ax.streamplot(lons, lats, u, v, density=3, linewidth=1, color='k', transform=ccrs.PlateCarree(), zorder=2)    

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(LONf, LONi, -8))
    gl.xlabels_top = False
    gl.xlines = False

    gl.yformatter = LATITUDE_FORMATTER
    gl.ylocator = mticker.FixedLocator(np.arange(LATi, LATf, -8))
    gl.ylabels_right = False
    gl.ylines = False
    gl.ylabel_style = {'rotation':'vertical'}

    nome_arq="SAT10_"+mes.replace("-","").replace(" ","").replace(":","")[:-6]
    plt.tight_layout()
    plt.savefig(nome_arq+".png", dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()


def ponto(ax, lat, lon, texto, color, mk='x'):
    ax.plot(lon,lat, mk, markerfacecolor=color, markeredgecolor=color, markersize=5)
    ax.text(lon-0.8,lat+0.5, texto, color=color)

def legenda(ax, lat, lon):
    ax.text(lon, lat,
    "Exclusive Economic Zone (EEZ)\n"+
    "Exclusion of areas with:\n"+
    "    Bathymetry < 1 km\n"+
    "    Annual average wind speed < 7 m/s\n"+
    "    Annual average SSRD < 113 W/m²\n"+
    "    Distante from the coastline > 18 km\n"+
    "    Environmental protection\n",
    zorder=6, fontsize=10, ha="left", va="center", family='fantasy',
    bbox=dict(boxstyle="square", ec="black", fc="white")
    )
    # ax.plot(lon+.3, lat+3.98, "_", markerfacecolor='#f818ea', markeredgecolor='k', markersize=8, zorder=7)
    # ax.plot(lon+.3, lat+3, "s", markerfacecolor='#f818ea', markeredgecolor='k', markersize=8, zorder=7)
    # ax.plot(lon+.3, lat+2, "s", markerfacecolor='red', markeredgecolor='k', markersize=8, zorder=7)
    # ax.plot(lon+.3, lat-0.71, "s", markerfacecolor='#e06969', markeredgecolor='k', markersize=8, zorder=7)
    # ax.plot(lon+.3, lat-1.59, "s", markerfacecolor='#741b47', markeredgecolor='k', markersize=8, zorder=7)
    # ax.plot(lon+.3, lat-2.48, "s", markerfacecolor='#6fa8dc', markeredgecolor='k', markersize=8, zorder=7)

def plota_mapa(var, nome_arq, bounds=None):
    plt.rcParams["figure.figsize"]=[10,10]    
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from matplotlib.colors import ListedColormap
    from matplotlib.colors import BoundaryNorm
    from cartopy.io.shapereader import Reader
    import matplotlib.ticker as mticker

    u, v = [], []
    norm = ""
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    if bounds:       
        ax.set_extent(*bounds, crs=ccrs.PlateCarree())

    SHP_AS=Reader("shps/AS/america_do_sul.shp").geometries()
    SHP_BRASIL=Reader("shps/Brasil/estados_2010.shp").geometries()
    SHP_1000=Reader("shps/Batimetria_1000m/Batimetria_1000m.shp").geometries()
    SHP_ZEE=Reader("shps/ZEE/ZEE.shp").geometries()

    ax.add_geometries(SHP_AS, ccrs.PlateCarree(), linewidth=0.4, facecolor="silver", edgecolor="grey", zorder=3.1)
    ax.add_geometries(SHP_BRASIL, ccrs.PlateCarree(), linewidth=0.4, facecolor="silver", edgecolor="grey", zorder=3.1)
    ax.add_geometries(SHP_ZEE, ccrs.PlateCarree(), linewidth=0.5, facecolor="none", edgecolor="black", zorder=3)
    ax.add_geometries(SHP_1000, ccrs.PlateCarree(), linewidth=0.5, facecolor="none", edgecolor="black", zorder=3)

    if "Capacity" in str(nome_arq) or "Disponibility" in str(nome_arq) or "Complements" in str(nome_arq) or "Synergy" in str(nome_arq):
        titulo = nome_arq.replace("_"," ")
        z = var
        levels = np.arange(0, 100.1, 10)
        cmap, label = "inferno_r", "%"

    elif "Wind_Speed" in str(nome_arq):
        titulo = nome_arq.replace("_"," ")
        z = var
        levels = np.arange(1,11.1,0.5)
        cmap, label = "afmhot_r", "m/s"

    elif "WPD" in str(nome_arq) or "SPD" in str(nome_arq):
        titulo = nome_arq.replace("_"," ")
        z= var
        levels = np.arange(0,601,50)
        cmap, label = "hot_r", "W/m²"

    ax.set_title(titulo, loc='left', fontsize=16)

    cf = ax.contourf(lons, lats, z, levels=levels, cmap=cmap, transform=ccrs.PlateCarree(), zorder=2)

    axins = inset_axes(ax, width="100%", height="100%", loc='lower left', borderpad=0, bbox_to_anchor=(0.03, 0.35, 0.06, 0.4),
                        bbox_transform=ax.transAxes, axes_kwargs= {'title': label})   
    fig.colorbar(cf, cax=axins)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(-20, -51, -5))
    gl.xlabels_top = False
    gl.xlines = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylocator = mticker.FixedLocator(np.arange(5, -36, -5))
    gl.ylabels_right = False
    gl.ylines = False

    plt.savefig("../SAIDAS/"+nome_arq+".png", bbox_inches='tight')
    #plt.show()              

def siglas(ax, r):
    SIGLAS = { 
                "CE": (-5.3067450, -39.843010), "AP": (1.258245, -51.609507),
                "MA": (-4.498199, -45.338160),
                "PA": (-2.679555, -47.934728),
                "PB": (-7.2663620, -36.642835), "PE": (-7.755057, -35.272919),
                "PI": (-7.5557450, -42.802614), "RN": (-5.7803050, -36.786438),
                "TO": (-6.466770, -48.178640), }

    SIGLAS2 = {  
                "ES": (-19.615438, -40.765840), "GO": (-17.138677, -50.355021),
                "MS": (-19.792475, -52.458537), "BA": (-16.992282, -39.690487),
                "MG": (-19.587865, -44.699646), "PR": (-25.108172, -49.621746),
                "RJ": (-22.170342, -42.644215), "RS": (-29.354136, -53.336608),
                "SC": (-27.207366, -50.449029), "SP": (-22.117968, -48.911867),
                "MT": (-16.599580, -53.698997),}

    if r=="NNE":
        for k, v in SIGLAS.items():
            ax.text(v[1], v[0], k, zorder=5)
            ax.text(-40.2, 5,"Atlantic Ocean", color='#000000', fontstyle='italic', fontsize=20, zorder=5)
    if r=="SSE":
        for k, v in SIGLAS2.items():
            ax.text(v[1], v[0], k, zorder=5)
            ax.text(-43, -34.5, "Atlantic Ocean", color='#000000', fontstyle='italic', fontsize=20, zorder=5)

def plota_mapa2(var, nome_arq):
    plt.rcParams["figure.figsize"]=[10,10]    
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from cartopy.io.shapereader import Reader
    import matplotlib.ticker as mticker

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    if "NNE" in nome_arq:
        siglas(ax, "NNE")
        LATi, LONi=6, -52
        LATf, LONf=-8, -34
    if "SSE" in nome_arq:
        siglas(ax, "SSE")
        LATi, LONi=-16, -54
        LATf, LONf=-35, -35.5

    bounds = [(LONi, LONf, LATi, LATf)]
    ax.set_extent(*bounds, crs=ccrs.PlateCarree())

    SHP_AS=Reader("shps/AS/america_do_sul.shp").geometries()
    SHP_BRASIL=Reader("shps/Brasil/estados_2010.shp").geometries()
    SHP_1000=Reader("shps/Batimetria_1000m/Batimetria_1000m.shp").geometries()
    SHP_ZEE=Reader("shps/ZEE/ZEE.shp").geometries()
    SHP_DIF = Reader("shps/Oceano_corte_ZEE/Oceano_ZEE.shp").geometries()
   
    ax.add_geometries(SHP_AS, ccrs.PlateCarree(), linewidth=0.4, facecolor="silver", edgecolor="black", zorder=3.1)
    ax.add_geometries(SHP_BRASIL, ccrs.PlateCarree(), linewidth=0.4, facecolor="white", edgecolor="black", zorder=3.1)
    ax.add_geometries(SHP_ZEE, ccrs.PlateCarree(), linewidth=2, facecolor="#bee8ff", edgecolor="white", zorder=2.9)
    ax.add_geometries(SHP_DIF, ccrs.PlateCarree(), facecolor="#bee8ff", edgecolor="none", zorder=3, linewidth=0.5)
    ax.add_geometries(SHP_1000, ccrs.PlateCarree(), linewidth=1, facecolor="none", edgecolor="blue", zorder=3)
    
    if "Wind_Speed" in str(nome_arq):
        titulo = nome_arq.replace("_"," ")
        z = var
        levels = np.arange(1,11.1,0.5)
        cmap, label = "afmhot_r", "m/s"

    elif "Capacity" in str(nome_arq) or "Disponibility" in str(nome_arq) or "Complements" in str(nome_arq) or "Synergy" in str(nome_arq):
        titulo = nome_arq.replace("_"," ")
        z = var
        levels = np.arange(0, 100.1, 10)
        cmap, label = "hot_r", "%"

    elif "Scale" in str(nome_arq):
        titulo = nome_arq.replace("_"," ")
        z = var
        levels = np.arange(0, 15, 1)
        cmap, label = "hot_r", "m/s"

    elif "Shape" in str(nome_arq):
        titulo = nome_arq.replace("_"," ")
        z = var
        levels = np.arange(0, 22, 2)
        cmap, label = "hot_r", ""

    elif "Total Annual Energy Production" in str(nome_arq):
        titulo = nome_arq.replace("_"," ")
        z= var
        levels = np.arange(0,601,50)
        cmap, label = "hot_r", "GWh"#"MWh"

    elif "WPD" in str(nome_arq) or "SPD" in str(nome_arq):
        titulo = nome_arq.replace("_"," ")
        z= var
        levels = np.arange(0,601,50)
        cmap, label = "afmhot_r", "W/m²"

    ax.set_title(titulo, loc='left', fontsize=16)

    # cf = ax.contourf(lons, lats, z, levels=levels, cmap=cmap, transform=ccrs.PlateCarree(), zorder=3)
    cf = ax.pcolormesh(lons, lats, z, cmap=cmap,  transform=ccrs.PlateCarree(), zorder=3)
    

    if "NNE" in nome_arq:
        axins = inset_axes(ax, width="100%", height="100%", loc='lower left',
                       borderpad=0, bbox_to_anchor=(0.03, 0.05, 0.06, 0.4),
                       bbox_transform=ax.transAxes, axes_kwargs= {'title': label})
    if "SSE" in nome_arq:
        axins = inset_axes(ax, width="100%", height="100%", loc='lower left',
                       borderpad=0, bbox_to_anchor=(0.03, 0.35, 0.06, 0.4),
                       bbox_transform=ax.transAxes, axes_kwargs= {'title': label})
    fig.colorbar(cf, cax=axins)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(-20, -51, -5))
    gl.xlabels_top = False
    gl.xlines = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylocator = mticker.FixedLocator(np.arange(5, -36, -5))
    gl.ylabels_right = False
    gl.ylines = False

    plt.savefig("../SAIDAS/"+nome_arq+".png", bbox_inches='tight')
    #plt.show()     

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.io.shapereader import Reader
import matplotlib.ticker as mticker
import matplotlib.lines as mlines

def plota_grid():
    plt.rcParams["figure.figsize"]=[10,10]
    
    LATi, LONi=-30, -30
    LATf, LONf=-30.75, -30.75
    bounds = [(LONi, LONf, LATi, LATf)]

    SHP_DIF = Reader("shps/Oceano_corte_ZEE/Oceano_ZEE.shp").geometries()
    area_turbina = (7*240)/(111.12*1000)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(*bounds, crs=ccrs.PlateCarree())
    ax.add_geometries(SHP_DIF, ccrs.PlateCarree(), facecolor="#bee8ff", edgecolor="none", linewidth=0.5)
    
    for ln in np.arange(LONi-0.25, LONf, -0.25):
        for lt in np.arange(LATi-0.25, LATf, -0.25):
            pt = plt.plot(ln, lt, 'o', color='red', markeredgecolor='black', markerfacecolor='red', markersize=10)

    for ln in np.arange(-30-0.125-area_turbina, -30.5+0.125, -area_turbina):
        for lt in np.arange(-30-0.125-area_turbina, -30.5+0.125, -area_turbina):
            t1 = plt.plot(ln, lt, marker="1", markeredgecolor='blue', markersize=6)

    for ln in np.arange(-30.25-0.125-area_turbina, -30.75+0.125, -area_turbina):
        for lt in np.arange(-30.25-0.125-area_turbina, -30.75+0.125, -area_turbina):
            t2 = plt.plot(ln, lt, marker="1", markeredgecolor='purple', markersize=6)

    pt = mlines.Line2D([], [], color='red', marker='o', markersize=10, markeredgecolor='black', label="Grid points")
    t1 = mlines.Line2D([], [], color='blue', marker='1', markersize=15, label="Example Turbines 1")
    t2 = mlines.Line2D([], [], color='purple', marker='1', markersize=15, label="Example Turbines 2")
    
    plt.legend(handles=[pt, t1, t2])

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=1, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(LONi, LONf, -0.25))
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylocator = mticker.FixedLocator(np.arange(LATi, LATf, -0.25))

    plt.savefig("../SAIDAS/grade.png", bbox_inches='tight')

def plota_gridzoom():
    import pandas as pd
    import matplotlib.patches as mpatches
    plt.rcParams["figure.figsize"]=[10,10]
    LATi, LONi=-30.38012, -30.38012
    LATf, LONf=-30.625, -30.625
    bounds = [(LONi, LONf, LATi, LATf)]
    SHP_DIF = Reader("shps/Oceano_corte_ZEE/Oceano_ZEE.shp").geometries()
    area_turbina = (7*240)/(111.12*1000)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(*bounds, crs=ccrs.PlateCarree())
    ax.add_geometries(SHP_DIF, ccrs.PlateCarree(), facecolor="#bee8ff", edgecolor="none", linewidth=0.5)
    lts, lns = [], []
    for ln in np.arange(-30.25-0.125-area_turbina, -30.75+0.125, -area_turbina):
        for lt in np.arange(-30.25-0.125-area_turbina, -30.75+0.125, -area_turbina):
            t2 = plt.plot(ln, lt, marker="1", markeredgecolor='purple', markersize=6)
            # ax.add_patch(mpatches.Circle(xy=[ln, lt],
            #             linewidth=0.5, radius=area_turbina*0.5, color='blue',
            #             fill=False, transform=ccrs.PlateCarree(), zorder=3))
            lns.append(ln)
            lts.append(lt)
    df = pd.DataFrame({'lon': lns,'lat': lts})
    df.to_csv("../SAIDAS/grade_zoom.csv")
    # pt = mlines.Line2D([], [], color='red', marker='o', markersize=10, markeredgecolor='black', label="Grid points")
    # t1 = mlines.Line2D([], [], color='blue', marker='1', markersize=15, label="Example Turbines 1")
    # t2 = mlines.Line2D([], [], color='purple', marker='1', markersize=15, label="Example Turbines 2")
    # plt.legend(handles=[pt, t1, t2])
    plt.savefig("../SAIDAS/grade_zoom.png", bbox_inches='tight')

# plota_grid()
plota_gridzoom()

def plota_mapa_24h(var, tipo):
    plt.rcParams["figure.figsize"]=[10,10]    
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from matplotlib.colors import ListedColormap
    from matplotlib.colors import BoundaryNorm
    from cartopy.io.shapereader import Reader
    import matplotlib.ticker as mticker

    def configura(ax, z, hora, h, v):
        SHP_AS=Reader("shps/AS/america_do_sul.shp").geometries()
        SHP_BRASIL=Reader("shps/Brasil/estados_2010.shp").geometries()
        SHP_ZEE=Reader("shps/ZEE/ZEE.shp").geometries()
        ax.add_geometries(SHP_AS, ccrs.PlateCarree(), linewidth=0.4, facecolor="silver", edgecolor="grey", zorder=3.1)
        ax.add_geometries(SHP_BRASIL, ccrs.PlateCarree(), linewidth=0.4, facecolor="silver", edgecolor="grey", zorder=3.1)
        ax.add_geometries(SHP_ZEE, ccrs.PlateCarree(), linewidth=0.5, facecolor="none", edgecolor="black", zorder=3)
        ax.set_title(hora)
        
        levels = np.arange(0, 100.1, 10)
        cmap = "inferno"
        cf = ax.contourf(lons, lats, z, levels=levels, cmap=cmap, transform=ccrs.PlateCarree(), zorder=2)
         
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
        gl.xformatter = LONGITUDE_FORMATTER
        gl.xlocator = mticker.FixedLocator(np.arange(-50, -22, 4))
        gl.yformatter = LATITUDE_FORMATTER
        gl.ylocator = mticker.FixedLocator(np.arange(5, -36, -5))
        gl.xlines = False
        gl.ylines = False
        gl.xlabel_style = {'size': 8, 'color': 'black'}
        gl.ylabel_style = {'size': 8, 'color': 'black'}
            
        if v and h:
            gl.xlabels_top = False
            gl.ylabels_right = False
        elif v and not h:
            gl.xlabels_top = False
            gl.xlabels_bottom = False
            gl.ylabels_right = False
        elif not v and h:
            gl.xlabels_top = False
            gl.ylabels_right = False
            gl.ylabels_left = False
        elif not v and not h:
            gl.xlabels_top = False
            gl.xlabels_bottom = False
            gl.ylabels_right = False
            gl.ylabels_left = False
         
        return cf

    fig = plt.figure()
    gs = fig.add_gridspec(4, 6, hspace=0.15, wspace=0.05)
    ax1 = fig.add_subplot(gs[0,0], projection=ccrs.PlateCarree())
    ax2 = fig.add_subplot(gs[0,1], projection=ccrs.PlateCarree())
    ax3 = fig.add_subplot(gs[0,2], projection=ccrs.PlateCarree())
    ax4 = fig.add_subplot(gs[0,3], projection=ccrs.PlateCarree())
    ax5 = fig.add_subplot(gs[0,4], projection=ccrs.PlateCarree())
    ax6 = fig.add_subplot(gs[0,5], projection=ccrs.PlateCarree())    
    ax7 = fig.add_subplot(gs[1,0], projection=ccrs.PlateCarree())
    ax8 = fig.add_subplot(gs[1,1], projection=ccrs.PlateCarree())
    ax9 = fig.add_subplot(gs[1,2], projection=ccrs.PlateCarree())
    ax10 = fig.add_subplot(gs[1,3], projection=ccrs.PlateCarree())
    ax11 = fig.add_subplot(gs[1,4], projection=ccrs.PlateCarree())
    ax12 = fig.add_subplot(gs[1,5], projection=ccrs.PlateCarree())
    ax13 = fig.add_subplot(gs[2,0], projection=ccrs.PlateCarree())
    ax14 = fig.add_subplot(gs[2,1], projection=ccrs.PlateCarree())
    ax15 = fig.add_subplot(gs[2,2], projection=ccrs.PlateCarree())
    ax16 = fig.add_subplot(gs[2,3], projection=ccrs.PlateCarree())
    ax17 = fig.add_subplot(gs[2,4], projection=ccrs.PlateCarree())
    ax18 = fig.add_subplot(gs[2,5], projection=ccrs.PlateCarree())
    ax19 = fig.add_subplot(gs[3,0], projection=ccrs.PlateCarree())
    ax20 = fig.add_subplot(gs[3,1], projection=ccrs.PlateCarree())
    ax21 = fig.add_subplot(gs[3,2], projection=ccrs.PlateCarree())
    ax22 = fig.add_subplot(gs[3,3], projection=ccrs.PlateCarree())
    ax23 = fig.add_subplot(gs[3,4], projection=ccrs.PlateCarree())
    ax24 = fig.add_subplot(gs[3,5], projection=ccrs.PlateCarree())
    
    configura(ax1, var['00'], "00 UTC", False, True)
    configura(ax2, var['01'], "01 UTC", False, False)
    configura(ax3, var['02'], "02 UTC", False, False)
    configura(ax4, var['03'], "03 UTC", False, False)
    configura(ax5, var['04'], "04 UTC", False, False)
    configura(ax6, var['05'], "05 UTC", False, False)
    configura(ax7, var['06'], "06 UTC", False, True)
    configura(ax8, var['07'], "07 UTC", False, False)
    configura(ax9, var['08'], "08 UTC", False, False)
    configura(ax10, var['09'], "09 UTC", False, False)
    configura(ax11, var['10'], "10 UTC", False, False)
    configura(ax12, var['11'], "11 UTC", False, False)
    configura(ax13, var['12'], "12 UTC", False, True)
    configura(ax14, var['13'], "13 UTC", False, False)
    configura(ax15, var['14'], "14 UTC", False, False)
    configura(ax16, var['15'], "15 UTC", False, False)
    configura(ax17, var['16'], "16 UTC", False, False)
    configura(ax18, var['17'], "17 UTC", False, False)
    configura(ax19, var['18'], "18 UTC", True, True)
    configura(ax20, var['19'], "19 UTC", True, False)
    configura(ax21, var['20'], "20 UTC", True, False)
    configura(ax22, var['21'], "21 UTC", True, False)
    configura(ax23, var['22'], "22 UTC", True, False)
    cf = configura(ax24, var['23'], "23 UTC", True, False)

    axins = inset_axes(ax24, width="100%", height="100%", 
                        loc='lower left', borderpad=0,
                        bbox_to_anchor=(-5.6, -0.5, 6.5, 0.2),
                        bbox_transform=ax24.transAxes, axes_kwargs= {'title': "%"})
    fig.colorbar(cf, cax=axins, orientation='horizontal')

    plt.savefig("../SAIDAS/"+tipo+".png", bbox_inches='tight')


def plota_mapa_tecnico_regiao(pontos_dentro, bar_dentro, bounds=None):
    plt.rcParams["figure.figsize"]=[15,7]    
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from matplotlib.colors import ListedColormap
    from matplotlib.colors import BoundaryNorm
    from cartopy.io.shapereader import Reader
    import matplotlib.ticker as mticker

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
 
    if bounds:
        ax.set_extent(*bounds, crs=ccrs.PlateCarree())

    SHP_AS = Reader("shps/AS/america_do_sul.shp").geometries()
    SHP_BRASIL = Reader("shps/Brasil/estados_2010.shp").geometries()
    SHP_BATIMETRIA1000 = Reader("shps/Batimetria_1000m/Batimetria_1000m.shp").geometries()
    SHP_UC = Reader("shps/UC/UC_NE.shp").geometries()
    SHP_18KM = Reader("shps/18km/distancia-18km.shp").geometries()
    SHP_DIF = Reader("shps/Oceano_corte_ZEE/Oceano_ZEE.shp").geometries()
    SHP_ZEE = Reader("shps/ZEE/ZEE.shp").geometries()
   
    # ax.add_geometries(SHP_AS, ccrs.PlateCarree(), facecolor="silver", edgecolor="gray", zorder=3.1, linewidth=0.4)
    ax.add_geometries(SHP_BRASIL, ccrs.PlateCarree(), facecolor="none", edgecolor="gray", zorder=3.1, linewidth=0.4)
    # ax.add_geometries(SHP_BATIMETRIA1000, ccrs.PlateCarree(), facecolor="#483D8B", edgecolor="black", zorder=2.9, linewidth=0.5)
    # ax.add_geometries(SHP_UC, ccrs.PlateCarree(), facecolor="#32CD32", edgecolor="#32CD32", zorder=3, linewidth=0.5)
    # ax.add_geometries(SHP_DIF, ccrs.PlateCarree(), facecolor="#90EE90", edgecolor="black", zorder=3, linewidth=0.5)
    ax.add_geometries(SHP_ZEE, ccrs.PlateCarree(), facecolor="none", edgecolor="silver", zorder=2.8, linewidth=0.5)

    ax.text(-36.0, -34.4, "Atlantic Ocean", style='italic', weight='bold', size=12, c='white', zorder=4)

    #legenda(ax, -28.6,-38.1)    
    # legenda(ax, -9,-32.0)
    # levels2 = np.arange(-1,1)
    # cmap2 = ListedColormap(["#e06969","#e06969"])
    # ax.contourf(lons, lats, ws7, levels=levels2, cmap=cmap2, transform=ccrs.PlateCarree(), zorder=3) 

    # cmap3 = ListedColormap(["red","red"])
    # ax.contourf(lons, lats, ROCI113, levels=levels2, cmap=cmap3, transform=ccrs.PlateCarree(), zorder=3)  

    # ax.add_geometries(SHP_18KM, ccrs.PlateCarree(), facecolor="#4c1130", edgecolor="black", linewidth=0.5, zorder=3)

    
    for lon, lat in pontos_dentro[:]:
        ax.plot(lon, lat, '.', color='red', markersize=1, transform=ccrs.PlateCarree(), zorder=5)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(-20, -51, -5))
    gl.xlabels_top = False
    gl.xlines = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylocator = mticker.FixedLocator(np.arange(5, -36, -5))
    gl.ylabels_right = False
    gl.ylines = False

    plt.savefig("../SAIDAS/Mapa_Tecnico.png", bbox_inches='tight')
    
def plota_power_curve(dados, f, k, c):
    plt.rcParams["figure.figsize"]=[8, 5]
    
    plt.title("Weibull", size=18)
    plt.xlabel('wind speed (m/s)', size=14)
    plt.ylabel('Probability', size=14)
    
    xvalues = np.arange(0, 28, .1)
    xlim = [0, 27]
    plt.xlim(xlim)
    plt.xticks(np.arange(0, 28))
    ylim = [0, np.max(f)]
    plt.ylim(ylim)
    yvalues = np.arange(0, np.max(f)+0.1, .05)
    yticks = [str(round(pd*100))+"%" for pd in yvalues]
    plt.yticks(yvalues, yticks)
    
    plt.plot(np.arange(0, 28, .1), f, linewidth=1.5, color="red" , label="k = "+str(round(k,1))+", c = "+str(round(c,1)))
    plt.legend(fontsize=16)
    plt.savefig("../SAIDAS/WEIBULL.png", bbox_inches='tight')
    plt.show()
    
    plt.rcParams["figure.figsize"]=[8, 4]
    plt.title("Power Curve", size=18)
    plt.ylabel('Turbine Power (MW)', size=14)
    plt.xlabel('wind speed (m/s)', size=14)
    
    plt.ylim(dados['Power (MW)'].min(), 16)
    plt.xlim([0, 27])
    plt.yticks(np.arange(1,17,2))
    plt.xticks(np.arange(0,28,1))
    
    plt.plot(dados['Wind (m/s)'], dados['Power (MW)'], label='IEC Class IB', color='blue')
    plt.legend(loc='best')
    plt.savefig("../SAIDAS/POWER_CURVE.png", bbox_inches='tight')
    plt.show()
