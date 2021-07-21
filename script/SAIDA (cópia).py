#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 22:40:43 2020

@author: administrador
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr

import pandas as pd
import shapefile as shp
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

PATH_HD="/media/administrador/6B8177A501318231/"
meses=["Jan","Fev","Mar","Abr","Mai","Jun","Jul","Ago","Set","Out","Nov","Dez",]
mesoes=["JAN","FEV","MAR","ABR","MAI","JUN","JUL","AGO","SET","OUT","NOV","DEZ",]

lons=np.load(PATH_HD+"DISSERTACAO/ARRAY/clons.npy")
lats=np.load(PATH_HD+"DISSERTACAO/ARRAY/clats.npy")

# PONTOS=np.load(PATH_HD+"DISSERTACAO/ARRAY/PONTOS.npy")
# PONTOS_COMPLEMENTARES=np.load(PATH_HD+"DISSERTACAO/ARRAY/PONTOS_COMPLEMENTARES.npy")
# sf=shp.Reader(PATH_HD+"SHAPES/AS - Brasil2/america_do_sul.shp")
# sf=shp.Reader("/home/administrador/Documentos/shapes/estados_br/BRUFE250GC_SIR.shp")

# fields = [x[0] for x in sf.fields][1:]
# records = [list(i) for i in sf.records()]
# shps = [s.points for s in sf.shapes()]
# df = pd.DataFrame(columns=fields, data=records)
# df = df.assign(coords=shps)
# cont=0
# xxx, yyy=[],[]
# for xy in PONTOS:
#     dentro=False
#     for index in range(0, len(df)):
#         polygon = Polygon(df.at[index,'coords'])
#         x=xy[0]
#         y=xy[1]
#         if polygon.contains(Point(x,y)):
#             dentro=True
#     if not dentro:
#         xxx.append(x)
#         yyy.append(y)    
#         plt.plot(x, y, ".", zorder=3)
#         cont=cont+1

# salvar=np.zeros((cont,2))
# cont=0
# for x,y in zip(xxx, yyy):
#     salvar[cont][0]=x
#     salvar[cont][1]=y
#     cont=cont+1

def plota_sat(mes, vel, u, v, lats, lons):
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    from matplotlib.colors import ListedColormap
    from cartopy.io.shapereader import Reader
    import matplotlib.ticker as mticker
    import math

    plt.rcParams["figure.figsize"]=[10,10]    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    LATi, LONi=0, -49
    LATf, LONf=-18, -30

    LATi, LONi=30, -42
    LATf, LONf=-50, -11

    bounds = [(LONi, LONf, LATi, LATf)]
    ax.set_extent(*bounds, crs=ccrs.PlateCarree())

    SHP_BRASIL=Reader(PATH_HD+"SHAPES/Brasil/estados_2010.shp").geometries()
    SHP_ZEE=Reader(PATH_HD+"SHAPES/ZEE/ZEE.shp").geometries()

    ax.add_geometries(SHP_ZEE, ccrs.PlateCarree(), linewidth=1, facecolor="none", edgecolor="w", zorder=3)
    ax.add_geometries(SHP_BRASIL, ccrs.PlateCarree(), linewidth=1, facecolor="w", edgecolor="gray", zorder=3)
    #ax.text(-47, -16.5, mes,
    #    zorder=10, fontsize=22, ha="center", va="center", family='monospace', weight='bold',
    #    bbox=dict(boxstyle="square", ec="black", fc="gray")
    #)
    
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
    plt.savefig(PATH_HD+"DISSERTACAO/SAIDAS/SATELITE/"+nome_arq+".png", dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()


def ponto(ax, lat, lon, texto, color, mk='x'):
    ax.plot(lon,lat, mk, markerfacecolor=color, markeredgecolor=color, markersize=5)
    ax.text(lon-0.8,lat+0.5, texto, color=color)

def legenda(ax, lat, lon):
    ax.text(lon, lat,
    "      Área Utilizável\n"+
    "      Estados brasileiros\n"+
    "Áreas excluídas dentro da ZEE:\n"+
    "      Área de Proteção Marinha\n"+
    "      Batimetria > 1000 m\n"+
    "      18 km da costa\n"+
    "      Velocidade do vento < 7 m/s",
    zorder=6, fontsize=12, ha="left", va="center", family='fantasy',
    bbox=dict(boxstyle="square", ec="black", fc="white")
    )
    ax.plot(lon+.3, lat+1.45, "s", markerfacecolor='orange', markeredgecolor='k', markersize=10, zorder=7)
    ax.plot(lon+.3, lat+1, "s", markerfacecolor='w', markeredgecolor='k', markersize=10, zorder=7)
    ax.plot(lon+.3, lat, "s", markerfacecolor='yellowgreen', markeredgecolor='k', markersize=10, zorder=7)
    ax.plot(lon+.3, lat-.45, "x", markerfacecolor='blue', markeredgecolor='b', markersize=10, zorder=7)
    ax.plot(lon+.3, lat-.9, "s", markerfacecolor='gray', markeredgecolor='k', markersize=10, zorder=7)
    ax.plot(lon+.3, lat-1.35, "s", markerfacecolor='#fbea51', markeredgecolor='k', markersize=10, zorder=7)

def plota_mapa(var, mes, tipo):
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    from matplotlib.colors import ListedColormap
    from matplotlib.colors import BoundaryNorm
    from cartopy.io.shapereader import Reader
    import matplotlib.ticker as mticker
    import math
    import geopandas as gpd

    u, v = [], []
    norm = ""
    plt.rcParams["figure.figsize"]=[10,10]    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    
    LONi, LONf = -53.5, -23.0
    LATi, LATf = 7.0, -36.75

    bounds = [(LONi, LONf, LATi, LATf)]
    bounds = [LONi, LONf, LATi, LATf]
    ax.set_extent(bounds, crs=ccrs.PlateCarree())

    SHP_AS=Reader(PATH_HD+"SHAPES/AS/america_do_sul.shp").geometries()
    SHP_BRASIL=Reader(PATH_HD+"SHAPES/Brasil/estados_2010.shp").geometries()
    SHP_ZEE=Reader(PATH_HD+"SHAPES/ZEE/ZEE.shp").geometries()
    ax.add_geometries(SHP_AS, ccrs.PlateCarree(), linewidth=0.4, facecolor="silver", edgecolor="grey", zorder=3.1)
    ax.add_geometries(SHP_BRASIL, ccrs.PlateCarree(), linewidth=0.4, facecolor="silver", edgecolor="grey", zorder=3.1)
    ax.add_geometries(SHP_ZEE, ccrs.PlateCarree(), linewidth=0.5, facecolor="none", edgecolor="black", zorder=3)

    # c=gpd.read_file(PATH_HD+"SHAPES/Brasil/estados_2010.shp")
    # c['coords'] = c['geometry'].apply(lambda x: x.representative_point().coords[:])
    # c['coords'] = [coords[0] for coords in c['coords']]
    # for idx, row in c.iterrows():
    #     if row['sigla'] not in ['AC','AM','RO','RR','MS','MT']:
    #         llon = row['coords'][0]
    #         llat = row['coords'][1]
    #         if row['sigla']=='RS':
    #             llon += 1.5
    #         ax.text(llon, llat, row['sigla'],
    #             fontsize=10, weight="bold", color='#696969', horizontalalignment='center', zorder=7)

    if tipo in ["AVA_DPE", "AVA_ROCI", "WCS", "SCW", "WSS"]:
        titulo = mes
        z = var
        levels = np.arange(0, 101, 5)
        cmap, label, cor_ZEE = "Blues", "%", "red"
        nome_arq = tipo

    if tipo in ["Median_DPE", "Median_ROCI"]:
        titulo = mes
        lmax, z = var
        levels = np.arange(0, lmax, 100)
        cmap, label, cor_ZEE = "Blues", "W/m²", "red"
        nome_arq = tipo
    if tipo in ["RCOV_DPE", "RCOV_ROCI"]:
        titulo = mes
        lmax, z = var
        levels = np.arange(0, lmax, 0.1)
        cmap, label, cor_ZEE = "Blues", "RCOV", "red"
        nome_arq = tipo
    if tipo in ["IQR_DPE", "IQR_ROCI"]:
        titulo = mes
        lmax, z = var
        levels = np.arange(0, lmax, 100)
        cmap, label, cor_ZEE = "Blues", "W/m²", "red"
        nome_arq = tipo

    if tipo==1:
        titulo = 'Área de estudo - Costa do Brasil'
        z = np.zeros((len(lats), len(lons)))
        levels = np.arange(0,2)
        cores = ["#B2DFEE","#B2DFEE"]
        cmap, label, cor_ZEE = ListedColormap(cores), False, "black"
        # legenda(ax, -48.3,-16.0)
        nome_arq = "MAPA"

    if tipo in [2, 3]:
        titulo = "WS a 1"+str("0")*(tipo-1)+" m e direção (setas)"
        u, v, z = var
        levels = np.arange(0,13,0.5)
        cmap, label, cor_ZEE = "jet", "m/s", "red"
        nome_arq = "WSD1"+str("0")*(tipo-1)+"_"+str(mes)

    if tipo == 4:
        titulo = 'Densidade de Potencial Eólico' # Densidade de Potencial Eólico
        z = var        
        colors = ["#FFFE9E","#FAE27B","#F7C95C","#F4B046","#F1963C","#EB7D48","#DB685E","#C9566A","#B5456F","#9F3870","#882E6C","#702664","#572159","#401C4B","#29173B","#151128","#040404"]
        from matplotlib.colors import ListedColormap
        cmap = ListedColormap(colors)
        levels1 = [i for i in range(0,200,20)]
        levels2 = [i for i in range(200,500,100)]
        levels3 = [i for i in range(500,1251,250)]
        levels = levels1+levels2+levels3
        norm = BoundaryNorm(levels, cmap.N)
        label, cor_ZEE = "W/m²", "black"
        nome_arq="DPE_"+str(mes+1)

    if tipo == 5:
        titulo = 'Radiação de Onda Curta Incidente' #Radiação de Onda Curta Incidente
        z = var
        colors = ["#FFFE9E","#FAE27B","#F7C95C","#F4B046","#F1963C","#EB7D48","#DB685E","#C9566A","#B5456F","#9F3870","#882E6C","#702664","#572159","#401C4B","#29173B","#151128","#040404"]
        from matplotlib.colors import ListedColormap
        cmap = ListedColormap(colors)
        levels1 = [i for i in range(60,100,10)]
        levels2 = [i for i in range(100,250,30)]
        levels3 = [i for i in range(250,330,10)]
        levels = levels1+levels2+levels3
        label, cor_ZEE = "W/m²", "black"
        nome_arq="ROCI_"+str(mes+1)

    if tipo in [20, 30]:
        titulo = "Média Sazonal de WS a 1"+str("0")*int((tipo/10)-1)+" m e direção (setas)"
        u, v, z= var
        levels = np.arange(0,13,0.5)
        cmap, label, cor_ZEE = "jet", "m/s", "red"
        nome_arq="SZN_WSD1"+str("0")*int((tipo/10)-1)+"_"+mes

    if tipo==40:
        titulo = "Média Sazonal de DPE"
        z = var
        colors = ["#FFFE9E","#FAE27B","#F7C95C","#F4B046","#F1963C","#EB7D48","#DB685E","#C9566A","#B5456F","#9F3870","#882E6C","#702664","#572159","#401C4B","#29173B","#151128","#040404"]
        from matplotlib.colors import ListedColormap
        cmap = ListedColormap(colors)
        levels1 = [i for i in range(0,200,20)]
        levels2 = [i for i in range(200,500,100)]
        levels3 = [i for i in range(500,1251,250)]
        levels = levels1+levels2+levels3
        norm = BoundaryNorm(levels, cmap.N)
        label, cor_ZEE = "W/m²", "black"
        nome_arq="SZN_DPE_"+mes

    if tipo==50:
        titulo = "Média Sazonal de ROCI"
        z = var
        colors = ["#FFFE9E","#FAE27B","#F7C95C","#F4B046","#F1963C","#EB7D48","#DB685E","#C9566A","#B5456F","#9F3870","#882E6C","#702664","#572159","#401C4B","#29173B","#151128","#040404"]
        from matplotlib.colors import ListedColormap
        cmap = ListedColormap(colors)
        levels = [i for i in range(80,250,10)]
        label, cor_ZEE = "W/m²", "black"
        nome_arq="SZN_ROCI_"+mes

    if tipo in [200, 300]:
        titulo = "Média Anual de WS a 1"+str("0")*int((tipo/100)-1)+" m e direção (setas)"
        u, v, z= var
        levels = np.arange(0,13,0.5)
        cmap, label, cor_ZEE = "jet", "m/s", "red"
        nome_arq="ANUAL_WSD1"+str("0")*int((tipo/100)-1)+"_"+mes

    if tipo==400:
        titulo = "Média Anual de DPE"
        z= var
        colors = ["#FFFE9E","#FAE27B","#F7C95C","#F4B046","#F1963C","#EB7D48","#DB685E","#C9566A","#B5456F","#9F3870","#882E6C","#702664","#572159","#401C4B","#29173B","#151128","#040404"]
        from matplotlib.colors import ListedColormap
        cmap = ListedColormap(colors)
        levels1 = [i for i in range(0,200,20)]
        levels2 = [i for i in range(200,600,100)]
        levels3 = [i for i in range(600,1001,200)]
        levels = levels1+levels2+levels3
        norm = BoundaryNorm(levels, cmap.N)
        label, cor_ZEE = "W/m²", "black"
        nome_arq="ANUAL_DPE_"+mes

    if tipo==500:
        titulo = "Média Anual de ROCI"
        z= var
        colors = ["#FFFE9E","#FAE27B","#F7C95C","#F4B046","#F1963C","#EB7D48","#DB685E","#C9566A","#B5456F","#9F3870","#882E6C","#702664","#572159","#401C4B","#29173B","#151128","#040404"]
        from matplotlib.colors import ListedColormap
        cmap = ListedColormap(colors)
        levels = [i for i in range(150,320,10)]
        label, cor_ZEE = "W/m²", "black"
        nome_arq="ANUAL_ROCI_"+mes

    if tipo in [2, 3, 4, 5]:
        # subtitulo = meses[mes]
        subtitulo = mes
    if tipo in [20, 30, 40, 50]:
        subtitulo = meses[mes]
    if tipo in [200, 300, 400, 500]:
        subtitulo = mes

    if tipo in np.arange(2, 501):
        ax.text(-28, -35, subtitulo,
            zorder=10, fontsize=16, ha="center", va="center", family='monospace', weight='bold',
            bbox=dict(boxstyle="square", ec="black", fc="gray")
        )

    ax.add_geometries(SHP_ZEE, ccrs.PlateCarree(), linewidth=1, facecolor="none", edgecolor=cor_ZEE, zorder=3)
    ax.set_title(titulo, loc='left')

    if norm!="":
        cf = ax.contourf(lons, lats, z, levels=levels, cmap=cmap, norm=norm, transform=ccrs.PlateCarree(), zorder=2)
    else:
        cf = ax.contourf(lons, lats, z, levels=levels, cmap=cmap, transform=ccrs.PlateCarree(), zorder=2)    

    if len(u)>0 and len(v)>0:
        ax.streamplot(lons, lats, u, v, density=3, linewidth=1, color='k', transform=ccrs.PlateCarree(), zorder=2)

    if label:
        cb = fig.colorbar(cf, orientation='horizontal', extend='max', aspect=65, shrink=0.5, pad=0.05, extendrect='True')
        cb.set_label(label)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(-25, -56, -5))
    gl.xlabels_top = False
    gl.xlines = False

    gl.yformatter = LATITUDE_FORMATTER
    gl.ylocator = mticker.FixedLocator(np.arange(5, -40, -5))
    gl.ylabels_right = False
    gl.ylines = False
    gl.ylabel_style = {'rotation':'vertical'}

    # plt.tight_layout()
    plt.savefig(PATH_HD+"DISSERTACAO/ARTIGO/"+nome_arq+".png", dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()

def plota_linear_pontos(V10, V100, rad, szn_V10, szn_V100, szn_RAD, PONTOS_VENTOMAIOR7MS):
    plt.rcParams["figure.figsize"]=[10,5]
    tam=len(PONTOS_VENTOMAIOR7MS)
    cont=0    
    global PONTOS_COMPLEMENTARES
    PONTOS_COMPLEMENTARES=np.zeros((60,2))
    for xy in PONTOS_VENTOMAIOR7MS:
        x=xy[0]
        y=xy[1]
        ponto_V10=V10[:,   list(lats.data).index(y),    list(lons.data).index(x)]
        ponto_V100=V100[:, list(lats.data).index(y),    list(lons.data).index(x)]
        ponto_RAD=rad[:,   list(lats.data).index(y),    list(lons.data).index(x)]

        # V100_VERAO=(ponto_V100[11]+ponto_V100[0]+ponto_V100[1])/3
        # V100_OUTONO=(ponto_V100[2]+ponto_V100[3]+ponto_V100[4])/3
        # V100_INVERNO=(ponto_V100[5]+ponto_V100[6]+ponto_V100[7])/3
        # V100_PRIMAVERA=(ponto_V100[8]+ponto_V100[9]+ponto_V100[10])/3
        # RAD_VERAO=(ponto_RAD[11]+ponto_RAD[0]+ponto_RAD[1])/3
        # RAD_OUTONO=(ponto_RAD[2]+ponto_RAD[3]+ponto_RAD[4])/3
        # RAD_INVERNO=(ponto_RAD[5]+ponto_RAD[6]+ponto_RAD[7])/3
        # RAD_PRIMAVERA=(ponto_RAD[8]+ponto_RAD[9]+ponto_RAD[10])/3

        V100_VERAO=(ponto_V100[11]+ponto_V100[0]+ponto_V100[1])/3
        V100_OUTONO=(ponto_V100[2]+ponto_V100[3]+ponto_V100[4])/3
        V100_INVERNO=(ponto_V100[5]+ponto_V100[6]+ponto_V100[7])/3
        V100_PRIMAVERA=(ponto_V100[8]+ponto_V100[9]+ponto_V100[10])/3
        RAD_VERAO=(ponto_RAD[11]+ponto_RAD[0]+ponto_RAD[1])/3
        RAD_OUTONO=(ponto_RAD[2]+ponto_RAD[3]+ponto_RAD[4])/3
        RAD_INVERNO=(ponto_RAD[5]+ponto_RAD[6]+ponto_RAD[7])/3
        RAD_PRIMAVERA=(ponto_RAD[8]+ponto_RAD[9]+ponto_RAD[10])/3

        # if  V100_INVERNO>V100_VERAO and RAD_VERAO>RAD_INVERNO:
            # PONTOS_COMPLEMENTARES[cont][0]=x
            # PONTOS_COMPLEMENTARES[cont][1]=y
        # if  V100_INVERNO>V100_VERAO and V100_INVERNO>V100_OUTONO and V100_INVERNO>V100_PRIMAVERA and RAD_VERAO>RAD_OUTONO and RAD_VERAO>RAD_INVERNO and RAD_VERAO>RAD_PRIMAVERA:
            # fig, ax=plt.subplots()
            # ax2=ax.twinx()    
            # datevar=np.arange(1,13)    
            # ax.set_title("Vel. Vento a 10m/100m (m/s) e Radiação (W/m²) p/ LAT:"+str(y)+" LON:"+str(x),loc='left')
    
            # ax.plot(datevar, ponto_V10, linewidth=0.7, color='black')
            # ax.plot(datevar, ponto_V100, linewidth=0.7, color='blue')
            # ax2.plot(datevar, ponto_RAD, linewidth=0.7, color='red')
    
            # ax.set_ylabel('Velocidade (m/s)', multialignment='center', fontsize=12)    
            # ax.set_ylim(0, 12)        
            # ax2.set_ylabel('Radiação (W/m²)', multialignment='center', fontsize=12)
            # ax2.set_ylim(100, 305)
        
            # ax.set_xlim(datevar[0], datevar[-1])
            # ax.set_xticks(datevar)
            # ax.set_xticklabels(meses, fontsize=12)        
    
            # ax.text(9, 2, "      vento 10 m\n      vento 100 m\n      radiação",
            #     zorder=6, fontsize=12, ha="left", va="center", family='fantasy',
            #     bbox=dict(boxstyle="square", ec="black", fc="white")
            # )
            # ax.plot(9.2, 2.65, "s", markerfacecolor='black', markeredgecolor='k', markersize=10, zorder=7)
            # ax.plot(9.2, 2.05, "s", markerfacecolor='blue', markeredgecolor='k', markersize=10, zorder=7)
            # ax.plot(9.2, 1.4, "s", markerfacecolor='red', markeredgecolor='k', markersize=10, zorder=7)
            # plt.savefig(PATH_HD+"DISSERTACAO/SAIDAS/137/Perfis_"+str(cont)+".png", dpi=300, bbox_inches='tight', pad_inches=0)
            # plt.show()
            # cont+=1