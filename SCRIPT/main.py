import os

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
import math
from constants import path_wnd, path_rad, wpdmin, spdmin
# import output

start = datetime.now()


def ln(x):
    from math import log
    return log(x)


def weib(v,k,c):
    return (k/c) * (v/c)**(k-1) * np.exp(-(v/c)**k)


def DPE(v):
    return 0.5*1.225*v**3


def extrapola(ws):
    return ws * ln(100/0.0002)/ln(10/0.0002)

w_file_list = [file_name for file_name in sorted(os.listdir(path_wnd))
                         if file_name.endswith('.nc')]
r_file_list = [file_name for file_name in sorted(os.listdir(path_rad))
                         if file_name.endswith('.nc')]


for wnd_file, rad_file in zip(w_file_list, r_file_list):
    year = wnd_file[7:11]
    print(year)
    wnd_ds = xr.open_dataset(path_wnd + wnd_file)
    rad_ds = xr.open_dataset(path_rad + rad_file)

    time = wnd_ds['time'].data
    tnh = len(time)

    r_data = rad_ds['ssrd'][:,:,:].data
    spd_data = r_data/3600

    u_data = wnd_ds['u10'][:,:,:].data
    v_data = wnd_ds['v10'][:,:,:].data
    wnd10_data = (u_data**2 + v_data**2)**0.5
    wnd100_data = extrapola(wnd10_data)

    nh_spd_available = spd_data.copy()
    nh_spd_available[nh_spd_available < spdmin] = 0
    nh_spd_available[nh_spd_available >= spdmin] = 1
    nh_spd_available = np.sum(nh_spd_available, axis=0)
    solar_availability = (nh_spd_available/tnh)*100

    nh_wpd_available = wnd100_data.copy()
    nh_wpd_available[nh_wpd_available < wpdmin] = 0
    nh_wpd_available[nh_wpd_available >= wpdmin] = 1
    nh_wpd_available = np.sum(nh_wpd_available, axis=0)
    wind_availability = (nh_wpd_available/tnh)*100
    
    wnd_ds.close()
    rad_ds.close()

    # output.plota_mapa(wind_availability, "Annual_WPD_Availability_" + year)
    # output.plota_mapa(solar_availability, "Annual_SPD_Availability_" + year)

    break



RPD_WS = np.zeros(r_data.shape)
RPD_R = np.zeros((ym, xm))
RPD_R[R >= limiar_r] = 1            

RPD_WS_NR = np.zeros((ym, xm))
RPD_WS_NR[(WS >= limiar_v) & (R < limiar_r)] = 1

RPD_NWS_R = np.zeros((ym, xm))
RPD_NWS_R[(WS < limiar_v) & (R >= limiar_r)] = 1

RPD_WS_R = np.zeros((ym, xm))
RPD_WS_R[(((WS >= limiar_v) & (R < limiar_r)) | ((WS < limiar_v) & (R >= limiar_r)))] = 1



#%%                       Carrega dados horários ERA5                         #

"""

NHRPD_WS,
NHRPD_R
NHRPD_WS_NR,
NHRPD_NWS_R,
NHRPD_WS_R 


for t in range(0, tam):
    WS, R = ws100[t], r[t]
    hora = str(tt[t]).split("T")[1][:2]

    RPD_WS = np.zeros((ym, xm))
    RPD_WS[WS >= limiar_v] = 1            

    RPD_R = np.zeros((ym, xm))
    RPD_R[R >= limiar_r] = 1            

    RPD_WS_NR = np.zeros((ym, xm))
    RPD_WS_NR[(WS >= limiar_v) & (R < limiar_r)] = 1

    RPD_NWS_R = np.zeros((ym, xm))
    RPD_NWS_R[(WS < limiar_v) & (R >= limiar_r)] = 1

    RPD_WS_R = np.zeros((ym, xm))
    RPD_WS_R[(((WS >= limiar_v) & (R < limiar_r)) | ((WS < limiar_v) & (R >= limiar_r)))] = 1

    if hora not in NHRPD_WS_R:
        NHRPD_WS[hora] = RPD_WS
        NHRPD_R[hora] = RPD_R
        NHRPD_WS_NR[hora] = RPD_WS_NR
        NHRPD_NWS_R[hora] = RPD_NWS_R
        NHRPD_WS_R[hora] = RPD_WS_R
    else:
        NHRPD_WS[hora] += RPD_WS
        NHRPD_R[hora] += RPD_R
        NHRPD_WS_NR[hora] += RPD_WS_NR
        NHRPD_NWS_R[hora] += RPD_NWS_R
        NHRPD_WS_R[hora] += RPD_WS_R
    if hora == "23":
        ttotal += 1



DISP_WS, DISP_R, WCS, SCW, WSS = {}, {}, {}, {}, {}
for h in ['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23']:
    if a==1:
        DISP_WS[h] = (NHRPD_WS[h]/ttotal)*100
        DISP_R[h] = (NHRPD_R[h]/ttotal)*100
        WCS[h] = (NHRPD_WS_NR[h]/ttotal)*100
        SCW[h] = (NHRPD_NWS_R[h]/ttotal)*100
        WSS[h] = (NHRPD_WS_R[h]/ttotal)*100

        np.save('../DADOS/NPY/24H/DISP_WS'+h+'.npy', DISP_WS[h])
        np.save('../DADOS/NPY/24H/DISP_R'+h+'.npy', DISP_R[h])
        np.save('../DADOS/NPY/24H/WCS_'+h+'.npy', WCS[h])
        np.save('../DADOS/NPY/24H/SCW_'+h+'.npy', SCW[h])
        np.save('../DADOS/NPY/24H/WSS_'+h+'.npy', WSS[h])
    else:
        # DISP_WS[h] = np.load('../DADOS/NPY/24H/DISP_WS'+h+'.npy')
        # DISP_R[h] = np.load('../DADOS/NPY/24H/DISP_R'+h+'.npy')
        WCS[h] = np.load('../DADOS/NPY/24H/WCS_'+h+'.npy')
        SCW[h] = np.load('../DADOS/NPY/24H/SCW_'+h+'.npy')
        WSS[h] = np.load('../DADOS/NPY/24H/WSS_'+h+'.npy')

import SAIDA as OUT
OUT.plota_mapa_24h(DISP_WS, 'Wind_Disponibility_24H')
OUT.plota_mapa_24h(DISP_R, 'Solar_Disponibility_24H')
OUT.plota_mapa_24h(WCS, 'Wind_Complements_Solar_24H')
OUT.plota_mapa_24h(SCW, 'Solar_Complements_Wind_24H')
OUT.plota_mapa_24h(WSS, 'Wind_and_Solar_Synergy_24H')
"""

#%%                        Número de horas com RPD                            #

"""
limiar_r = 113*2
limiar_v = 7
ttotal = 0
# anoi, anof=1990, 1991
anoi, anof=1990, 2020
a = 12
if a==1:
    for ano in range(anoi, anof):
        print(ano)
        netcdf_r = xr.open_dataset("../DADOS/RAD/ERA5_RAD_"+str(ano)+".nc")
        netcdf_v = xr.open_dataset("../DADOS/VENTO/ERA_UV_"+str(ano)+".nc")
        tt = netcdf_r['time'].data
        tam = len(tt)
        r = np.array(netcdf_r.variables['ssrd'][:tam,:,:]/3600)
        u = np.array(netcdf_v.variables['u10'][:tam,:,:])
        v = np.array(netcdf_v.variables['v10'][:tam,:,:])
        ws = np.sqrt(u*u+v*v)
        ws100 = ws * (ln(100/0.0002)/ln(10/0.0002))

        for t in range(0, tam):
            WS, R = ws100[t], r[t]

            # RPD_WS = np.zeros((ym, xm))
            # RPD_WS[WS >= limiar_v] = 1            

            RPD_R = np.zeros((ym, xm))
            RPD_R[R >= limiar_r] = 1            

            RPD_WS_NR = np.zeros((ym, xm))
            RPD_WS_NR[(WS >= limiar_v) & (R < limiar_r)] = 1

            RPD_NWS_R = np.zeros((ym, xm))
            RPD_NWS_R[(WS < limiar_v) & (R >= limiar_r)] = 1

            RPD_WS_R = np.zeros((ym, xm))
            RPD_WS_R[(((WS >= limiar_v) & (R < limiar_r)) | ((WS < limiar_v) & (R >= limiar_r)))] = 1

            if ttotal==0:
                # NHRPD_WS = RPD_WS
                NHRPD_R = RPD_R
                NHRPD_WS_NR = RPD_WS_NR
                NHRPD_NWS_R = RPD_NWS_R
                NHRPD_WS_R = RPD_WS_R
            else:
                # NHRPD_WS += RPD_WS
                NHRPD_R += RPD_R
                NHRPD_WS_NR += RPD_WS_NR
                NHRPD_NWS_R += RPD_NWS_R
                NHRPD_WS_R += RPD_WS_R
            ttotal += 1
        time.sleep(10)

    # DISP_WS = (NHRPD_WS/ttotal)*100
    DISP_R = (NHRPD_R/ttotal)*100
    WCS = (NHRPD_WS_NR/ttotal)*100
    SCW = (NHRPD_NWS_R/ttotal)*100
    WSS = (NHRPD_WS_R/ttotal)*100

    # np.save('../DADOS/NPY/30ANOS/NHRPD_WS.npy', NHRPD_WS)
    # np.save('../DADOS/NPY/30ANOS/NHRPD_R.npy', NHRPD_R)
    # np.save('../DADOS/NPY/30ANOS/NHRPD_WS_NR.npy', NHRPD_WS_NR)
    # np.save('../DADOS/NPY/30ANOS/NHRPD_NWS_R.npy', NHRPD_NWS_R)
    # np.save('../DADOS/NPY/30ANOS/NHRPD_WS_R.npy', NHRPD_WS_R)
    # np.save('../DADOS/NPY/30ANOS/ttotal.npy', np.array(ttotal))
    
    # np.save('../DADOS/NPY/30ANOS/DISP_WPD.npy', DISP_WS)
    np.save('../DADOS/NPY/30ANOS/DISP_SPD2.npy', DISP_R)
    np.save('../DADOS/NPY/30ANOS/WCS2.npy', WCS)
    np.save('../DADOS/NPY/30ANOS/SCW2.npy', SCW)
    np.save('../DADOS/NPY/30ANOS/WSS2.npy', WSS)
else:
    DISP_WS = np.load('../DADOS/NPY/30ANOS/DISP_WPD.npy')
    DISP_R = np.load('../DADOS/NPY/30ANOS/DISP_SPD2.npy')
    WCS = np.load('../DADOS/NPY/30ANOS/WCS2.npy')
    SCW = np.load('../DADOS/NPY/30ANOS/SCW2.npy')
    WSS = np.load('../DADOS/NPY/30ANOS/WSS2.npy')

import SAIDA as OUT
OUT.plota_mapa(DISP_WS, "Annual_WPD_Disponibility")
OUT.plota_mapa(DISP_R, "Annual_SPD_Disponibility")
OUT.plota_mapa(WCS, "Annual_Wind_Complements_Solar")
OUT.plota_mapa(SCW, "Annual_Solar_Complements_Wind")
OUT.plota_mapa(WSS, "Annual_Wind_and_Solar_Synergy")
"""

#%%                        Número de horas com RPD                            #

"""
limiar_r = 113
limiar_v = 7
ttotal = 0
anoi, anof=1990, 1991
# anoi, anof=1990, 2020
a = 1 
for ano in range(anoi, anof):
    print(ano)
    netcdf_r = xr.open_dataset("../DADOS/RAD/ERA5_RAD_"+str(ano)+".nc")
    netcdf_v = xr.open_dataset("../DADOS/VENTO/ERA_UV_"+str(ano)+".nc")
    tt = netcdf_r['time'].data
    if anof-anoi == 2020-1990:
        tam = len(tt)
    if anof-anoi==1991-1990:
        tam = 31*24
    if a==1:
        r = np.array(netcdf_r.variables['ssrd'][:tam,:,:]/3600)
        u = np.array(netcdf_v.variables['u10'][:tam,:,:])
        v = np.array(netcdf_v.variables['v10'][:tam,:,:])
        ws = np.sqrt(u*u+v*v)
        ws100 = ws * (ln(100/0.0002)/ln(10/0.0002))

        NHRPD_WS = np.zeros((ym, xm))
        NHRPD_R = np.zeros((ym, xm))
        NHRPD_WS_NR = np.zeros((ym, xm))
        NHRPD_NWS_R = np.zeros((ym, xm))
        NHRPD_WS_R = np.zeros((ym, xm))

        for t in range(u.shape[0]):
            WS, R = ws100[t], r[t]
    
            RPD_WS = np.zeros((ym, xm))
            RPD_WS[WS >= limiar_v] = 1
            NHRPD_WS = NHRPD_WS + RPD_WS

            RPD_R = np.zeros((ym, xm))
            RPD_R[R >= limiar_r] = 1
            NHRPD_R = NHRPD_R + RPD_R
       
            RPD_WS_NR = np.zeros((ym, xm))
            RPD_WS_NR[(WS >= limiar_v) & (R < limiar_r)] = 1
            NHRPD_WS_NR = NHRPD_WS_NR + RPD_WS_NR

            RPD_NWS_R = np.zeros((ym, xm))
            RPD_NWS_R[(WS < limiar_v) & (R >= limiar_r)] = 1        
            NHRPD_NWS_R = NHRPD_NWS_R + RPD_NWS_R

            RPD_WS_R = np.zeros((ym, xm))    
            RPD_WS_R[(((WS >= limiar_v) & (R < limiar_r)) | ((WS < limiar_v) & (R >= limiar_r)))] = 1
            #RPD_WS_R[(WS >= limiar_v) & (R >= limiar_r)] = 1
            NHRPD_WS_R = NHRPD_WS_R + RPD_WS_R
    else:
        NHRPD_WS = np.load('../DADOS/NPY/30ANOS/NHRPD_WS.npy')
        NHRPD_R = np.load('../DADOS/NPY/30ANOS/NHRPD_R.npy')
        NHRPD_WS_NR = np.load('../DADOS/NPY/30ANOS/NHRPD_WS_NR.npy')
        NHRPD_NWS_R = np.load('../DADOS/NPY/30ANOS/NHRPD_NWS_R.npy')
        NHRPD_WS_R = np.load('../DADOS/NPY/30ANOS/NHRPD_WS_R.npy')
    ttotal += tam

if a==1:
    DISP_WS = (NHRPD_WS/ttotal)*100
    DISP_R = (NHRPD_R/ttotal)*100
    WCS = (NHRPD_WS_NR/ttotal)*100
    SCW = (NHRPD_NWS_R/ttotal)*100
    WSS = (NHRPD_WS_R/ttotal)*100

    np.save('../DADOS/NPY/30ANOS/DISP_WS.npy', DISP_WS)
    np.save('../DADOS/NPY/30ANOS/DISP_R.npy', DISP_R)
    np.save('../DADOS/NPY/30ANOS/WCS.npy', WCS)
    np.save('../DADOS/NPY/30ANOS/SCW.npy', SCW)
    np.save('../DADOS/NPY/30ANOS/WSS.npy', WSS)
else:
    DISP_WS = np.load('../DADOS/NPY/30ANOS/DISP_WS.npy')
    DISP_R = np.load('../DADOS/NPY/30ANOS/DISP_R.npy')
    WCS = np.load('../DADOS/NPY/30ANOS/WCS.npy')
    SCW = np.load('../DADOS/NPY/30ANOS/SCW.npy')
    WSS = np.load('../DADOS/NPY/30ANOS/WSS.npy')

import SAIDA as OUT
OUT.plota_mapa(DISP_WS, "Annual_DPE_Disponibility")
OUT.plota_mapa(DISP_R, "Annual_ROCI_Disponibility")
OUT.plota_mapa(WCS, "Annual_Wind_Complements_Solar")
OUT.plota_mapa(SCW, "Annual_Solar_Complements_Wind")
OUT.plota_mapa(WSS, "Annual_Wind_and_Solar_Synergy")
"""

#%%                       Carrega dados horários ERA5                         #

"""
dpe24 = {}
anoi, anof=1990, 2020
c = 0
for ano in range(anoi, anof):
    print(ano)

    netcdf_v = xr.open_dataset("../DADOS/VENTO/ERA_UV_"+str(ano)+".nc")
    u = np.array(netcdf_v.variables['u10'][:,:,:])
    v = np.array(netcdf_v.variables['v10'][:,:,:])
    ws = np.sqrt(u*u+v*v)
    ws100 = ws * (ln(100/0.0002)/ln(10/0.0002))
    dpe =  0.5 * 1.225 * 1 * (ws100**3)
    tt = netcdf_v['time'].data
  
    ttotal = len(tt)
    for t in range(0,ttotal):
        hora = int(str(tt[t]).split("T")[1][:2])
        if int(hora) not in dpe24:
            dpe24[hora] = dpe[t]
        else:
            dpe24[hora] = dpe24[hora] + dpe[t]
        if hora == 0:
            c+=1
import SAIDA as OUT
for h in range(0, 24):    
    np.save("/home/administrador/Documentos/DISSERTACAO/SAIDAS/dpe_24h/mh_"+str(h)+".npy", dpe24[h]/c)
    OUT.plota_mapa(dpe24[h]/c, 'Media_Horaria_de_DPE_('+str(h)+'H)')
"""
"""
rad24 = {}
anoi, anof=1990, 2020
c = 0
for ano in range(anoi, anof):
    print(ano)

    netcdf_v = xr.open_dataset("../DADOS/RAD/ERA5_RAD_"+str(ano)+".nc")
    r = np.array(netcdf_v.variables['ssrd'][:,:,:]/3600)
    tt = netcdf_v['time'].data
  
    ttotal = len(tt)
    for t in range(0,ttotal):
        hora = int(str(tt[t]).split("T")[1][:2])
        if int(hora) not in rad24:
            rad24[hora] = r[t]
        else:
            rad24[hora] = rad24[hora] + r[t]
        if hora == 0:
            c+=1
AA = []
import SAIDA as OUT
for h in range(0, 24):
    #np.save("/home/administrador/Documentos/DISSERTACAO/SAIDAS/roci_24h/mh_"+str(h)+".npy", rad24[h]/c)
    AA.append(rad24[h]/c)
    campo = rad24[h]/c
    campo[campo < 0] = 0
    OUT.plota_mapa(campo, 'Media_Horaria_de_ROCI_('+str(h)+'H)')
"""

#%%                       Carrega dados horários ERA5                         #

"""
c = 0
anoi, anof=1990, 2020
for ano in range(anoi, anof):
    print(ano)
    netcdf_v = xr.open_dataset("../DADOS/VENTO/ERA_UV_"+str(ano)+".nc")
    u = np.array(netcdf_v.variables['u10'][:,:,:])
    v = np.array(netcdf_v.variables['v10'][:,:,:])
    ws = np.sqrt(u*u+v*v)
    
    wsm = np.mean(ws, axis=0)
    um = np.mean(u, axis=0)
    vm = np.mean(v, axis=0)

    if c==0:
        soma_ws = wsm
        soma_u = um
        soma_v = vm
    else:
        soma_ws = soma_ws + wsm
        soma_u = soma_u + um
        soma_v = soma_v + vm
    c = c+1

ws_medio_geral = soma_ws/c
u_medio_geral = soma_u/c
v_medio_geral = soma_v/c
"""
"""
c = 0
anoi, anof=1990, 2020
for ano in range(anoi, anof):
    print(ano)

    netcdf_r = xr.open_dataset("../DADOS/RAD/ERA5_RAD_"+str(ano)+".nc")
    tt = netcdf_r['time'].data
    ttotal = len(tt)
    
    for t in range(0,ttotal):
        hora = int(str(tt[t]).split("T")[1][:2])
        if int(hora)>=10 and int(hora)<=20:
            r = np.array(netcdf_r.variables['ssrd'][t,:,:]/3600)

            if c==0:
                soma_r = r
            else:
                soma_r = soma_r + r
            c = c+1

ROCI_medio_geral = soma_r/c
"""

#%%                        Plotagem da Área de Estudo                         #

# ws100_medio_geral = np.load(PATH_DADOS + "NPY/30ANOS/ws100_medio_geral.npy")
# DPE_medio_geral = np.load(PATH_DADOS + "NPY/30ANOS/DPE_medio_geral.npy")
# ROCI_medio_geral = np.load(PATH_DADOS + "NPY/30ANOS/ROCI_medio_geral.npy")

# import SAIDA as OUT
# WS7 = ws100_medio_geral.copy()
# WS7[ws100_medio_geral < 7] = -1
# ROCI113 = ROCI_medio_geral.copy()
# ROCI113[ROCI_medio_geral < 113] = -1

# import SAIDA as OUT
# OUT.plota_mapa_tecnico_regiao(WS7, ROCI113)

# ------------------------------------------------------------------------- #
#                          Contagem de DPE na ZEE                           #
# ------------------------------------------------------------------------- #

ZEE_pontos = pd.read_csv("ZEE_pontos/ZEE1000_TEC.csv")
NNE_pontos = pd.read_csv("ZEE_pontos/NNE1000_TEC.csv")
SSE_pontos = pd.read_csv("ZEE_pontos/SSE1000_TEC.csv")

# anoi, anof=1990, 2020
# arqs = ["../DADOS/VENTO/ERA_UV_"+str(ano)+".nc" for ano in range(anoi, anof)]
# netcdf_v = xr.open_mfdataset(arqs)

# # -- NNE
# print('NNE')
# NNE_MED_WS100 = np.ones((177,125)) * np.nan
# NNE_STD_WS100 = np.ones((177,125)) * np.nan
# for ii in range(0, NNE_pontos.shape[0]):
#     lat, lon = NNE_pontos['lat'].iloc[ii], NNE_pontos['lon'].iloc[ii]
#     j, i = list(lats).index(lat), list(lons).index(lon)
    
#     u = np.array(netcdf_v.variables['u10'][:,j,i])
#     v = np.array(netcdf_v.variables['v10'][:,j,i])
#     ws = np.sqrt(u*u+v*v)
#     ws100 = ws * (ln(100/0.0002)/ln(10/0.0002))

#     NNE_STD_WS100[j, i] = np.nanstd(ws100)
#     NNE_MED_WS100[j, i] = np.nanmean(ws100)

# # -- SSE
# print('SSE')
# SSE_MED_WS100 = np.ones((177,125)) * np.nan
# SSE_STD_WS100 = np.ones((177,125)) * np.nan
# for ii in range(0, SSE_pontos.shape[0]):
#     lat, lon = SSE_pontos['lat'].iloc[ii], SSE_pontos['lon'].iloc[ii]
#     j, i = list(lats).index(lat), list(lons).index(lon)
    
#     u = np.array(netcdf_v.variables['u10'][:,j,i])
#     v = np.array(netcdf_v.variables['v10'][:,j,i])
#     ws = np.sqrt(u*u+v*v)
#     ws100 = ws * (ln(100/0.0002)/ln(10/0.0002))

#     SSE_MED_WS100[j, i] = np.nanmean(ws100)
#     SSE_STD_WS100[j, i] = np.nanstd(ws100)

NNE_MED_WS100 = np.load(PATH_DADOS + "K_C/NNE_MED_WS100.npy")
NNE_STD_WS100 = np.load(PATH_DADOS + "K_C/NNE_STD_WS100.npy")
SSE_MED_WS100 = np.load(PATH_DADOS + "K_C/SSE_MED_WS100.npy")
SSE_STD_WS100 = np.load(PATH_DADOS + "K_C/SSE_STD_WS100.npy")

#%% Numero de turbinas

nTurbinas = 169/1000
rated_speed = 10.88
power_rating = 15
dados = pd.read_csv('../DADOS/power15M.csv', sep=";")
Pt = (dict(zip(dados['Wind (m/s)'], dados['Power (MW)'])))

NNE_k = (NNE_STD_WS100/NNE_MED_WS100)**-1.086
NNE_c = np.ones(NNE_k.shape) * np.nan
for j in range(0, NNE_c.shape[0]):
    for i in range(0, NNE_c.shape[1]):
        NNE_c[j,i] = float(NNE_MED_WS100[j,i]/math.gamma(1+(1/NNE_k[j,i])))

NNE_TAEP = np.ones(NNE_k.shape) * np.nan
for j in range(0, NNE_c.shape[0]):
    for i in range(0, NNE_c.shape[1]):
        NNE_PE = [Pt[V] * weib(V, NNE_k[j,i], NNE_c[j,i]) * 8760/1000 for V in np.arange(0, 28)]
        NNE_TAEP[j,i] = np.nansum(NNE_PE) * nTurbinas
NNE_AE = np.ones(NNE_k.shape) * nTurbinas * power_rating * 8760/1000 #MWh
NNE_CF = NNE_TAEP/NNE_AE * 100

for ii in range(0, NNE_pontos.shape[0]):
    lat, lon = NNE_pontos['lat'].iloc[ii], NNE_pontos['lon'].iloc[ii]
    j, i = list(lats).index(lat), list(lons).index(lon)    
    NNE_AE[NNE_TAEP==0] = np.nan
    NNE_TAEP[NNE_TAEP==0] = np.nan
    NNE_CF[NNE_CF==0] = np.nan



SSE_k = (SSE_STD_WS100/SSE_MED_WS100)**-1.086
SSE_c = np.ones(SSE_k.shape) * np.nan
for j in range(0, SSE_c.shape[0]):
    for i in range(0, SSE_c.shape[1]):
        SSE_c[j,i] = float(SSE_MED_WS100[j,i]/math.gamma(1+(1/SSE_k[j,i])))

SSE_TAEP = np.ones(SSE_k.shape) * np.nan
for j in range(0, SSE_c.shape[0]):
    for i in range(0, SSE_c.shape[1]):
        SSE_PE = [Pt[V] * weib(V, SSE_k[j,i], SSE_c[j,i]) * 8760/1000 for V in np.arange(0, 28)]
        SSE_TAEP[j,i] = np.nansum(SSE_PE) * nTurbinas
SSE_AE = np.ones(SSE_k.shape) * nTurbinas * power_rating * 8760/1000 #MWh
SSE_CF = SSE_TAEP/SSE_AE * 100

for ii in range(0, SSE_pontos.shape[0]):
    lat, lon = SSE_pontos['lat'].iloc[ii], SSE_pontos['lon'].iloc[ii]
    j, i = list(lats).index(lat), list(lons).index(lon)    
    SSE_AE[SSE_TAEP==0] = np.nan
    SSE_TAEP[SSE_TAEP==0] = np.nan
    SSE_CF[SSE_CF==0] = np.nan

print("Area", len(NNE_pontos)*961)
print("N Turbines", len(NNE_pontos)*169)
print("NNE_AE: {:6.2f} GWh".format(np.nansum(NNE_AE)), "{:6.2f} TWh".format(np.nansum(NNE_AE)/1000), "TWh")
print("NNE_TAEP: {:6.2f} GWh".format(np.nansum(NNE_TAEP)), "{:6.2f} TWh".format(np.nansum(NNE_TAEP)/1000), "TWh")
print("NNE_CF: {:6.2f}%".format(np.nanmean(NNE_CF)))
print("")
print("Area", len(SSE_pontos)*961)
print("N Turbines", len(SSE_pontos)*169)
print("SSE_AE: {:6.2f} GWh".format(np.nansum(SSE_AE)), "{:6.2f} TWh".format(np.nansum(SSE_AE)/1000), "TWh")
print("SSE_TAEP: {:6.2f} GWh".format(np.nansum(SSE_TAEP)), "{:6.2f} TWh".format(np.nansum(SSE_TAEP)/1000), "TWh")
print("SSE_CF: {:6.2f}%".format(np.nanmean(SSE_CF)))

#%%
print('Plotando')
OUT.plota_mapa2(NNE_c, 'NNE_-_Scale Parameter')
OUT.plota_mapa2(NNE_k, 'NNE_-_Shape Parameter')
OUT.plota_mapa2(NNE_TAEP, 'NNE_-_Total Annual Energy Production')
OUT.plota_mapa2(NNE_CF, 'NNE_-_Factor Capacity')

OUT.plota_mapa2(SSE_c, 'SSE_-_Scale Parameter')
OUT.plota_mapa2(SSE_k, 'SSE_-_Shape Parameter')
OUT.plota_mapa2(SSE_TAEP, 'SSE_-_Total Annual Energy Production')
OUT.plota_mapa2(SSE_CF, 'SSE_-_Factor Capacity')

#%% Seila
# with open("./latlon.csv", 'w') as arquivo:
#     arquivo.writelines('lat, lon'+'\n')
#     for lon, lat in pontos_dentro_ZEE:
#         arquivo.writelines(str(lat)+','+str(lon)+'\n')
#     arquivo.close()

# #%%
# print("1")
# # Carrega o shape e transforma ele em um dataframe, 
# # e pega os valores de lat e lon que formam o shape
# SHP_ZEE = shp.Reader("shps/ZEE/ZEE.shp")
# coords = [s.points for s in SHP_ZEE.shapes()]
# #Cada linha do dataframe tem uma celula de coords,
# #que são todos os pontos limites do shape
# i = 0  #indice do poligono que pode ser um bairro ou uma divisão da bacia
# polygon = Polygon(coords[i])
# pontos_dentro_ZEE = []
# polygon = Polygon(coords[i])
# for lat in lats:
#     for lon in lons:
#         p = Point(lon, lat)
#         if p.within(polygon):
#             pontos_dentro_ZEE.append((lon, lat))

# print("2")
# SHP_BRASIL = shp.Reader("shps/Brasil/estados_2010.shp")
# coords = [s.points for s in SHP_BRASIL.shapes()]
# pontos_dentro_Brasil = []
# for i in range(0, 27):
#     print(i)
#     polygon = Polygon(coords[i])
#     for lat in lats:
#         for lon in lons:
#             p = Point(lon, lat)
#             if p.within(polygon):
#                 pontos_dentro_Brasil.append((lon, lat))

# print("3")
# ZEE_pontos_plotar = []
# for ponto in pontos_dentro_ZEE:
#     if ponto not in pontos_dentro_Brasil:
#         ZEE_pontos_plotar.append(ponto)

# ROCI_medio_geral = np.load(PATH_DADOS + "NPY/30ANOS/ROCI_medio_geral.npy")
# numero_de_pontos = 0
# masked_ROCI_medio_geral = np.zeros((ROCI_medio_geral.shape))
# for lon, lat in ZEE_pontos_plotar:
#     j = list(lats).index(lat)
#     i = list(lons).index(lon)
#     masked_ROCI_medio_geral[j, i] = ROCI_medio_geral[j, i].copy()
#     numero_de_pontos += 1

# # -- NNE
# LATi, LONi=1, -45
# LATf, LONf=-8, -34
# NNE_pontos_plotar = []
# for lon, lat in ZEE_pontos_plotar:
#     if lat<=LATi and lat>=LATf and lon>=LONi and lon<=LONf:
#         NNE_pontos_plotar.append((lon, lat))

# # -- SSE
# LATi, LONi=-20, -54
# LATf, LONf=-35, -35
# SSE_pontos_plotar = []
# for lon, lat in ZEE_pontos_plotar:
#     if lat<=LATi and lat>=LATf and lon>=LONi and lon<=LONf:
#         SSE_pontos_plotar.append((lon, lat))


# DPE_medio_geral = np.load(PATH_DADOS + "NPY/30ANOS/DPE_medio_geral.npy")

# numero_de_pontos = 0
# ZEE_masked_DPE_medio_geral = np.zeros((DPE_medio_geral.shape))
# ZEE_masked_ROCI_medio_geral = np.zeros((ROCI_medio_geral.shape))
# for lon, lat in ZEE_pontos_plotar:
#     j = list(lats).index(lat)
#     i = list(lons).index(lon)
#     ZEE_masked_DPE_medio_geral[j, i] = DPE_medio_geral[j, i].copy()
#     ZEE_masked_ROCI_medio_geral[j, i] = ROCI_medio_geral[j, i].copy()
#     numero_de_pontos += 1
# ZEE_soma_dpe = np.sum(ZEE_masked_DPE_medio_geral)
# ZEE_media_dpe = round(ZEE_soma_dpe/numero_de_pontos, 2)
# ZEE_soma_dpe = round(ZEE_soma_dpe*0.000001, 2)
# ZEE_soma_dps = np.sum(ZEE_masked_ROCI_medio_geral)
# ZEE_media_dps = round(ZEE_soma_dps/numero_de_pontos, 2)
# ZEE_soma_dps = round(ZEE_soma_dps*0.000001, 2)

# numero_de_pontos = 0
# NNE_masked_DPE_medio_geral = np.zeros((DPE_medio_geral.shape))
# NNE_masked_ROCI_medio_geral = np.zeros((ROCI_medio_geral.shape))
# for lon, lat in NNE_pontos_plotar:
#     j = list(lats).index(lat)
#     i = list(lons).index(lon)
#     NNE_masked_DPE_medio_geral[j, i] = DPE_medio_geral[j, i].copy()
#     NNE_masked_ROCI_medio_geral[j, i] = ROCI_medio_geral[j, i].copy()
#     numero_de_pontos += 1
# NNE_soma_dpe = np.sum(NNE_masked_DPE_medio_geral)
# NNE_media_dpe = round(NNE_soma_dpe/numero_de_pontos, 2)
# NNE_soma_dpe = round(NNE_soma_dpe*0.000001, 2)
# NNE_soma_dps = np.sum(NNE_masked_ROCI_medio_geral)
# NNE_media_dps = round(NNE_soma_dps/numero_de_pontos, 2)
# NNE_soma_dps = round(NNE_soma_dps*0.000001, 2)

# numero_de_pontos = 0
# SSE_masked_DPE_medio_geral = np.zeros((DPE_medio_geral.shape))
# SSE_masked_ROCI_medio_geral = np.zeros((ROCI_medio_geral.shape))
# for lon, lat in SSE_pontos_plotar:
#     j = list(lats).index(lat)
#     i = list(lons).index(lon)
#     SSE_masked_DPE_medio_geral[j, i] = DPE_medio_geral[j, i].copy()
#     SSE_masked_ROCI_medio_geral[j, i] = ROCI_medio_geral[j, i].copy()
#     numero_de_pontos += 1

# SSE_soma_dpe = np.sum(SSE_masked_DPE_medio_geral)
# SSE_media_dpe = round(SSE_soma_dpe/numero_de_pontos, 2)
# SSE_soma_dpe = round(SSE_soma_dpe*0.000001, 2)
# SSE_soma_dps = np.sum(SSE_masked_ROCI_medio_geral)
# SSE_media_dps = round(SSE_soma_dps/numero_de_pontos, 2)
# SSE_soma_dps = round(SSE_soma_dps*0.000001, 2)

# import SAIDA as OUT
# LATi, LONi=7, -54
# LATf, LONf=-35, -23
# bounds = [(LONi, LONf, LATi, LATf)]
# OUT.plota_mapa(ZEE_masked_DPE_medio_geral, "ZEE_Annual_Average_WPD", bounds)
# OUT.plota_mapa(ZEE_masked_ROCI_medio_geral, "ZEE_Annual_Average_SPD", bounds)

# LATi, LONi=1, -45
# LATf, LONf=-8, -34
# bounds = [(LONi, LONf, LATi, LATf)]
# OUT.plota_mapa(NNE_masked_DPE_medio_geral, "NNE_Annual_Average_WPD", bounds)
# OUT.plota_mapa(NNE_masked_ROCI_medio_geral, "NNE_Annual_Average_SPD", bounds)

# LATi, LONi=-20, -54
# LATf, LONf=-35, -35
# bounds = [(LONi, LONf, LATi, LATf)]
# OUT.plota_mapa(SSE_masked_DPE_medio_geral, "SSE_Annual_Average_WPD", bounds)
# OUT.plota_mapa(SSE_masked_ROCI_medio_geral, "SSE_Annual_Average_SPD", bounds)




# ws100_medio_geral = np.load(PATH_DADOS + "NPY/30ANOS/ws100_medio_geral.npy")
# WS7 = np.zeros(ws100_medio_geral.shape)
# WS7[ws100_medio_geral<7] = 1

# ncout = Dataset('/home/administrador/Documentos/DISSERTACAO/SCRIPT/myfile.nc','w','NETCDF3'); # using netCDF3 for output format 
# ncout.createDimension('lon',125);
# ncout.createDimension('lat',177);
# lonvar = ncout.createVariable('lon','float32',('lon'));
# lonvar[:] = lons;
# latvar = ncout.createVariable('lat','float32',('lat'));
# latvar[:] = lats;
# myvar = ncout.createVariable('ws100','float32',('lat','lon'));
# myvar.setncattr('units','m/s');
# myvar[:] = WS7;

print ('- finished! Tempo:', datetime.now() - start)