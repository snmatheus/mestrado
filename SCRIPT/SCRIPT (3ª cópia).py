import SAIDA as OUT
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"]=[10,4]

lats = OUT.lats
lons = OUT.lons

ym, xm = lats.shape[0], lons.shape[0]

def ln(x):
    from math import log
    return log(x)

# ------------------------------------------------------------------------- #
#                         Histograma de radiação                            #
# ------------------------------------------------------------------------- #
"""
medias = {}
frequencias = {}
anoi, anof=1990, 2020
anoi, anof=1990, 2020
for ano in range(anoi, anof):
    print(ano)
    if 1==2:
        netcdf_r = xr.open_dataset("../DADOS/RAD/ERA5_RAD_"+str(ano)+".nc")
           
        tt = netcdf_r['time'].data
        ttotal = len(tt)
        
        c = 0
        mes_ant='-'
        for t in range(0,ttotal):
            mes = str(tt[t]).split("-")[1]
            r = np.array(netcdf_r.variables['ssrd'][t,:,:]/3600)
            if c==0:
                soma_r = r
            else:
                soma_r = soma_r + r
            c = c+1

            if mes != mes_ant:
                #print(mes)
                medias[str(ano)+mes] = np.reshape(soma_r/c, soma_r.shape[0]*soma_r.shape[1])
                np.save('../DADOS/NPY/DIST/MEDIA_RAD_DIA_'+str(ano)+mes+'.npy', medias[str(ano)+mes])
                c = 0
            mes_ant = mes

    else:
        for mes in ['01','02','03','04','05','06','07','08','09','10','11','12']:
            medias[str(ano)+mes] = np.load('../DADOS/NPY/DIST/MEDIA_RAD_DIA_'+str(ano)+mes+'.npy')

radiacoes = []
for campo in medias.values():
    for rad in campo:
        if rad>10:
            radiacoes.append(rad)
#print(radiacoes[0],radiacoes[1])
perc = np.percentile(radiacoes, 5)

bins = np.arange(0, np.max(radiacoes), np.max(radiacoes)/100)

plt.xlim([min(radiacoes)-5, max(radiacoes)+5])
plt.ylim([0, 250000])

plt.hist(radiacoes, bins=bins, color='#0504aa', alpha=0.7, rwidth=0.85)
plt.plot([perc, perc], [0, 250000], linewidth=1.5, c='magenta')

plt.xlabel('Radiação (W/m²)')
plt.ylabel('Ocorrências')
ax = plt.gca()

text = 'Média = '+str(round(np.mean(radiacoes),1))+'\nMáximo = '+str(round(np.max(radiacoes),1))+'\nMínimo = '+str(round(np.min(radiacoes),1))+'\nPercentil 5% = '+str(round(perc))
plt.text(0.95, 0.95, text, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
plt.savefig('hist_rad.png')
"""
# ------------------------------------------------------------------------- #
#                              Plotagem Anual                               #
# ------------------------------------------------------------------------- #
"""
ws_medio_geral = np.load("/home/administrador/Documentos/DISSERTACAO/DADOS/NPY/30ANOS/ws_medio_geral.npy")
ws100_medio_geral = np.load("/home/administrador/Documentos/DISSERTACAO/DADOS/NPY/30ANOS/ws100_medio_geral.npy")
DPE_medio_geral = np.load("/home/administrador/Documentos/DISSERTACAO/DADOS/NPY/30ANOS/DPE_medio_geral.npy")
ROCI_medio_geral = np.load("/home/administrador/Documentos/DISSERTACAO/DADOS/NPY/30ANOS/ROCI_medio_geral.npy")

import SAIDA as OUT
OUT.plota_mapa(ws_medio_geral, "Annual_Average_Wind_Speed_at_10_m")
OUT.plota_mapa(ws100_medio_geral, "Annual_Average_Wind_Speed_at_100_m")
OUT.plota_mapa(DPE_medio_geral, "Annual_Average_WPD")
OUT.plota_mapa(ROCI_medio_geral, "Annual_Average_SPD")
"""
# ------------------------------------------------------------------------- #
#                       Carrega dados horários ERA5                         #
# ------------------------------------------------------------------------- #

anoi, anof=1990, 1991
anoi, anof=1990, 2020
NHRPD_WS, NHRPD_R, NHRPD_WS_NR, NHRPD_NWS_R, NHRPD_WS_R = {}, {}, {}, {}, {}
limiar_r = 113
limiar_v = 7
a=1
ttotal = 0
xxx = 31*24
xxx = 0
for ano in range(anoi, anof):
    print(ano)
    netcdf_r = xr.open_dataset("../DADOS/RAD/ERA5_RAD_"+str(ano)+".nc")
    netcdf_v = xr.open_dataset("../DADOS/VENTO/ERA_UV_"+str(ano)+".nc")
    tt = netcdf_r['time'].data
    if xxx!=0:
        tam = xxx
    else:
        tam = len(tt)
    if a==1:
        r = np.array(netcdf_r.variables['ssrd'][:tam,:,:]/3600)
        u = np.array(netcdf_v.variables['u10'][:tam,:,:])
        v = np.array(netcdf_v.variables['v10'][:tam,:,:])
        ws = np.sqrt(u*u+v*v)
        ws100 = ws * (ln(100/0.0002)/ln(10/0.0002))

        tt = netcdf_v['time'].data
        #for t in range(0, len(tt)):
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
    else:
        for t in range(0, tam):
            hora = str(tt[t]).split("T")[1][:2]
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

# ------------------------------------------------------------------------- #
#                        Número de horas com RPD                            #
# ------------------------------------------------------------------------- #
""" 
limiar_r = 113
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
        if anof-anoi == 2020-1990:
            tam = len(tt)
        if anof-anoi==1991-1990:
            tam = 31*24
            tam = len(tt)
        r = np.array(netcdf_r.variables['ssrd'][:tam,:,:]/3600)
        u = np.array(netcdf_v.variables['u10'][:tam,:,:])
        v = np.array(netcdf_v.variables['v10'][:tam,:,:])
        ws = np.sqrt(u*u+v*v)
        ws100 = ws * (ln(100/0.0002)/ln(10/0.0002))

        RPD_WS = np.zeros(v.shape)
        RPD_R = np.zeros(r.shape)
        RPD_WS_NR = np.zeros(v.shape)
        RPD_NWS_R = np.zeros(v.shape)
        RPD_WS_R = np.zeros(v.shape)

        RPD_WS[ws100 >= limiar_v] = 1
        RPD_R[r >= limiar_r] = 1
        RPD_WS_NR[(ws100 >= limiar_v) & (r < limiar_r)] = 1
        RPD_NWS_R[(ws100 < limiar_v) & (r >= limiar_r)] = 1        
        RPD_WS_R[(((ws100 >= limiar_v) & (r < limiar_r)) | ((ws100 < limiar_v) & (r >= limiar_r)))] = 1

        if ttotal==0:
            NHRPD_WS = np.sum(RPD_WS, axis=0)
            NHRPD_R = np.sum(RPD_R, axis=0)
            NHRPD_WS_NR = np.sum(RPD_WS_NR, axis=0)
            NHRPD_NWS_R = np.sum(RPD_NWS_R, axis=0)
            NHRPD_WS_R = np.sum(RPD_WS_R, axis=0)
        else:
            NHRPD_WS += np.sum(RPD_WS, axis=0)
            NHRPD_R += np.sum(RPD_R, axis=0)
            NHRPD_WS_NR += np.sum(RPD_WS_NR, axis=0)
            NHRPD_NWS_R += np.sum(RPD_NWS_R, axis=0)
            NHRPD_WS_R += np.sum(RPD_WS_R, axis=0)
        ttotal += r.shape[0]

    DISP_WS = (NHRPD_WS/ttotal)*100
    DISP_R = (NHRPD_R/ttotal)*100
    WCS = (NHRPD_WS_NR/ttotal)*100
    SCW = (NHRPD_NWS_R/ttotal)*100
    WSS = (NHRPD_WS_R/ttotal)*100

    np.save('../DADOS/NPY/30ANOS/NHRPD_WS.npy', NHRPD_WS)
    np.save('../DADOS/NPY/30ANOS/NHRPD_R.npy', NHRPD_R)
    np.save('../DADOS/NPY/30ANOS/NHRPD_WS_NR.npy', NHRPD_WS_NR)
    np.save('../DADOS/NPY/30ANOS/NHRPD_NWS_R.npy', NHRPD_NWS_R)
    np.save('../DADOS/NPY/30ANOS/NHRPD_WS_R.npy', NHRPD_WS_R)
    np.save('../DADOS/NPY/30ANOS/ttotal.npy', np.array(ttotal))
    
    np.save('../DADOS/NPY/30ANOS/DISP_WPD.npy', DISP_WS)
    np.save('../DADOS/NPY/30ANOS/DISP_SPD.npy', DISP_R)
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
OUT.plota_mapa(DISP_WS, "Annual_WPD_Disponibility")
OUT.plota_mapa(DISP_R, "Annual_SPD_Disponibility")
OUT.plota_mapa(WCS, "Annual_Wind_Complements_Solar")
OUT.plota_mapa(SCW, "Annual_Solar_Complements_Wind")
OUT.plota_mapa(WSS, "Annual_Wind_and_Solar_Synergy")
"""
 
# ------------------------------------------------------------------------- #
#                        Número de horas com RPD                            #
# ------------------------------------------------------------------------- #

# limiar_r = 113
# limiar_v = 7
# ttotal = 0
# anoi, anof=1990, 1991
# # anoi, anof=1990, 2020
# a = 1 
# for ano in range(anoi, anof):
#     print(ano)
#     netcdf_r = xr.open_dataset("../DADOS/RAD/ERA5_RAD_"+str(ano)+".nc")
#     netcdf_v = xr.open_dataset("../DADOS/VENTO/ERA_UV_"+str(ano)+".nc")
#     tt = netcdf_r['time'].data
#     if anof-anoi == 2020-1990:
#         tam = len(tt)
#     if anof-anoi==1991-1990:
#         tam = 31*24
#     if a==1:
#         r = np.array(netcdf_r.variables['ssrd'][:tam,:,:]/3600)
#         u = np.array(netcdf_v.variables['u10'][:tam,:,:])
#         v = np.array(netcdf_v.variables['v10'][:tam,:,:])
#         ws = np.sqrt(u*u+v*v)
#         ws100 = ws * (ln(100/0.0002)/ln(10/0.0002))

#         NHRPD_WS = np.zeros((ym, xm))
#         NHRPD_R = np.zeros((ym, xm))
#         NHRPD_WS_NR = np.zeros((ym, xm))
#         NHRPD_NWS_R = np.zeros((ym, xm))
#         NHRPD_WS_R = np.zeros((ym, xm))

#         for t in range(u.shape[0]):
#             WS, R = ws100[t], r[t]
    
#             RPD_WS = np.zeros((ym, xm))
#             RPD_WS[WS >= limiar_v] = 1
#             NHRPD_WS = NHRPD_WS + RPD_WS

#             RPD_R = np.zeros((ym, xm))
#             RPD_R[R >= limiar_r] = 1
#             NHRPD_R = NHRPD_R + RPD_R
       
#             RPD_WS_NR = np.zeros((ym, xm))
#             RPD_WS_NR[(WS >= limiar_v) & (R < limiar_r)] = 1
#             NHRPD_WS_NR = NHRPD_WS_NR + RPD_WS_NR

#             RPD_NWS_R = np.zeros((ym, xm))
#             RPD_NWS_R[(WS < limiar_v) & (R >= limiar_r)] = 1        
#             NHRPD_NWS_R = NHRPD_NWS_R + RPD_NWS_R

#             RPD_WS_R = np.zeros((ym, xm))    
#             RPD_WS_R[(((WS >= limiar_v) & (R < limiar_r)) | ((WS < limiar_v) & (R >= limiar_r)))] = 1
#             #RPD_WS_R[(WS >= limiar_v) & (R >= limiar_r)] = 1
#             NHRPD_WS_R = NHRPD_WS_R + RPD_WS_R
    #else:
        #NHRPD_WS = np.load('../DADOS/NPY/30ANOS/NHRPD_WS.npy')
        #NHRPD_R = np.load('../DADOS/NPY/30ANOS/NHRPD_R.npy')
        #NHRPD_WS_NR = np.load('../DADOS/NPY/30ANOS/NHRPD_WS_NR.npy')
        #NHRPD_NWS_R = np.load('../DADOS/NPY/30ANOS/NHRPD_NWS_R.npy')
        #NHRPD_WS_R = np.load('../DADOS/NPY/30ANOS/NHRPD_WS_R.npy')
    # ttotal += tam

# if a==1:
#     DISP_WS = (NHRPD_WS/ttotal)*100
#     DISP_R = (NHRPD_R/ttotal)*100
#     WCS = (NHRPD_WS_NR/ttotal)*100
#     SCW = (NHRPD_NWS_R/ttotal)*100
#     WSS = (NHRPD_WS_R/ttotal)*100

#     np.save('../DADOS/NPY/30ANOS/DISP_WS.npy', DISP_WS)
#     np.save('../DADOS/NPY/30ANOS/DISP_R.npy', DISP_R)
#     np.save('../DADOS/NPY/30ANOS/WCS.npy', WCS)
#     np.save('../DADOS/NPY/30ANOS/SCW.npy', SCW)
#     np.save('../DADOS/NPY/30ANOS/WSS.npy', WSS)
# else:
#     DISP_WS = np.load('../DADOS/NPY/30ANOS/DISP_WS.npy')
#     DISP_R = np.load('../DADOS/NPY/30ANOS/DISP_R.npy')
#     WCS = np.load('../DADOS/NPY/30ANOS/WCS.npy')
#     SCW = np.load('../DADOS/NPY/30ANOS/SCW.npy')
#     WSS = np.load('../DADOS/NPY/30ANOS/WSS.npy')

# import SAIDA as OUT
# OUT.plota_mapa(DISP_WS, "Annual_DPE_Disponibility")
# OUT.plota_mapa(DISP_R, "Annual_ROCI_Disponibility")
# OUT.plota_mapa(WCS, "Annual_Wind_Complements_Solar")
# OUT.plota_mapa(SCW, "Annual_Solar_Complements_Wind")
# OUT.plota_mapa(WSS, "Annual_Wind_and_Solar_Synergy")

# ------------------------------------------------------------------------- #
#                       Carrega dados horários ERA5                         #
# ------------------------------------------------------------------------- #

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
# ------------------------------------------------------------------------- #
#                       Carrega dados horários ERA5                         #
# ------------------------------------------------------------------------- #
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
# ------------------------------------------------------------------------- #
#                        Calcula Mediana, IQR, RCOV                         #
# ------------------------------------------------------------------------- #
"""
from scipy.stats import iqr

ym, xm = lats.shape[0], lons.shape[0]

MEDIAN_DPE = np.zeros((ym, xm))
RCOV_DPE = np.zeros((ym, xm))
IQR_DPE = np.zeros((ym, xm))
del ym, xm

anoi, anof=2010, 2020
arqs = ["../DADOS/VENTO/ERA_UV_"+str(ano)+".nc" for ano in range(anoi, anof)]
netcdf_v = xr.open_mfdataset(arqs)
del anoi, anof, arqs

for lat in range(0, lats.shape[0]):
    for lon in range(0, lons.shape[0]):
        u = np.array(netcdf_v.variables['u10'][:, lat, lon])
        v = np.array(netcdf_v.variables['v10'][:, lat, lon])
        ws = np.sqrt(u*u+v*v)

        ws100 = ws * (ln(100/0.0002)/ln(10/0.0002))
        DPE = 0.5 * 1.225 * 1 * (ws100**3)

        MEDIAN_DPE[lat, lon] = np.median(ws100)
        RCOV_DPE[lat, lon] = np.std(DPE, ddof=1)/np.median(DPE)
        IQR_DPE[lat, lon] = iqr(DPE)
    print(lat,lats.shape[0],"|",lon,lons.shape[0])


import SAIDA as OUT

MEDIAN_ROCI = [0,0]
RCOV_ROCI = [0,0]
IQR_ROCI = [0,0]

lmax = round(np.max([np.max(MEDIAN_DPE), np.max(MEDIAN_ROCI)]))
OUT.plota_mapa([lmax, MEDIAN_DPE], "Mediana_de_DPE")
#OUT.plota_mapa([lmax, MEDIAN_ROCI],  "Mediana_de_ROCI")

lmax = round(np.max([np.max(RCOV_DPE), np.max(RCOV_ROCI)]))
OUT.plota_mapa([2, RCOV_DPE], "RCOV_DPE")
#OUT.plota_mapa([1.1, RCOV_ROCI], "RCOV_ROCI")

lmax = round(np.max([np.max(IQR_DPE), np.max(IQR_ROCI)]))
OUT.plota_mapa([lmax, IQR_DPE], "IQR_DPE")
#OUT.plota_mapa([lmax, IQR_ROCI], "IQR_ROCI")
"""
# ------------------------------------------------------------------------- #
#                        Plotagem da Área de Estudo                         #
# ------------------------------------------------------------------------- #

# ws100_medio_geral = np.load("/home/administrador/Documentos/DISSERTACAO/DADOS/NPY/30ANOS/ws100_medio_geral.npy")
# DPE_medio_geral = np.load("/home/administrador/Documentos/DISSERTACAO/DADOS/NPY/30ANOS/DPE_medio_geral.npy")
# ROCI_medio_geral = np.load("/home/administrador/Documentos/DISSERTACAO/DADOS/NPY/30ANOS/ROCI_medio_geral.npy")

# import SAIDA as OUT
# WS7 = ws100_medio_geral.copy()
# WS7[ws100_medio_geral < 7] = -1
# ROCI113 = ROCI_medio_geral.copy()
# ROCI113[ROCI_medio_geral < 113] = -1

# import SAIDA as OUT
# OUT.plota_mapa_tecnico_regiao(WS7, ROCI113)

