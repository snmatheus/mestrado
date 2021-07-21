import SAIDA as OUT
import numpy as np
import xarray as xr

lats = OUT.lats
lons = OUT.lons
ym, xm = lats.shape[0], lons.shape[0]

def ln(x):
    from math import log
    return log(x)

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
    #np.save("/home/administrador/Documentos/DISSERTACAO/SAIDAS/dpe_24h/mh_"+str(h)+".npy", dpe24[h]/c)
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

c = 0
anoi, anof=1990, 2020
for ano in range(anoi, anof):
    print(ano)

    netcdf_r = xr.open_dataset("../DADOS/RAD/ERA5_RAD_"+str(ano)+".nc")
    r = np.array(netcdf_r.variables['ssrd'][:,:,:]/3600)    
    rm = np.mean(r, axis=0)

    if c==0:
        soma_r = rm
    else:
        soma_r = soma_r + rm
    c = c+1

ROCI_medio_geral = soma_r/c
# ------------------------------------------------------------------------- #
#                              Plotagem Anual                               #
# ------------------------------------------------------------------------- #

# ws100_medio_geral = ws_medio_geral * (ln(100/0.0002)/ln(10/0.0002))
# DPE_medio_geral = 0.5 * 1.225 * 1 * (ws100_medio_geral**3)

# DPE_medio_geral = np.load("/home/administrador/Documentos/DISSERTACAO/DADOS/NPY/DPE_medio_geral.npy")
# ROCI_medio_geral = np.load("/home/administrador/Documentos/DISSERTACAO/DADOS/NPY/ROCI_medio_geral.npy")
import SAIDA as OUT
# OUT.plota_mapa([u_medio_geral, v_medio_geral, ws_medio_geral], 2)
# OUT.plota_mapa([u_medio_geral, v_medio_geral, ws100_medio_geral], 3)
# OUT.plota_mapa(DPE_medio_geral, 4)
OUT.plota_mapa(ROCI_medio_geral, 5)


# ------------------------------------------------------------------------- #
#                        Número de horas com RPD                            #
# ------------------------------------------------------------------------- #

# NHRPD_WS = np.zeros((ym, xm))
# NHRPD_R = np.zeros((ym, xm))

# NHRPD_WS_NR = np.zeros((ym, xm))
# NHRPD_NWS_R = np.zeros((ym, xm))
# NHRPD_WS_R = np.zeros((ym, xm))

# ttotal = 0
# anoi, anof=1990, 2020
# for ano in range(anoi, anof):
#     print(ano)

#     netcdf_r = xr.open_dataset("../DADOS/RAD/ERA5_RAD_"+str(ano)+".nc")
    # r = np.array(netcdf_r.variables['ssrd'][:,:,:]/3600)

    # netcdf_v = xr.open_dataset("../DADOS/VENTO/ERA_UV_"+str(ano)+".nc")
    # u = np.array(netcdf_v.variables['u10'][:,:,:])
    # v = np.array(netcdf_v.variables['v10'][:,:,:])
    # ws = np.sqrt(u*u+v*v)
    
    # for t in range(u.shape[0]):
    #     WS, R = ws[t], r[t]
    
    #     RPD_WS = np.zeros((ym, xm))
    #     RPD_R = np.zeros((ym, xm))
    #     RPD_WS[WS >= 7] = 1
    #     RPD_R[R >= 170] = 1
    #     NHRPD_WS = NHRPD_WS + RPD_WS
    #     NHRPD_R = NHRPD_R + RPD_R
        
    #     RPD_WS_NR = np.zeros((ym, xm))
    #     RPD_NWS_R = np.zeros((ym, xm))
    #     RPD_WS_R = np.zeros((ym, xm))    
    #     RPD_WS_NR[(WS >= 7) & (R < 170)] = 1
    #     RPD_NWS_R[(WS < 7) & (R >= 170)] = 1        
    #     RPD_WS_R[(((WS >= 7) & (R < 170)) | ((WS < 7) & (R >= 170)))] = 1
    #     NHRPD_WS_NR = NHRPD_WS_NR + RPD_WS_NR
    #     NHRPD_NWS_R = NHRPD_NWS_R + RPD_NWS_R
    #     NHRPD_WS_R = NHRPD_WS_R + RPD_WS_R        
    # ttotal_WS += u.shape[0]
    # ttotal_R += r.shape[0]
#     ttotal += netcdf_r.variables['ssrd'].shape[0]

# NHRPD_WS = np.load('NHRPD_WS.npy')
# NHRPD_R = np.load('NHRPD_R.npy')
# NHRPD_WS_NR = np.load('NHRPD_WS_NR.npy')
# NHRPD_NWS_R = np.load('NHRPD_NWS_R.npy')
# NHRPD_WS_R = np.load('NHRPD_WS_R.npy')

# Availability_WS = (NHRPD_WS/ttotal)*100
# Availability_R = (NHRPD_R/ttotal)*100

# WCS = (NHRPD_WS_NR/ttotal)*100
# SCW = (NHRPD_NWS_R/ttotal)*100
# WSS = (NHRPD_WS_R/ttotal)*100

# ------------------------------------------------------------------------- #
#                           Plotagem Sinergia                               #
# ------------------------------------------------------------------------- #
"""
import SAIDA as OUT
OUT.plota_mapa(Availability_WS, "Avaliabilidade_de_DPE")
OUT.plota_mapa(Availability_R, "Avaliabilidade_de_ROCI")
OUT.plota_mapa(WCS, "Wind_Complements_Solar")
OUT.plota_mapa(SCW, "Solar_Complements_Wind")
OUT.plota_mapa(WSS, "Wind_and_Solar_Synergy")
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
"""
import SAIDA as OUT
OUT.plota_mapa('', -9, 1)
"""
# ------------------------------------------------------------------------- #
#                           Plotagem de Satélite                            #
# ------------------------------------------------------------------------- #
"""
import SAIDA as OUT
QS_nc=xr.open_dataset(PATH_HD+"DISSERTACAO/DADOS/QUIKSCAT/WIND_L3-QUIKSCAT_SEAWINDS_25_ASC.nc")
QS_wd=QS_nc["wind_to_dir"][:]
QS_ws=QS_nc["wind_speed"][:]


QS_nc=xr.open_dataset(PATH_HD+"DISSERTACAO/DADOS/QUIKSCAT/KNMI-GLO-WIND_L3-REP-OBS_QUIKSCAT_SEAWINDS_25_DES_1595484912636.nc")
QS_time=[str(dt) for dt in QS_nc["time"].data]
QS_lats=QS_nc["lat"][:].data
QS_lons=QS_nc["lon"][:].data-360

for i in range(1, len(QS_time)):
    t=QS_time[i]        
    QS_u=QS_nc["eastward_wind"][i].data
    QS_v=QS_nc["northward_wind"][i].data
    QS_vu=(QS_v**2+QS_u**2)**0.5
    plota_sat(t, QS_vu, QS_u, QS_v, QS_lats, QS_lons)
    i=i+1
"""