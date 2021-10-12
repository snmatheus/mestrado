import datetime
import numpy as np
from netCDF4 import Dataset, num2date
import itertools
#import xarray as xr

def plota(lat, lon, u, v,intensi, mes):
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    from cartopy.io.shapereader import Reader

    meses=["Jan","Fev","Mar","Abr","Mai","Jun","Jul","Aug","Set","Out","Nov","Dec",]
    plt.rcParams["figure.figsize"]=[6,6]

    # Coordinates to limit map areas
    LATi, LONi=26, -76
    LATf, LONf=-37, -11
    
    LATi, LONi=6, -60
    LATf, LONf=-37, -30
    bounds = [(min(lon), max(lon), min(lat), max(lat))] 
    
    fig = plt.figure(figsize=(17., 12.))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(*bounds, crs=ccrs.PlateCarree())

    cf = ax.contourf(lon, lat, intensi[:, :], levels=np.arange(0,50,0.25), cmap="viridis", transform=ccrs.PlateCarree()) #plt.cm.gist_earth_r
    cb = fig.colorbar(cf, orientation='horizontal', extend='max', aspect=65, shrink=0.5, pad=0.05, extendrect='True')
    cb.set_label('Intensidade do vento', size='x-large')
    
#    ax.streamplot(lon, lat, u, v, density=5, linewidth=1, color='k', transform=ccrs.PlateCarree())

#    def plota_ponto(lon, lat, texto):
#        ax.plot(lon,lat, '*', markerfacecolor='r', markeredgecolor='r', markersize=10)
#        ax.text(lon-0.8,lat+1, texto, color='r')
#
#    plota_ponto(-49.00, -33.00, 'S1') #sul
#    plota_ponto(-46.00, -28.00, 'S2') #sul
#    plota_ponto(-42.50, -24.50, 'SE1') #sudeste
#    plota_ponto(-38.25, -21.25, 'SE2') #sudeste
#    plota_ponto(-35.50, -12.00, 'NE1') #Nordeste
#    plota_ponto(-36.50, -2.00,  'NE2')  #Nordeste
#    plota_ponto(-46.00,  2.50,  'NE3')  #Nordeste

    ax.add_geometries(Reader("../../SHAPES/Brasil/estados_2010.shp").geometries(), ccrs.PlateCarree(), facecolor="silver", edgecolor="black")
    ax.add_geometries(Reader("../../SHAPES/AS/america_do_sul.shp").geometries(), ccrs.PlateCarree(), facecolor="silver", edgecolor="black")
    ax.set_title('Intensidade (m/s) e direção do vento ('+meses[mes-1]+')',loc='left')
    
    fig.tight_layout() 
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlines = False
    gl.ylines = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylabel_style = {'rotation':'vertical'}
    plt.savefig("../SAIDAS/TESTE.png", dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()

def plota2(lat, lon, var):
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    from cartopy.io.shapereader import Reader

    plt.rcParams["figure.figsize"]=[6,6]    
    # Set up our projection
    crs = ccrs.PlateCarree()    
    # Coordinates to limit map areas
    LATi, LONi=26, -76
    LATf, LONf=-37, -11
    
    LATi, LONi=6, -60
    LATf, LONf=-37, -30
    bounds = [(LONi, LONf, LATi, LATf)] 
    
    fig = plt.figure(figsize=(17., 12.))
    ax = fig.add_subplot(1, 1, 1, projection=crs)
    ax.set_extent(*bounds, crs=ccrs.PlateCarree())

    cf = ax.contourf(lon, lat, var[:, :], cmap="viridis", transform=ccrs.PlateCarree()) # levels=np.arange(0,18,0.25), 
    cb = fig.colorbar(cf, orientation='horizontal', extend='max', aspect=65, shrink=0.5, pad=0.05, extendrect='True')
    cb.set_label('Intensidade do vento', size='x-large')
    
    #sp=ax.streamplot(lon, lat, var, density=5, linewidth=1, color='k', transform=ccrs.PlateCarree())

    ax.add_geometries(Reader("../../Shapes/Brasil/estados_2010.shp").geometries(), ccrs.PlateCarree(), facecolor="silver", edgecolor="black")
    #ax.add_geometries(Reader("../../Shapes/Oceano/Oceano_br2.shp").geometries(), ccrs.PlateCarree(), linewidth=0.2, facecolor="paleturquoise", edgecolor="none")
    ax.add_geometries(Reader("../../Shapes/AS/america_do_sul.shp").geometries(), ccrs.PlateCarree(), facecolor="silver", edgecolor="black")
    ax.set_title('Intensidade (m/s) e direção do vento',loc='left')
    
    fig.tight_layout() 
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlines = False
    gl.ylines = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylabel_style = {'rotation':'vertical'}
    plt.show()

def plota_linear_pontos(ws, ponto):
    import matplotlib.pyplot as plt
    global datevar

    plt.rcParams["figure.figsize"]=[10,3]
    fig, ax = plt.subplots()
    fig.tight_layout()
    
    ax.set_title(ponto,loc='left')
    ax.plot(datevar, ws, label='Velocidade do Vento', color='blue')
    ax.fill_between(datevar, ws, 0, facecolor='blue', alpha=0.5)
    ax.set_xlim(datevar[0], datevar[-1])
    plot_range = [0, max(ws)+1, 1]
    ax.set_ylabel('Velocidade do Vento(m/s)', multialignment='center')
    ax.set_ylim(plot_range[0], plot_range[1], plot_range[2])    
    ax.grid(b=True, which='major', axis='y', color='k', linestyle='--', linewidth=0.5)    
    plt.show()

def nc_media_mensal(nc, campo, datevar, soma, cont):
    for idx, data in enumerate(datevar):
        var=nc.variables[campo][idx,:,:]
        for mes in range(1,13):
            if data.month==mes:
                if data.month==1 and data.day==1 and data.hour==0 and data.minute==0:
                    print("Janeiro de",data.year)
                if soma[mes]=="":
                    soma[mes]=var
                else:
                    soma[mes]=soma[mes]+var
                cont[mes]=cont[mes]+1
                break
    return soma, cont


#t_unit = "hours since 1900-01-01 00:00:00.0"
#t_cal = "gregorian"
#
#nc1_u = Dataset("ERA5_U_1981-1990.nc", mode='r')
#nc2_u = Dataset("ERA5_U_1991-2000.nc", mode='r')
#nc3_u = Dataset("ERA5_U_2001-2010.nc", mode='r')
#
#nc1_v = Dataset("ERA5_V_1981-1990.nc", mode='r')
#nc2_v = Dataset("ERA5_V_1991-2000.nc", mode='r')
#nc3_v = Dataset("ERA5_V_2001-2010.nc", mode='r')
#
#lat=nc1_u.variables["latitude"][:]
#lon=nc1_u.variables["longitude"][:]
#
#dt1=num2date(nc1_u.variables["time"][:],units = t_unit,calendar = t_cal)
#dt2=num2date(nc2_u.variables["time"][:],units = t_unit,calendar = t_cal)
#dt3=num2date(nc3_u.variables["time"][:],units = t_unit,calendar = t_cal)
#
#somau={1:"", 2:"", 3:"", 4:"", 5:"", 6:"", 7:"", 8:"", 9:"", 10:"", 11:"", 12:""}
#somav={1:"", 2:"", 3:"", 4:"", 5:"", 6:"", 7:"", 8:"", 9:"", 10:"", 11:"", 12:""}
#cont={1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0}
#mediau={}
#mediav={}
#
#somau, cont=nc_media_mensal(nc1_u, "u10", dt1, somav, cont)
#somau, cont=nc_media_mensal(nc2_u, "u10", dt2, somav, cont)
#somau, cont=nc_media_mensal(nc3_u, "u10", dt3, somav, cont)
#
#somav, cont=nc_media_mensal(nc1_v, "v10", dt1, somau, cont)
#somav, cont=nc_media_mensal(nc2_v, "v10", dt2, somau, cont)
#somav, cont=nc_media_mensal(nc3_v, "v10", dt3, somau, cont)


intensidade={}
for mes in range(1,13):
    mediau[mes]=somau[mes]/(cont[mes]+1)
    mediav[mes]=somav[mes]/(cont[mes]+1)

#    np.save("matrizes_mensais/media_u_"+str(mes)+".npy", mediau[mes])
#    np.save("matrizes_mensais/media_v_"+str(mes)+".npy", mediav[mes])
# 
#    mediau[mes]=np.load("matrizes_mensais/media_u_"+str(mes)+".npy")
#    mediav[mes]=np.load("matrizes_mensais/media_v_"+str(mes)+".npy")

    intensidade[mes]=(mediau[mes].data**2+mediav[mes].data**2)**0.5
    plota(lat, lon, mediau[mes].data, mediav[mes].data, intensidade[mes], mes)
#
#r=1
#u1=nc1_u.variables["u10"][:r, :, :]
#v1=nc1_v.variables["v10"][:r, :, :]
#a=""
#for uu in u1:
#    if a=="":
#        a=uu
#    else:
#        a=a+uu
#aa1=a/r
#
#aa2=np.sum(u1, axis=0)/r
##
#inte=(u1.data**2+v1.data**2)**0.5
#plota(lat, lon, u1.data, v1.data, inte, mes)

#plota_linear_pontos(intensi[:,252,108], 'S1 - LAT: -33.00, LON: -49.00')
#plota_linear_pontos(intensi[:,232,120], 'S2 - LAT: -28.00, LON: -46.00')
#plota_linear_pontos(intensi[:,218,134], 'SE1 - LAT: -24.50, LON: -42.50')
#plota_linear_pontos(intensi[:,205,151], 'SE2 - LAT: -21.25, LON: -38.25')
#plota_linear_pontos(intensi[:,168,162], 'NE1 - LAT: -12.00, LON: -35.50')
#plota_linear_pontos(intensi[:,128,158], 'NE2 - LAT: -2.00, LON: -36.50')
#plota_linear_pontos(intensi[:,110,120], 'NE3 - LAT: 2.50, LON: -46.00')


#X = [(14,2, (1,14)),(14,2,2), (14,2,4), (14,2,6), (14,2,8), (14,2,10), (14,2,12), (14,2,14)]
##X = [  (4,2,1),(4,2,2), (4,2,3), (4,2,5), (4,2,(4,6))]
#plt.subplots_adjust(bottom=0, left=0, top = 0.975, right=1)
#for nrows, ncols, plot_number in X:
#    plt.subplot(nrows, ncols, plot_number)
#    plt.xticks([])
#    plt.yticks([])
#
#plt.show()