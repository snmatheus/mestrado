import SAIDA as OUT
import numpy as np
import xarray as xr

lats = OUT.lats
lons = OUT.lons
print(np.min(lats), np.max(lats))
print(np.min(lons), np.max(lons))
exit
ym, xm = lats.shape[0], lons.shape[0]

def ln(x):
    from math import log
    return log(x)

# ------------------------------------------------------------------------- #
#                         Histograma de radiação                            #
# ------------------------------------------------------------------------- #

frequencias = {}
anoi, anof=1990, 2020
anoi, anof=1990, 1991
c = 0
for ano in range(anoi, anof):
    print(ano)
    netcdf_r = xr.open_dataset("../DADOS/RAD/ERA5_RAD_"+str(ano)+".nc")
    tt = netcdf_r['time'].data
    for t in range(0, len(tt)):
        print(t)
        campo = netcdf_r.variables['ssrd'][t,:,:].data/3600
        if campo.min()!=0 and campo.max()!=0:
            for lat in range(0, lats.shape[0]):
                for lon in range(0, lons.shape[0]):
                    r_pontual = int(round(campo[lat,lon]))
                    chave = int(str(r_pontual)[:-1]+'0')

                    if chave > 0:
                        if chave not in frequencias:
                            frequencias[chave] = 1
                        else:
                            frequencias[chave] = frequencias[chave]+1

import pandas as pd

df_freq = pd.DataFrame({'radiacao': list(frequencias.keys()), 'frequencia': list(frequencias.values())})
ftotal = np.sum(list(frequencias.values()))
df_freq['porcentagem'] = (df_freq['frequencia']/ftotal)*100
df_freq.to_csv('df_freq.csv', index=False)

# s = 0
# for i, x in enumerate(df_freq_p['porcentagem']):
#     s += x
#     if s>5:
#         break
# radiacao95 = df_freq_p['radiacao'].iloc[i]

#%%
import matplotlib.pyplot as plt
x = list(df_freq['radiacao'])
y = list(df_freq['porcentagem'])
#y = list(df_freq['frequencia'])
plt.xlim(0, np.max(x))
plt.ylim(0, 101)
plt.bar(list(frequencias.keys()), list(frequencias.values()))
plt.savefig('./teste.png')
plt.show()