#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 21:57:44 2020

@author: administrador
"""

import requests
import pandas as pd


urls = [
    "https://www.marinha.mil.br/chm/sites/www.marinha.mil.br.chm/files/historico_santos.txt",
    "https://www.marinha.mil.br/chm/sites/www.marinha.mil.br.chm/files/historico_itajai_0.txt",
    "https://www.marinha.mil.br/chm/sites/www.marinha.mil.br.chm/files/historico_riogrande_0.txt",
    "https://www.marinha.mil.br/chm/sites/www.marinha.mil.br.chm/files/historico_minuano.txt",
    "https://www.marinha.mil.br/chm/sites/www.marinha.mil.br.chm/files/historico_itaguai_1.txt",
    "https://www.marinha.mil.br/chm/sites/www.marinha.mil.br.chm/files/historico_itaoca.txt",
    "https://www.marinha.mil.br/chm/sites/www.marinha.mil.br.chm/files/historico_vitoria_0.txt",
    "https://www.marinha.mil.br/chm/sites/www.marinha.mil.br.chm/files/historico_portoseguro_0.txt",
    "https://www.marinha.mil.br/chm/sites/www.marinha.mil.br.chm/files/historico_recife_0.txt",
    "https://www.marinha.mil.br/chm/sites/www.marinha.mil.br.chm/files/historico_fortaleza.txt",
    "https://www.marinha.mil.br/chm/sites/www.marinha.mil.br.chm/files/historico_cabofrio2_0.txt",
]
        
boias = ["santos", "itajai", "riogrande", "minuano", "itaguai", "itaoca", "vitoria",
         "portoseguro", "recife", "fortaleza", "cabofrio2"]        

for i in range(0, len(urls)):
    url = urls[i]
    boia = boias[i]

    html = requests.get(url)
    
    conteudo = html.text.split("\n")
    lista = []
    
    conteudo.pop(0)
    conteudo.pop(-1)
    for l in conteudo:    
        lista.append(l.split(","))
    
    columns = ["Datetime","Lat","Lon","Battery","bHead","Wspd","Wdir","Gust","Atmp","Pres","Dewp","Humi",
              "Wtmp","Cvel1","Cdir1","Cvel2","Cdir2","Cvel3","Cdir3","Wvht","Wmax","Dpd","Mwd","Spread"]
    DF_BOIA = pd.DataFrame(lista, columns = columns)
    
    DF = DF_BOIA[["Datetime","Lat","Lon","Wspd","Wdir"]]
    DF.to_csv("../DADOS/PNBOIA/"+boia+".csv", index=False)
    
    ws = []
    c = 0
    tam = len(DF["Wspd"])
    for i in range(0, tam):
        l = DF["Wspd"].iloc[i]
        ws.append(float(l))
        if float(l) == -9999:
            c+=1
    
    print(boia, round((c/tam)*100,2),'% dos dados inválidos')

