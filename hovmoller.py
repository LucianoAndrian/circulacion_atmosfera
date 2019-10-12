#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 21:01:02 2019

@author: auri
"""
import matplotlib.pyplot as plt 
import numpy as np
from numpy import empty

def hovmoller(psi_c, psi, lon, lat, estado_basico, perturbacion,nombre):
    anomalia_10dias = empty([10,128,256]) 
    for i in np.arange(0,10):
        anomalia_10dias[i,:,:]=psi_c - psi[49+i,:,:]
    
    dias = np.arange(1,11)
    
    plt.figure()
    plt.contourf(lon[58:256],dias, anomalia_10dias[:,43,58:256], cmap = 'viridis') 
    plt.colorbar()
    plt.contour(lon[58:256],dias, anomalia_10dias[:,43,58:256], levels = 0, colors = "w" )
    plt.axvline(220, color = "r")
    plt.xlabel('longitud')
    plt.ylabel('Dias')
    plt.title("Hovmoller zonal $\Psi$ Estado Basico" + estado_basico + "Perturbacion" + perturbacion) 
    plt.axvline(x = 220,ymin = 0, ymax = 10 , color = "r")
    plt.savefig("Hovmoller_zonal_" + nombre, dpi = 200)
    
    plt.figure()
    plt.contourf(lat[0:86], dias, anomalia_10dias[:,0:86,220], cmap = 'viridis') 
    plt.colorbar()
    plt.contour(lat[0:86], dias, anomalia_10dias[:,0:86,220], levels = 0, colors = "w" )
    plt.axvline(-28, color = "r")
    plt.xlabel('longitud')
    plt.ylabel('Dias')
    plt.title("Hovmoller meridional $\Psi$ Estado Basico" + estado_basico + "Perturbacion" + perturbacion) 
    plt.savefig("Hovmoller_meridonal_" + nombre, dpi = 200)
