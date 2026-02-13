# -*- coding: utf-8 -*-
"""
Probabilistic mapper of mineral resource
INPUTS:
    source (str): path to the geology survey data
    rock* (dict): substrings present in the respective data column to select rock type *
    xmin (float): minimum x [m]
    xmax (float): maximum x [m]
    ymin (float): minimum y [m]
    ymax (float): maximum y [m]
    dx (float>0): grid resolution in x [m]
    dy (float>0): grid resolution in y [m]
    fall_off (float>0): fall-off distance of probability [m]

OUTPUTS:
    Figures of probability maps
    Computation time
"""

import os
cd=os.getcwd()
import numpy as np
from shapely.ops import unary_union
import matplotlib
from matplotlib import pyplot as plt
import geopandas as gpd
from utils import probability, distance_map
import time
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['savefig.dpi'] = 300
plt.close('all')

#%% Inputs
source=os.path.join(cd,'data','proprietary.gpkg')#source file

#rock types selection
rock1={'rock_class':'ultramafic','rock_type':'serpentinite'}
rock2={'rock_type':'granodioritic'}

#grid
xmin=480000
xmax=640000
ymin=5520000
ymax=5680000
dx=1000
dy=1000

fall_off=10000#fall-off distance

#%% Initialization

#read data
Data = gpd.read_file(source)
n=len(Data.index)

#data selection
sel1=0
for c in rock1:
    sel1+=np.array([rock1[c] in Data[c].iloc[i] for i in range(n)])
sel1=np.where(sel1>0)[0]

sel2=0
for c in rock2:
    sel2+=np.array([rock2[c] in Data[c].iloc[i] for i in range(n)])
sel2=np.where(sel2>0)[0]

#grid definition
x=np.arange(xmin,xmax+dx/1,dx)
y=np.arange(ymin,ymax+dy/1,dy)
X,Y=np.meshgrid(x,y,indexing='ij')

#%% Main

#exctract releveant polygons
poly1=Data['geometry'].iloc[sel1]
poly2=Data['geometry'].iloc[sel2]

#probability map using distance from nearest interface
t0=time.time()
distance1=distance_map(unary_union(poly1),X,Y)
distance2=distance_map(unary_union(poly2),X,Y)
prob_near=probability(distance1,distance2,fall_off)
print(f'Elapsed time for nearest probability calculation: {time.time()-t0:.02f} s.')

#probability map superposing individual pair contributions
t0=time.time()
prob_inv=1
ctr=0
for p1 in poly1:
    distance1=distance_map(p1,X,Y)
    for p2 in poly2:
        distance2=distance_map(p2,X,Y)
        prob_inv*=1-probability(distance1,distance2,fall_off)
    ctr+=1
    print(f'{ctr}/{len(sel1)} pairwise maps done.')
    
prob_cum=1-prob_inv
print(f'Elapsed time for cumulated probability calculation: {time.time()-t0:.02f} s.')

#%% Plots

#nearest probabilit heat map
fig = plt.figure(figsize=(18,10))
plt.pcolor(x/1000,y/1000,prob_near.T,cmap=r'Greys',vmin=0,vmax=1)
for p1 in poly1:
    x1, y1 = p1.exterior.xy
    plt.plot(np.array(x1)/1000, np.array(y1)/1000, color='b')
    plt.fill(np.array(x1)/1000, np.array(y1)/1000, color='b',alpha=0.2)

for p2 in poly2:
    x2, y2 = p2.exterior.xy
    plt.plot(np.array(x2)/1000, np.array(y2)/1000, color='orange')
    plt.fill(np.array(x2)/1000, np.array(y2)/1000, color='orange',alpha=0.2)
    
plt.grid()
plt.xlim([xmin/1000,xmax/1000])
plt.ylim([ymin/1000,ymax/1000])
plt.xticks(np.arange(xmin/1000,xmax/1000,fall_off/1000))
plt.yticks(np.arange(ymin/1000,ymax/1000,fall_off/1000))
plt.gca().set_aspect("equal")
plt.xlabel(r'$x$ [km]')
plt.ylabel(r'$y$ [km]')
plt.colorbar(label='Probability (nearest)')

#cumulated probability heat map
fig = plt.figure(figsize=(18,10))
plt.pcolor(x/1000,y/1000,prob_cum.T,cmap=r'Greys',vmin=0,vmax=1)
for p1 in poly1:
    x1, y1 = p1.exterior.xy
    plt.plot(np.array(x1)/1000, np.array(y1)/1000, color='b')
    plt.fill(np.array(x1)/1000, np.array(y1)/1000, color='b',alpha=0.2)

for p2 in poly2:
    x2, y2 = p2.exterior.xy
    plt.plot(np.array(x2)/1000, np.array(y2)/1000, color='orange')
    plt.fill(np.array(x2)/1000, np.array(y2)/1000, color='orange',alpha=0.2)
    
plt.grid()
plt.xlim([xmin/1000,xmax/1000])
plt.ylim([ymin/1000,ymax/1000])
plt.xticks(np.arange(xmin/1000,xmax/1000,fall_off/1000))
plt.yticks(np.arange(ymin/1000,ymax/1000,fall_off/1000))
plt.gca().set_aspect("equal")
plt.xlabel(r'$x$ [km]')
plt.ylabel(r'$y$ [km]')
plt.colorbar(label='Probability (cumulated)')