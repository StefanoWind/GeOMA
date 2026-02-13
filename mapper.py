# -*- coding: utf-8 -*-
"""
Probabilistic mapper of resource
"""
import os
cd=os.getcwd()
import numpy as np
import shapely
from shapely.ops import unary_union
import matplotlib
from matplotlib import pyplot as plt
import geopandas as gpd
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['savefig.dpi'] = 300
plt.style.use("seaborn-v0_8-bright")
plt.close('all')

#%% Inputs
source=os.path.join(cd,'data','proprietary.gpkg')

rock1={'rock_class':'ultramafic','rock_type':'serpentinite'}
rock2={'rock_type':'granodioritic'}

xmin=480000
xmax=640000
ymin=5520000
ymax=5680000
dx=1000
dy=1000
fall_off=10000

#%%Functions
def probability(distance1,distance2,fall_off,p0=1,tolerance=0.05):
    from scipy.special import erfinv
    '''
    Gaussian probability shape function
    
    INPUTS:
    distance: Euclidean distance from interface
    p0: maximum probability
    cut_out: fall-off distance
    tolerance: fraction of probability integral outside fall-off (double-sided)
    '''
    sigma=fall_off/(2**0.5*erfinv(1-tolerance))/2
    prob=p0*np.exp(-distance1**2/(2*sigma**2))*np.exp(-distance2**2/(2*sigma**2))
    
    return prob

def distance_map(points,X,Y):
    grid = shapely.points(X.ravel(), Y.ravel())
    distance =  shapely.distance(points, grid).reshape(X.shape)
    
    return distance

#%% Initialization
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

x=np.arange(xmin,xmax+dx/1,dx)
y=np.arange(ymin,ymax+dy/1,dy)
X,Y=np.meshgrid(x,y,indexing='ij')

#%% Main

poly1=Data['geometry'].iloc[sel1]
poly2=Data['geometry'].iloc[sel2]

#calculate distance grid
distance1=distance_map(unary_union(poly1),X,Y)
distance2=distance_map(unary_union(poly2),X,Y)

#probability map using distance from nearest interface
prob_near=probability(distance1,distance2,fall_off)

#probability map superposing individual pair contributions
prob_inv=1
ctr=0
for p1 in poly1:
    distance1=distance_map(p1,X,Y)
    for p2 in poly2:
        distance2=distance_map(p2,X,Y)
        prob_inv*=1-probability(distance1,distance2,fall_off)
    ctr+=1
    print(f'{ctr}/{len(sel1)} polygons done')
    
prob_cum=1-prob_inv

#%% Plots
#heat map
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