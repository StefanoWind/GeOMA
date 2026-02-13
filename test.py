# -*- coding: utf-8 -*-
"""
Test mapping operations

INPUTS:
    poly1 (list of shapely polygons): list of polygons of type 1 [arbitrary units]
    poly2 (list of shapely polygons): list of polygons of type 2 [arbitrary units]
    xmin (float): minimum x [arbitrary units]
    xmax (float): maximum x [arbitrary units]
    ymin (float): minimum y [arbitrary units]
    ymax (float): maximum y [arbitrary units]
    dx (float>0): grid resolution in x [arbitrary units]
    dy (float>0): grid resolution in y [arbitrary units]
    fall_off (float>0): fall-off distance of probability [arbitrary units]

OUTPUTS:
    Figures of sample probability distribution and maps
"""

from shapely.geometry import Polygon
from shapely.ops import unary_union
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from utils import probability, distance_map

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['savefig.dpi'] = 300

plt.close('all')

#%% Inputs

#define lists of polygons
poly1 = [Polygon([(-1,0),(0,0.5),(0,0.5),(0,1),(-1,1),(-1,0)]),
      Polygon([(-2.5,-1.5),(-1.25,-1.5),(-1.25,1),(-2.5,-1)])]
poly2 = [Polygon([(0,-1),(0,1),(1,1),(1,-1)]),
      Polygon([(-1,-2),(-1,0),(0,0),(-1,-2)])] 

#grid
xmin=-3
xmax=2
ymin=-2
ymax=1.5
dx=0.01
dy=0.01

fall_off=0.5#fall-off distance

#%% Initialization

#grid definition
x=np.arange(xmin,xmax+dx/1,dx)
y=np.arange(ymin,ymax+dy/1,dy)
X,Y=np.meshgrid(x,y,indexing='ij')

#%% Main

#calculate distance grid
distance1=distance_map(unary_union(poly1),X,Y)
distance2=distance_map(unary_union(poly2),X,Y)

#probability map using distance from nearest interface
prob_near=probability(distance1,distance2,fall_off)

#probability map superposing individual pair contributions
prob_inv=1
for p1 in poly1:
    distance1=distance_map(p1,X,Y)
    for p2 in poly2:
        distance2=distance_map(p2,X,Y)
        prob_inv*=1-probability(distance1,distance2,fall_off)

prob_cum=1-prob_inv

#%% Plots

#probability function
cmap = plt.cm.jet
dist_plot=np.arange(0,1.1,0.1)
colors = cmap(np.linspace(0, 1, len(dist_plot)))
plt.figure()
x_plot=np.arange(-1,1.1,0.001)
ctr=0
for d in dist_plot:
    prob_plot=probability(x_plot-d/2,x_plot+d/2, fall_off)
    plt.plot(x_plot,prob_plot,color=colors[ctr],label=f'Distance = {d:.2f}')
    plt.plot(-d/2,1,'.',color=colors[ctr])
    plt.plot(d/2,1,'.',color=colors[ctr])

    ctr+=1
plt.title(f'Fall off distance = {fall_off}')
plt.xlabel('$x$')
plt.ylabel('$P$')
plt.grid()
plt.legend()

#heat maps
fig = plt.figure(figsize=(18,3))
ax=plt.subplot(1,3,1)
plt.pcolor(x,y,prob_near.T,cmap='Greys',vmin=0,vmax=1)
for p1 in poly1:
    x1, y1 = p1.exterior.xy
    plt.plot(np.array(x1), np.array(y1), color='b')
    plt.fill(np.array(x1), np.array(y1), color='b',alpha=0.2)

for p2 in poly2:
    x2, y2 = p2.exterior.xy
    plt.plot(np.array(x2), np.array(y2), color='orange')
    plt.fill(np.array(x2), np.array(y2), color='orange',alpha=0.2)
    
plt.grid()
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks(np.arange(xmin,xmax,fall_off),rotation=90)
plt.yticks(np.arange(ymin,ymax,fall_off))
ax.set_aspect("equal")
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.colorbar(label='Probability (nearest)')

ax=plt.subplot(1,3,2)
plt.pcolor(x,y,prob_cum.T,cmap='Greys',vmin=0,vmax=1)
for p1 in poly1:
    x1, y1 = p1.exterior.xy
    plt.plot(np.array(x1), np.array(y1), color='b')
    plt.fill(np.array(x1), np.array(y1), color='b',alpha=0.2)

for p2 in poly2:
    x2, y2 = p2.exterior.xy
    plt.plot(np.array(x2), np.array(y2), color='orange')
    plt.fill(np.array(x2), np.array(y2), color='orange',alpha=0.2)

plt.grid()
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks(np.arange(xmin,xmax,fall_off),rotation=90)
plt.yticks(np.arange(ymin,ymax,fall_off))
ax.set_aspect("equal")
plt.xlabel(r'$x$')
plt.colorbar(label='Probability (cumulated)')

ax=plt.subplot(1,3,3)
plt.pcolor(x,y,prob_cum.T-prob_near.T,cmap='seismic',vmin=-0.2,vmax=0.2)
for p1 in poly1:
    x1, y1 = p1.exterior.xy
    plt.plot(np.array(x1), np.array(y1), color='b')
    plt.fill(np.array(x1), np.array(y1), color='b',alpha=0.2)

for p2 in poly2:
    x2, y2 = p2.exterior.xy
    plt.plot(np.array(x2), np.array(y2), color='orange')
    plt.fill(np.array(x2), np.array(y2), color='orange',alpha=0.2)

plt.grid()
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xticks(np.arange(xmin,xmax,fall_off),rotation=90)
plt.yticks(np.arange(ymin,ymax,fall_off))
ax.set_aspect("equal")
plt.xlabel(r'$x$')
plt.colorbar(label='Probability difference \n (cumulated-nearest)')
plt.tight_layout()