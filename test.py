# -*- coding: utf-8 -*-
"""
Test mapping operations
"""

from shapely.geometry import Polygon
from shapely.ops import unary_union
import shapely
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['savefig.dpi'] = 300

plt.close('all')

#%% Inputs
poly1 = [Polygon([(-1,0),(0,0.5),(0,0.5),(0,1),(-1,1),(-1,0)]),
      Polygon([(-2.5,-1.5),(-1.25,-1.5),(-1.25,1),(-2.5,-1)])]
poly2 = [Polygon([(0,-1),(0,1),(1,1),(1,-1)]),
      Polygon([(-1,-2),(-1,0),(0,0),(-1,-2)])] 

xmin=-3
xmax=2
ymin=-2
ymax=1.5
dx=0.01
dy=0.01

fall_off=0.5

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

#heat map
fig = plt.figure(figsize=(18,4))
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
plt.xticks(np.arange(xmin,xmax,fall_off))
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
plt.xticks(np.arange(xmin,xmax,fall_off))
plt.yticks(np.arange(ymin,ymax,fall_off))
ax.set_aspect("equal")
plt.xlabel(r'$x$')
plt.colorbar(label='Probability difference \n (cumulated-nearest)')
plt.tight_layout()