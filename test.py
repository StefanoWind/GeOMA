# -*- coding: utf-8 -*-
"""
Test mapping operations
"""

from shapely.geometry import Polygon
import shapely
import numpy as np
from matplotlib import pyplot as plt
plt.close('all')

#%% Inputs
p1 = Polygon([(0,0), (1,0), (1,1), (0,1)])
p2 = Polygon([(1,0), (2,0), (2,1), (1,1)]) 
xmin=-2
xmax=2
ymin=-1.5
ymax=1.5
dx=0.1
dy=0.1

sigma=0.25

#%% Initialization
x=np.arange(xmin,xmax+dx/1,dx)
y=np.arange(ymin,ymax+dy/1,dy)
X,Y=np.meshgrid(x,y,indexing='ij')

#%% Main
p_int=p1.intersection(p2)

mask = shapely.contains_xy(p_int.buffer((dx**2+dy**2)**0.5/2), X, Y)

prob_inv=np.ones_like(X)
for i,j in zip(*np.where(mask)):
    R=((X-x[i])**2+(Y-y[j])**2)**0.5
    prob_inv*=(1-np.exp(-R**2/(2*sigma**2)))
    
prob=1-prob_inv

#%% Plots
fig, ax = plt.subplots()
plt.pcolor(x,y,prob.T,cmap='RdYlGn',edgecolor=(0,0,0,0.1))
# polygon 1
x1, y1 = p1.exterior.xy
ax.fill(x1, y1, alpha=0.5, label="Polygon 1")

# polygon 2
x2, y2 = p2.exterior.xy
ax.fill(x2, y2, alpha=0.5, label="Polygon 2")

# intersection (if any)
if not p_int.is_empty:
    if p_int.geom_type == "Polygon":
        x_int, y_int = p_int.exterior.xy
        ax.fill(x_int, y_int, color="red", alpha=0.8, label="Intersection")
    else:
        ax.plot(*p_int.xy, color="k", linewidth=2, label="Intersection")

ax.set_aspect("equal")
ax.legend()
plt.show

