# -*- coding: utf-8 -*-
"""
Utilities for GeOMA
"""
from scipy.special import erfinv
import shapely
import numpy as np

def probability(distance1,distance2,fall_off,p0=1,tolerance=0.05):
    
    '''
    Calculate probability using double Gaussian shape function
    
    INPUTS:
        distance1 (np.array): map of Euclidean distance from nearest type 1 polygon [arbitrary units]
        distance2 (np.array): map of Euclidean distance from nearest type 2 polygon [arbitrary units]
        fall_off (float): fall-off distance [arbitrary units]
        p0 (0<float<1): maximum probability
        tolerance (0<float<1): fraction of integral of each Gaussian leaking outside of half the fall-off distance (double-sided)
    
    OUTPUTS:
        prob: probability map 
    '''
    
    sigma=fall_off/(2**0.5*erfinv(1-tolerance))/2#calculate Gaussian spread for given tolerance/fall-off
    prob=p0*np.exp(-(distance1**2+distance2**2)/(2*sigma**2))#probability function
    
    return prob

def distance_map(points,X,Y):
    '''
    Build of distance
    INPUTS:
        points (shapely points): points from which minimum distance is calculated
        X,Y (np.array): grid [arbitrary units]
    OUTPUTS:
        distance (np.array): map of minimum distance [arbitrary units]
    '''
    grid = shapely.points(X.ravel(), Y.ravel())
    distance =  shapely.distance(points, grid).reshape(X.shape)
    
    return distance