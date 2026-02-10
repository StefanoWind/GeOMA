# -*- coding: utf-8 -*-
"""
Probabilistic mapper of resource
"""
import os
cd=os.getcwd()
import numpy as np
from matplotlib import pyplot as plt
import warnings
import pandas as pd
import scipy as sp
import matplotlib
import re
import matplotlib.gridspec as gridspec
import geopandas as gpd
import xarray as xr
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['savefig.dpi'] = 300
plt.style.use("seaborn-v0_8-bright")
plt.close('all')
warnings.filterwarnings('ignore')


#%% Inputs
source=os.path.join(cd,'data','proprietary.gpkg')

#%% Initialization
Data = gpd.read_file(source)