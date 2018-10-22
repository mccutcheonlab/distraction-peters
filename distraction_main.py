#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 15:50:31 2018

@author: u1490431
"""
# Moved these 2 imports into jmf 
import JM_general_functions as jmf
import JM_custom_figs as jmfig
from distraction_assemble import *
from distraction_sessionfigs import *

import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import scipy.io as sio

# Extracts data from metafiles and sets up ppp_sessions dictionary

# Where to deposit the data 
picklefolder = '/Volumes/KP_HARD_DRI/distraction_paper/output/'
xls_file = '/Volumes/KP_HARD_DRI/distraction_paper/THPH1_THPH2.xlsx'
metafile_out = '/Volumes/KP_HARD_DRI/distraction_paper/THPH1_THPH2_metafile'
datafolder = '/Volumes/KP_HARD_DRI/distraction_paper/THPH matfiles/'
# outputfolder = picklefolder
sheetname = 'THPH1&2Metafile' 

distraction_sessions = metafile2sessions(xls_file, metafile_out, datafolder, picklefolder, sheetname)

#
## Code to indictae which files to assemble and whether to save and/or make figures
#assemble_sacc = False
#assemble_cond1 = False
#assemble_cond2 = True
#assemble_pref = False
assemble_single = True 
#
savefile=False
makefigs=False

# Code to run for single rat
if assemble_single:
    sessions_to_add = assemble_sessions(distraction_sessions,
                  rats_to_include = ['thph1-1'],
                  rats_to_exclude = ['thph1-2'],
                  sessions_to_include = ['lick1'], # come up with new name for session / day etc.
                  outputfile=picklefolder + 'ppp_test.pickle',
                  savefile=savefile,
                  makefigs=makefigs)
