#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 10:49:16 2018

@author: u1490431
"""

import numpy as np
import csv

def medfilereader(filename, varsToExtract = 'all',
                  sessionToExtract = 1,
                  verbose = False,
                  remove_var_header = False):
    if varsToExtract == 'all':
        numVarsToExtract = np.arange(0,26)
    else:
        numVarsToExtract = [ord(x)-97 for x in varsToExtract]
    
    f = open(filename, 'r')
    f.seek(0)
    filerows = f.readlines()[8:]
    datarows = [isnumeric(x) for x in filerows]
    matches = [i for i,x in enumerate(datarows) if x == 0.3]
    if sessionToExtract > len(matches):
        print('Session ' + str(sessionToExtract) + ' does not exist.')
    if verbose == True:
        print('There are ' + str(len(matches)) + ' sessions in ' + filename)
        print('Analyzing session ' + str(sessionToExtract))
    
    varstart = matches[sessionToExtract - 1]
    medvars = [[] for n in range(26)]
    
    k = int(varstart + 27)
    for i in range(26):
        medvarsN = int(datarows[varstart + i + 1])
        
        medvars[i] = datarows[k:k + int(medvarsN)]
        k = k + medvarsN
        
    if remove_var_header == True:
        varsToReturn = [medvars[i][1:] for i in numVarsToExtract]
    else:
        varsToReturn = [medvars[i] for i in numVarsToExtract]

    if np.shape(varsToReturn)[0] == 1:
        varsToReturn = varsToReturn[0]
    return varsToReturn

def metafilereader(filename):
    
    f = open(filename, 'r')
    f.seek(0)
    header = f.readlines()[0]
    f.seek(0)
    filerows = f.readlines()[1:]
    
    tablerows = []
    
    for i in filerows:
        tablerows.append(i.split('\t'))
        
    header = header.split('\t')
    # need to find a way to strip end of line \n from last column - work-around is to add extra dummy column at end of metafile
    return tablerows, header

def isnumeric(s):
    try:
        x = float(s)
        return x
    except ValueError:
        return float('nan')
    
def writemed2csv(medfile, varsToExtract):
    varsout = []
    for v in varsToExtract:
        varsout.append(medfilereader(medfile,
                                     varsToExtract=v,
                                     remove_var_header = True))
    print(varsout)
    
    with open(medfile+'.csv', 'w', newline="") as f:
        c = csv.writer(f)
        for v in varsout:
            c.writerows(zip(v[0], v[1]))

# Values etc to change



metafile = '/Volumes/KP_HARD_DRI/Michael_DPCP/DPCP_Metafile_Michael.txt'
medfolder = '/Volumes/KP_HARD_DRI/Michael_DPCP/DPCP_medfiles/'

# If only one bottle then a single list WITHIN a list e.g. [['b', 'c']]
varsToExtract = [['e', 'f']]


metadata, header = metafilereader(metafile)

for row in metadata:
    medfile = medfolder + row[0]
    writemed2csv(medfile, varsToExtract=varsToExtract)