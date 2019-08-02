#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 09:16:54 2018

@author: u1490431
"""
'''
   All imports and all functions needed for Distraction Paper 
   Run this script before anything else 
   
   Contains:
       
   loadmatfile, distractedOrNot, remcheck, distractionCalc2...
'''

# Import modules --------------------------------------------------------

import numpy as np
import scipy.io as sio

import matplotlib.pyplot as plt
import pandas as pd
import os
import matplotlib as mpl
import itertools
import matplotlib.mlab as mlab
#import seaborn as sb
import statistics as stats
# Set plot parameters and styles

#sb.set_context("paper")
#sb.set_style("white")

# Plot settings, font / size / styles
Calibri = {'fontname':'Calibri'}
Size = {'fontsize': 20}
label_size = 14
plt.rcParams['xtick.labelsize'] = label_size 
plt.rcParams['ytick.labelsize'] = label_size 

# Functions -------------------------------------------------------------
def mapTTLs(matdict):
    for x in ['LiA_', 'La2_']:
        try:
            licks = getattr(matdict['output'], x)
        except:
            print('File has no ' + x)
            
    lickson = licks.onset
    licksoff = licks.offset 
            
    return lickson, licksoff 


'''
    Loads a matlab converted TDT file, produces output session 
    dictionary of data from Synapse tankfiles. Photometry data
    and TTLs for licks, distractors etc.

'''

def loadmatfile(file):
    a = sio.loadmat(file, squeeze_me=True, struct_as_record=False)
    print(type(a))
    sessiondict = {}
    sessiondict['blue'] = a['output'].blue
    sessiondict['uv'] = a['output'].uv
    sessiondict['fs'] = a['output'].fs   
    
    try:
        
        sessiondict['licks'] = a['output'].licks.onset
        sessiondict['licks_off'] = a['output'].licks.offset
    except:  ## find which error it is
    
       sessiondict['licks'], sessiondict['licks_off'] = mapTTLs(a)
        
        

    sessiondict['distractors'] = distractionCalc2(sessiondict['licks'])

#   #write distracted or not to produce 2 lists of times, distracted and notdistracted
    #distracted, notdistracted= distractedOrNot(sessiondict['distractors'], sessiondict['licks'])
#    
    sessiondict['distracted'], sessiondict['notdistracted'] = distractedOrNot(sessiondict['distractors'], sessiondict['licks'])
  #  sessiondict['notdistracted'] = notdistracted
   
 # ''' sessiondict['lickRuns'] = lickRunCalc(sessiondict['licks']) ''' 
    
    return sessiondict
# -----------------------------------------------------------------

def distractedOrNot(distractors, licks):
    distracted = []
    notdistracted = []
    lickList = []
    for l in licks:
        lickList.append(l)
    

    for index, distractor in enumerate(distractors):
        if distractor in licks:

            ind = lickList.index(distractor)
            try:
                if (licks[ind+1] - licks[ind]) > 1:
                    distracted.append(licks[ind])
                else:
                    if (licks[ind+1] - licks[ind]) < 1:
                        notdistracted.append(licks[ind])
            except IndexError:
                print('last lick was a distractor!!!')
                distracted.append(licks[ind])

    return(distracted, notdistracted)


def remcheck(val, range1, range2):
    # function checks whether value is within range of two decimels
    if (range1 < range2):
        if (val > range1) and (val < range2):
            return True
        else:
            return False
    else:
        if (val > range1) or (val < range2):
            return True
        else:
            return False


def distractionCalc2(licks, pre=1, post=1):
    licks = np.insert(licks, 0, 0)
    b = 0.001
    d = []
    idx = 3
    
    while idx < len(licks):
        if licks[idx]-licks[idx-2] < 1 and remcheck(b, licks[idx-2] % 1, licks[idx] % 1) == False:
                d.append(licks[idx])
                b = licks[idx] % 1
                idx += 1
                try:
                    while licks[idx]-licks[idx-1] < 1:
                        b = licks[idx] % 1
                        idx += 1
                except IndexError:
                    pass
        else:
            idx +=1
#    print(len(d))
    
#    print(d[-1])
    if d[-1] > 3599:
        d = d[:-1]
        
#    print(len(d))
    
    return d


# LickCalc ============================================================
# Looking at function from Murphy et al (2017)

"""
This function will calculate data for bursts from a train of licks. The threshold
for bursts and clusters can be set. It returns all data as a dictionary.
"""
def lickCalc(licks, offset = [], burstThreshold = 0.25, runThreshold = 10, 
             binsize=60, histDensity = False):
    
    # makes dictionary of data relating to licks and bursts
    if type(licks) != np.ndarray or type(offset) != np.ndarray:
        try:
            licks = np.array(licks)
            offset = np.array(offset)
        except:
            print('Licks and offsets need to be arrays and unable to easily convert.')
            return

    lickData = {}
    
    if len(offset) > 0:
        lickData['licklength'] = offset - licks
        lickData['longlicks'] = [x for x in lickData['licklength'] if x > 0.3]
    else:
        lickData['licklength'] = []
        lickData['longlicks'] = []

    lickData['licks'] = np.concatenate([[0], licks])
    lickData['ilis'] = np.diff(lickData['licks'])
    lickData['freq'] = 1/np.mean([x for x in lickData['ilis'] if x < burstThreshold])
    lickData['total'] = len(licks)
    
    # Calculates start, end, number of licks and time for each BURST 
    lickData['bStart'] = [val for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > burstThreshold)]  
    lickData['bInd'] = [i for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > burstThreshold)]
    lickData['bEnd'] = [lickData['licks'][i-1] for i in lickData['bInd'][1:]]
    lickData['bEnd'].append(lickData['licks'][-1])

    lickData['bLicks'] = np.diff(lickData['bInd'] + [len(lickData['licks'])])    
    lickData['bTime'] = np.subtract(lickData['bEnd'], lickData['bStart'])
    lickData['bNum'] = len(lickData['bStart'])
    if lickData['bNum'] > 0:
        lickData['bMean'] = np.nanmean(lickData['bLicks'])
    else:
        lickData['bMean'] = 0
    
    lickData['bILIs'] = [x for x in lickData['ilis'] if x > burstThreshold]
    
    lickData['bILIs'] = [x for x in lickData['ilis'] if x > burstThreshold]

    # Calculates start, end, number of licks and time for each RUN
    lickData['rStart'] = [val for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > runThreshold)]  
    lickData['rInd'] = [i for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > runThreshold)]
    lickData['rEnd'] = [lickData['licks'][i-1] for i in lickData['rInd'][1:]]
    lickData['rEnd'].append(lickData['licks'][-1])

    lickData['rLicks'] = np.diff(lickData['rInd'] + [len(lickData['licks'])])    
    lickData['rTime'] = np.subtract(lickData['rEnd'], lickData['rStart'])
    lickData['rNum'] = len(lickData['rStart'])

    lickData['rILIs'] = [x for x in lickData['ilis'] if x > runThreshold]
    try:
        lickData['hist'] = np.histogram(lickData['licks'][1:], 
                                    range=(0, 3600), bins=int((3600/binsize)),
                                          density=histDensity)[0]
    except TypeError:
        print('Problem making histograms of lick data')
        
    return lickData  

def asnumeric(s):
    try:
        x = float(s)
        return x
    except ValueError:
        return float('nan')
    
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
    datarows = [asnumeric(x) for x in filerows]
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


def MetaExtractor (metafile):
    f = open(metafile, 'r')
    f.seek(0)
    Metafilerows = f.readlines()[1:]
    tablerows = []

    for row in Metafilerows: 
        items = row.split(',')
        tablerows.append(items)

    MedFilenames, RatID, Date, Day, Session, Drug, TotLicks, Distractions, \
    NonDistractions, PercentDistracted = [], [], [], [], [], [], [], [], [], []

    for i, lst in enumerate(tablerows):
       MedFilenames = MedFilenames + [lst[1]]
       RatID = RatID + [lst[2]]
       Date = Date + [lst[3]]
       #Day = Day + [lst[3]]
       Session = Session + [lst[4]]
     #  Drug = Drug + [lst[5]]
       TotLicks = TotLicks + [lst[6]]
       Distractions = Distractions + [lst[8]] 
       #NonDistractions = NonDistractions + [lst[8]]
       PercentDistracted = PercentDistracted + [lst[9]]
 
    return ({'MedFilenames':MedFilenames, 'RatID':RatID, 'Date':Date, 'Session':Session, \
             'TotLicks':TotLicks, 'Distractions':Distractions, \
             'PercentDistracted':PercentDistracted})
    
def time2samples(self):
    tick = self.output.Tick.onset
    maxsamples = len(tick)*int(self.fs)
    if (len(self.data) - maxsamples) > 2*int(self.fs):
        print('Something may be wrong with conversion from time to samples')
        print(str(len(self.data) - maxsamples) + ' samples left over. This is more than double fs.')
    
    self.t2sMap = np.linspace(min(tick), max(tick), maxsamples)
    
def snipper(data, timelock, fs = 1, t2sMap = [], preTrial=10, trialLength=30,
                 adjustBaseline = True,
                 bins = 0):

    if len(timelock) == 0:
        print('No events to analyse! Quitting function.')
        raise Exception('no events')
    nSnips = len(timelock)
    pps = int(fs) # points per sample
    pre = int(preTrial*pps) 
#    preABS = preTrial
    length = int(trialLength*pps)
# converts events into sample numbers
    event=[]
    if len(t2sMap) > 1:
        for x in timelock:
            event.append(np.searchsorted(t2sMap, x, side="left"))
    else:
        event = [x*fs for x in timelock]

    avgBaseline = []
    snips = np.empty([nSnips,length])

    for i, x in enumerate(event):
        start = int(x) - pre
        avgBaseline.append(np.mean(data[start : start + pre]))
#        print(x)
        try:
            snips[i] = data[start : start+length]
        except: # Deals with recording arrays that do not have a full final trial
            snips = snips[:-1]
            avgBaseline = avgBaseline[:-1]
            nSnips = nSnips-1

    if adjustBaseline == True:
        snips = np.subtract(snips.transpose(), avgBaseline).transpose()
        snips = np.divide(snips.transpose(), avgBaseline).transpose()

    if bins > 0:
        if length % bins != 0:
            snips = snips[:,:-(length % bins)]
        totaltime = snips.shape[1] / int(fs)
        snips = np.mean(snips.reshape(nSnips,bins,-1), axis=2)
        pps = bins/totaltime
              
    return snips, pps

def trialsFig(ax, trials1, trials2, pps=1, preTrial=10, scale=5, noiseindex = [],
              plotnoise=True,
              eventText='event', 
              ylabel=''):

    if len(noiseindex) > 0:
        trialsNoise = np.array([i for (i,v) in zip(trials1, noiseindex) if v])
        trials1 = np.array([i for (i,v) in zip(trials1, noiseindex) if not v])
        if plotnoise == True:
            ax.plot(trialsNoise.transpose(), c='red', alpha=0.4)
        
    ax.plot(trials1.transpose(), c='lightblue', alpha=0.6)
    
    
    ax.plot(trials2.transpose(), c='thistle', alpha=0.6)
    ax.plot(np.mean(trials2, axis=0), c='purple', linewidth=2)
    ax.plot(np.mean(trials1,axis=0), c='blue', linewidth=2)
    ax.set(ylabel = ylabel)
    ax.xaxis.set_visible(False)
            
    scalebar = scale * pps

    yrange = ax.get_ylim()[1] - ax.get_ylim()[0]
    scalebary = (yrange / 10) + ax.get_ylim()[0]
    scalebarx = [ax.get_xlim()[1] - scalebar, ax.get_xlim()[1]]
    
    ax.plot(scalebarx, [scalebary, scalebary], c='k', linewidth=2)
    ax.text((scalebarx[0] + (scalebar/2)), scalebary-(yrange/50), str(scale) +' s', ha='center',va='top', **Calibri, **Size)
 
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    xevent = pps * preTrial  
    ax.plot([xevent, xevent],[ax.get_ylim()[0], ax.get_ylim()[1] - yrange/20],'--')
    ax.text(xevent, ax.get_ylim()[1], eventText, ha='center',va='bottom', **Calibri, **Size)
    
    return ax




def med_abs_dev(data, b=1.4826):
    median = np.median(data)
    devs = [abs(i-median) for i in data]
    mad = np.median(devs)*b
                   
    return mad

def findnoise(data, background, t2sMap = [], fs = 1, bins=0, method='sd'):
    
    bgSnips, _ = snipper(data, background, t2sMap=t2sMap, fs=fs, bins=bins)
    
    if method == 'sum':
        bgSum = [np.sum(abs(i)) for i in bgSnips]
        bgMAD = med_abs_dev(bgSum)
        bgMean = np.mean(bgSum)
    elif method == 'sd':
        bgSD = [np.std(i) for i in bgSnips]
        bgMAD = med_abs_dev(bgSD)
        bgMean = np.mean(bgSD)
   
    return bgMAD, bgMean


def makerandomevents(minTime, maxTime, spacing = 77, n=100):
    events = []
    total = maxTime-minTime
    start = 0
    for i in np.arange(0,n):
        if start > total:
            start = start - total
        events.append(start)
        start = start + spacing
    events = [i+minTime for i in events]
    return events


def makephotoTrials(self, bins, events, threshold=10):
    bgMAD = findnoise(self.data, self.randomevents,
                          t2sMap = self.t2sMap, fs = self.fs, bins=bins,
                          method='sum')          
    blueTrials, self.pps = snipper(self.data, events,
                                        t2sMap = self.t2sMap, fs = self.fs, bins=bins)        
    UVTrials, self.pps = snipper(self.dataUV, events,
                                        t2sMap = self.t2sMap, fs = self.fs, bins=bins)
    sigSum = [np.sum(abs(i)) for i in blueTrials]
    sigSD = [np.std(i) for i in blueTrials]
    noiseindex = [i > bgMAD*threshold for i in sigSum]

    return blueTrials, UVTrials, noiseindex


def removenoise(snipsIn, noiseindex):
    snipsOut = np.array([x for (x,v) in zip(snipsIn, noiseindex) if not v])   
    return snipsOut

def trialsMultShadedFig(ax, trials, pps = 1, scale = 5, preTrial = 10,
                      eventText = 'event', ylabel = '',
                      linecolor=['purple', 'blue'], errorcolor=['thistle', 'lightblue'],
                        title=''):
    
    for i in [0, 1]:
        yerror = [np.std(i)/np.sqrt(len(i)) for i in trials[i].T]
        y = np.mean(trials[i],axis=0)
        x = np.arange(0,len(y))
    
        ax.plot(x, y, c=linecolor[i], linewidth=2)

        errorpatch = ax.fill_between(x, y-yerror, y+yerror, color=errorcolor[i], alpha=0.8)
    
    ax.set(ylabel = ylabel)
    ax.xaxis.set_visible(False)
            
    scalebar = scale * pps

    yrange = ax.get_ylim()[1] - ax.get_ylim()[0]
    scalebary = (yrange / 10) + ax.get_ylim()[0]
    scalebarx = [ax.get_xlim()[1] - scalebar, ax.get_xlim()[1]]
    
    ax.plot(scalebarx, [scalebary, scalebary], c='k', linewidth=2) # below in '' = 5
    ax.text((scalebarx[0] + (scalebar/2)), scalebary-(yrange/50), '5 s', ha='center',va='top', **Calibri, **Size)
 
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    xevent = pps * preTrial
    ax.plot([xevent, xevent],[ax.get_ylim()[0], ax.get_ylim()[1] - yrange/20],'--')
    ax.text(xevent, ax.get_ylim()[1], eventText, ha='center',va='bottom', **Calibri, **Size)
    ax.set_title(title, fontsize=14)
    
    return ax, errorpatch

def nearestevents(timelock, events, preTrial=10, trialLength=30):
#    try:
#        nTrials = len(timelock)
#    except TypeError:
#        nTrials = 1
    data = []
    start = [x - preTrial for x in timelock]
    end = [x + trialLength - preTrial for x in start]
    for start, end in zip(start, end):
        data.append([x for x in events if (x > start) & (x < end)])
    for i, x in enumerate(data):
        data[i] = x - timelock[i]      
    
    return data


