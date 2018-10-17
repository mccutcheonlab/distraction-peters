#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 09:08:16 2018

@author: u1490431
"""

"""
Chapter 4 - Distraction and photometry in VTA 


"""

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
# Distracted or not peaks
# Manuall add the lists here 
# Start with THPH1 and 2 lick days
# Then THPH1 and 2 distraction 
# Then THPH1 and 2 habituation 

# WHICH RATS DID NOT HAVE SIGNAL?
# THPH1 AND 2
# Lick day 
TDTfiles_thph_lick = ['thph1-1_lick6', 'thph1-2_lick6', 'thph1-3_lick6', 'thph1-4_lick6', 'thph1-5_lick6',\
                'thph1-6_lick6', 'thph2-1_lick3', 'thph2-2_lick3','thph2-3_lick3','thph2-4_lick3', \
                'thph2-5_lick3','thph2-6_lick3', 'thph2-7_lick6', 'thph2-8_lick6']

# Modelled distractors change variable names for this script or script section 
# Distraction day 
TDTfiles_thph_dis = ['thph1-3_distraction1', \
                'thph1-4_distraction1','thph1-5_distraction1','thph1-6_distraction1', \
                'thph2-1_distraction', 'thph2-2_distraction', 'thph2-3_distraction', \
                'thph2-4_distraction', 'thph2-5_distraction', 'thph2-6_distraction', \
                'thph2-7_distraction', 'thph2-8_distraction'] #'thph1-1_distraction1', 'thph1-2_distraction1'

# Habituation day 
TDTfiles_thph_hab = ['thph1-3_distraction2',\
                'thph1-4_distraction2', 'thph1-5_distraction2','thph1-6_distraction2',\
                'thph2-1_habituation', 'thph2-2_habituation', 'thph2-3_habituation', \
                'thph2-4_habituation', 'thph2-5_habituation', 'thph2-6_habituation', \
                'thph2-7_habituation'] #['thph1-1_distraction2','thph1-2_distraction2',


TDTfilepath = '/Volumes/KP_HARD_DRI/All_Matlab_Converts/BIG CONVERSION 14 AUG 2018/THPH matfiles/'

savepath = '/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures'

# LICKING ANALYSIS **************************************************************************
# Loop through files and calculate burst and run lengths
# Assign empty lists for storing arrays of burst/run lengths

# Lengths, all burst or run lengths here NOT
allBursts = []
allBurstsTimes = []
allRuns = []
allRunTimes = []

allRunIndices = []
allrILIs = []
allbILIs = []
allRatBlue = []
allRatUV = []
allRatFS = []
allRatLicks = []

for filename in TDTfiles_thph_lick:
    
    file = TDTfilepath + filename
    ratdata = loadmatfile(file)
    allRatBlue.append(ratdata['blue'])
    allRatUV.append(ratdata['uv'])
    allRatFS.append(ratdata['fs'])
    allRatLicks.append(ratdata['licks'])
    burstanalysis = lickCalc(ratdata['licks'], offset=ratdata['licks_off'])
    burstList = burstanalysis['bLicks'] # n licks per burst 
    runList = burstanalysis['rLicks'] # n licks per run
    burstListTimes = burstanalysis['bStart'] # Actual times of start of runs  
    runListTimes = burstanalysis['rStart'] # Actual times of start of bursts 
    allBursts.append(burstList)
    allRuns.append(runList)
    allRunTimes.append(runListTimes)
    allBurstsTimes.append(burstListTimes)
#    allRunIndices.append(indexRunList)
#    allrILIs.append(runILIs)
#    allbILIs.append(burstILIs)
    
# Make the list of lists into one long list for histogram 
MergedBurstList = list(itertools.chain.from_iterable(allBursts)) 
MergedRunList = list(itertools.chain.from_iterable(allRuns)) 
    
# Descriptives - aggregated data
meanburstlength = round(np.mean(MergedBurstList))
medburstlen = round(np.median(MergedBurstList))
meanrunlength = round(np.mean(MergedRunList))
medrunlen = round(np.median(MergedRunList))

##  Make the snips 


blueMeansBurst = []
uvMeansBurst = []
blueMeansRuns = []
uvMeansRuns = []
allbluesnips = []
alluvsnips = []

for i, val in enumerate(allRunTimes):
    
    # make a blue and uv snip for all 14, and noise remover / index
    blueSnips, ppsBlue = snipper(allRatBlue[i], allRunTimes[i], fs=allRatFS[i], bins=300)
    uvSnips, ppsUV = snipper(allRatUV[i], allRunTimes[i], fs=allRatFS[i], bins=300)

    randevents = makerandomevents(allRatBlue[i][300], allRatBlue[i][-300])
    bgMad, bgMean = findnoise(allRatBlue[i], randevents, fs=allRatFS[i], method='sum', bins=300)
    threshold = 1
    sigSum = [np.sum(abs(i)) for i in blueSnips]
    noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]

#     Might not need the noise index, this is just for trials fig 
    
#    fig = plt.figure()
#    ax = plt.subplot(1,1,1)
#    ax.set_ylim([-0.03, 0.03])
#    #ax.set_ylim([-0.05, 0.05])
#    trialsMultShadedFig(ax, [uvSnips,blueSnips], ppsBlue, eventText='First Lick in Run')
#    plt.text(250,0.03, '{}'.format(len(allRunTimes[i])) + ' Runs' )
#    
#    fig2 = plt.figure()
#    ax2 = plt.subplot(1,1,1)
#    ax2.set_ylim([-0.2, 0.2])
#    trialsFig(ax2, blueSnips, uvSnips, ppsBlue, eventText='First Lick in Run') #noiseindex=noiseindex) #, )
#    plt.text(250,0.2, '{}'.format(len(allRunTimes[i])) + ' Runs' )
#
#    filepath ='/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/'
#    ratname = str(i+1) +'.pdf'
#
#    fig2.savefig(filepath+ratname)

 # these four lines used later to define means plot (made after runs)
    blueMean = np.mean(blueSnips, axis=0)
    blueMeansRuns.append(blueMean)
    uvMean = np.mean(uvSnips, axis=0)
    uvMeansRuns.append(uvMean)
    allbluesnips.append(blueSnips)
    alluvsnips.append(uvSnips)

# All runs and all bursts (representative rat, etc.)
# Then segregate by long and short (for all rats quartiles not for each rat)


#########################################################################################
# Average of all runs, all rats, all trials 
## Mean of ALL runs and ALL rats on multishaded figure

#linecolor=['purple', 'blue'], errorcolor=['thistle', 'lightblue']


fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.03, 0.03])
ax.set_ylim([-0.05, 0.05])
trialsMultShadedFig(ax, [np.asarray(uvMeansRuns[2:]),np.asarray(blueMeansRuns)], ppsBlue, eventText='First Lick in Run', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'])
plt.text(250,0.03, '{}'.format(len(MergedRunList)) + ' Runs' ) ## Edit this to be all
#fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/All_Runs_All_Rats.pdf')


## Shows every single trial for each rat for runs - to choose representative sample
#for index, sniplist in enumerate(allbluesnips):
#    for ind, lis in enumerate(sniplist):
#        fig = plt.figure(figsize=(6,3))
#        ax = plt.subplot(1,1,1)
#        ax = plt.plot(allbluesnips[index][ind])
#        plt.text(250,0, '{}'.format([index,ind]))
        
#########################################################################################


# Individual trial 1 - 13,1
f = plt.figure(figsize=(6,2))
ax = plt.subplot(111)
ax.plot(allbluesnips[13][1] , color='blue')
ax.plot(alluvsnips[13][1] , color='purple')
triallicks = nearestevents(allRunTimes[13],allRatLicks[13])# allRatLicks[13], allRatLicks[13]) 
xvals1 = [(x+10)*10 for x in triallicks[1]] 
yvals1 = [ax.get_ylim()[1]] * len(xvals1)
ax.scatter(xvals1, yvals1, marker='|', s=90, c='k')  
ax.set_ylim([-0.2, 0.2])    
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
#f.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/SingleTrial_THPH2.8_Run2.pdf')


# Individual trial 1 - 13,6
f = plt.figure(figsize=(6,2))
ax = plt.subplot(111)
ax.plot(allbluesnips[13][6] , color='blue')
ax.plot(alluvsnips[13][6] , color='purple')
triallicks = nearestevents(allRunTimes[13],allRatLicks[13])# allRatLicks[13], allRatLicks[13]) 
xvals1 = [(x+10)*10 for x in triallicks[6]] 
yvals1 = [ax.get_ylim()[1]] * len(xvals1)
ax.scatter(xvals1, yvals1, marker='|', s=90, c='k')  
ax.set_ylim([-0.2, 0.2])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
#f.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/SingleTrial_THPH2.8_Run7.pdf',  bbox_inches="tight")


# Individual trial 1 - 12, 11
f = plt.figure(figsize=(6,2))
ax = plt.subplot(111)
ax.plot(allbluesnips[12][11] , color='blue')
ax.plot(alluvsnips[12][11] , color='purple')
triallicks = nearestevents(allRunTimes[12],allRatLicks[12])# allRatLicks[13], allRatLicks[13]) 
xvals1 = [(x+10)*10 for x in triallicks[11]] 
yvals1 = [ax.get_ylim()[1]] * len(xvals1)
ax.scatter(xvals1, yvals1, marker='|', s=90,c='k')  
ax.set_ylim([-0.2, 0.2])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
#f.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/SingleTrial_THPH2.7_Run12.pdf',  bbox_inches="tight")

# Individual trial 1 - 11, 10
f = plt.figure(figsize=(6,2))
ax = plt.subplot(111)
ax.plot(allbluesnips[11][10] , color='blue')
ax.plot(alluvsnips[11][10] , color='purple')
triallicks = nearestevents(allRunTimes[11],allRatLicks[11])# allRatLicks[13], allRatLicks[13]) 
xvals1 = [(x+10)*10 for x in triallicks[10]] 
yvals1 = [ax.get_ylim()[1]] * len(xvals1)
ax.scatter(xvals1, yvals1, marker='|', s=90,c='k')
ax.set_ylim([-0.2, 0.2])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
#f.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/SingleTrial_THPH2.6_Run11.pdf',  bbox_inches="tight")

# Individual trial 1 - 4, 33
f = plt.figure(figsize=(6,2))
ax = plt.subplot(111)
ax.plot(allbluesnips[4][33] , color='blue')
ax.plot(alluvsnips[4][33] , color='purple')
triallicks = nearestevents(allRunTimes[4],allRatLicks[4])# allRatLicks[13], allRatLicks[13]) 
xvals1 = [(x+10)*10 for x in triallicks[33]] 
yvals1 = [ax.get_ylim()[1]] * len(xvals1)
ax.scatter(xvals1, yvals1, marker='|', s=90,c='k')
ax.set_ylim([-0.2, 0.2])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
#f.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/SingleTrial_THPH1.5_Run34.pdf',  bbox_inches="tight")

# Individual trial 1 - 3, 13
f = plt.figure(figsize=(6,2))
ax = plt.subplot(111)
ax.plot(allbluesnips[3][13] , color='blue')
ax.plot(alluvsnips[3][13] , color='purple')
triallicks = nearestevents(allRunTimes[3],allRatLicks[3])# allRatLicks[13], allRatLicks[13]) 
xvals1 = [(x+10)*10 for x in triallicks[13]] 
yvals1 = [ax.get_ylim()[1]] * len(xvals1)
ax.scatter(xvals1, yvals1, marker='|', s=90,c='k')
ax.set_ylim([-0.2, 0.2])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
#f.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/SingleTrial_THPH1.4_Run14.pdf',  bbox_inches="tight")

# Individual trial 1 - 3, 13 ----- TESTER FOR SCALE 
f = plt.figure(figsize=(6,2))
ax = plt.subplot(111)
ax.plot(allbluesnips[3][13] , color='black')
ax.plot(alluvsnips[3][13] , color='black')
triallicks = nearestevents(allRunTimes[3],allRatLicks[3])# allRatLicks[13], allRatLicks[13]) 
xvals1 = [(x+10)*10 for x in triallicks[13]] 
yvals1 = [ax.get_ylim()[1]] * len(xvals1)
ax.scatter(xvals1, yvals1, marker='|', s=90,c='k')
ax.set_ylim([-0.2, 0.2])
# Making x scale
scale = 5
scalebar = scale * ppsBlue
yrange = ax.get_ylim()[1] - ax.get_ylim()[0]
scalebary = (yrange / 10) + ax.get_ylim()[0]
scalebarx = [ax.get_xlim()[1] - scalebar, ax.get_xlim()[1]]
ax.plot(scalebarx, [scalebary, scalebary], c='k', linewidth=2)
ax.text((scalebarx[0] + (scalebar/2)), scalebary-(yrange/50), str(scale) +' s', ha='center',va='top', **Calibri, **Size)
#f.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/SingleTrial_SCALE.pdf',  bbox_inches="tight")


###########################################################################################

# Long and short runs, 25th and 75th percentiles 

## Find short and long run lengths 
aggregateLowerQuart = np.percentile(MergedRunList,25)
aggregateUpperQuart = np.percentile(MergedRunList,75) 
    
allLogIndRuns = []
for runLicksList in allRuns:
    logIndRuns = [] # so 14 empty arrays exist and then don't (add to larger)
    for item in runLicksList:
        if item < aggregateLowerQuart:
            logIndRuns.append('LOWER')
        else:
            if item > aggregateUpperQuart:
                logIndRuns.append('UPPER')
            
            else:
                logIndRuns.append('MIDDLE')
    allLogIndRuns.append(logIndRuns)
            
# 14 lists of logical indices for whether the n licks in that run was L,M,U


 
# Final lists of run times sorted by u.m.l and stored by rat (14 lists)
lowerqRunTimes = []
uppqRunTimes = []
mid50RunTimes = []
lowerLengths = []
upperLengths = []

# Might need the index and the value for both lists ??
for i, listofindices in enumerate(allLogIndRuns):
    templower = []
    tempupper = []
    tempmid = []
    for j, runIndex in enumerate(listofindices):
        
        #print(i,j) i is the list index and j is the item index
        
        if runIndex == 'LOWER':
            lowrun = allRunTimes[i][j]
            templower.append(lowrun)
            lowerLengths.append(allRuns[i][j])
               
        if runIndex == 'UPPER':
            upprun = allRunTimes[i][j]
            tempupper.append(upprun)
            upperLengths.append(allRuns[i][j])
#                
        if runIndex == 'MIDDLE':
            midrun = allRunTimes[i][j]
            tempmid.append(midrun)
            
    lowerqRunTimes.append(templower)
    uppqRunTimes.append(tempupper)
    mid50RunTimes.append(tempmid)
 
    ## Figure how to also get the actual lengths? Must have done this. Line 603
    ## runLicksList inside allRuns 
# _________________________________________________________________________

         
#allign to lowerqRunTimes
#allign to uppqRunTimes

uvMeans_short_run = []
blueMeans_short_run = []
uvMeans_long_run = []
blueMeans_long_run = []
# Makes tonnes of individual plots             
# Individual rats might look odd as low Ns, but mean of mean will be better 
# Repeat this with bursts if looks like might be useful 

for i, val in enumerate(lowerqRunTimes):
    try:
        # make a blue and uv snip for all 14, and noise remover / index
        blueSnips, ppsBlue = snipper(allRatBlue[i], lowerqRunTimes[i], fs=allRatFS[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUV[i], lowerqRunTimes[i], fs=allRatFS[i], bins=300)
    
        randevents = makerandomevents(allRatBlue[i][300], allRatBlue[i][-300])
        bgMad, bgMean = findnoise(allRatBlue[i], randevents, fs=allRatFS[i], method='sum', bins=300)
        threshold = 1
        sigSum = [np.sum(abs(i)) for i in blueSnips]
        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
        # Might not need the noise index, this is just for trials fig 
    except: 
        pass
#    
#    fig12 = plt.figure()
#    ax10 = plt.subplot(1,1,1)
#    ax10.set_ylim([-0.03, 0.03])
#    #ax.set_ylim([-0.05, 0.05])
#    trialsMultShadedFig(ax10, [uvSnips,blueSnips], ppsBlue, eventText='First Lick Short Run')
#    plt.text(250,0.03, '{}'.format(len(lowerqRunTimes[i])) + ' short runs' )
##    
#    fig13 = plt.figure()
#    ax11 = plt.subplot(1,1,1)
#    ax11.set_ylim([-0.2, 0.2])
#    trialsFig(ax11, blueSnips, uvSnips, ppsBlue, eventText='First Lick in Short Run', noiseindex=noiseindex) #, )
#    plt.text(250,0.2, '{}'.format(len(lowerqRunTimes[i])) + ' short runs' )

# # these four lines used later to define means plot (made after runs)
    blueMeanSHORT = np.mean(blueSnips, axis=0)
    blueMeans_short_run.append(blueMeanSHORT)
    uvMeanSHORT = np.mean(uvSnips, axis=0)
    uvMeans_short_run.append(uvMeanSHORT)

MergedRunList_Short = list(itertools.chain.from_iterable(lowerqRunTimes)) 
# Average of all SHORT runs, all rats, all trials 
## Mean of ALL SHORT runs and ALL rats on multishaded figure

#linecolor=['purple', 'blue'], errorcolor=['thistle', 'lightblue']
fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.03, 0.03])
ax.set_ylim([-0.05, 0.05])
trialsMultShadedFig(ax, [np.asarray(uvMeans_short_run),np.asarray(blueMeans_short_run)], ppsBlue, eventText='First Lick in Short Run', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'])
plt.text(250,0.03, '{}'.format(len(MergedRunList_Short)) + ' Short Runs' ) ## Edit this to be all
fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Short_Runs_All_Rats.pdf', bbox_inches="tight")


## ===============================================================

# Photometry figures (individual) for LONG runs 

for i, val in enumerate(uppqRunTimes):
    try:
        # make a blue and uv snip for all 14, and noise remover / index
        blueSnips, ppsBlue = snipper(allRatBlue[i], uppqRunTimes[i], fs=allRatFS[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUV[i], uppqRunTimes[i], fs=allRatFS[i], bins=300)
    
        randevents = makerandomevents(allRatBlue[i][300], allRatBlue[i][-300])
        bgMad, bgMean = findnoise(allRatBlue[i], randevents, fs=allRatFS[i], method='sum', bins=300)
        threshold = 1
        sigSum = [np.sum(abs(i)) for i in blueSnips]
        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
        # Might not need the noise index, this is just for trials fig 
    except: 
        pass
    
#    fig14 = plt.figure()
#    ax12 = plt.subplot(1,1,1)
#    ax12.set_ylim([-0.03, 0.03])
#    #ax.set_ylim([-0.05, 0.05])
#    trialsMultShadedFig(ax12, [uvSnips,blueSnips], ppsBlue, eventText='First Lick Long Run')
#    plt.text(250,0.03, '{}'.format(len(uppqRunTimes[i])) + ' long runs' )
#    
#    fig14 = plt.figure()
#    ax13 = plt.subplot(1,1,1)
#    ax13.set_ylim([-0.2, 0.2])
#    trialsFig(ax13, blueSnips, uvSnips, ppsBlue, eventText='First Lick in Long Run', noiseindex=noiseindex) #, )
#    plt.text(250,0.2, '{}'.format(len(uppqRunTimes[i])) + ' long runs' )

# # these four lines used later to define means plot (made after runs) 
    blueMeanLONG = np.mean(blueSnips, axis=0)
    blueMeans_long_run.append(blueMeanLONG)
    uvMeanLONG = np.mean(uvSnips, axis=0)
    uvMeans_long_run.append(uvMeanLONG)

MergedRunList_Long = list(itertools.chain.from_iterable(uppqRunTimes)) 
# Average of all SHORT runs, all rats, all trials 
## Mean of ALL SHORT runs and ALL rats on multishaded figure

# Not sure how to turn the scale off here (removed ppsBlue)
#linecolor=['purple', 'blue'], errorcolor=['thistle', 'lightblue']
fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.03, 0.03])
ax.set_ylim([-0.05, 0.05])
trialsMultShadedFig(ax, [np.asarray(uvMeans_long_run),np.asarray(blueMeans_long_run)], ppsBlue, eventText='First Lick in Long Run', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'], scale=0)
plt.text(250,0.03, '{}'.format(len(MergedRunList_Long)) + ' Long Runs' ) ## Edit this to be all
fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Long_Runs_All_Rats.pdf', bbox_inches="tight")


fig = plt.figure(figsize=(6,3))
ax= plt.subplot(1,1,1)
ax.set_ylim([-0.05, 0.05])
LONG_SHORTrunMultFig = trialsMultShadedFig(ax, [np.asarray(blueMeans_short_run),np.asarray(blueMeans_long_run)], ppsBlue, eventText=''
                                                  , linecolor=['k', 'firebrick'], errorcolor=['darkgrey', 'darkorange'], scale=0)
ax.set(ylabel = chr(916) + 'F')
ax.yaxis.label.set_size(14)


# Simple figure legend
import matplotlib.patches as mpatches

orange_patch = mpatches.Patch(color='darkorange', label='Long Runs')
grey_patch = mpatches.Patch(color='darkgrey', label='Short Runs')
plt.legend(handles=[orange_patch, grey_patch], fontsize=14)
plt.show()
fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/ShortANDLong_Runs_All_Rats.pdf', bbox_inches="tight")


######
##############################################################################################################
##############################################################################################################
##############################################################################################################

# DISTRACTION FILES (minus the first 2) - this was run with all included 
### Distractors, distracted and not distracted, licks and blue / uv signals 


allRatBlue = []
allRatUV = []
allRatFS = []
allRatLicks = []
allRatDistractors = []
allRatDistracted = []
allRatNotDistracted = []
blueMeans_distractor = []
uvMeans_distractor = [] 
blueMeans_distracted = []
uvMeans_distracted = []
blueMeans_notdistracted = []
uvMeans_notdistracted = [] 
allbluesnips = []
alluvsnips = []

for filename in TDTfiles_thph_dis:
    
    file = TDTfilepath + filename
    ratdata = loadmatfile(file)
    allRatBlue.append(ratdata['blue'])
    allRatUV.append(ratdata['uv'])
    allRatFS.append(ratdata['fs'])
    allRatLicks.append(ratdata['licks'])
    allRatDistractors.append(ratdata['distractors'])
    allRatDistracted.append(ratdata['distracted'])
    allRatNotDistracted.append(ratdata['notdistracted'])


for i, val in enumerate(allRatDistractors):
    try:
        # make a blue and uv snip for all 14, and noise remover / index
        blueSnips, ppsBlue = snipper(allRatBlue[i], allRatDistractors[i], fs=allRatFS[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUV[i], allRatDistractors[i], fs=allRatFS[i], bins=300)
    
        randevents = makerandomevents(allRatBlue[i][300], allRatBlue[i][-300])
        bgMad, bgMean = findnoise(allRatBlue[i], randevents, fs=allRatFS[i], method='sum', bins=300)
        threshold = 1
        sigSum = [np.sum(abs(i)) for i in blueSnips]
        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
        # Might not need the noise index, this is just for trials fig 
    except: 
        pass
# Individual plots to choose a representative rat 
    
#    fig14 = plt.figure()
#    ax13 = plt.subplot(1,1,1)
#    ax13.set_ylim([-0.15, 0.15])
#    trialsFig(ax13, blueSnips, uvSnips, ppsBlue, eventText='Distractor') #, noiseindex=noiseindex) #, )
#    plt.text(250,0.2, '{}'.format(len(allRatDistractors[i])) + ' distractors' )
#    fig14.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Distractors_' + str(i) + '.pdf', bbox_inches="tight")

    
    blueMeanDISTRACTOR = np.mean(blueSnips, axis=0)
    blueMeans_distractor.append(blueMeanDISTRACTOR)
    uvMeanDISTRACTOR = np.mean(uvSnips, axis=0)
    uvMeans_distractor.append(uvMeanDISTRACTOR)


# Means for distractORS trials here MULT SHADED FIG 
fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.04, 0.04])
trialsMultShadedFig(ax, [np.asarray(uvMeans_distractor),np.asarray(blueMeans_distractor)], ppsBlue, eventText='Distractor', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'], scale=0)
# EDIT THIS TEXT TO SHOW NUMBER OF TOTAL DISTRACTORS OR TRIALS ON THE AVERAGED PLOT 
#plt.text(250,0.03, '{}'.format(len(MergedRunList_Long)) + ' Long Runs' ) ## Edit this to be all
#fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Distractors_All_Rats.pdf', bbox_inches="tight")



for i, val in enumerate(allRatDistracted):
    try:
        # make a blue and uv snip for all 14, and noise remover / index
        blueSnips, ppsBlue = snipper(allRatBlue[i], allRatDistracted[i], fs=allRatFS[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUV[i], allRatDistracted[i], fs=allRatFS[i], bins=300)
    
        randevents = makerandomevents(allRatBlue[i][300], allRatBlue[i][-300])
        bgMad, bgMean = findnoise(allRatBlue[i], randevents, fs=allRatFS[i], method='sum', bins=300)
        threshold = 1
        sigSum = [np.sum(abs(i)) for i in blueSnips]
        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
        # Might not need the noise index, this is just for trials fig 
    except: 
        pass
# Individual plots to choose a representative rat 
#    fig14 = plt.figure()
#    ax13 = plt.subplot(1,1,1)
#    ax13.set_ylim([-0.15, 0.15])
#    trialsFig(ax13, blueSnips, uvSnips, ppsBlue, eventText='Distracted') #, noiseindex=noiseindex) #, )
#    plt.text(250,0.2, '{}'.format(len(allRatDistracted[i])) + ' distracted' )
#    fig14.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Distracted_' + str(i) + '.pdf', bbox_inches="tight")
#

    blueMeanDISTRACTED = np.mean(blueSnips, axis=0)
    blueMeans_distracted.append(blueMeanDISTRACTED)
    uvMeanDISTRACTED = np.mean(uvSnips, axis=0)
    uvMeans_distracted.append(uvMeanDISTRACTED)
    allbluesnips.append(blueSnips)
    alluvsnips.append(uvSnips)
# Means for distracted trials here MULT SHADED FIG 
fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.04, 0.04])
trialsMultShadedFig(ax, [np.asarray(uvMeans_distracted),np.asarray(blueMeans_distracted)], ppsBlue, eventText='Distracted trial', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'], scale=0)
# EDIT THIS TEXT TO SHOW NUMBER OF TOTAL DISTRACTORS OR TRIALS ON THE AVERAGED PLOT 
#plt.text(250,0.03, '{}'.format(len(MergedRunList_Long)) + ' Long Runs' ) ## Edit this to be all
#fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Distracted_All_Rats.pdf', bbox_inches="tight")




for i, val in enumerate(allRatNotDistracted):
    try:
        # make a blue and uv snip for all 14, and noise remover / index
        blueSnips, ppsBlue = snipper(allRatBlue[i], allRatNotDistracted[i], fs=allRatFS[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUV[i], allRatNotDistracted[i], fs=allRatFS[i], bins=300)
    
        randevents = makerandomevents(allRatBlue[i][300], allRatBlue[i][-300])
        bgMad, bgMean = findnoise(allRatBlue[i], randevents, fs=allRatFS[i], method='sum', bins=300)
        threshold = 1
        sigSum = [np.sum(abs(i)) for i in blueSnips]
        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
        # Might not need the noise index, this is just for trials fig 
    except: 
        pass
    
# Individual plots to choose a representative rat 
    fig14 = plt.figure()
    ax13 = plt.subplot(1,1,1)
    ax13.set_ylim([-0.15, 0.15])
    trialsFig(ax13, blueSnips, uvSnips, ppsBlue, eventText='Not Distracted') #, noiseindex=noiseindex) #, )
    plt.text(250,0.2, '{}'.format(len(allRatNotDistracted[i])) + ' not distracted' )
    fig14.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/NotDistracted_' + str(i) + '.pdf', bbox_inches="tight")

# Means for not distracted trials here MULT SHADED FIG 

    blueMeanNOTDISTRACTED = np.mean(blueSnips, axis=0)
    blueMeans_notdistracted.append(blueMeanNOTDISTRACTED)
    uvMeanNOTDISTRACTED = np.mean(uvSnips, axis=0)
    uvMeans_notdistracted.append(uvMeanNOTDISTRACTED)
# Means for distracted trials here MULT SHADED FIG 
fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.04, 0.04])
trialsMultShadedFig(ax, [np.asarray(uvMeans_notdistracted),np.asarray(blueMeans_notdistracted)], ppsBlue, eventText='Not Distracted trial', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'], scale=0)
# EDIT THIS TEXT TO SHOW NUMBER OF TOTAL DISTRACTORS OR TRIALS ON THE AVERAGED PLOT 
#plt.text(250,0.03, '{}'.format(len(MergedRunList_Long)) + ' Long Runs' ) ## Edit this to be all
#fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/NotDistracted_All_Rats.pdf', bbox_inches="tight")


########################################################################################################


# (1) Make plots of each trial or pick ones that look good 

# Individual trial Run this with [0] --> [13] for each rat all trials saved 
# Manually chose trials to run single trials for

# Indices all moved by +2
# 2.7 - 3, 6, 9
# 2.6 - 3, 9, 14
# 2.5 - 5
# 2.1 - 0
# 1.5 - 0

#for ind, snip in enumerate(allbluesnips[13]):
#    print(ind)
#    f = plt.figure(figsize=(6,3))
#    ax = plt.subplot(111)
#    ax.plot(allbluesnips[13][ind] , color='blue')
#    ax.plot(alluvsnips[13][ind] , color='purple')
#
#
#    triallicks = nearestevents(allRatDistractors[13],allRatLicks[13])
#    trialdistractors = nearestevents(allRatDistractors[13],allRatDistractors[13])
#    trialdistracted = nearestevents(allRatDistractors[13],allRatDistracted[13])
#    trialnotdistracted = nearestevents(allRatDistractors[13],allRatNotDistracted[13])
#     
#    xvals1 = [(x+10)*10 for x in triallicks[ind]] 
#    xvals1 = [(x+10)*10 for x in triallicks[ind]]
#    #xvals2 = [(x+10)*10 for x in trialdistractors[trial]]
#    xvals3 = [(x+10)*10 for x in trialdistracted[ind]]
#    xvals4 = [(x+10)*10 for x in trialnotdistracted[ind]]
#    yvals1 = [ax.get_ylim()[1]] * len(xvals1)
#    #yvals2 = [ax.get_ylim()[1] + 0.005] * len(xvals2)
#    yvals3 = [ax.get_ylim()[1] + 0.02] * len(xvals3) 
#    yvals4 = [ax.get_ylim()[1] + 0.02] * len(xvals4)
#    #ax.scatter(xvals, yvals)
#    ax.scatter(xvals1, yvals1, marker='|', s=90, c='k')
#    #ax.scatter(xvals2, yvals2, marker='*')
#    ax.scatter(xvals3, yvals3, marker='o', facecolors= 'k', edgecolors='k', linewidth=2, s=60)
#    ax.scatter(xvals4, yvals4, marker='o', facecolors= 'none', edgecolors='k', linewidth=2, s=60)
#    
#  #  ax.set_ylim([-0.0, 0.05])
#    ax.spines['right'].set_visible(False)
#    ax.spines['top'].set_visible(False)
#    ax.spines['left'].set_visible(False)
#    ax.spines['bottom'].set_visible(False)
#    ax.xaxis.set_visible(False)
#    ax.yaxis.set_visible(False)
#    f.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/THPH2.8_trial_' + str(ind) +'.pdf',  bbox_inches="tight")

### REMEMBER the axes have been removed and the scales are different 

### Are these just alligned to distracted or both??? Check which list you use and how many trials 

## Smaller time period before ? Less than a second so you know no distractors are present???? 

'''

#### Representative rat 


## repeat on LICKING DAY - with modelled distractors and modelled distraction / or not 



### Then look at white noise vs non white noise 
## Then white noise on habituation day AND all distractors on habituation day (not distracted trials as very few, maybe)
## THEN compare the peaks with stats and export to SPSS for t-tests?

'''
################################################################################################
################################################################################################

## Subsetting data, to get peaks 

## Variables : stored blue and uv snips for different conditions (means for each rat, lists of 12 - not 14 1.1 and 1.2 removed)
## Do this all on a uv subtracted baseline?

# UV is essentially zero but check this 

# (1) UV signal subtracted from the blue 
# (2) Derive the peaks by following rules:
# (3) Compare using ANOVA or t-tests

## Expects list of 12 rats with mean snips in each field 
def uvSubtractor(rat_snip_means_list_blue, uv):
    subtractedSignal = []
    for ind, rat in enumerate(rat_snip_means_list_blue):
        subtractedSignal.append(rat_snip_means_list_blue[ind] - uv[ind])
        
    return subtractedSignal

bkgnd_sub_Runs = uvSubtractor(blueMeansRuns, uvMeansRuns)
bkgnd_sub_Short_Runs = uvSubtractor(blueMeans_short_run, uvMeans_short_run)
bkgnd_sub_Long_Runs = uvSubtractor(blueMeans_long_run, uvMeans_long_run)
bkgnd_sub_Distractor = uvSubtractor(blueMeans_distractor, uvMeans_distractor)
bkgnd_sub_Distracted = uvSubtractor(blueMeans_distracted, uvMeans_distracted)
bkgnd_sub_Notdistracted = uvSubtractor(blueMeans_notdistracted, uvMeans_notdistracted)

## Expects list of 12 rats with mean snips in each field 
## Give it the background subtracted snips or just the blue (as list of lists with eachlist a rat)

def PhotoPeaksCalc(snips_all_rats):
    
    allRat_peak = []
    allRat_t = []
    allRat_pre = []
    allRat_post = []
    allRat_base = []
    
    for rat in snips_all_rats:
        pre_event = np.mean(rat[0:50]) # Average for 5 seconds, 10 seconds before event 
        peak = np.max(rat[100:300]) ## Minus the average of the first 5 seconds and after 100 points (slice)
        peak_range = rat[100:300]
        a = peak_range.tolist()
        peak_index = a.index(peak)
        t = peak_index / 10
        pre_event = np.mean(rat[50:100])
        post_event = np.mean(rat[100:300])
        baseline = np.mean(rat[0:50])
        
        allRat_peak.append(peak)
        allRat_t.append(t)
        allRat_pre.append(pre_event)
        allRat_post.append(post_event)
        allRat_base.append(baseline) 
        
    
    return allRat_peak, allRat_t, allRat_pre, allRat_post, allRat_base


### Photometry peak variables - remember these lists will be unequal distractors 12 rats not 14?
### Should maybe take out the first 2 of the licks too? As they had no signal ???
## Manually removed the first 2 rats in SPSS -- should also remove them here (check later about the photo plots)

# All runs
peak_runs, t_runs, pre_runs, post_runs, baseline_runs = PhotoPeaksCalc(bkgnd_sub_Runs[2:])
# Short runs
peak_short_runs, t_short_runs, pre_short_runs, post_short_runs, baseline_short_runs = PhotoPeaksCalc(bkgnd_sub_Short_Runs[2:])
# Long runs
peak_long_runs, t_long_runs, pre_long_runs, post_long_runs, baseline_long_runs = PhotoPeaksCalc(bkgnd_sub_Long_Runs[2:])
# Distractors
peak_distractor, t_distractor, pre_distractor, post_distractor, baseline_distractor = PhotoPeaksCalc(bkgnd_sub_Distractor)
# Distracted
peak_distracted, t_distracted, pre_distracted, post_distracted, baseline_distracted = PhotoPeaksCalc(bkgnd_sub_Distracted)
# Not distracted 
peak_notdistracted, t_notdistracted, pre_notdistracted, post_notdistracted, baseline_notdistracted = PhotoPeaksCalc(bkgnd_sub_Notdistracted)


## Run the distraction code on the lick day data (minus rats 1 and 2) to get modelled distractor
## peaks 
##############################################################################################################
##############################################################################################################
##############################################################################################################

# DISTRACTION FILES (minus the first 2) - this was run with all included 
### Distractors, distracted and not distracted, licks and blue / uv signals 


allRatBlueMOD = []
allRatUVMOD = []
allRatFSMOD = []
allRatLicksMOD = []
allRatDistractorsMOD = []
allRatDistractedMOD = []
allRatNotDistractedMOD = []
blueMeans_distractorMOD = []
uvMeans_distractorMOD = [] 
blueMeans_distractedMOD = []
uvMeans_distractedMOD = []
blueMeans_notdistractedMOD = []
uvMeans_notdistractedMOD = [] 
allbluesnipsMOD = []
alluvsnipsMOD = []

for filename in TDTfiles_thph_lick[2:]:
    
    file = TDTfilepath + filename
    ratdata = loadmatfile(file)
    allRatBlueMOD.append(ratdata['blue'])
    allRatUVMOD.append(ratdata['uv'])
    allRatFSMOD.append(ratdata['fs'])
    allRatLicksMOD.append(ratdata['licks'])
    allRatDistractorsMOD.append(ratdata['distractors'])
    allRatDistractedMOD.append(ratdata['distracted'])
    allRatNotDistractedMOD.append(ratdata['notdistracted'])


for i, val in enumerate(allRatDistractorsMOD):
    try:
        # make a blue and uv snip for all 14, and noise remover / index
        blueSnips, ppsBlue = snipper(allRatBlueMOD[i], allRatDistractorsMOD[i], fs=allRatFSMOD[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUVMOD[i], allRatDistractorsMOD[i], fs=allRatFSMOD[i], bins=300)
    
        randevents = makerandomevents(allRatBlueMOD[i][300], allRatBlueMOD[i][-300])
        bgMad, bgMean = findnoise(allRatBlueMOD[i], randevents, fs=allRatFSMOD[i], method='sum', bins=300)
        threshold = 1
        sigSum = [np.sum(abs(i)) for i in blueSnips]
        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
        # Might not need the noise index, this is just for trials fig 
    except: 
        pass
# Individual plots to choose a representative rat 
    
#    fig14 = plt.figure()
#    ax13 = plt.subplot(1,1,1)
#    ax13.set_ylim([-0.15, 0.15])
#    trialsFig(ax13, blueSnips, uvSnips, ppsBlue, eventText='Distractor') #, noiseindex=noiseindex) #, )
#    plt.text(250,0.2, '{}'.format(len(allRatDistractors[i])) + ' distractors' )
#    fig14.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Distractors_' + str(i) + '.pdf', bbox_inches="tight")

    
    blueMeanDISTRACTOR = np.mean(blueSnips, axis=0)
    blueMeans_distractorMOD.append(blueMeanDISTRACTOR)
    uvMeanDISTRACTOR = np.mean(uvSnips, axis=0)
    uvMeans_distractorMOD.append(uvMeanDISTRACTOR)


# Means for distractORS trials here MULT SHADED FIG 
fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.04, 0.04])
trialsMultShadedFig(ax, [np.asarray(uvMeans_distractorMOD),np.asarray(blueMeans_distractorMOD)], ppsBlue, eventText='Modelled Distractor', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'], scale=0)
# EDIT THIS TEXT TO SHOW NUMBER OF TOTAL DISTRACTORS OR TRIALS ON THE AVERAGED PLOT 
#plt.text(250,0.03, '{}'.format(len(MergedRunList_Long)) + ' Long Runs' ) ## Edit this to be all
fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Modelled_Distractors_All_Rats.pdf', bbox_inches="tight")



for i, val in enumerate(allRatDistractedMOD):
    try:
        # make a blue and uv snip for all 14, and noise remover / index
        blueSnips, ppsBlue = snipper(allRatBlueMOD[i], allRatDistractedMOD[i], fs=allRatFSMOD[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUVMOD[i], allRatDistractedMOD[i], fs=allRatFSMOD[i], bins=300)
    
        randevents = makerandomevents(allRatBlueMOD[i][300], allRatBlueMOD[i][-300])
        bgMad, bgMean = findnoise(allRatBlueMOD[i], randevents, fs=allRatFSMOD[i], method='sum', bins=300)
        threshold = 1
        sigSum = [np.sum(abs(i)) for i in blueSnips]
        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
        # Might not need the noise index, this is just for trials fig 
    except: 
        pass
# Individual plots to choose a representative rat 
#    fig14 = plt.figure()
#    ax13 = plt.subplot(1,1,1)
#    ax13.set_ylim([-0.15, 0.15])
#    trialsFig(ax13, blueSnips, uvSnips, ppsBlue, eventText='Distracted') #, noiseindex=noiseindex) #, )
#    plt.text(250,0.2, '{}'.format(len(allRatDistracted[i])) + ' distracted' )
#    fig14.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Distracted_' + str(i) + '.pdf', bbox_inches="tight")
#

    blueMeanDISTRACTED = np.mean(blueSnips, axis=0)
    blueMeans_distractedMOD.append(blueMeanDISTRACTED)
    uvMeanDISTRACTED = np.mean(uvSnips, axis=0)
    uvMeans_distractedMOD.append(uvMeanDISTRACTED)
    allbluesnipsMOD.append(blueSnips)
    alluvsnipsMOD.append(uvSnips)
# Means for distracted trials here MULT SHADED FIG 
fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.04, 0.04])
trialsMultShadedFig(ax, [np.asarray(uvMeans_distractedMOD),np.asarray(blueMeans_distractedMOD)], ppsBlue, eventText='Distracted trial MOD', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'], scale=0)
# EDIT THIS TEXT TO SHOW NUMBER OF TOTAL DISTRACTORS OR TRIALS ON THE AVERAGED PLOT 
#plt.text(250,0.03, '{}'.format(len(MergedRunList_Long)) + ' Long Runs' ) ## Edit this to be all
#fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Distracted_All_Rats.pdf', bbox_inches="tight")

for i, val in enumerate(allRatNotDistractedMOD):
    try:
        # make a blue and uv snip for all 14, and noise remover / index
        blueSnips, ppsBlue = snipper(allRatBlueMOD[i], allRatNotDistractedMOD[i], fs=allRatFSMOD[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUVMOD[i], allRatNotDistractedMOD[i], fs=allRatFSMOD[i], bins=300)
    
        randevents = makerandomevents(allRatBlueMOD[i][300], allRatBlueMOD[i][-300])
        bgMad, bgMean = findnoise(allRatBlueMOD[i], randevents, fs=allRatFSMOD[i], method='sum', bins=300)
        threshold = 1
        sigSum = [np.sum(abs(i)) for i in blueSnips]
        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
        # Might not need the noise index, this is just for trials fig 
    except: 
        pass
    
## Individual plots to choose a representative rat 
#    fig14 = plt.figure()
#    ax13 = plt.subplot(1,1,1)
#    ax13.set_ylim([-0.15, 0.15])
#    trialsFig(ax13, blueSnips, uvSnips, ppsBlue, eventText='Not Distracted') #, noiseindex=noiseindex) #, )
#    plt.text(250,0.2, '{}'.format(len(allRatNotDistractedMOD[i])) + ' not distracted' )
##    fig14.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/NotDistracted_' + str(i) + '.pdf', bbox_inches="tight")

# Means for not distracted trials here MULT SHADED FIG 

    blueMeanNOTDISTRACTED = np.mean(blueSnips, axis=0)
    blueMeans_notdistractedMOD.append(blueMeanNOTDISTRACTED)
    uvMeanNOTDISTRACTED = np.mean(uvSnips, axis=0)
    uvMeans_notdistractedMOD.append(uvMeanNOTDISTRACTED)
# Means for distracted trials here MULT SHADED FIG 
fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.04, 0.04])
trialsMultShadedFig(ax, [np.asarray(uvMeans_notdistractedMOD),np.asarray(blueMeans_notdistractedMOD)], ppsBlue, eventText='Not Distracted trial MOD', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'], scale=0)
# EDIT THIS TEXT TO SHOW NUMBER OF TOTAL DISTRACTORS OR TRIALS ON THE AVERAGED PLOT 
#plt.text(250,0.03, '{}'.format(len(MergedRunList_Long)) + ' Long Runs' ) ## Edit this to be all
#fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/NotDistracted_All_Rats.pdf', bbox_inches="tight")

bkgnd_sub_Distractor_MOD = uvSubtractor(blueMeans_distractorMOD, uvMeans_distractorMOD)
bkgnd_sub_Distracted_MOD = uvSubtractor(blueMeans_distractedMOD, uvMeans_distractedMOD)
bkgnd_sub_Notdistracted_MOD = uvSubtractor(blueMeans_notdistractedMOD, uvMeans_notdistractedMOD)
# Distractors (in this case only these peaks in SPSS, may decide to use DIS and NOTDIS later though)
peak_distractorMOD, t_distractorMOD, pre_distractorMOD, post_distractorMOD, baseline_distractorMOD = PhotoPeaksCalc(bkgnd_sub_Distractor_MOD)
# Distracted
peak_distractedMOD, t_distractedMOD, pre_distractedMOD, post_distractedMOD, baseline_distractedMOD = PhotoPeaksCalc(bkgnd_sub_Distracted_MOD)
# Not distracted 
peak_notdistractedMOD, t_notdistractedMOD, pre_notdistractedMOD, post_notdistractedMOD, baseline_notdistractedMOD = PhotoPeaksCalc(bkgnd_sub_Notdistracted_MOD)






## NOW - make all photo peaks again from HABITUATION day not lick day or distraction day 

# Repeat everything (not the licking just distraction) and make all of the snips and peaks 
# Then transfer these to SPSS 

## Run the distraction code on the HABITUATION DAY data (minus rats 1 and 2) 
##############################################################################################################
##############################################################################################################
##############################################################################################################

# HABITUATION FILES (minus the first 2) 
### Distractors, distracted and not distracted, licks and blue / uv signals 

allRatBlueHAB = []
allRatUVHAB = []
allRatFSHAB = []
allRatLicksHAB = []
allRatDistractorsHAB = []
allRatDistractedHAB = []
allRatNotDistractedHAB = []
blueMeans_distractorHAB = []
uvMeans_distractorHAB = [] 
blueMeans_distractedHAB = []
uvMeans_distractedHAB = []
blueMeans_notdistractedHAB = []
uvMeans_notdistractedHAB = [] 
allbluesnipsHAB = []
alluvsnipsHAB = []

for filename in TDTfiles_thph_hab:
    
    file = TDTfilepath + filename
    ratdata = loadmatfile(file)
    allRatBlueHAB.append(ratdata['blue'])
    allRatUVHAB.append(ratdata['uv'])
    allRatFSHAB.append(ratdata['fs'])
    allRatLicksHAB.append(ratdata['licks'])
    allRatDistractorsHAB.append(ratdata['distractors'])
    allRatDistractedHAB.append(ratdata['distracted'])
    allRatNotDistractedHAB.append(ratdata['notdistracted'])


for i, val in enumerate(allRatDistractorsHAB):
    try:
        # make a blue and uv snip for all 14, and noise remover / index
        blueSnips, ppsBlue = snipper(allRatBlueHAB[i], allRatDistractorsHAB[i], fs=allRatFSHAB[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUVHAB[i], allRatDistractorsHAB[i], fs=allRatFSHAB[i], bins=300)
    
        randevents = makerandomevents(allRatBlueHAB[i][300], allRatBlueHAB[i][-300])
        bgMad, bgMean = findnoise(allRatBlueHAB[i], randevents, fs=allRatFSHAB[i], method='sum', bins=300)
        threshold = 1
        sigSum = [np.sum(abs(i)) for i in blueSnips]
        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
        # Might not need the noise index, this is just for trials fig 
    except: 
        pass
# Individual plots to choose a representative rat 
    
#    fig14 = plt.figure()
#    ax13 = plt.subplot(1,1,1)
#    ax13.set_ylim([-0.15, 0.15])
#    trialsFig(ax13, blueSnips, uvSnips, ppsBlue, eventText='Distractor') #, noiseindex=noiseindex) #, )
#    plt.text(250,0.2, '{}'.format(len(allRatDistractors[i])) + ' distractors' )
#    fig14.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Distractors_' + str(i) + '.pdf', bbox_inches="tight")

    
    blueMeanDISTRACTOR = np.mean(blueSnips, axis=0)
    blueMeans_distractorHAB.append(blueMeanDISTRACTOR)
    uvMeanDISTRACTOR = np.mean(uvSnips, axis=0)
    uvMeans_distractorHAB.append(uvMeanDISTRACTOR)


# Means for distractORS trials here MULT SHADED FIG 
fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.04, 0.04])
trialsMultShadedFig(ax, [np.asarray(uvMeans_distractorHAB),np.asarray(blueMeans_distractorHAB)], ppsBlue, eventText='Distractor', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'], scale=0)
# EDIT THIS TEXT TO SHOW NUMBER OF TOTAL DISTRACTORS OR TRIALS ON THE AVERAGED PLOT 
#plt.text(250,0.03, '{}'.format(len(MergedRunList_Long)) + ' Long Runs' ) ## Edit this to be all
#fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Habituation_Distractors_All_Rats.pdf', bbox_inches="tight")



for i, val in enumerate(allRatDistractedHAB):
    try:
        # make a blue and uv snip for all 14, and noise remover / index
        blueSnips, ppsBlue = snipper(allRatBlueHAB[i], allRatDistractedHAB[i], fs=allRatFSHAB[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUVHAB[i], allRatDistractedHAB[i], fs=allRatFSHAB[i], bins=300)
    
        randevents = makerandomevents(allRatBlueHAB[i][300], allRatBlueHAB[i][-300])
        bgMad, bgMean = findnoise(allRatBlueHAB[i], randevents, fs=allRatFSHAB[i], method='sum', bins=300)
        threshold = 1
        sigSum = [np.sum(abs(i)) for i in blueSnips]
        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
        # Might not need the noise index, this is just for trials fig 
    except: 
        pass
# Individual plots to choose a representative rat 
#    fig14 = plt.figure()
#    ax13 = plt.subplot(1,1,1)
#    ax13.set_ylim([-0.15, 0.15])
#    trialsFig(ax13, blueSnips, uvSnips, ppsBlue, eventText='Distracted') #, noiseindex=noiseindex) #, )
#    plt.text(250,0.2, '{}'.format(len(allRatDistracted[i])) + ' distracted' )
#    fig14.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Distracted_' + str(i) + '.pdf', bbox_inches="tight")
#

    blueMeanDISTRACTED = np.mean(blueSnips, axis=0)
    blueMeans_distractedHAB.append(blueMeanDISTRACTED)
    uvMeanDISTRACTED = np.mean(uvSnips, axis=0)
    uvMeans_distractedHAB.append(uvMeanDISTRACTED)
    allbluesnipsHAB.append(blueSnips)
    alluvsnipsHAB.append(uvSnips)
# Means for distracted trials here MULT SHADED FIG 
fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.04, 0.04])
trialsMultShadedFig(ax, [np.asarray(uvMeans_distractedHAB),np.asarray(blueMeans_distractedHAB)], ppsBlue, eventText='Distracted trial', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'], scale=0)
# EDIT THIS TEXT TO SHOW NUMBER OF TOTAL DISTRACTORS OR TRIALS ON THE AVERAGED PLOT 
#plt.text(250,0.03, '{}'.format(len(MergedRunList_Long)) + ' Long Runs' ) ## Edit this to be all
#fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Distracted_All_Rats.pdf', bbox_inches="tight")

for i, val in enumerate(allRatNotDistractedHAB):
    try:
        # make a blue and uv snip for all 14, and noise remover / index
        blueSnips, ppsBlue = snipper(allRatBlueHAB[i], allRatNotDistractedHAB[i], fs=allRatFSHAB[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUVHAB[i], allRatNotDistractedHAB[i], fs=allRatFSHAB[i], bins=300)
    
        randevents = makerandomevents(allRatBlueHAB[i][300], allRatBlueHAB[i][-300])
        bgMad, bgMean = findnoise(allRatBlueHAB[i], randevents, fs=allRatFSHAB[i], method='sum', bins=300)
        threshold = 1
        sigSum = [np.sum(abs(i)) for i in blueSnips]
        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
        # Might not need the noise index, this is just for trials fig 
    except: 
        pass
    
## Individual plots to choose a representative rat 
#    fig14 = plt.figure()
#    ax13 = plt.subplot(1,1,1)
#    ax13.set_ylim([-0.15, 0.15])
#    trialsFig(ax13, blueSnips, uvSnips, ppsBlue, eventText='Not Distracted') #, noiseindex=noiseindex) #, )
#    plt.text(250,0.2, '{}'.format(len(allRatNotDistractedMOD[i])) + ' not distracted' )
##    fig14.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/NotDistracted_' + str(i) + '.pdf', bbox_inches="tight")

# Means for not distracted trials here MULT SHADED FIG 

    blueMeanNOTDISTRACTED = np.mean(blueSnips, axis=0)
    blueMeans_notdistractedHAB.append(blueMeanNOTDISTRACTED)
    uvMeanNOTDISTRACTED = np.mean(uvSnips, axis=0)
    uvMeans_notdistractedHAB.append(uvMeanNOTDISTRACTED)
# Means for distracted trials here MULT SHADED FIG 
fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.04, 0.04])
trialsMultShadedFig(ax, [np.asarray(uvMeans_notdistractedHAB),np.asarray(blueMeans_notdistractedHAB)], ppsBlue, eventText='Not Distracted trial', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'], scale=0)
# EDIT THIS TEXT TO SHOW NUMBER OF TOTAL DISTRACTORS OR TRIALS ON THE AVERAGED PLOT 
#plt.text(250,0.03, '{}'.format(len(MergedRunList_Long)) + ' Long Runs' ) ## Edit this to be all
#fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/NotDistracted_All_Rats.pdf', bbox_inches="tight")

bkgnd_sub_Distractor_HAB = uvSubtractor(blueMeans_distractorHAB, uvMeans_distractorHAB)
bkgnd_sub_Distracted_HAB = uvSubtractor(blueMeans_distractedHAB, uvMeans_distractedHAB)
bkgnd_sub_Notdistracted_HAB = uvSubtractor(blueMeans_notdistractedHAB, uvMeans_notdistractedHAB)
# Distractors (in this case only these peaks in SPSS, may decide to use DIS and NOTDIS later though)
peak_distractorHAB, t_distractorHAB, pre_distractorHAB, post_distractorHAB, baseline_distractorHAB = PhotoPeaksCalc(bkgnd_sub_Distractor_HAB)
# Distracted
peak_distractedHAB, t_distractedHAB, pre_distractedHAB, post_distractedHAB, baseline_distractedHAB = PhotoPeaksCalc(bkgnd_sub_Distracted_HAB)
# Not distracted 
peak_notdistractedHAB, t_notdistractedHAB, pre_notdistractedHAB, post_notdistractedHAB, baseline_notdistractedHAB = PhotoPeaksCalc(bkgnd_sub_Notdistracted_HAB)

