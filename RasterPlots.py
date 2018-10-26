#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 21:17:27 2018

@author: u1490431
"""

# Post distraction pause rater plot 

def distractionrasterFig(ax, timelock, events,
                         pre = 1, post = 1,
                         sortevents=None, sortdirection='ascending'):

    if sortevents != None:
        if len(timelock) != len(sortevents):
            print('Length of sort events does not match timelock events; no sorting')
            
        ## Rats (like 6) where the last distractor is the last lick (so no pdp at the end)
        ## Repeated code here which could be turned into a "sortevents" function if wanted
        if len(timelock) == (len(sortevents) + 1):
            sortevents.append(0)
            
            if sortdirection == 'ascending':
                sortOrder = np.argsort(sortevents)
            else:
                sortOrder = np.argsort(sortevents)[::-1]
                
            timelock = [timelock[i] for i in sortOrder]
    

            
        else:
            if sortdirection == 'ascending':
                sortOrder = np.argsort(sortevents)
            else:
                sortOrder = np.argsort(sortevents)[::-1]
                
            timelock = [timelock[i] for i in sortOrder]
    
    rasterData = [[] for i in timelock]
    
    for i,x in enumerate(timelock):
        rasterData[i] = [j-x for j in events if (j > x-pre) & (j < x+post)]
 

# Works out how many of the trials are distracted (pauses longer than 1 sec in sortevents)            
    count = 0
    for a, b in enumerate(sortevents):
        if b > 1:
            count += 1

    
    for ith, trial in enumerate(rasterData):
        xvals = [x for x in trial] 
        yvals = [1+ith] * len(xvals)

# if sorted events != None 
        if ith < count:  # 26 if ascending --> MANUALLY CHECKED? --> could check the len of PDP if not ordered
            ax.scatter(xvals, yvals, marker=',', s=1, color='r')
     
        else:
            ax.scatter(xvals, yvals, marker=',', s=1, color='k')

# else:
    

# Taken from Ch4 analysis_licking - currently contains a lot of redundant code 

# On licking day (modelled distractors --> all rats)

# On distraction day (real distractors)
# On habituation day (real distractors but less distracted)

# Ordered then not ordered (check how you will change the colours)
 
### SORTED PLOTS HERE - MODELLED DAY, DISTRACTION DAY, HABITUATION DAY 

## lick day (modelled distractors) 

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

    allRatDistractorsMOD.append(ratdata['distractors'])
    allRatDistractedMOD.append(ratdata['distracted'])
    allRatNotDistractedMOD.append(ratdata['notdistracted'])
  #  print(ratdata['distractors'])

    indices1 = []       
    for index, value in enumerate(ratdata['distractors']):
        a = np.where(ratdata['licks'] == value) 
        indices1.append(a)       
    
    pdps = []
    for tupl in indices1:
        i = tupl[0][0]
        if i+1 < len(ratdata['licks']):
            
            pdp = (ratdata['licks'][i+1] - ratdata['licks'][i])
            pdps.append(pdp)    
    
   # pdps.append(0) ## Add in a check, if PDPs len does not equal len(distractors)
    figure12 = plt.figure(figsize=(6,6))
    ax6 = plt.subplot(111)
    ax6.spines['right'].set_visible(False)
    ax6.xaxis.set_visible(False)
    ax6.spines['top'].set_visible(False)
    ax6.spines['bottom'].set_visible(False)
    ax6.set(ylabel = 'Trials')
    ax6.yaxis.label.set_size(14)
    
    scale = 1
    scalebar = 1
    yrange = ax6.get_ylim()[1] - ax6.get_ylim()[0]
    scalebary = (yrange / 10) + ax6.get_ylim()[0]
    scalebarx = [ax6.get_xlim()[1] - scalebar, ax6.get_xlim()[1]]
    ax6.plot(scalebarx, [scalebary, scalebary], c='k', linewidth=2)
    ax6.text((scalebarx[0] + (scalebar/2)), scalebary-(yrange/50), str(scale) +' s', ha='center',va='top', **Calibri, **Size) 
    
    rasterPlot = distractionrasterFig(ax6, ratdata['distractors'], ratdata['licks'], pre=1, post=10, sortevents=pdps, sortdirection='dec')

#    fig12.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Distractors_' + str(i) + '.pdf', bbox_inches="tight")


#
#Make plots for every rat on lick day, distraction day and habituation day
#Choose the best lookng representative rat for the ordered, non ordered (by time)
#Make the dots where it is distracted a different colour (this was added into a previous script)
#Sort out the spacing issues (talk to Jaime perhaps) - replace with a different marker or something 
#    Could make the plots as large as possible 
#    Think about issues of different scales on (1) Different days and (2) Different rats 
#    Consider using ',' rather than '.' (a pixel and not a dot)
#    

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
    indices1 = []       
    
    for index, value in enumerate(ratdata['distractors']):
        a = np.where(ratdata['licks'] == value) 
        indices1.append(a)       
    
    pdps = []
    for tupl in indices1:
        i = tupl[0][0]
        if i+1 < len(ratdata['licks']):
            
            pdp = (ratdata['licks'][i+1] - ratdata['licks'][i])
            pdps.append(pdp)    
    
   # pdps.append(0) ## Add in a check, if PDPs len does not equal len(distractors)
    figure12 = plt.figure(figsize=(6,6))
    ax6 = plt.subplot(111)
    ax6.spines['right'].set_visible(False)
    ax6.xaxis.set_visible(False)
    ax6.spines['top'].set_visible(False)
    ax6.spines['bottom'].set_visible(False)
    ax6.set(ylabel = 'Trials')
    ax6.yaxis.label.set_size(14)
    
    scale = 1
    scalebar = 1
    yrange = ax6.get_ylim()[1] - ax6.get_ylim()[0]
    scalebary = (yrange / 10) + ax6.get_ylim()[0]
    scalebarx = [ax6.get_xlim()[1] - scalebar, ax6.get_xlim()[1]]
    ax6.plot(scalebarx, [scalebary, scalebary], c='k', linewidth=2)
    ax6.text((scalebarx[0] + (scalebar/2)), scalebary-(yrange/50), str(scale) +' s', ha='center',va='top', **Calibri, **Size) 
    
    rasterPlot = distractionrasterFig(ax6, ratdata['distractors'], ratdata['licks'], pre=1, post=10, sortevents=pdps, sortdirection='dec')

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

for filename in TDTfiles_thph_hab:
    
    file = TDTfilepath + filename
    ratdata = loadmatfile(file)
    allRatBlue.append(ratdata['blue'])
    allRatUV.append(ratdata['uv'])
    allRatFS.append(ratdata['fs'])
    allRatLicks.append(ratdata['licks'])
    allRatDistractors.append(ratdata['distractors'])
    allRatDistracted.append(ratdata['distracted'])
    allRatNotDistracted.append(ratdata['notdistracted'])
    indices1 = []       
    
    for index, value in enumerate(ratdata['distractors']):
        a = np.where(ratdata['licks'] == value) 
        indices1.append(a)       
    
    pdps = []
    for tupl in indices1:
        i = tupl[0][0]
        if i+1 < len(ratdata['licks']):
            
            pdp = (ratdata['licks'][i+1] - ratdata['licks'][i])
            pdps.append(pdp)    
    
   # pdps.append(0) ## Add in a check, if PDPs len does not equal len(distractors)
    figure12 = plt.figure(figsize=(6,6))
    ax6 = plt.subplot(111)
    ax6.spines['right'].set_visible(False)
    ax6.xaxis.set_visible(False)
    ax6.spines['top'].set_visible(False)
    ax6.spines['bottom'].set_visible(False)
    ax6.set(ylabel = 'Trials')
    ax6.yaxis.label.set_size(14)
    
    scale = 1
    scalebar = 1
    yrange = ax6.get_ylim()[1] - ax6.get_ylim()[0]
    scalebary = (yrange / 10) + ax6.get_ylim()[0]
    scalebarx = [ax6.get_xlim()[1] - scalebar, ax6.get_xlim()[1]]
    ax6.plot(scalebarx, [scalebary, scalebary], c='k', linewidth=2)
    ax6.text((scalebarx[0] + (scalebar/2)), scalebary-(yrange/50), str(scale) +' s', ha='center',va='top', **Calibri, **Size) 
    
    rasterPlot = distractionrasterFig(ax6, ratdata['distractors'], ratdata['licks'], pre=1, post=10, sortevents=pdps, sortdirection='dec')


