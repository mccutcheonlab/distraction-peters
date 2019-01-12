#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 17:37:57 2019

@author: kate
"""

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
    figure12 = plt.figure(figsize=(6,3))
    ax6 = plt.subplot(111)
    ax6.spines['right'].set_visible(False)
    ax6.xaxis.set_visible(False)
    ax6.spines['top'].set_visible(False)
    ax6.spines['bottom'].set_visible(False)
    ax6.set(ylabel = 'Trials')
    ax6.yaxis.label.set_size(14)
#    distractionrasterFig(ax6, ratdata['distractors'], ratdata['licks'], pre=1, post=10, sortevents='yes', sortdirection='dec')
    

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



# RAT3 - 10 minutes  (value zero is rat3 here as 1 and 2 deleted in distraciton)
fig9 = plt.figure(figsize=(12,2))
ax7 = plt.subplot(1,1,1)
plt.plot(allRatBlue[10], color='royalblue')
plt.plot(allRatUV[10] + 150, color='darkorchid') ### OFFSET THE UV 
ax7.set_xticks([0,(10*60*allRatFS[0]),(20*60*allRatFS[0]),(30*60*allRatFS[0]),(40*60*allRatFS[0]),(50*60*allRatFS[0]),(60*60*allRatFS[0])] )
ax7.set_xticklabels([0,10,20,30,40,50,60])
ax7.set_xlabel('Mins', fontsize=14)
#ax7.set_xlim([500000,700000]) # looks really nice scale wise, approx 3 mins

## FIRSRT 10 MINS
#ax7.set_xlim([122070.31494140625,732421.8896484375]) # 2 mins to 12 mins, a 10 min snip without noise at start
# Second 10 mins 
ax7.set_xlim([122070.31494140625,1220703.1494140625]) # 2 mins to 12 mins, a 10 min snip without noise at start

ax7.set_ylim([500,900])

multipliedLicks = []
for element in allRatLicks[10]:
    multElement = element*allRatFS[0]
    multipliedLicks.append([multElement])
    
multipliedDistractors = []
for element in allRatDistractors[10]:
    multElement = element*allRatFS[0]
    multipliedDistractors.append([multElement])  
    
multipliedDistracted = []
for element in allRatDistracted[10]:
    multElement = element*allRatFS[0]
    multipliedDistracted.append([multElement])  
    
xvals = multipliedLicks
yvals = [ax7.get_ylim()[1] - 10] * len(xvals)
ax7.scatter(xvals, yvals, marker='|', color='k', linewidth=0.2)
ax7.set(ylabel = '∆F')
ax7.yaxis.label.set_size(14)
ax7.xaxis.set_visible(False)


xvals = multipliedDistractors
yvals = [ax7.get_ylim()[1] - 40] * len(xvals)
ax7.scatter(xvals, yvals, marker='o', color='k', linewidth=0.2)
ax7.set(ylabel = '∆F')
ax7.yaxis.label.set_size(14)
ax7.xaxis.set_visible(False)

xvals = multipliedDistracted
yvals = [ax7.get_ylim()[1] - 40] * len(xvals)
ax7.scatter(xvals, yvals, marker='o', color='r', linewidth=0.2)
ax7.set(ylabel = '∆F')
ax7.yaxis.label.set_size(14)
ax7.xaxis.set_visible(False)


            
scalebar = 1*allRatFS[0]*60 # 1 minute

yrange = ax7.get_ylim()[1] - ax7.get_ylim()[0]
scalebary = (yrange / 10) + ax7.get_ylim()[0]
scalebarx = [ax7.get_xlim()[1] - scalebar, ax7.get_xlim()[1]]
ax7.plot(scalebarx, [scalebary, scalebary], c='k', linewidth=2)
ax7.text((scalebarx[0] + (scalebar/2)), scalebary-(yrange/50), '1 Min', ha='center',va='top', **Calibri, **Size)
ax7.spines['right'].set_visible(False)
ax7.spines['top'].set_visible(False)
ax7.spines['bottom'].set_visible(False)
#fig9.savefig('/Volumes/KPMSB352/PHOTOMETRY MMIN18/PDF figures/LongTimeCourse.pdf', bbox_inches="tight") 
fig9.savefig('/Users/kate/Desktop/Peters, McCutcheon & Young, 2019/Draft 1/LongTimeCourseDISRat10_20min.pdf', bbox_inches="tight") 
       


### Short time course 

# RAT3 - 10 minutes  (value zero is rat3 here as 1 and 2 deleted in distraciton)
fig9 = plt.figure(figsize=(4,2))
ax7 = plt.subplot(1,1,1)
plt.plot(allRatBlue[10], color='royalblue')
plt.plot(allRatUV[10] + 160, color='darkorchid') ### OFFSET THE UV 
ax7.set_xticks([0,(10*60*allRatFS[0]),(20*60*allRatFS[0]),(30*60*allRatFS[0]),(40*60*allRatFS[0]),(50*60*allRatFS[0]),(60*60*allRatFS[0])] )
ax7.set_xticklabels([0,10,20,30,40,50,60])
ax7.set_xlabel('Mins', fontsize=14)
#ax7.set_xlim([500000,700000]) # looks really nice scale wise, approx 3 mins

## FIRSRT 10 MINS
#ax7.set_xlim([122070.31494140625,732421.8896484375]) # 2 mins to 12 mins, a 10 min snip without noise at start
# Second 10 mins 
ax7.set_xlim([396728.5235595703,518798.83850097656]) # 2 mins to 12 mins, a 10 min snip without noise at start

ax7.set_ylim([500,800])

multipliedLicks = []
for element in allRatLicks[10]:
    multElement = element*allRatFS[0]
    multipliedLicks.append([multElement])
    
multipliedDistractors = []
for element in allRatDistractors[10]:
    multElement = element*allRatFS[0]
    multipliedDistractors.append([multElement])  
    
multipliedDistracted = []
for element in allRatDistracted[10]:
    multElement = element*allRatFS[0]
    multipliedDistracted.append([multElement])  
    
xvals = multipliedLicks
yvals = [ax7.get_ylim()[1] - 10] * len(xvals)
ax7.scatter(xvals, yvals, marker='|', color='k', linewidth=0.2)
ax7.set(ylabel = '∆F')
ax7.yaxis.label.set_size(14)
ax7.xaxis.set_visible(False)


xvals = multipliedDistractors
yvals = [ax7.get_ylim()[1] - 40] * len(xvals)
ax7.scatter(xvals, yvals, marker='o', color='k', linewidth=0.2)
ax7.set(ylabel = '∆F')
ax7.yaxis.label.set_size(14)
ax7.xaxis.set_visible(False)

xvals = multipliedDistracted
yvals = [ax7.get_ylim()[1] - 40] * len(xvals)
ax7.scatter(xvals, yvals, marker='o', color='r', linewidth=0.2)
ax7.set(ylabel = '∆F')
ax7.yaxis.label.set_size(14)
ax7.xaxis.set_visible(False)


            
scalebar = 1*allRatFS[0]*60 # 1 minute

yrange = ax7.get_ylim()[1] - ax7.get_ylim()[0]
scalebary = (yrange / 10) + ax7.get_ylim()[0]
scalebarx = [ax7.get_xlim()[1] - scalebar, ax7.get_xlim()[1]]
ax7.plot(scalebarx, [scalebary, scalebary], c='k', linewidth=2)
ax7.text((scalebarx[0] + (scalebar/2)), scalebary-(yrange/50), '1 Min', ha='center',va='top', **Calibri, **Size)
ax7.spines['right'].set_visible(False)
ax7.spines['top'].set_visible(False)
ax7.spines['bottom'].set_visible(False)
#fig9.savefig('/Volumes/KPMSB352/PHOTOMETRY MMIN18/PDF figures/LongTimeCourse.pdf', bbox_inches="tight") 
fig9.savefig('/Users/kate/Desktop/Peters, McCutcheon & Young, 2019/Draft 1/LongTimeCourseDISRat10_2min.pdf', bbox_inches="tight") 
       


## then individual trials 


## Converts snips to z-scores based on baseline period before event 
# Default is 10 seconds as snips set up to be 30 seconds with 10 min 
# before the event 
    
def zscore(snips, baseline_points=100):

   BL_range = range(baseline_points)
   z_snips = []
   for i in snips:
       mean = np.mean(i[BL_range])
       sd = np.std(i[BL_range])
       z_snips.append([(x-mean)/sd for x in i])

   return z_snips    





