#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 19:09:07 2019

@author: kate
"""

## Functions and code to creat the long time course figures for distraction paper

# Step 1 - run the all functions file 
# Step 2 - run the Ch4_analysis_licking and Ch4_analysis_distraction files 
# Step 3 - change the snipper function to make longer snips and add this to this file
# Step 4 - copy over only the code that is needed from the licking and distraciton file 
#        - to make the long time course figures 

#! ! ! Decide if you are going to use data from one representative rats or all of the animals 
        # Check back to MMIN18 and see whether you used group or individual data 
        # Maybe only need the blue signal or maybe also violet is useful 

# Step 5 - make and save long timecourse figures 
        
        
# Step 6 - subset the required data for the barscatter plots, decide which info is important 
        
# Step 7 - calculate peaks within time region and look at these statistically 
        
        
# Step 8 - check which figures and statistics from thesis / and corrections (time bins)
#        - will be used

# Step 9 - make barscatter of 'probability' of distraction in 10 min time bins 
#        - averaged for all of the rats (make the percent for each and then average)
#        - show statistically the lack of habituaiton effect within the session 
        
        

## Step 10 - is it possible to make barscatters for average post distraction pauses 
#          - or other measure median then bin these into 10 mins as with percent distracted?
        
        
### Could plot all the rats as 1hr then as 30 mins to see which rat to use for long time course 
        
# RAT3

fig9 = plt.figure(figsize=(12,2))
ax7 = plt.subplot(1,1,1)
plt.plot(allRatBlue[3], color='royalblue')
plt.plot(allRatUV[3], color='darkorchid')
ax7.set_xticks([0,(10*60*allRatFS[0]),(20*60*allRatFS[0]),(30*60*allRatFS[0]),(40*60*allRatFS[0]),(50*60*allRatFS[0]),(60*60*allRatFS[0])] )
ax7.set_xticklabels([0,10,20,30,40,50,60])
ax7.set_xlabel('Mins', fontsize=14)
#ax7.set_xlim([500000,700000]) # looks really nice scale wise, approx 3 mins
ax7.set_xlim([122070.31494140625,732421.8896484375]) # 2 mins to 12 mins, a 10 min snip without noise at start
ax7.set_ylim([400,800])



# Adding the scatter to long time course plot of photo signals
#allRatLicks.append(ratdata['licks'])
multipliedLicks = []
for element in allRatLicks[3]:
    multElement = element*allRatFS[0]
    multipliedLicks.append([multElement])
    
xvals = multipliedLicks
yvals = [ax7.get_ylim()[1] - 100] * len(xvals)
ax7.scatter(xvals, yvals, marker='|', color='k', linewidth=0.2)

# Get rid of the spines and add labels and ticks to plot 
# Add a 1 minute scale bar OR tick labels for mins 
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
        

# RAT 12
fig9 = plt.figure(figsize=(12,2))
ax7 = plt.subplot(1,1,1)
plt.plot(allRatBlue[12], color='royalblue')
plt.plot(allRatUV[12], color='darkorchid')
ax7.set_xticks([0,(10*60*allRatFS[0]),(20*60*allRatFS[0]),(30*60*allRatFS[0]),(40*60*allRatFS[0]),(50*60*allRatFS[0]),(60*60*allRatFS[0])] )
ax7.set_xticklabels([0,10,20,30,40,50,60])
ax7.set_xlabel('Mins', fontsize=14)
#ax7.set_xlim([500000,700000]) # looks really nice scale wise, approx 3 mins
ax7.set_xlim([122070.31494140625,732421.8896484375]) # 2 mins to 12 mins, a 10 min snip without noise at start
ax7.set_ylim([300,800])



# Adding the scatter to long time course plot of photo signals
#allRatLicks.append(ratdata['licks'])
multipliedLicks = []
for element in allRatLicks[12]:
    multElement = element*allRatFS[0]
    multipliedLicks.append([multElement])
    
xvals = multipliedLicks
yvals = [ax7.get_ylim()[1] - 50] * len(xvals)
ax7.scatter(xvals, yvals, marker='|', color='k', linewidth=0.2)

# Get rid of the spines and add labels and ticks to plot 
# Add a 1 minute scale bar OR tick labels for mins 
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
        
