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
        else:
            if sortdirection == 'ascending':
                sortOrder = np.argsort(sortevents)
            else:
                sortOrder = np.argsort(sortevents)[::-1]
                
            timelock = [timelock[i] for i in sortOrder]
    
    rasterData = [[] for i in timelock]
    
    for i,x in enumerate(timelock):
        rasterData[i] = [j-x for j in events if (j > x-pre) & (j < x+post)]

#    for ith, trial in enumerate(rasterData):
#        if ith < 26:
# 
#            ax.vlines(trial, ith + .5, ith + 1.5)
#        else:
#            ax.vlines(trial, ith + .5, ith + 1.5, color='blue')
#            
#            
    for ith, trial in enumerate(rasterData): 
        xvals = [x for x in trial]  
        yvals = [1+ith] * len(xvals)
        
        #if ith<40:  # 26 if ascending
            #ax.scatter(xvals, yvals, marker='.', color='b')
            
        #else:
         #   ax.scatter(xvals, yvals, marker='.', color='k')
        ax.scatter(xvals, yvals, marker=',',s=2,color='k')
        
        
       
       
# produces the index in the lick data where the distractor was (indices1)
# now use these indices to add one and subtract the VALUE at index+1 from the VALUE at index
indices1 = []       
for index, value in enumerate(examplerat['distractors']):
    a = np.where(examplerat['licks'] == value) 
    indices1.append(a)       

pdps = []
for tupl in indices1:
    i = tupl[0][0]
    if i+1 < len(examplerat['licks']):
        
        pdp = (examplerat['licks'][i+1] - examplerat['licks'][i])
        pdps.append(pdp)    

# Check the PDPs first one is very long?Yes it is 
#
pdps.append(0)        
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

rasterPlot = distractionrasterFig(ax6, examplerat['distractors'], examplerat['licks'], pre=1, post=10, sortevents=pdps, sortdirection='dec')
#figure12.savefig('/Volumes/KPMSB352/PHOTOMETRY MMIN18/PDF figures/RasterLickDay2.3.pdf') 

#
#
#Make plots for every rat on lick day, distraction day and habituation day
#Choose the best lookng representative rat for the ordered, non ordered (by time)
#Make the dots where it is distracted a different colour (this was added into a previous script)
#Sort out the spacing issues (talk to Jaime perhaps) - replace with a different marker or something 
#    Could make the plots as large as possible 
#    Think about issues of different scales on (1) Different days and (2) Different rats 
#    Consider using ',' rather than '.' (a pixel and not a dot)
#    
#   s = 1  
    
    