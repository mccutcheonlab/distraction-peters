#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 17:29:57 2018

@author: u1490431
"""



## these are created in RATSER PLOT 
##pdp_day_list = [pdps_lickday, pdps_disday, pdps_habday]

## PDPS on licking day 

fig = plt.figure()
#plt.title('Lickday SAL_M', **Calibri, **Size)
ax = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

for index, licklist in enumerate(pdps_lickday):
    plot = cumulativelickFig(ax, pdps_lickday[index], normed=True, color='lightgrey', log=True)
avg = [item for rat in pdps_lickday for item in rat] 
cumulativelickFig(ax, avg, normed=True, color='dimgrey', log=True)
ax.set(ylabel = 'Probability')
ax.yaxis.label.set_size(16)
ax.set(xlabel = 'post-distraction pause log(s)')
ax.xaxis.label.set_size(16)
fig.savefig('/Volumes/KP_HARD_DRI/distraction_paper/figs/CumulativeLickDay.pdf', bbox_inches='tight') 


## PDPS on distraction day 

fig = plt.figure()
#plt.title('Lickday SAL_M', **Calibri, **Size)
ax = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

for index, licklist in enumerate(pdps_disday):
    plot = cumulativelickFig(ax, pdps_disday[index], normed=True, color='lightgrey', log=True)
avg = [item for rat in pdps_disday for item in rat] 
cumulativelickFig(ax, avg, normed=True, color='darkblue', log=True)
ax.set(ylabel = 'Probability')
ax.yaxis.label.set_size(16)
ax.set(xlabel = 'post-distraction pause log(s)')
ax.xaxis.label.set_size(16)
fig.savefig('/Volumes/KP_HARD_DRI/distraction_paper/figs/CumulativeDisDay.pdf', bbox_inches='tight') 



## PDPS on habituation day 

fig = plt.figure()
#plt.title('Lickday SAL_M', **Calibri, **Size)
ax = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

for index, licklist in enumerate(pdps_habday):
    plot = cumulativelickFig(ax, pdps_habday[index], normed=True, color='lightgrey', log=True)
avg = [item for rat in pdps_habday for item in rat] 
cumulativelickFig(ax, avg, normed=True, color='green', log=True)
ax.set(ylabel = 'Probability')
ax.yaxis.label.set_size(16)
ax.set(xlabel = 'post-distraction pause log(s)')
ax.xaxis.label.set_size(16)
fig.savefig('/Volumes/KP_HARD_DRI/distraction_paper/figs/CumulativeHabDay.pdf', bbox_inches='tight') 


## Plot all three lines on one plot (NOT AVERAGES, but the cPDPs for rat 1.4, index 3)
fig = plt.figure()
#plt.title('Lickday SAL_M', **Calibri, **Size)
ax = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

cumulativelickFig(ax, pdps_lickday[3], normed=True, color='dimgrey', log=True)
cumulativelickFig(ax, pdps_disday[1], normed=True, color='green', log=True)
cumulativelickFig(ax, pdps_habday[1], normed=True, color='black', log=True)
ax.set(ylabel = 'Probability')
ax.yaxis.label.set_size(16)
ax.set(xlabel = 'post-distraction pause log(s)')
ax.xaxis.label.set_size(16)
fig.savefig('/Volumes/KP_HARD_DRI/distraction_paper/figs/Cumulative1.4Rep.pdf', bbox_inches='tight') 

make the long course potometry plots with the zoomed 

