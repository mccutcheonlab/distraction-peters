#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 11:03:30 2019

@author: kate
"""

def PhotoPeaksCalc(snips_all_rats):
    
    allRat_peak = []
    allRat_t = []
    allRat_pre = []
    allRat_post = []
    allRat_base = []
    
    allRatten_percent_baseline= []
    
    for rat in snips_all_rats:
        pre_event = np.mean(rat[0:50]) # Average for 5 seconds, 10 seconds before event 
        absolutepeak = np.max(rat[100:300]) ## Added this - finds highest (the t is first)
        peak = np.max(rat[100:130]) # to 300 in origianl (can this work with both?) ## Minus the average of the first 5 seconds and after 100 points (slice)
        peak_range = rat[100:130]
        a = peak_range.tolist()
        peak_index = a.index(peak) 
        t = peak_index / 10
        pre_event = np.mean(rat[50:100])
        post_event = np.mean(rat[100:300])
        baseline = np.mean(rat[0:50])
        
        allRat_peak.append(absolutepeak)
        allRat_t.append(t)
        allRat_pre.append(pre_event)
        allRat_post.append(post_event)
        allRat_base.append(baseline)   
        
# Calculate the slope (time to decay to 10% of baseline)
        
        total_points = 0
        pointindices = []
        for index, point in enumerate(rat[peak_index:300]):
            if point < baseline*0.9:
                pointindices.append(index)
                total_points += 1
        allRatten_percent_baseline.append(pointindices[0])  # the first time there is a point below 10% 
    
    
    return allRat_peak, allRat_t, allRat_pre, allRat_post, allRat_base, allRatten_percent_baseline


# Distractors
peak_distractor, t_distractor, pre_distractor, post_distractor, baseline_distractor, slope_distractor = PhotoPeaksCalc(bkgnd_sub_Distractor)
peak_distracted, t_distracted, pre_distracted, post_distracted, baseline_distracted, slope_distracted = PhotoPeaksCalc(bkgnd_sub_Distracted)
# Not distracted 
peak_notdistracted, t_notdistracted, pre_notdistracted, post_notdistracted, baseline_notdistracted, slope_notdistracted = PhotoPeaksCalc(bkgnd_sub_Notdistracted)


print(np.mean(slope_distractor))
print(np.mean(slope_distracted))
print(np.mean(slope_notdistracted))
