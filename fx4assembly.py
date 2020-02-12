# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 12:39:35 2020

@author: admin
"""
import numpy as np
import scipy.signal as sig

def correctforbaseline(blue, uv):
    pt = len(blue)
    X = np.fft.rfft(uv, pt)
    Y = np.fft.rfft(blue, pt)
    Ynet = Y-X

    datafilt = np.fft.irfft(Ynet)

    datafilt = sig.detrend(datafilt)

    b, a = sig.butter(9, 0.012, 'low', analog=True)
    datafilt = sig.filtfilt(b, a, datafilt)
    
    return datafilt

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
            
    if d[-1] > 3599:
        d = d[:-1]
    
    return d

def distracted_or_not(distractors, licks, delay=1):   
    pdp = [] # post-distraction pause
    distractedArray = []

    for d in distractors:               
        try:
            pdp.append([i-d for i in licks if (i > d)][0])
        except IndexError:
            pdp.append(3600-d) # designates end of session as max pdp

    distracted_boolean_array = np.array([i>delay for i in pdp], dtype=bool)
    
    distracted = [d for d,l in zip(distractors, distracted_boolean_array) if l]
    notdistracted = [d for d,l in zip(distractors, distracted_boolean_array) if not l] 
    
    if np.isnan(pdp)[-1] == 1: 
        distracted_boolean_array[-1] = True
    
    return [distracted, notdistracted], distracted_boolean_array, pdp