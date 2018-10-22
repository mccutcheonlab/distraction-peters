#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 16:08:19 2018

@author: u1490431
"""

import JM_general_functions as jmf
import JM_custom_figs as jmfig
import distraction_sessionfigs as sessionfigs
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import scipy.io as sio
import dill


## Install DILL --> pip (source activate)

''' Makes a txt file from the excel file
    Reads in as row and header
'''

class Session(object):
# each attribute for each "session" --> one row 
    
    def __init__(self, sessionID, metafiledata, hrows, datafolder, outputfolder):
        self.sessionID = sessionID
#        self.sessiontype = metafiledata[hrows['stype']]
        self.medfile = metafiledata[hrows['medfile']]
        self.rat = metafiledata[hrows['rat']].replace('.', '-')
        self.session = metafiledata[hrows['session']]
        self.box = metafiledata[hrows['box']] 
#        self.bottleL = metafiledata[hrows['bottleL']]
#        self.bottleR = metafiledata[hrows['bottleR']]
        
        self.ttl_licks = metafiledata[hrows['ttl-licks']]
        self.ttl_distractors = metafiledata[hrows['ttl-distractors']]
        self.ttl_distracted = metafiledata[hrows['ttl-distracted']]

        self.matlabfile = datafolder + self.rat + '_' + self.session + '.mat'   
        self.outputfolder = outputfolder
        
    def loadmatfile(self):
        a = sio.loadmat(self.matlabfile, squeeze_me=True, struct_as_record=False) 
        self.output = a['output']
        self.fs = self.output.fs
        self.data = self.output.blue
        self.dataUV = self.output.uv




        
def metafile2sessions(xlfile, metafile, datafolder, outputfolder, sheetname='metafile'):
    jmf.metafilemaker(xlfile, metafile, sheetname=sheetname, fileformat='txt')
    rows, header = jmf.metafilereader(metafile + '.txt')
    
    hrows = {}
    for idx, field in enumerate(header):
        hrows[field] = idx
    
    sessions = {}

# Turns each row into a session object - contains all the information 
    
    for row in rows:
        # Doesnt like dots 
        sessionID = row[hrows['rat']].replace('.','-') + '_' + row[hrows['session']]
        sessions[sessionID] = Session(sessionID, row, hrows, datafolder, outputfolder)
    
    return sessions   



def assemble_sessions(sessions,
                      rats_to_include=[],
                      rats_to_exclude=[],
                      sessions_to_include=[],
                      outputfile=[],
                      savefile=False,
                      makefigs=False):
    
    rats = []
    for session in sessions:
        x = sessions[session]
        if x.rat not in rats:
            rats.append(x.rat)

    if len(rats_to_include) > 0:
        print('Overriding values in rats_to_exclude because of entry in rats_to_include.')
        rats_to_exclude = list(rats)
        for rat in rats_to_include:
            rats_to_exclude.remove(rat)
    
    sessions_to_remove = []
    
    for session in sessions:
          
        x = sessions[session]
        
        if x.rat not in rats_to_exclude and x.session in sessions_to_include: 

             print('hello')    
             x.loadmatfile()
             try:
                x.loadmatfile()
        
                print('\nAnalysing rat ' + x.rat + ' in session ' + x.session)


             except:
                print('Could not extract data from ' + x.sessionID)
            
                
#                # Load in data from .mat file (convert from Tank first using Matlab script)
#                x.loadmatfile()
#                # Work out time to samples
#                x.setticks()
#                x.time2samples()       
#                # Find out which bottles have TTLs/Licks associated with them     
#                x.check4events()
#    
#                x.setbottlecolors()
#                try:
#                    x.left['lickdata'] = jmf.lickCalc(x.left['licks'],
#                                      offset = x.left['licks_off'],
#                                      burstThreshold = 0.50)
#                except IndexError:
#                    x.left['lickdata'] = 'none'
#                    
#                try:
#                    x.right['lickdata'] = jmf.lickCalc(x.right['licks'],
#                              offset = x.right['licks_off'],
#                              burstThreshold = 0.50)
#                except IndexError:
#                    x.right['lickdata'] = 'none'
#                
#                bins = 300
#                
#                x.randomevents = jmf.makerandomevents(120, max(x.tick)-120)
#                x.bgTrials, x.pps = jmf.snipper(x.data, x.randomevents,
#                                                t2sMap = x.t2sMap, fs = x.fs, bins=bins)
#                
#                for side in [x.left, x.right]:   
#                    if side['exist'] == True:
#                        side['snips_sipper'] = jmf.mastersnipper(x, side['sipper'])
#                        side['snips_licks'] = jmf.mastersnipper(x, side['lickdata']['rStart'])
#                        try:
#                            timelock_events = [licks for licks in side['lickdata']['rStart'] if licks in side['licks-forced']]
#                            latency_events = side['sipper']
#                            side['snips_licks_forced'] = jmf.mastersnipper(x, timelock_events, latency_events=latency_events, latency_direction='pre')
#                        except KeyError:
#                            pass
#                        try:
#                            side['lats'] = jmf.latencyCalc(side['lickdata']['licks'], side['sipper'], cueoff=side['sipper_off'], lag=0)
#                        except TypeError:
#                            print('Cannot work out latencies as there are lick and/or sipper values missing.')
#                            side['lats'] = []
#                x.side2subs()
#                 
#                if makefigs == True:
#                    pdf_pages = PdfPages(x.outputfolder + session + '.pdf')
#                    sessionfigs.makeBehavFigs(x, pdf_pages)
#                    sessionfigs.makePhotoFigs(x, pdf_pages)
#            except:
#                print('Could not extract data from ' + x.sessionID)
#            
#            try:
#                pdf_pages.close()
#                plt.close('all')
#            except:
#                print('Nothing to close')
#        else:
#            sessions_to_remove.append(session)
#    
#    for session in sessions_to_remove:
#        sessions.pop(session)
#        
#    for rat in rats_to_exclude:
#        idx = rats.index(rat)
#        del rats[idx]
#    
#    if savefile == True:
#        pickle_out = open(outputfile, 'wb')
#        dill.dump([sessions, rats], pickle_out)
#        pickle_out.close()
        
    return sessions

