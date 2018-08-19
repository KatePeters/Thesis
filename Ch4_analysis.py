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
import seaborn as sb
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
TDTfiles_thph_dis = ['thph1-1_distraction1', 'thph1-2_distraction1','thph1-3_distraction1', \
                'thph1-4_distraction1','thph1-5_distraction1','thph1-6_distraction1', \
                'thph2-1_distraction', 'thph2-2_distraction', 'thph2-3_distraction', \
                'thph2-4_distraction', 'thph2-5_distraction', 'thph2-6_distraction', \
                'thph2-7_distraction', 'thph2-8_distraction']

# Habituation day 
TDTfiles_thph_hab = ['thph1-1_distraction2','thph1-2_distraction2','thph1-3_distraction2',\
                'thph1-4_distraction2', 'thph1-5_distraction2','thph1-6_distraction2',\
                'thph2-1_habituation', 'thph2-2_habituation', 'thph2-3_habituation', \
                'thph2-4_habituation', 'thph2-5_habituation', 'thph2-6_habituation', \
                'thph2-7_habituation']


TDTfilepath = '/Volumes/KP_HARD_DRI/All_Matlab_Converts/BIG CONVERSION 14 AUG 2018/THPH matfiles/'


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
    
    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    ax.set_ylim([-0.03, 0.03])
    #ax.set_ylim([-0.05, 0.05])
    trialsMultShadedFig(ax, [uvSnips,blueSnips], ppsBlue, eventText='First Lick in Run')
    plt.text(250,0.03, '{}'.format(len(allRunTimes[i])) + ' Runs' )
    
    fig2 = plt.figure()
    ax2 = plt.subplot(1,1,1)
    ax2.set_ylim([-0.2, 0.2])
    trialsFig(ax2, blueSnips, uvSnips, ppsBlue, eventText='First Lick in Run', noiseindex=noiseindex) #, )
    plt.text(250,0.2, '{}'.format(len(allRunTimes[i])) + ' Runs' )

 # these four lines used later to define means plot (made after runs)
    blueMean = np.mean(blueSnips, axis=0)
    blueMeansRuns.append(blueMean)
    uvMean = np.mean(uvSnips, axis=0)
    uvMeansRuns.append(uvMean)
    
# All runs and all bursts (representative rat, etc.)
# Then segregate by long and short (for all rats quartiles not for each rat)

# Get the correct allignment from the snipper, that uses the times taken from where? 
uvMeans_all_run = []
blueMeans_all_run = []

## Mean of ALL runs and ALL rats on multishaded figure

#linecolor=['purple', 'blue'], errorcolor=['thistle', 'lightblue']

fig = plt.figure()
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.03, 0.03])
#ax.set_ylim([-0.05, 0.05])
trialsMultShadedFig(ax, [np.asarray(uvMeansRuns),np.asarray(blueMeansRuns)], ppsBlue, eventText='First Lick in Run', linecolor = ['black','green'], errorcolor = ['lightgray','lightgreen'])
plt.text(250,0.03, '{}'.format(len(allRunTimes[i])) + ' Runs' ) ## Edit this to be all


'''

# Makes tonnes of individual plots             
# Individual rats might look odd as low Ns, but mean of mean will be better 
# Repeat this with bursts if looks like might be useful 

for i, val in enumerate(allRuns):
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
#    
#    fig13 = plt.figure()
#    ax11 = plt.subplot(1,1,1)
#    ax11.set_ylim([-0.2, 0.2])
#    trialsFig(ax11, blueSnips, uvSnips, ppsBlue, eventText='First Lick in Short Run', noiseindex=noiseindex) #, )
#    plt.text(250,0.2, '{}'.format(len(lowerqRunTimes[i])) + ' short runs' )

# # these four lines used later to define means plot (made after runs)
    blueMeanALL = np.mean(blueSnips, axis=0)
    blueMeans_all_run.append(blueMeanSHORT)
    uvMeanALL = np.mean(uvSnips, axis=0)
    uvMeans_all_run.append(uvMeanSHORT)
















## Find short and long run lengths 
aggregateLowerQuart = np.percentile(MergedRunList,25)
aggregateUpperQuart = np.percentile(MergedRunList,75) 
    
allLogIndRuns = []
for runLicksList in allRuns:
    logIndRuns = [] # so 14 empty arrays exist and then don't (add to larger)
    for item in runLicksList:
        if item < aggregateLowerQuart:
            logIndRuns.append('UPPER')
        else:
            if item > aggregateUpperQuart:
                logIndRuns.append('LOWER')
            
            else:
                logIndRuns.append('MIDDLE')
    allLogIndRuns.append(logIndRuns)
            
# 14 lists of logical indices for whether the n licks in that run was L,M,U

 
# Final lists of run times sorted by u.m.l and stored by rat (14 lists)
lowerqRunTimes = []
uppqRunTimes = []
mid50RunTimes = []

# Might need the index and the value for both lists ??
for i, listofindices in enumerate(allLogIndRuns):
    templower = []
    tempupper = []
    tempmid = []
    for j, runIndex in enumerate(listofindices):
        
        #print(i,j) i is the list index and j is the item index
        
        if runIndex == 'LOWER':
            lowrun = allRunsTimes[i][j]
            templower.append(lowrun)
               
        if runIndex == 'UPPER':
            upprun = allRunsTimes[i][j]
            tempupper.append(upprun)
#                
        if runIndex == 'MIDDLE':
            midrun = allRunsTimes[i][j]
            tempmid.append(midrun)
            
    lowerqRunTimes.append(templower)
    uppqRunTimes.append(tempupper)
    mid50RunTimes.append(tempmid)
    
# _________________________________________________________________________

    
''' 
         
#allign to lowerqRunTimes
#allign to uppqRunTimes

# Then, if not too complex, add BOTH BLUE SIGNALS high and low burst numbers to the 
    # same photometry plot 

'''

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
#    
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


## ===============================================================

# Photometry figures (individual) for LONG bursts 

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
'''


'''



# Called means because give 14 arrays which are the MEAN for each rat 

# Probably want to store the blueSnips and uvSnips so that representative RAT can be plotted
allblueMeans = []
alluvMeans = []
for i, val in enumerate(allDistracted):

        # make a blue and uv snip for all 14
        blueSnips, ppsBlue = snipper(allRatBlue[i], allDistracted[i], fs=allRatFS[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUV[i], allDistracted[i], fs=allRatFS[i], bins=300)
# # these four lines used later to define means plot (made after runs) 
        # Makes a mean for each rat's snips
        blueMean = np.mean(blueSnips, axis=0)
        allblueMeans.append(blueMean)
        uvMean = np.mean(uvSnips, axis=0)
        alluvMeans.append(uvMean)
        
        

'''        