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
allBursts = []
allRuns = []
allRunIndices = []
allrILIs = []
allbILIs = []

for filename in TDTfiles_thph_lick:
    
    file = TDTfilepath + filename
    ratdata = loadmatfile(file)
    burstanalysis = lickCalc(ratdata['licks'], offset=ratdata['licks_off'])
    burstList = burstanalysis['bLicks'] # type, array 
    runList = burstanalysis['rLicks'] # type array
#    indexRunList = burstanalysis['rInd'] 
    runILIs = burstanalysis['rILIs']
    burstILIs = burstanalysis['bILIs']
    
    allBursts.append(burstList)
    allRuns.append(runList)
#    allRunIndices.append(indexRunList)
    allrILIs.append(runILIs)
    allbILIs.append(burstILIs)
    
# Make the list of lists into one long list for histogram 
MergedBurstList = list(itertools.chain.from_iterable(allBursts)) 
MergedRunList = list(itertools.chain.from_iterable(allRuns)) 
    
# Descriptives - aggregated data
meanburstlength = round(np.mean(MergedBurstList))
medburstlen = round(np.median(MergedBurstList))
meanrunlength = round(np.mean(MergedRunList))
medrunlen = round(np.median(MergedRunList))

##  Make the snips 


# All runs and all bursts (representative rat, etc.)
# Then segregate by long and short (for all rats quartiles not for each rat)















'''

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
         
allign to lowerqRunTimes
allign to uppqRunTimes

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

# Assign empty lists for storing arrays of burst/run lengths? where is this code
allDistracted = []
allNotDistracted = []
allDistractors = []
allRatLicks = []
allRatBlue = []
allRatUV = []
allRatFS = []


# Loop through files and extracts all info needed for snips 
for filename in TDTfiles_thph_dis:
    
    file = TDTfilepath + filename
    ratdata = loadmatfile(file)
    distracted = ratdata['distracted']
    notdistracted = ratdata['notdistracted']
    distractors = ratdata['distractors']
        
    allRatBlue.append(ratdata['blue'])
    allRatUV.append(ratdata['uv'])
    allRatFS.append(ratdata['fs'])
    
    allRatLicks.append(ratdata['licks'])
    allDistracted.append(distracted)
    allNotDistracted.append(notdistracted)
    allDistractors.append(distractors)






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
        
        
# Make the snips 
'''
# Should run different analysis for licking and distraction data 

# Find notes of which figures to have for each photomoetry chapter 

# Calculating different peaks, and means based on on DISTRACTED and NOT DISTRACTED
# To find peaks and averaged activity across time 

# Called means because give 14 arrays which are the MEAN for each rat 
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

#values for averaged activity for the 2 seconds preceeding (14 values, 1 per rat) =

#TwoSecBEFOREactivity = (np.mean(allblueMeans, axis=1)) # axis 1 gives 14 values, axis 0 gives 1 list of 20 (the averaged average snip)

# Want to produce slices of the lists already here
# I have 14 lists of 300. 0 to 100 are the 10 seconds preding distraction
# 100 to 300 are 20 seconds after 

# Split into (1) 80 - 100 (the 2 seconds before distraction)
           # (2) 110 - 300 (the 19 seconds after - to see if supression)

# JUST BLUE ACTIVITY IGNORING THE UV (NO SUBTRACTION AS PRETTY MUCH ZERO)



all2SecBefore = []
all20SecAfter = []
for eachlist in allblueMeans:
    slice1 = eachlist[80:100]
    #print(len(slice1))
    
    slice2 = eachlist[100:300]
    #print(len(slice2))
    
    all2SecBefore.append(slice1)
    all20SecAfter.append(slice2) 
    
    #then out of the loop mean these (axis 1 not 0) = gives 14 points for the barscatter

      
MeanOf2SecDISTRACTED = (np.mean(all2SecBefore, axis=1))
MeanOf20SecDISTRACTED = (np.mean(all20SecAfter, axis=1))

# Repeat for not distracted - reassignes all values 

allblueMeans = []
alluvMeans = []
for i, val in enumerate(allNotDistracted):

        # make a blue and uv snip for all 14
        blueSnips, ppsBlue = snipper(allRatBlue[i], allNotDistracted[i], fs=allRatFS[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUV[i], allNotDistracted[i], fs=allRatFS[i], bins=300)
# # these four lines used later to define means plot (made after runs) 
        # Makes a mean for each rat's snips
        blueMean = np.mean(blueSnips, axis=0)
        allblueMeans.append(blueMean)
        uvMean = np.mean(uvSnips, axis=0)
        alluvMeans.append(uvMean)


all2SecBefore = []
all20SecAfter = []
for eachlist in allblueMeans:
    slice1 = eachlist[80:100]
    #print(len(slice1))
    
    slice2 = eachlist[100:300]
    #print(len(slice2))
    
    all2SecBefore.append(slice1)
    all20SecAfter.append(slice2) 
    
MeanOf2SecNOT_DISTRACTED = (np.mean(all2SecBefore, axis=1))
MeanOf20SecNOT_DISTRACTED = (np.mean(all20SecAfter, axis=1))

PEAK2SEC = np.empty((2,), dtype=np.object)
PEAK2SEC[0] = MeanOf2SecDISTRACTED 
PEAK2SEC[1] = MeanOf2SecNOT_DISTRACTED

peak20sec = np.empty((2,), dtype=np.object)
peak20sec[0] = MeanOf20SecDISTRACTED
peak20sec[1] = MeanOf20SecNOT_DISTRACTED


'''        
        
        
        