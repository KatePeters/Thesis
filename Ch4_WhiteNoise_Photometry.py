#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 09:57:41 2018

@author: u1490431
"""
### NEED TO RUN PARTS OF CHAPTER 4 FIRST (analysis_licking) to read meta 
## Need JEM_Figures and FUnctions imports too 

# NECESSARY IMPORTS AND FUNCTIONS 

# Imports
import numpy as np
from scipy import stats
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib as mpl 

#import seaborn as sn
def uvSubtractor(rat_snip_means_list_blue, uv):
    subtractedSignal = []
    for ind, rat in enumerate(rat_snip_means_list_blue):
        subtractedSignal.append(rat_snip_means_list_blue[ind] - uv[ind])
        
    return subtractedSignal

def PhotoPeaksCalc(snips_all_rats):
    
    allRat_peak = []
    allRat_t = []
    allRat_pre = []
    allRat_post = []
    allRat_base = []
    
    for rat in snips_all_rats:
        pre_event = np.mean(rat[0:50]) # Average for 5 seconds, 10 seconds before event 
        peak = np.max(rat[100:130]) ## Minus the average of the first 5 seconds and after 100 points (slice)
        peak_range = rat[100:130]
        a = peak_range.tolist()
        peak_index = a.index(peak)
        t = peak_index / 10
        pre_event = np.mean(rat[50:100])
        post_event = np.mean(rat[100:300])
        baseline = np.mean(rat[0:50])
        
        allRat_peak.append(peak)
        allRat_t.append(t)
        allRat_pre.append(pre_event)
        allRat_post.append(post_event)
        allRat_base.append(baseline) 
        
    
    return allRat_peak, allRat_t, allRat_pre, allRat_post, allRat_base

def MetaExtractor (metafile):
    f = open(metafile, 'r')
    f.seek(0)
    Metafilerows = f.readlines()[1:]
    tablerows = []

    for row in Metafilerows: 
        items = row.split(',')
        tablerows.append(items)

    MedFilenames, RatID, Date, Day, Session, Drug, TotLicks, Distractions, \
    NonDistractions, PercentDistracted = [], [], [], [], [], [], [], [], [], []

    for i, lst in enumerate(tablerows):
       MedFilenames = MedFilenames + [lst[0]]
       RatID = RatID + [lst[1]]
       Date = Date + [lst[2]]
       Session = Session + [lst[3]]
       Drug = Drug + [lst[4]]
       TotLicks = TotLicks + [lst[5]]
#       Distractions = Distractions + [lst[7]] 
#       NonDistractions = NonDistractions + [lst[8]]
       PercentDistracted = PercentDistracted + [lst[8]]
    return ({'MedFilenames':MedFilenames, 'RatID':RatID, 'Date':Date, 'Session':Session, \
             'Drug':Drug, 'TotLicks':TotLicks})



def subsetter(dictionary, dates, drug, dis=False, verbose=False):
    '''
    SUBSETTER KP
    # Subsets data according to date, reads in dictionnary produced from metafile
    # and subsets into variable based on date(s) and drug condition 
    # if distraction day argument is given as True adds the distractor type 
    # to the output lists for later processing 
    
    '''
    subset = []
    for ind, filename in enumerate(dictionary['MedFilenames']):
        path = medfolder + filename
        onsets, offsets, med_dis_times, dis_type = medfilereader(path, ['b', 'c', 'i', 'j'], remove_var_header = True)  # e onset, f offset

        if dis == True:
            if dictionary['Date'][ind] in dates and dictionary['Drug'][ind] == drug:
                subset.append([onsets, offsets, dis_type, dictionary['RatID'][ind]])
                
        elif dis==False:   
            if dictionary['Date'][ind] in dates and dictionary['Drug'][ind] == drug:
                subset.append([onsets, offsets, dictionary['RatID'][ind]])
            
        if verbose: #assumes true
            print('filename, or comment ...') 
    return subset

    
def medfilereader(filename, varsToExtract = 'all',
                  sessionToExtract = 1,
                  verbose = False,
                  remove_var_header = False):
    if varsToExtract == 'all':
        numVarsToExtract = np.arange(0,26)
    else:
        numVarsToExtract = [ord(x)-97 for x in varsToExtract]
    
    f = open(filename, 'r')
    f.seek(0)
    filerows = f.readlines()[8:]
    datarows = [asnumeric(x) for x in filerows]
    matches = [i for i,x in enumerate(datarows) if x == 0.3]
    if sessionToExtract > len(matches):
        print('Session ' + str(sessionToExtract) + ' does not exist.')
    if verbose == True:
        print('There are ' + str(len(matches)) + ' sessions in ' + filename)
        print('Analyzing session ' + str(sessionToExtract))
    
    varstart = matches[sessionToExtract - 1]
    medvars = [[] for n in range(26)]
    
    k = int(varstart + 27)
    for i in range(26):
        medvarsN = int(datarows[varstart + i + 1])
        
        medvars[i] = datarows[k:k + int(medvarsN)]
        k = k + medvarsN
        
    if remove_var_header == True:
        varsToReturn = [medvars[i][1:] for i in numVarsToExtract]
    else:
        varsToReturn = [medvars[i] for i in numVarsToExtract]

    if np.shape(varsToReturn)[0] == 1:
        varsToReturn = varsToReturn[0]
    return varsToReturn

def asnumeric(s):
    try:
        x = float(s)
        return x
    except ValueError:
        return float('nan')



#now get the function to give the indices or the time stamps of the 
#distractors that are white noise and non whitenoise 
def discalc_modalities(dictionary, modalitykey):
    ''' Calculates distractors, distracted and modalities for dictionary of 
    rats by group for just distraction day only 
    
    Calculates a grouped percentage of each modality and how distracted 
    by that modality rats are on average (by group)
    
    '''
    percent_dis_whitenoise_group = []
    percent_dis_tone_group = []
    percent_dis_combined_group = []
    percent_dis_all_non_WN_group = []
    discalcgroup = []
    
    ## SAL MALES - DISTRACTION DAY ONLY - DISTRACTOR TYPE ANALYSIS INCLUDED
# Finds distracted or not (corrects for med slipping issue)
    for rat in dictionary:
        
        discalc = distractionCalc2(rat[0])
        distracted, notdistracted = distractedOrNot(discalc, rat[0])
      # too simplistic, here there are distracted, not dis and then ALL DISTRACTORS TYPES
       
    
        dis_numeric = []
        ndis_numeric = []
    # Modality analysis - calculates which distractors contain different features (whitenoise, tone or combination)
    # Then works out on how many of these trials rats are distracted (individual) before creating a mean 
        for d in distracted:
            dis_numeric.append([rat[2][idx] for idx, val in enumerate(discalc) if val == d][0])
            
        for nd in notdistracted:
            ndis_numeric.append([rat[2][idx] for idx, val in enumerate(discalc) if val == nd][0])   
      
        
        discalcgroup.append([distracted, dis_numeric, notdistracted, ndis_numeric])
        # Makes the distracted trial types into integers 
        
 ### 4 lists, distracted time stamp then identity of cue, not distracted timestamp and identity       

        
        dis_numeric = [int(d) for d in dis_numeric]
        print(dis_numeric)
        # Counts to work out percentages after finding how many are each modality 
        d_whitenoise_count = 0
        d_tone_count = 0
        d_combined_count = 0 
        
        dis_type_text = [] #labels the distypes with text labels and adds to the counts
        for d in dis_numeric:
            if d in modalitykey['whitenoise']:
                dis_type_text.append('whitenoise')
                d_whitenoise_count += 1
            elif d in modalitykey['tone']:
                dis_type_text.append('tone')
                d_tone_count += 1
            elif d in modalitykey['combined3']:
                dis_type_text.append('combined3')
                d_combined_count += 1 
        
        #print(d_whitenoise_count, d_tone_count, d_combined_count)
        

            
        # Non-distracted trials by modality 
        ndis_numeric = [int(d) for d in ndis_numeric]
        nd_whitenoise_count = 0
        nd_tone_count = 0
        nd_combined_count = 0 
        
        ndis_type_text = []
        for d in ndis_numeric:
            if d in modalitykey['whitenoise']:
                ndis_type_text.append('whitenoise')
                nd_whitenoise_count += 1
            elif d in modalitykey['tone']:
                ndis_type_text.append('tone')
                nd_tone_count += 1
            elif d in modalitykey['combined3']:
                ndis_type_text.append('combined3')
                nd_combined_count += 1 



### Added 20/07 to fix division by zero 
#        if d_whitenoise_count != 0:
##            
#            d_percent_white_noise = (d_whitenoise_count / (d_whitenoise_count+nd_whitenoise_count))*100
#            d_percent_tone = (d_tone_count / (d_tone_count+nd_tone_count))*100
#            d_percent_combined = (d_combined_count / (d_combined_count+nd_combined_count))*100
#      
#        else:
#            d_percent_white_noise = 0 
#            d_percent_tone = 0
#            d_percent_combined = 0                
#
#        nd_percent_white_noise = (nd_whitenoise_count /(d_whitenoise_count+nd_whitenoise_count))*100
#        nd_percent_tone = (nd_tone_count /(d_tone_count+nd_tone_count))*100
#        nd_percent_combined =  (nd_combined_count /(d_combined_count+nd_combined_count))*100
#        
       # print(d_percent_white_noise, nd_percent_white_noise)
        
        percent_distracted_whitenoise = d_whitenoise_count / (d_whitenoise_count + nd_whitenoise_count) *100
        percent_distracted_tone = d_tone_count / (d_tone_count + nd_tone_count) *100
        percent_distracted_combined = d_combined_count / (d_combined_count + nd_combined_count) *100  
        percent_distracted_all_non_WN = (d_tone_count + d_combined_count) / ((d_tone_count + d_combined_count )+ (nd_tone_count + nd_combined_count)) *100  
        
        percent_dis_whitenoise_group.append(percent_distracted_whitenoise)
        percent_dis_tone_group.append(percent_distracted_tone)
        percent_dis_combined_group.append(percent_distracted_combined)
        percent_dis_all_non_WN_group.append(percent_distracted_all_non_WN)
    
     
    mean_percent_WHITENOISE = np.mean(percent_dis_whitenoise_group) # the average percentage of JUST whitenoise trials that rats are distracted on 
    mean_percent_TONE = np.mean(percent_dis_tone_group)
    mean_percent_COMBINED = np.mean(percent_dis_combined_group)
    mean_percent_ALL_NON_WN = np.mean(percent_dis_all_non_WN_group)
    
    return discalcgroup, percent_dis_whitenoise_group, percent_dis_tone_group, \
            percent_dis_combined_group, percent_dis_all_non_WN_group, mean_percent_WHITENOISE, mean_percent_TONE, \
            mean_percent_COMBINED, mean_percent_ALL_NON_WN


## Add the outputs to discalc group (maybe the list of dis and not dis by WN and NWN)





## Isolate the code that  makes the TrialsMultFig for distractors,
## distracted and not distracted (chapter 4 - analysis_licking)


## Isolate the code that calculates TIMES of the whitenoise containing 
## and non-whitenoise containing distractors (chapter 3_ modalities)


## May need to add a time indexing seciton (as it currently just uses counts
## to calculate percent distracted rather than the stamps of where 
"""
(1) Reads in metafiles from filepath, contains information on filenames and basic descriptives
(2) Extracts lick data from MED files. Uses the names files in metafile and stores licking 
    data into variables by chosen day (ie. last lick day / distraction day). Stores lick data
    by group as large list of lists. Eg. last_lick_sal_M contains 16 lists of lick onsets for 
    just the last day of licking, each list is a rat, with onsets and offsets stored.
(3) Lick analysis - 
(4) Distraction analysis -  
(5) Post distractor pauses and predistractor pauses - calculates the time periods 
    for both distracted and non distracted trials separately for all groups

"""
################################################################################
# LICK ANALYSIS THPH (for med not tdt)
################################################################################
# (1) Find and extract info from the metafile(s) 
metafileTHPH1_2 = '/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Ch4_whitenoise/metafile.csv'
extractTHPH1_2 = MetaExtractor(metafileTHPH1_2)
# Folder with all medfiles (THPH1, THPH2)
medfolder = '/Volumes/KP_HARD_DRI/kp259/THPH1AND2/med/' # sometimes need space1 after DRI

# (2) - Get all lick onset/offset data for all rats/sessions (make this a function later)
''' 
  Info THPH1 and THPH2:
      
  THPH1 |THPH2 |CONDITION
  ------|------|-----------
  170418|171007|dis

'''
# Subsetting all data by day / drug (drug = 1 is like an include col here)

distraction = subsetter(extractTHPH1_2, ['170621','170810', '170815'], '1', dis=True)

modalitykey = {'whitenoise':[1,4], 'tone':[2,5], 'combined3':[3,6]}

discalc, percent_dis_whitenoise, percent_dis_tone,\
percent_dis_combined, percent_dis_all_non_WN, mean_percent_WHITENOISE, mean_percent_TONE,\
mean_percent_COMBINED, mean_percent_ALL_NON_WN = discalc_modalities(distraction, modalitykey)


#1,3,4,6 = white noise 
#2,5 = non white noise 



## Simply exclude rat (11 / whichever rat it is )
## HERE WRITE THE PEAKS FOR DISTRACTED FROM WN AND NWN 


## Code here using (1) List of white noise (2) Non white noise

## need rat data (use distracted or distractors as the blue and uv)
## if run the analysis programme first should work 

# DISTRACTION FILES (minus the first 2) - this was run with all included 
### Distractors, distracted and not distracted, licks and blue / uv signals 
ratdis2 = []
ratdis3 = []

allRatBlue = []
allRatUV = []
allRatFS = []
allRatLicks = []
allRatDistractors = []
allRatDistracted = []
allRatNotDistracted = []

blueMeans_whitenoise = []
uvMeans_whitenoise = [] 

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



### index out of range issue, need to get rid of distractors where it is the last
## distractor / last lick     

ratdis2 = []    
for index, rat in enumerate(allRatDistracted):
    ratdislistWN = []
    for index2, item in enumerate(rat): ## NO - if index2 in discalc index[0]
        if discalc[index][1][index2] == 1 or discalc[index][1][index2]== 4 or discalc[index][1][index2] == 3 or discalc[index][1][index2] == 4 or discalc[index][1][index2] == 6:
            a = allRatDistracted[index][index2]
            ratdislistWN.append(a)
            
    ratdis2.append(ratdislistWN)


ratdis3 = []
for index, rat in enumerate(allRatDistracted):
    ratdislistNWN = []
    for index2, item in enumerate(rat): ## NO - if index2 in discalc index[0]
#        print(index2)
        if discalc[index][1][index2] == 2 or discalc[index][1][index2] == 5:
            a = allRatDistracted[index][index2]
            ratdislistNWN.append(a)
            
    ratdis3.append(ratdislistNWN)    

ratdis4 = []
for index, element in enumerate(ratdis2):
    ratdis4.append(ratdis2[index] + ratdis3[index])    

for i, val in enumerate(ratdis4):

    try:
        # make a blue and uv snip for all 14, and noise remover / index
        blueSnips, ppsBlue = snipper(allRatBlue[i], ratdis4[i], fs=allRatFS[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUV[i], ratdis4[i], fs=allRatFS[i], bins=300)
    
        randevents = makerandomevents(allRatBlue[i][300], allRatBlue[i][-300])
        bgMad, bgMean = findnoise(allRatBlue[i], randevents, fs=allRatFS[i], method='sum', bins=300)
        threshold = 1
        sigSum = [np.sum(abs(i)) for i in blueSnips]
        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
        # Might not need the noise index, this is just for trials fig 
    except: 
        pass

    blueMeanDISTRACTED = np.mean(blueSnips, axis=0)
    blueMeans_distracted.append(blueMeanDISTRACTED)
    uvMeanDISTRACTED = np.mean(uvSnips, axis=0)
    uvMeans_distracted.append(uvMeanDISTRACTED)
    allbluesnips.append(blueSnips)
    alluvsnips.append(uvSnips)
# Means for distracted trials here MULT SHADED FIG 
fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.06, 0.06])
trialsMultShadedFig(ax, [np.asarray(uvMeans_distracted),np.asarray(blueMeans_distracted)], ppsBlue, eventText='Distracted trial', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'], scale=0)
# EDIT THIS TEXT TO SHOW NUMBER OF TOTAL DISTRACTORS OR TRIALS ON THE AVERAGED PLOT 
#plt.text(250,0.03, '{}'.format(len(MergedRunList_Long)) + ' Long Runs' ) ## Edit this to be all
#fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Distracted_All_Rats.pdf', bbox_inches="tight")



################################################################
################################################################
################################################################
################################################################
################################################################ 
ratdis2 = []
ratdis3 = []

allRatBlue = []
allRatUV = []
allRatFS = []
allRatLicks = []
allRatDistractors = []
allRatDistracted = []
allRatNotDistracted = []

blueMeans_whitenoise = []
uvMeans_whitenoise = [] 

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

    

ratdis2 = []    
for index, rat in enumerate(allRatDistracted):
    ratdislistWN = []
    for index2, item in enumerate(rat): ## NO - if index2 in discalc index[0]
        if discalc[index][1][index2] == 1 or discalc[index][1][index2]== 4 or discalc[index][1][index2] == 3 or discalc[index][1][index2] == 4 or discalc[index][1][index2] == 6:
            a = allRatDistracted[index][index2]
            ratdislistWN.append(a)
            
    ratdis2.append(ratdislistWN)


ratdis3 = []
for index, rat in enumerate(allRatDistracted):
    ratdislistNWN = []
    for index2, item in enumerate(rat): ## NO - if index2 in discalc index[0]
#        print(index2)
        if discalc[index][1][index2] == 2 or discalc[index][1][index2] == 5:
            a = allRatDistracted[index][index2]
            ratdislistNWN.append(a)
            
    ratdis3.append(ratdislistNWN)    

ratdis4 = []
for index, element in enumerate(ratdis2):
    ratdis4.append(ratdis2[index] + ratdis3[index])    

for i, val in enumerate(ratdis2):

    try:
        # make a blue and uv snip for all 14, and noise remover / index
        blueSnips, ppsBlue = snipper(allRatBlue[i], ratdis2[i], fs=allRatFS[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUV[i], ratdis2[i], fs=allRatFS[i], bins=300)
    
        randevents = makerandomevents(allRatBlue[i][300], allRatBlue[i][-300])
        bgMad, bgMean = findnoise(allRatBlue[i], randevents, fs=allRatFS[i], method='sum', bins=300)
        threshold = 1
        sigSum = [np.sum(abs(i)) for i in blueSnips]
        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
        # Might not need the noise index, this is just for trials fig 
    except: 
        pass


    blueMeanDISTRACTED = np.mean(blueSnips, axis=0)
    blueMeans_distracted.append(blueMeanDISTRACTED)
    uvMeanDISTRACTED = np.mean(uvSnips, axis=0)
    uvMeans_distracted.append(uvMeanDISTRACTED)
    allbluesnips.append(blueSnips)
    alluvsnips.append(uvSnips)
# Means for distracted trials here MULT SHADED FIG 
fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.06, 0.06])
trialsMultShadedFig(ax, [np.asarray(uvMeans_distracted),np.asarray(blueMeans_distracted)], ppsBlue, eventText='White noise distractor', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'], scale=0)
fig.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/WhiteNoiseDistractedPhoto.pdf', bbox_inches="tight")


bkgnd_sub_DistractorWN = uvSubtractor(blueMeans_distracted, uvMeans_distracted)




ratdis2 = []
ratdis3 = []

allRatBlue = []
allRatUV = []
allRatFS = []
allRatLicks = []
allRatDistractors = []
allRatDistracted = []
allRatNotDistracted = []

blueMeans_whitenoise = []
uvMeans_whitenoise = [] 

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

  

ratdis2 = []    
for index, rat in enumerate(allRatDistracted):
    ratdislistWN = []
    for index2, item in enumerate(rat): ## NO - if index2 in discalc index[0]
        if discalc[index][1][index2] == 1 or discalc[index][1][index2]== 4 or discalc[index][1][index2] == 3 or discalc[index][1][index2] == 4 or discalc[index][1][index2] == 6:
            a = allRatDistracted[index][index2]
            ratdislistWN.append(a)
            
    ratdis2.append(ratdislistWN)


ratdis3 = []
for index, rat in enumerate(allRatDistracted):
    ratdislistNWN = []
    for index2, item in enumerate(rat): ## NO - if index2 in discalc index[0]
#        print(index2)
        if discalc[index][1][index2] == 2 or discalc[index][1][index2] == 5:
            a = allRatDistracted[index][index2]
            ratdislistNWN.append(a)
            
    ratdis3.append(ratdislistNWN)    

ratdis4 = []
for index, element in enumerate(ratdis2):
    ratdis4.append(ratdis2[index] + ratdis3[index])    

for i, val in enumerate(ratdis3):

    try:
        # make a blue and uv snip for all 14, and noise remover / index
        blueSnips, ppsBlue = snipper(allRatBlue[i], ratdis3[i], fs=allRatFS[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUV[i], ratdis3[i], fs=allRatFS[i], bins=300)
    
        randevents = makerandomevents(allRatBlue[i][300], allRatBlue[i][-300])
        bgMad, bgMean = findnoise(allRatBlue[i], randevents, fs=allRatFS[i], method='sum', bins=300)
        threshold = 1
        sigSum = [np.sum(abs(i)) for i in blueSnips]
        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
        # Might not need the noise index, this is just for trials fig 
    except: 
        pass

    blueMeanDISTRACTED = np.mean(blueSnips, axis=0)
    blueMeans_distracted.append(blueMeanDISTRACTED)
    uvMeanDISTRACTED = np.mean(uvSnips, axis=0)
    uvMeans_distracted.append(uvMeanDISTRACTED)
    allbluesnips.append(blueSnips)
    alluvsnips.append(uvSnips)
# Means for distracted trials here MULT SHADED FIG 
fig = plt.figure(figsize=(6,3))
ax = plt.subplot(1,1,1)
ax.set_ylim([-0.06, 0.06])
trialsMultShadedFig(ax, [np.asarray(uvMeans_distracted),np.asarray(blueMeans_distracted)], ppsBlue, eventText='Non white noise distractor', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'], scale=0)
fig.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/NonWhiteNoiseDistractedPhoto.pdf', bbox_inches="tight")

bkgnd_sub_DistractorNWN = uvSubtractor(blueMeans_distracted, uvMeans_distracted)


peak_distractorWN, t_distractorWN, pre_distractorWN, post_distractorWN, baseline_distractorWN = PhotoPeaksCalc(bkgnd_sub_DistractorWN)
peak_distractorNWN, t_distractorNWN, pre_distractorNWN, post_distractorNWN, baseline_distractorNWN = PhotoPeaksCalc(bkgnd_sub_DistractorNWN)

np.mean(peak_distractorWN)
np.mean(peak_distractorNWN)

## Add in barscatter here too for all 4 measures 

mpl.rcParams['figure.subplot.wspace'] = 0.6
mpl.rcParams['figure.subplot.right'] = 1
mpl.rcParams['font.size'] = 14
def MultBy100(list):
    output = [x*100 for x in list]
    
    return output
## Modelled versus distracion day presented distractors (GREY and GREEN)
wnVnwnPeak = [MultBy100(peak_distractorWN), MultBy100(peak_distractorNWN)]
wnVnwnt = [t_distractorWN, t_distractorNWN]
wnVnwnPre = [MultBy100(pre_distractorWN), MultBy100(pre_distractorNWN)]
wnVnwnPost = [MultBy100(post_distractorWN), MultBy100(post_distractorNWN)]

figureA, ax = plt.subplots(nrows=1, ncols=4, figsize=(10,4)) ### x,y 
figureA.tight_layout(pad=3, w_pad=3, h_pad=1.0)

labels = []
ax[0], barx, barlist, sclist = barscatter(wnVnwnPeak, ax=ax[0],transpose=False, paired=True, barfacecolor=['powderblue','gold'], barfacecoloroption='individual',  ylabel='Peak (%)', barlabels=labels, baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[1], barx, barlist, sclist = barscatter(wnVnwnt, ax=ax[1], transpose=False, paired=True, barfacecolor=['powderblue','gold'], barfacecoloroption='individual',  ylabel='t (s)', barlabels=labels, baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[2], barx, barlist, sclist = barscatter(wnVnwnPre, ax=ax[2],transpose=False, paired=True, barfacecolor=['powderblue','gold'], barfacecoloroption='individual',  ylabel='Pre-event period (mean %)', barlabels=labels,  baredgecolor=[''] )#,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[3], barx, barlist, sclist = barscatter(wnVnwnPost, ax=ax[3],transpose=False, paired=True, barfacecolor=['powderblue','gold'], barfacecoloroption='individual',  ylabel='Post-event period (mean %)', barlabels=labels, baredgecolor=[''] )#,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])

ax[0].set_ylabel('Peak (% ΔF)')
ax[1].set_ylabel('t (s)')
ax[2].set_ylabel('Pre-event period (mean % ΔF)')
ax[3].set_ylabel('Post-event period (mean % ΔF)')

ax[0].set_xticks([])
#ax[0].set_ylim([0,4000])
ax[1].set_xticks([])
#ax[1].set_ylim([0,25])
ax[2].set_xticks([])
#ax[2].set_ylim(0,1200)
ax[3].set_xticks([])
#ax[3].set_ylim(0,1200)

ax[0].spines['bottom'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[2].spines['bottom'].set_visible(False)
ax[3].spines['bottom'].set_visible(False)


figureA.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/WhitenoiseVNonWhitenoiseBarScatter.pdf', bbox_inches="tight")







# OLD SCRIPT BELOW:
    
    
    
    
    
##!/usr/bin/env python3
## -*- coding: utf-8 -*-
#"""
#Created on Mon Dec 10 09:57:41 2018
#
#@author: u1490431
#"""
## NECESSARY IMPORTS AND FUNCTIONS 
#
## Imports
#import numpy as np
#from scipy import stats
#from scipy.stats.stats import pearsonr
#import matplotlib.pyplot as plt
#import matplotlib as mpl 
#
##import seaborn as sn
#
#def uvSubtractor(rat_snip_means_list_blue, uv):
#    subtractedSignal = []
#    for ind, rat in enumerate(rat_snip_means_list_blue):
#        subtractedSignal.append(rat_snip_means_list_blue[ind] - uv[ind])
#        
#    return subtractedSignal
#def PhotoPeaksCalc(snips_all_rats):
#    
#    allRat_peak = []
#    allRat_t = []
#    allRat_pre = []
#    allRat_post = []
#    allRat_base = []
#    
#    for rat in snips_all_rats:
#        pre_event = np.mean(rat[0:50]) # Average for 5 seconds, 10 seconds before event 
#        peak = np.max(rat[100:300]) ## Minus the average of the first 5 seconds and after 100 points (slice)
#        peak_range = rat[100:300]
#        a = peak_range.tolist()
#        peak_index = a.index(peak)
#        t = peak_index / 10
#        pre_event = np.mean(rat[50:100])
#        post_event = np.mean(rat[100:300])
#        baseline = np.mean(rat[0:50])
#        
#        allRat_peak.append(peak)
#        allRat_t.append(t)
#        allRat_pre.append(pre_event)
#        allRat_post.append(post_event)
#        allRat_base.append(baseline) 
#        
#    
#    return allRat_peak, allRat_t, allRat_pre, allRat_post, allRat_base
#
#def MetaExtractor (metafile):
#    f = open(metafile, 'r')
#    f.seek(0)
#    Metafilerows = f.readlines()[1:]
#    tablerows = []
#
#    for row in Metafilerows: 
#        items = row.split(',')
#        tablerows.append(items)
#
#    MedFilenames, RatID, Date, Day, Session, Drug, TotLicks, Distractions, \
#    NonDistractions, PercentDistracted = [], [], [], [], [], [], [], [], [], []
#
#    for i, lst in enumerate(tablerows):
#       MedFilenames = MedFilenames + [lst[0]]
#       RatID = RatID + [lst[1]]
#       Date = Date + [lst[2]]
#       Session = Session + [lst[3]]
#       Drug = Drug + [lst[4]]
#       TotLicks = TotLicks + [lst[5]]
##       Distractions = Distractions + [lst[7]] 
##       NonDistractions = NonDistractions + [lst[8]]
#       PercentDistracted = PercentDistracted + [lst[8]]
#    return ({'MedFilenames':MedFilenames, 'RatID':RatID, 'Date':Date, 'Session':Session, \
#             'Drug':Drug, 'TotLicks':TotLicks})
#
#
#
#def subsetter(dictionary, dates, drug, dis=False, verbose=False):
#    '''
#    SUBSETTER KP
#    # Subsets data according to date, reads in dictionnary produced from metafile
#    # and subsets into variable based on date(s) and drug condition 
#    # if distraction day argument is given as True adds the distractor type 
#    # to the output lists for later processing 
#    
#    '''
#    subset = []
#    for ind, filename in enumerate(dictionary['MedFilenames']):
#        path = medfolder + filename
#        onsets, offsets, med_dis_times, dis_type = medfilereader(path, ['b', 'c', 'i', 'j'], remove_var_header = True)  # e onset, f offset
#
#        if dis == True:
#            if dictionary['Date'][ind] in dates and dictionary['Drug'][ind] == drug:
#                subset.append([onsets, offsets, dis_type, dictionary['RatID'][ind]])
#                
#        elif dis==False:   
#            if dictionary['Date'][ind] in dates and dictionary['Drug'][ind] == drug:
#                subset.append([onsets, offsets, dictionary['RatID'][ind]])
#            
#        if verbose: #assumes true
#            print('filename, or comment ...') 
#    return subset
#
#    
#def medfilereader(filename, varsToExtract = 'all',
#                  sessionToExtract = 1,
#                  verbose = False,
#                  remove_var_header = False):
#    if varsToExtract == 'all':
#        numVarsToExtract = np.arange(0,26)
#    else:
#        numVarsToExtract = [ord(x)-97 for x in varsToExtract]
#    
#    f = open(filename, 'r')
#    f.seek(0)
#    filerows = f.readlines()[8:]
#    datarows = [asnumeric(x) for x in filerows]
#    matches = [i for i,x in enumerate(datarows) if x == 0.3]
#    if sessionToExtract > len(matches):
#        print('Session ' + str(sessionToExtract) + ' does not exist.')
#    if verbose == True:
#        print('There are ' + str(len(matches)) + ' sessions in ' + filename)
#        print('Analyzing session ' + str(sessionToExtract))
#    
#    varstart = matches[sessionToExtract - 1]
#    medvars = [[] for n in range(26)]
#    
#    k = int(varstart + 27)
#    for i in range(26):
#        medvarsN = int(datarows[varstart + i + 1])
#        
#        medvars[i] = datarows[k:k + int(medvarsN)]
#        k = k + medvarsN
#        
#    if remove_var_header == True:
#        varsToReturn = [medvars[i][1:] for i in numVarsToExtract]
#    else:
#        varsToReturn = [medvars[i] for i in numVarsToExtract]
#
#    if np.shape(varsToReturn)[0] == 1:
#        varsToReturn = varsToReturn[0]
#    return varsToReturn
#
#def asnumeric(s):
#    try:
#        x = float(s)
#        return x
#    except ValueError:
#        return float('nan')
#
#
#
##now get the function to give the indices or the time stamps of the 
##distractors that are white noise and non whitenoise 
#def discalc_modalities(dictionary, modalitykey):
#    ''' Calculates distractors, distracted and modalities for dictionary of 
#    rats by group for just distraction day only 
#    
#    Calculates a grouped percentage of each modality and how distracted 
#    by that modality rats are on average (by group)
#    
#    '''
#    percent_dis_whitenoise_group = []
#    percent_dis_tone_group = []
#    percent_dis_combined_group = []
#    percent_dis_all_non_WN_group = []
#    discalcgroup = []
#    
#    ## SAL MALES - DISTRACTION DAY ONLY - DISTRACTOR TYPE ANALYSIS INCLUDED
## Finds distracted or not (corrects for med slipping issue)
#    for rat in dictionary:
#        
#        discalc = distractionCalc2(rat[0])
#        distracted, notdistracted = distractedOrNot(discalc, rat[0])
#      # too simplistic, here there are distracted, not dis and then ALL DISTRACTORS TYPES
#       
#    
#        dis_numeric = []
#        ndis_numeric = []
#    # Modality analysis - calculates which distractors contain different features (whitenoise, tone or combination)
#    # Then works out on how many of these trials rats are distracted (individual) before creating a mean 
#        for d in distracted:
#            dis_numeric.append([rat[2][idx] for idx, val in enumerate(discalc) if val == d][0])
#            
#        for nd in notdistracted:
#            ndis_numeric.append([rat[2][idx] for idx, val in enumerate(discalc) if val == nd][0])   
#      
#        
#        discalcgroup.append([distracted, dis_numeric, notdistracted, ndis_numeric])
#        # Makes the distracted trial types into integers 
#        
# ### 4 lists, distracted time stamp then identity of cue, not distracted timestamp and identity       
#
#        
#        dis_numeric = [int(d) for d in dis_numeric]
#        print(dis_numeric)
#        # Counts to work out percentages after finding how many are each modality 
#        d_whitenoise_count = 0
#        d_tone_count = 0
#        d_combined_count = 0 
#        
#        dis_type_text = [] #labels the distypes with text labels and adds to the counts
#        for d in dis_numeric:
#            if d in modalitykey['whitenoise']:
#                dis_type_text.append('whitenoise')
#                d_whitenoise_count += 1
#            elif d in modalitykey['tone']:
#                dis_type_text.append('tone')
#                d_tone_count += 1
#            elif d in modalitykey['combined3']:
#                dis_type_text.append('combined3')
#                d_combined_count += 1 
#        
#        #print(d_whitenoise_count, d_tone_count, d_combined_count)
#        
#
#            
#        # Non-distracted trials by modality 
#        ndis_numeric = [int(d) for d in ndis_numeric]
#        nd_whitenoise_count = 0
#        nd_tone_count = 0
#        nd_combined_count = 0 
#        
#        ndis_type_text = []
#        for d in ndis_numeric:
#            if d in modalitykey['whitenoise']:
#                ndis_type_text.append('whitenoise')
#                nd_whitenoise_count += 1
#            elif d in modalitykey['tone']:
#                ndis_type_text.append('tone')
#                nd_tone_count += 1
#            elif d in modalitykey['combined3']:
#                ndis_type_text.append('combined3')
#                nd_combined_count += 1 
#
#
#
#### Added 20/07 to fix division by zero 
##        if d_whitenoise_count != 0:
###            
##            d_percent_white_noise = (d_whitenoise_count / (d_whitenoise_count+nd_whitenoise_count))*100
##            d_percent_tone = (d_tone_count / (d_tone_count+nd_tone_count))*100
##            d_percent_combined = (d_combined_count / (d_combined_count+nd_combined_count))*100
##      
##        else:
##            d_percent_white_noise = 0 
##            d_percent_tone = 0
##            d_percent_combined = 0                
##
##        nd_percent_white_noise = (nd_whitenoise_count /(d_whitenoise_count+nd_whitenoise_count))*100
##        nd_percent_tone = (nd_tone_count /(d_tone_count+nd_tone_count))*100
##        nd_percent_combined =  (nd_combined_count /(d_combined_count+nd_combined_count))*100
##        
#       # print(d_percent_white_noise, nd_percent_white_noise)
#        
#        percent_distracted_whitenoise = d_whitenoise_count / (d_whitenoise_count + nd_whitenoise_count) *100
#        percent_distracted_tone = d_tone_count / (d_tone_count + nd_tone_count) *100
#        percent_distracted_combined = d_combined_count / (d_combined_count + nd_combined_count) *100  
#        percent_distracted_all_non_WN = (d_tone_count + d_combined_count) / ((d_tone_count + d_combined_count )+ (nd_tone_count + nd_combined_count)) *100  
#        
#        percent_dis_whitenoise_group.append(percent_distracted_whitenoise)
#        percent_dis_tone_group.append(percent_distracted_tone)
#        percent_dis_combined_group.append(percent_distracted_combined)
#        percent_dis_all_non_WN_group.append(percent_distracted_all_non_WN)
#    
#     
#    mean_percent_WHITENOISE = np.mean(percent_dis_whitenoise_group) # the average percentage of JUST whitenoise trials that rats are distracted on 
#    mean_percent_TONE = np.mean(percent_dis_tone_group)
#    mean_percent_COMBINED = np.mean(percent_dis_combined_group)
#    mean_percent_ALL_NON_WN = np.mean(percent_dis_all_non_WN_group)
#    
#    return discalcgroup, percent_dis_whitenoise_group, percent_dis_tone_group, \
#            percent_dis_combined_group, percent_dis_all_non_WN_group, mean_percent_WHITENOISE, mean_percent_TONE, \
#            mean_percent_COMBINED, mean_percent_ALL_NON_WN
#
#
### Add the outputs to discalc group (maybe the list of dis and not dis by WN and NWN)
#
#
#
#
#
### Isolate the code that  makes the TrialsMultFig for distractors,
### distracted and not distracted (chapter 4 - analysis_licking)
#
#
### Isolate the code that calculates TIMES of the whitenoise containing 
### and non-whitenoise containing distractors (chapter 3_ modalities)
#
#
### May need to add a time indexing seciton (as it currently just uses counts
### to calculate percent distracted rather than the stamps of where 
#"""
#(1) Reads in metafiles from filepath, contains information on filenames and basic descriptives
#(2) Extracts lick data from MED files. Uses the names files in metafile and stores licking 
#    data into variables by chosen day (ie. last lick day / distraction day). Stores lick data
#    by group as large list of lists. Eg. last_lick_sal_M contains 16 lists of lick onsets for 
#    just the last day of licking, each list is a rat, with onsets and offsets stored.
#(3) Lick analysis - 
#(4) Distraction analysis -  
#(5) Post distractor pauses and predistractor pauses - calculates the time periods 
#    for both distracted and non distracted trials separately for all groups
#
#"""
#################################################################################
## LICK ANALYSIS THPH (for med not tdt)
#################################################################################
## (1) Find and extract info from the metafile(s) 
#metafileTHPH1_2 = '/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Ch4_whitenoise/metafile.csv'
#extractTHPH1_2 = MetaExtractor(metafileTHPH1_2)
## Folder with all medfiles (THPH1, THPH2)
#medfolder = '/Volumes/KP_HARD_DRI/kp259/THPH1AND2/med/' # sometimes need space1 after DRI
#
## (2) - Get all lick onset/offset data for all rats/sessions (make this a function later)
#''' 
#  Info THPH1 and THPH2:
#      
#  THPH1 |THPH2 |CONDITION
#  ------|------|-----------
#  170418|171007|dis
#
#'''
## Subsetting all data by day / drug (drug = 1 is like an include col here)
#
#distraction = subsetter(extractTHPH1_2, ['170621','170810', '170815'], '1', dis=True)
#
#modalitykey = {'whitenoise':[1,4], 'tone':[2,5], 'combined3':[3,6]}
#
#discalc, percent_dis_whitenoise, percent_dis_tone,\
#percent_dis_combined, percent_dis_all_non_WN, mean_percent_WHITENOISE, mean_percent_TONE,\
#mean_percent_COMBINED, mean_percent_ALL_NON_WN = discalc_modalities(distraction, modalitykey)
#
#
##1,3,4,6 = white noise 
##2,5 = non white noise 
#
#'''
#ratdis2 = []  ### This is the list of timestamps for distracted trials 
#                ### that are WHITE NOISE CONTAINING 
#
#for index, rat in enumerate(discalc):
#    ratdislistWN = []
#    for index2, item in enumerate(rat[1]):
#        if item == 1 or item == 3 or item == 4 or item == 6:
#            a = discalc[index][0][index2]
#            ratdislistWN.append(a)
#        
#    ratdis2.append(ratdislistWN)
#    
#
#
#ratdis3 = []  ### This is the list of timestamps for distracted trials 
#                ### that are NON - WHITE NOISE CONTAINING 
#
#for index, rat in enumerate(discalc):
#    ratdislistNWN = []
#    for index2, item in enumerate(rat[1]):
#        if item == 2 or item == 5:
#            b = discalc[index][0][index2]
#            ratdislistNWN.append(b)
#        
#    ratdis3.append(ratdislistNWN)
#'''
#
### Simply exclude rat (11 / whichever rat it is )
### HERE WRITE THE PEAKS FOR DISTRACTED FROM WN AND NWN 
##peak_distractor, t_distractor, pre_distractor, post_distractor, baseline_distractor = PhotoPeaksCalc(bkgnd_sub_Distractor)
##peak_distractor, t_distractor, pre_distractor, post_distractor, baseline_distractor = PhotoPeaksCalc(bkgnd_sub_Distractor)
#
#
### Code here using (1) List of white noise (2) Non white noise
#
### need rat data (use distracted or distractors as the blue and uv)
### if run the analysis programme first should work 
#
## DISTRACTION FILES (minus the first 2) - this was run with all included 
#### Distractors, distracted and not distracted, licks and blue / uv signals 
#ratdis2 = []
#ratdis3 = []
#
#allRatBlue = []
#allRatUV = []
#allRatFS = []
#allRatLicks = []
#allRatDistractors = []
#allRatDistracted = []
#allRatNotDistracted = []
#
#blueMeans_whitenoise = []
#uvMeans_whitenoise = [] 
#
#blueMeans_distracted = []
#uvMeans_distracted = []
#blueMeans_notdistracted = []
#uvMeans_notdistracted = [] 
#allbluesnips = []
#alluvsnips = []
#
#for filename in TDTfiles_thph_dis:
#    
#    file = TDTfilepath + filename
#    ratdata = loadmatfile(file)
#    allRatBlue.append(ratdata['blue'])
#    allRatUV.append(ratdata['uv'])
#    allRatFS.append(ratdata['fs'])
#    allRatLicks.append(ratdata['licks'])
#    allRatDistractors.append(ratdata['distractors'])
#    allRatDistracted.append(ratdata['distracted'])
#
#
#
#### index out of range issue, need to get rid of distractors where it is the last
### distractor / last lick     
#
#ratdis2 = []    
#for index, rat in enumerate(allRatDistracted):
#    ratdislistWN = []
#    for index2, item in enumerate(rat): ## NO - if index2 in discalc index[0]
#        print(index2)
#        if discalc[index][1][index2] == 1 or discalc[index][1][index2] == 3 or discalc[index][1][index2] == 4 or discalc[index][1][index2] == 6:
#            a = allRatDistracted[index][index2]
#            ratdislistWN.append(a)
#            
#    ratdis2.append(ratdislistWN)
#
#
#ratdis3 = []
#for index, rat in enumerate(allRatDistracted):
#    ratdislistNWN = []
#    for index2, item in enumerate(rat): ## NO - if index2 in discalc index[0]
#        print(index2)
#        if discalc[index][1][index2] == 2 or discalc[index][1][index2] == 5:
#            a = allRatDistracted[index][index2]
#            ratdislistNWN.append(a)
#            
#    ratdis3.append(ratdislistNWN)    
#
#ratdis4 = []
#for index, element in enumerate(ratdis2):
#    ratdis4.append(ratdis2[index] + ratdis3[index])    
#
#for i, val in enumerate(ratdis4):
#
#    try:
#        # make a blue and uv snip for all 14, and noise remover / index
#        blueSnips, ppsBlue = snipper(allRatBlue[i], ratdis4[i], fs=allRatFS[i], bins=300)
#        uvSnips, ppsUV = snipper(allRatUV[i], ratdis4[i], fs=allRatFS[i], bins=300)
#    
#        randevents = makerandomevents(allRatBlue[i][300], allRatBlue[i][-300])
#        bgMad, bgMean = findnoise(allRatBlue[i], randevents, fs=allRatFS[i], method='sum', bins=300)
#        threshold = 1
#        sigSum = [np.sum(abs(i)) for i in blueSnips]
#        noiseindex = [i > bgMean + bgMad*threshold for i in sigSum]
#        # Might not need the noise index, this is just for trials fig 
#    except: 
#        pass
## Individual plots to choose a representative rat 
##    fig14 = plt.figure()
##    ax13 = plt.subplot(1,1,1)
##    ax13.set_ylim([-0.15, 0.15])
##    trialsFig(ax13, blueSnips, uvSnips, ppsBlue, eventText='Distracted') #, noiseindex=noiseindex) #, )
##    plt.text(250,0.2, '{}'.format(len(allRatDistracted[i])) + ' distracted' )
##    fig14.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Distracted_' + str(i) + '.pdf', bbox_inches="tight")
##
#
#    blueMeanDISTRACTED = np.mean(blueSnips, axis=0)
#    blueMeans_distracted.append(blueMeanDISTRACTED)
#    uvMeanDISTRACTED = np.mean(uvSnips, axis=0)
#    uvMeans_distracted.append(uvMeanDISTRACTED)
#    allbluesnips.append(blueSnips)
#    alluvsnips.append(uvSnips)
## Means for distracted trials here MULT SHADED FIG 
#fig = plt.figure(figsize=(6,3))
#ax = plt.subplot(1,1,1)
#ax.set_ylim([-0.06, 0.06])
#trialsMultShadedFig(ax, [np.asarray(uvMeans_distracted),np.asarray(blueMeans_distracted)], ppsBlue, eventText='Distracted trial', linecolor = ['purple','blue'], errorcolor = ['thistle','lightblue'], scale=0)
## EDIT THIS TEXT TO SHOW NUMBER OF TOTAL DISTRACTORS OR TRIALS ON THE AVERAGED PLOT 
##plt.text(250,0.03, '{}'.format(len(MergedRunList_Long)) + ' Long Runs' ) ## Edit this to be all
##fig.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Distracted_All_Rats.pdf', bbox_inches="tight")
#
## delete rat 11 (index -1) before making the photometry plots 
#
#
##deleted rats 1.1 and 1.2 from metafile for this analysis 
## weird -5 thing (possible with backup save file? header is different place?)
## Doesn't work, issue with rat 7, doesnt remove the header 
##rat7 = 1,4,5,6,2,3,3,1,4,5,6,2,3,6,5,1,4,2,4,5,6,1,2,3,1,3,4,5
##discalc[10][1] = rat7
## for now exclude, could add logical index of 1,0 for dis or not and 
## later categorise into list A and list B then use these.
#
#
