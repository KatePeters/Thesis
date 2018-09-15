#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 08:59:29 2018
@author: u1490431

#############################################################################

CHAPTER 2 - DISTRACTION PILOT

#############################################################################

This file contains the DATA PROCESSING and ANALYSIS script 
All functions and imports are at the beginning and no others are required
A separate PLOTTING script should be run second to produce chapter figures            

Imports:
    numpy=np, stats, pearsonr, plt, pml, sn  
Functions: 
    MetaExtractor(),subsetter(),medfilereader(),asnumeric(),lickanalysis()
    lickCalc(),grouped_lickanalysis(),discalc_modalities(),distractionCalc2()
    remcheck(),distractedOrNot(),disbygroup(),percentdisgroup()
"""

# Imports
import numpy as np
from scipy import stats
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib as mpl 
#import seaborn as sn
# Required functions in order of occurrence 
def MetaExtractor (metafile):
    f = open(metafile, 'r')
    f.seek(0)
    Metafilerows = f.readlines()[1:]
    tablerows = []

    for row in Metafilerows: 
        items = row.split(',')
        tablerows.append(items)

    MedFilenames, RatID, Date, Day, Session, Drug, TotLicks, Distractions, \
    NonDistractions, PercentDistracted, End = [], [], [], [], [], [], [], [], [], [], []

    for i, lst in enumerate(tablerows):
       MedFilenames = MedFilenames + [lst[0]]
       RatID = RatID + [lst[1]]
       Date = Date + [lst[2]]
       Day = Day + [lst[3]]
       Session = Session + [lst[4]]
       Drug = Drug + [lst[5]]
       TotLicks = TotLicks + [lst[6]]
       Distractions = Distractions + [lst[7]] 
       NonDistractions = NonDistractions + [lst[8]]
       PercentDistracted = PercentDistracted + [lst[9]]
       End = End + [lst[10]]
 
    return ({'MedFilenames':MedFilenames, 'RatID':RatID, 'Date':Date, 'Day':Day, 'Session':Session, \
             'Drug':Drug, 'TotLicks':TotLicks, 'Distractions':Distractions, \
             'PercentDistracted':PercentDistracted})



def subsetter2(dictionary, dates, dis=False, verbose=False):
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
        onsets,dis_type = medfilereader(path, ['d', 'j'], remove_var_header = True)  # e onset, f offset

        if dis == True:
            if dictionary['Date'][ind] in dates:
                subset.append([onsets, dis_type, dictionary['RatID'][ind]])
                
        elif dis==False:   
            if dictionary['Date'][ind] in dates:
                subset.append([onsets, dictionary['RatID'][ind]])
            
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
        
def lickanalysis2(lickdata, burstThreshold=0.25, runThreshold=10):
    '''
    LICKANALYSIS KP
    # Takes lickdata from subset lists/dictionary
    # produces 25 item dictionary for each rat / day 
    # lick information on bursts, clusters(runs) etc. 
    
    '''
    analysis = []
    for lists in lickdata:
        lickon = lists[0]

        lick_analysis = lickCalc(lickon, burstThreshold=burstThreshold, runThreshold=runThreshold)
        analysis.append(lick_analysis)
    return analysis



def lickCalc(licks, offset = [], burstThreshold = 0.25, runThreshold = 10, 
             binsize=60, histDensity = False):
    """
    This function will calculate data for bursts from a train of licks. The threshold
    for bursts and clusters can be set. It returns all data as a dictionary.
    """
    
    # makes dictionary of data relating to licks and bursts
    if type(licks) != np.ndarray or type(offset) != np.ndarray:
        try:
            licks = np.array(licks)
            offset = np.array(offset)
        except:
            print('Licks and offsets need to be arrays and unable to easily convert.')
            return

    lickData = {}

    if len(licks) == len(offset) + 1:
        licks = licks[0:-1]
    
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
        lickData['bMed'] = np.median(lickData['bLicks'])
    else:
        lickData['bMean'] = 0
        lickData['bMed'] = 0 
    
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
    
    if lickData['rNum'] > 0:
        lickData['rMean'] = np.nanmean(lickData['rLicks'])
        lickData['rMed'] = np.median(lickData['rLicks'])
    else:
        lickData['rMean'] = 0
        lickData['rMed'] = 0 

    lickData['rILIs'] = [x for x in lickData['ilis'] if x > runThreshold]
    try:
        lickData['hist'] = np.histogram(lickData['licks'][1:], 
                                    range=(0, 3600), bins=int((3600/binsize)),
                                          density=histDensity)[0]
    except TypeError:
        print('Problem making histograms of lick data')
        
    return lickData 


def grouped_lickanalysis(groupdicts):
    '''
    GROUPED_LICKANALYSIS 
    
    Takes list of dictionaries previously sorted by subsetter
    Finds lick analysis information on bursts, clusters, intervals 
    taken from individual lick analysis 
    Runs need to be defined somewhere????
    '''

    all_n_bursts, all_n_runs, all_mean_IBI, all_mean_burst_length, \
        all_mean_IRI, all_mean_run_length = [], [], [], [], [], []
    for dictionary in groupdicts:     
      
        n_bursts = dictionary['bNum']
        n_runs = dictionary['rNum']
        #Mean ILI for each burst for each rat then caclulate a mean of mean for the groups
        mean_inter_burst = np.mean(dictionary['bILIs']) 
        mean_burst_length = dictionary['bMean'] # bMean uses bLicks (n licks not ILIs)
        mean_inter_run = np.mean(dictionary['rILIs'])
        mean_run_length = dictionary['rMean']
        # median burst lengths, median inter-burst-intervals (all measures with medians)
        all_n_bursts.append(n_bursts)
        all_n_runs.append(n_runs)
        all_mean_IBI.append(mean_inter_burst)
        all_mean_burst_length.append(mean_burst_length) # rename this variable 
        all_mean_IRI.append(mean_inter_run)
        all_mean_run_length.append(mean_run_length)
    # Can use these means to make plots, use the full lists to do statistics 
        # comparing saline to pcp for each variable - is there a difference between 
        # the numbers of bursts, the IBIs the runs etc. in sal and pcp (m then f)    

    mean_n_bursts = np.mean(all_n_bursts)
    mean_n_runs = np.mean(all_n_runs)
    mean_mean_IBI = np.mean(all_mean_IBI)
    mean_mean_IRI = np.mean(all_mean_IRI)
    
    return mean_n_bursts, mean_n_runs, mean_mean_IBI, mean_mean_IRI,\
    all_n_bursts, all_n_runs, all_mean_IBI, all_mean_IRI, all_mean_burst_length, all_mean_run_length 
    

def discalc_modalities2(dictionary, modalitykey):
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
        print('a')
        
        discalc = distractionCalc2(rat[0])
        distracted, notdistracted = distractedOrNot(discalc, rat[0])
      #  print(len(distracted), len(notdistracted))
      #  work out percentage and add this too 
        discalcgroup.append([distracted, notdistracted])
    
        dis_numeric = []
        ndis_numeric = []
        
    # Modality analysis - calculates which distractors contain different features (whitenoise, tone or combination)
    # Then works out on how many of these trials rats are distracted (individual) before creating a mean 
        for d in distracted:
            dis_numeric.append([rat[1][idx] for idx, val in enumerate(discalc) if val == d][0])
        for nd in notdistracted:
            ndis_numeric.append([rat[1][idx] for idx, val in enumerate(discalc) if val == nd][0])   
        # Makes the distracted trial types into integers 
        dis_numeric = [int(d) for d in dis_numeric]
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
 
## DEBUGGING SECTION        

# Not an issue if there are zero distracted trials, because calculates 0%
# without problem, will be an issue if both 0 dis and 0 non-dis (because no 
# trials of that modality, need to exclude)
         
#        if d_whitenoise_count < 1:
#            print('y', d_whitenoise_count, nd_whitenoise_count)
#        if d_tone_count < 1:
#            print('x', d_tone_count, nd_tone_count)
#        if d_combined_count < 1:
#            print('z', d_combined_count, nd_combined_count)
 
# WHITENOISE        
        if nd_whitenoise_count > 0:            
            percent_distracted_whitenoise = d_whitenoise_count / (d_whitenoise_count + nd_whitenoise_count) *100
        elif d_whitenoise_count < 1:
            try:
                percent_distracted_whitenoise = d_whitenoise_count / (d_whitenoise_count + nd_whitenoise_count) *100
            except ZeroDivisionError:
                print('oops')
                percent_distracted_whitenoise = -1
        elif d_whitenoise_count > 1:
            percent_distracted_whitenoise = 100 
# TONE            
        ## what to do if both zero 
        if nd_tone_count > 0:
            percent_distracted_tone = d_tone_count / (d_tone_count + nd_tone_count) *100
        elif d_tone_count < 1:
            try:
                percent_distracted_tone = d_tone_count / (d_tone_count + nd_tone_count) *100
            except ZeroDivisionError:
                print('oopsT')
                percent_distracted_tone = -1
        elif d_tone_count > 1:
            percent_distracted_tone = 100 

# COMBINED        
        if nd_combined_count > 0:
            percent_distracted_combined = d_combined_count / (d_combined_count + nd_combined_count) *100  

        elif d_combined_count < 1:
            try:
                percent_distracted_combined = d_combined_count / (d_combined_count + nd_combined_count) *100  
            except ZeroDivisionError:
                print('oopsC')
                percent_distracted_combined = -1
        elif d_combined_count > 1:
            percent_distracted_combined = 100 

# NON WHITE NOISE

        if nd_combined_count > 0 or nd_tone_count > 0:
            percent_distracted_all_non_WN = (d_tone_count + d_combined_count) / ((d_tone_count + d_combined_count )+ (nd_tone_count + nd_combined_count)) *100  
                
        elif d_combined_count < 1 and d_tone_count < 1:
            print('all less')
            try:
                percent_distracted_all_non_WN = (d_tone_count + d_combined_count) / ((d_tone_count + d_combined_count )+ (nd_tone_count + nd_combined_count)) *100  
            except ZeroDivisionError:
                print('oopsNWN')
                percent_distracted_all_non_WN = -1     
        elif d_combined_count or d_tone_count > 1:
            percent_distracted_all_non_WN = 100
     
      
        percent_dis_whitenoise_group.append(percent_distracted_whitenoise)
        percent_dis_tone_group.append(percent_distracted_tone)
        percent_dis_combined_group.append(percent_distracted_combined)
        percent_dis_all_non_WN_group.append(percent_distracted_all_non_WN)
    
## Removing the occurrences of -1 which were added to account for zero distractors
    
    percent_dis_whitenoise_group = [x for x in percent_dis_whitenoise_group if x != -1]
    percent_dis_tone_group = [x for x in percent_dis_tone_group if x != -1]
    percent_dis_combined_group = [x for x in percent_dis_combined_group if x != -1]
    percent_dis_all_non_WN_group = [x for x in percent_dis_all_non_WN_group if x != -1]       

        
    mean_percent_WHITENOISE = np.mean(percent_dis_whitenoise_group) # the average percentage of JUST whitenoise trials that rats are distracted on 
    mean_percent_TONE = np.mean(percent_dis_tone_group)
    mean_percent_COMBINED = np.mean(percent_dis_combined_group)
    mean_percent_ALL_NON_WN = np.mean(percent_dis_all_non_WN_group)
    
    return discalcgroup, percent_dis_whitenoise_group, percent_dis_tone_group, \
            percent_dis_combined_group, percent_dis_all_non_WN_group, mean_percent_WHITENOISE, mean_percent_TONE, \
            mean_percent_COMBINED, mean_percent_ALL_NON_WN



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
    if len(d) > 1:
        if d[-1] > 3599:
            d = d[:-1]
        
#    print(len(d))
    
    return d

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


def disbygroup(dictionary):
    ''' Prodcues times of distracted and not distracted as 2 lists
    takes a dictionary of grouped rat data 
    '''

    dis = []
    for rat in dictionary:
        
        discalc = distractionCalc2(rat[0])         
        distracted, notdistracted = distractedOrNot(discalc, rat[0])
        dis.append([distracted, notdistracted])
        
    return dis

def pdpbygroup(distractiondict, groupdict):
    ''' Distraction dict = discalc for the chosen group 
        groupdict is the whole disctionary on distraction 
        day from subsetter
    
    '''
    pdps_dis_group, med_pdps_dis_group, preDPs_dis_group, \
    pdps_notdis_group, med_pdps_notdis_group, preDPs_notdis_group = [],[],[],[],[],[]
    
    for index, rat in enumerate(distractiondict):
        pdps_dis = []
        preDPs_dis = []
        for distractorlick in rat[0]:
            
            if distractorlick in groupdict[index][0] and distractorlick != groupdict[index][0][-1]:
                lick_index = groupdict[index][0].index(distractorlick) 
                lick_index_plus1 = lick_index+1
                lick_index_minus3 = lick_index-3
                distracted_PDP = groupdict[index][0][lick_index_plus1] - groupdict[index][0][lick_index]
                distracted_preDP = groupdict[index][0][lick_index] - groupdict[index][0][lick_index_minus3]
            
            pdps_dis.append(distracted_PDP)
            preDPs_dis.append(distracted_preDP)
        pdps_dis_group.append(pdps_dis)
        med_pdps_dis_group.append(np.median(pdps_dis))
        preDPs_dis_group.append(preDPs_dis)
    
    
    # Not distracted PDPs 
    
    for index, rat in enumerate(distractiondict):
        pdps_notdis = []
        preDPs_notdis = []
        for notdistractedlick in rat[1]:
            if notdistractedlick in groupdict[index][0] and notdistractedlick != groupdict[index][0][-1]:
                lick_index = groupdict[index][0].index(notdistractedlick) 
                lick_index_plus1 = lick_index+1
                lick_index_minus3 = lick_index-3
                notdistracted_PDP = groupdict[index][0][lick_index_plus1] - groupdict[index][0][lick_index]
                notdistracted_preDP = groupdict[index][0][lick_index] - groupdict[index][0][lick_index_minus3]
            
            pdps_notdis.append(notdistracted_PDP)
            preDPs_notdis.append(notdistracted_preDP)
        pdps_notdis_group.append(pdps_notdis)
        med_pdps_notdis_group.append(np.median(pdps_notdis))
        preDPs_notdis_group.append(preDPs_notdis)
    
    
    return pdps_dis_group, med_pdps_dis_group, preDPs_dis_group, \
        pdps_notdis_group, med_pdps_notdis_group, preDPs_notdis_group
      
''' Function calculates the median values for ALL pdps, 
    not separated by distracted or not distracted, but as one list
'''

def pdp_mean_calc(pdp_dis_list, pdp_not_dis_list): # mostly variables from pdps by group func
    
    all_pdps = []
    all_pdps_mean = []
    for index, val in enumerate(pdp_dis_list):
        ongoing_all_PDPs = pdp_dis_list[index] + pdp_not_dis_list[index]
        all_pdps.append(ongoing_all_PDPs)
        mean = np.mean(ongoing_all_PDPs)
        all_pdps_mean.append(mean)    
        
    return all_pdps, all_pdps_mean



def percentdisgroup(distractiondict):
    ''' Discalc_sal_M == distractiondict '''
    
    percent_dis_group = []
    for rat in distractiondict: 
        percentage = len(rat[0]) / (len(rat[0])+len(rat[1])) * 100
        percent_dis_group.append(percentage)
    return percent_dis_group

    
def nlicksgrouped(licklist1, licklist2, licklist3):
    """ Calculates n licks from lists already subset, for the 3 
        training days before distraction. Subset list expected to 
        contain n lists by group, n = rat
    """
    lick = []
    lickminus1 = []
    lickminus2 = []

    
    for rat in licklist1:
        lick.append(len(rat[0]))
        
    for rat in licklist2:
        lickminus1.append(len(rat[0]))
        
    for rat in licklist3:
        lickminus2.append(len(rat[0]))
        


    return lick, lickminus1, lickminus2   

def cumulativelickFig(ax, firstlick, normed=True, color='g', log=True, ylabel='none', xlabel='none', title='none'):
    sorted_data = np.sort(firstlick)
    yvals = np.arange(len(sorted_data)+1)
    
    if normed == True:
        nlicks = len(sorted_data)
        yvals =yvals/nlicks
        
    a = ax.step(np.concatenate([sorted_data, sorted_data[[-1]]]),
             yvals, color=color)
    
    if log == True:
        ax.set_xscale('log', basex=10) # option in funcion as True or False (log setting)
        ax.set_xlim([0.1, 1000])
        
        # Label axes
#    if ylabel != 'none':
#        plt.ylabel(ylabel, fontsize=14)
#    
#    if xlabel != 'none':
#        plt.xlabel(xlabel)
#        
#    if title != 'none':
#        plt.title(title, fontsize=14)
#        
    return ax, a


 
#############################################################################
#############################################################################
#############################################################################                
    
# DATA PROCESSING AND ANALYSIS    
    
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
# LICK ANALYSIS DIS1 (Pilot)
################################################################################
# (1) Find and extract info from the metafile
#
metafile_data = '/Volumes/KP_HARD_DRI/kp259/DIS1/DIS1Masterfile.csv'
extract_data = MetaExtractor(metafile_data)
# Folder with all medfiles 
medfolder = '/Volumes/KP_HARD_DRI/kp259/DIS1/' 

# (2) - Get all lick onset data for all rats/sessions - dis 1 have no offsets
''' 
  Info DIS1 cohort pilot experiment :
      
  DATE  |CONDITION
  ------|------------
  170225|last lick
  170226|dis
  170227|hab 1
  170228|hab 2
  170301|saline IP
  170302|sal hab
  170304|amphetamine
  170305|amph hab
  170307|nicotine
  
'''
# Subsetting all data by day / drug 
# MALES ***********************************************************************
# SALINE 

last_lick = subsetter2(extract_data, ['170225'])
lick_minus1 = subsetter2(extract_data, ['170224'])
lick_minus2 = subsetter2(extract_data, ['170222'])
distraction = subsetter2(extract_data, ['170226'], dis=True)
hab1 = subsetter2(extract_data, ['170227'], dis=True)
hab2 = subsetter2(extract_data, ['170228'], dis=True)
salineIP = subsetter2(extract_data, ['170301'], dis=True)
amph = subsetter2(extract_data, ['170304'], dis=True)

# (3) Lick calc for last lick day 

# assign empty variables to store outputs from lick calc (to find means/groups stats)
# lists where each item is a dictionary (25) derived from lickCalc for each rat / day
lick_analysis = lickanalysis2(last_lick)

# ***********************************************************************************!!!!!    
## LICK DAY ANALYSIS  
# Produce medians/means for individual rats and group means (just one group DIS1) 

mean_n_bursts, mean_n_runs, mean_mean_IBI, mean_mean_IRI,\
all_n_bursts, all_n_runs, all_mean_IBI, all_mean_IRI, \
all_mean_burst_length, all_mean_run_length = grouped_lickanalysis(lick_analysis)

################################################################################
# DISTRACTION ANALYSIS DIS1 - MODALITIES

################################################################################
# Distraction day analysis (including modalities)
modalitykey = {'whitenoise':[1,4], 'tone':[2,5], 'combined3':[3,6]}
discalc, percent_dis_whitenoise, percent_dis_tone,\
percent_dis_combined, percent_dis_all_non_WN, mean_percent_WHITENOISE, mean_percent_TONE,\
mean_percent_COMBINED, mean_percent_ALL_NON_WN = discalc_modalities2(distraction, modalitykey)


# Modelled distractors ˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚
# Not including modalities just where distractors occur and dis vs notdis
# Modelled distractors by group (last lick day)    
mod_dis = disbygroup(last_lick)
# Distraction day
dis_dis = disbygroup(distraction)
# Habituation days by group 
hab1_dis = disbygroup(hab1)
hab2_dis = disbygroup(hab2)
# Saline day
sal_dis = disbygroup(salineIP)
# Amphetamine days by group 
amph_dis = disbygroup(amph)


# POST DISTRACTION PAUSES 
# For both distracted and non-distracted trials 
###################################################################################
'''
Info: Structure of the calculated distractors lists 
discalc[0][0][0] # [rat][list][licktimestamp]
'''

# Distraction day 
pdps_dis, med_pdps_dis, preDPs_dis,\
pdps_notdis, med_pdps_notdis, preDPs_notdis,\
 = pdpbygroup(discalc, distraction) 

# Modelled last lick day (pdps_mod)
pdps_mod_dis, med_pdps_mod_dis, preDPs_mod_dis,\
pdps_mod_notdis, med_pdps_mod_notdis, preDPs_mod_notdis,\
= pdpbygroup(mod_dis, last_lick)
#pdps_hab1
pdps_hab1_dis, med_pdps_hab1_dis, preDPs_hab1_dis,\
pdps_hab1_notdis, med_pdps_hab1_notdis, preDPs_hab1_notdis,\
= pdpbygroup(hab1_dis, hab1) 
#pdps_hab2
pdps_hab2_dis, med_pdps_hab2_dis, preDPs_hab2_dis,\
pdps_hab2_notdis, med_pdps_hab2_notdis, preDPs_hab2_notdis,\
= pdpbygroup(hab2_dis, hab2)  

#pdps_salineIP
pdps_sal_dis, med_pdps_sal_dis, preDPs_sal_dis,\
pdps_sal_notdis, med_sal_hab2_notdis, preDPs_sal_notdis,\
= pdpbygroup(sal_dis, salineIP) 

#pdps_amph
pdps_amph_dis, med_pdps_amph_dis, preDPs_amph_dis,\
pdps_amph_notdis, med_pdps_amph_notdis, preDPs_amph_notdis,\
= pdpbygroup(amph_dis, amph)  


#### PDPs all - not separated by distracted and not distracted 
## Medians for all 

# SALINE MALES - DONE 
all_pdps_dis , all_pdps_mean = pdp_mean_calc(pdps_dis, pdps_notdis)
all_pdps_mod, all_pdps_mean_mod = pdp_mean_calc(pdps_mod_dis, pdps_mod_notdis)
all_pdps_hab1 , all_pdps_mean_hab1 = pdp_mean_calc(pdps_hab1_dis, pdps_hab1_notdis)
all_pdps_hab2 , all_pdps_mean_hab2 = pdp_mean_calc(pdps_hab2_dis, pdps_hab2_notdis)
all_pdps_sal , all_pdps_mean_sal = pdp_mean_calc(pdps_sal_dis, pdps_sal_notdis)
all_pdps_amph , all_pdps_mean_amph = pdp_mean_calc(pdps_amph_dis, pdps_amph_notdis)

##################################################################
# Percentage distracted    ##################################################################
# Saline Males
percent_dis_dis = percentdisgroup(discalc)
percent_dis_modelled = percentdisgroup(mod_dis)
percent_dis_hab1 = percentdisgroup(hab1_dis)
percent_dis_hab2 = percentdisgroup(hab2_dis)
percent_dis_sal = percentdisgroup(sal_dis)
percent_dis_amph = percentdisgroup(amph_dis)

# Using subset lists, finds N licks for the last 3 saccharin training days 
# Stores as list variables for later plotting     
nlicks, nlicks_minus1, nlicks_minus2 = nlicksgrouped(last_lick, lick_minus1, lick_minus2)

 
   
  
# ****************************************************************************
# ****************************************************************************

### Cumulative licking plots - Post distraction pauses for DIS1 animals 
## Code taken from the THPH1 and 2 analysis 

#import seaborn as sn
#sn.set_style("white")

# Plot settings, font / size / styles
Calibri = {'fontname':'Calibri'}
Size = {'fontsize': 22}
label_size = 18
plt.rcParams['xtick.labelsize'] = label_size 
plt.rcParams['ytick.labelsize'] = label_size 
plt.rcParams['lines.linewidth'] = 2

# CUMULATIVE PDPS plots 
## LAST LICK DAY AND DISTRACTION DAY PDPS ALL 
fig = plt.figure()
#plt.title('Lickday SAL_M', **Calibri, **Size)
ax = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Plots all for lick day with average
## MODELLED (adds all dis and not dis PDPs together so no zero index issues) 
all_pdps_mod = []
for index, pdplists in enumerate(pdps_mod_dis):    
    C = pdps_mod_dis[index] + pdps_mod_notdis[index]
    all_pdps_mod.append(C)
# Plots all for MODELLED day with average 
for index, licklist in enumerate(all_pdps_mod):
    plot = cumulativelickFig(ax, all_pdps_mod[index], normed=True, color='lightgrey', log=True)
avg = [item for rat in all_pdps_mod for item in rat] 
cumulativelickFig(ax, avg, normed=True, color='#0277BD', log=True)
ax.set(ylabel = 'Probability')
ax.yaxis.label.set_size(16)
ax.set(xlabel = 'post-distractor pause log(s)')
ax.xaxis.label.set_size(16)
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/Cumulative_allPDPs_mod.pdf', bbox_inches='tight')

## Distraction day  
# CUMULATIVE PDPS
fig = plt.figure()
#plt.title('Distraction SAL_M', **Calibri, **Size)
ax = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
  
all_pdps = []
for index, pdplists in enumerate(pdps_dis):    
    C = pdps_dis[index] + pdps_notdis[index]
    all_pdps.append(C)
# Plots all for MODELLED day with average 
for index, licklist in enumerate(all_pdps):
    plot = cumulativelickFig(ax, all_pdps[index], normed=True, color='lightgrey', log=True)
avg = [item for rat in all_pdps for item in rat] 
cumulativelickFig(ax, avg, normed=True, color='#880E4F', log=True)
#avg2 = [item for rat in pdps_dis_sal_M for item in rat] 
#cumulativelickFig(ax, avg2, normed=True, color='green', log=True)
#avg3 = [item for rat in pdps_notdis_sal_M for item in rat] 
#cumulativelickFig(ax, avg3, normed=True, color='blue', log=True)
ax.set(ylabel = 'Probability')
ax.yaxis.label.set_size(16)
ax.set(xlabel = 'post-distractor pause log(s)')
ax.xaxis.label.set_size(16)
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/Cumulative_allPDPs_dis.pdf', bbox_inches='tight')

### Both lines on one graph no individual rats at all 

fig = plt.figure()
#plt.title('PDPs by day and group', **Calibri, **Size)
ax = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
  

avg = [item for rat in all_pdps_mod for item in rat] 
avg1 = [item for rat in all_pdps for item in rat] 
cumulativelickFig(ax, avg, normed=True, color='#0277BD', log=True)
cumulativelickFig(ax, avg1, normed=True, color='#880E4F', log=True)

ax.set(ylabel = 'Probability')
ax.yaxis.label.set_size(16)
ax.set(xlabel = 'post-distractor pause log(s)')
ax.xaxis.label.set_size(16)

plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/Cumulative_allPDPs_modVSdis.pdf', bbox_inches='tight')



#############
#######################################################################

## Distraction day  - Distracted vs Not distracted
fig = plt.figure()
#plt.title('Dis vs Not dis M', **Calibri, **Size)
ax = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for index, licklist in enumerate(pdps_dis):
    plot = cumulativelickFig(ax, pdps_dis[index], normed=True, color='lightgrey', log=True)
avg = [item for rat in pdps_dis for item in rat] 
cumulativelickFig(ax, avg, normed=True, color='#E91E63', log=True)

for index, licklist in enumerate(pdps_notdis):
    plot = cumulativelickFig(ax, pdps_notdis[index], normed=True, color='lightgrey', log=True)
avg = [item for rat in pdps_notdis for item in rat] 
cumulativelickFig(ax, avg, normed=True, color='#4FC3F7', log=True)
ax.set(ylabel = 'Probability')
ax.yaxis.label.set_size(16)
ax.set(xlabel = 'post-distractor pause log(s)')
ax.xaxis.label.set_size(16)
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/Cumulative__disvsnot.pdf', bbox_inches='tight')


