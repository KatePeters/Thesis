#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 09:05:37 2018
@author: u1490431 - Kate Z Peters 

#############################################################################

CHAPTER 3 - DISTRACTION FROM ONGOING SACCHARIN CONSUMPTION: EXPERIMENTS
            IN SALINE AND PHENCYCLIDINE TREATED RATS

#############################################################################

This file contains the DATA PROCESSING and ANALYSIS script 
All functions and imports are at the beginning and no others are required
A separate PLOTTING script should be run second to produce chapter figures            

Imports:
    numpy as np  
Functions: 
    MetaExtractor(),subsetter(),medfilereader(),asnumeric(),lickanalysis()
    lickCalc(),grouped_lickanalysis(),discalc_modalities(),distractionCalc2()
    remcheck(),distractedOrNot(),disbygroup(),percentdisgroup()
"""

# Imports
import numpy as np
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
    NonDistractions, PercentDistracted = [], [], [], [], [], [], [], [], [], []

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
 
    return ({'MedFilenames':MedFilenames, 'RatID':RatID, 'Date':Date, 'Day':Day, 'Session':Session, \
             'Drug':Drug, 'TotLicks':TotLicks, 'Distractions':Distractions, \
             'PercentDistracted':PercentDistracted})



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
        onsets, offsets, med_dis_times, dis_type = medfilereader(path, ['e', 'f', 'i', 'j'], remove_var_header = True)  # e onset, f offset

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
        
def lickanalysis(lickdata, burstThreshold=0.25, runThreshold=10):
    '''
    LICKANALYSIS KP
    # Takes lickdata from subset lists/dictionary
    # produces 25 item dictionary for each rat / day 
    # lick information on bursts, clusters(runs) etc. 
    
    '''
    analysis = []
    for lists in lickdata:
        lickon = lists[0]
        offset = lists[1]
        lick_analysis = lickCalc(lickon, offset, burstThreshold=burstThreshold, runThreshold=runThreshold)
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
    
def discalc_modalities(dictionary, modalitykey):
    ''' Calculates distractors, distracted and modalities for dictionary of 
    rats by group for just distraction day only 
    
    Calculates a grouped percentage of each modality and how distracted 
    by that modality rats are on average (by group)
    
    '''
    percent_dis_whitenoise_group = []
    percent_dis_tone_group = []
    percent_dis_combined_group = []
    discalcgroup = []
    ## SAL MALES - DISTRACTION DAY ONLY - DISTRACTOR TYPE ANALYSIS INCLUDED
# Finds distracted or not (corrects for med slipping issue)
    for rat in dictionary:
        
        discalc = distractionCalc2(rat[0])
        distracted, notdistracted = distractedOrNot(discalc, rat[0])
        print(len(distracted), len(notdistracted))
      #  work out percentage and add this too 
        discalcgroup.append([distracted, notdistracted])
    
        dis_numeric = []
        ndis_numeric = []
    # Modality analysis - calculates which distractors contain different features (whitenoise, tone or combination)
    # Then works out on how many of these trials rats are distracted (individual) before creating a mean 
        for d in distracted:
            dis_numeric.append([rat[2][idx] for idx, val in enumerate(discalc) if val == d][0])
        for nd in notdistracted:
            ndis_numeric.append([rat[2][idx] for idx, val in enumerate(discalc) if val == nd][0])   
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
        
        #print(d_whitenoise_count, d_tone_count, d_combined_count)
        
### Added 20/07 to fix division by zero 
        if d_whitenoise_count != 0:
#            
            d_percent_white_noise = d_whitenoise_count / (len(dis_type_text))*100
            d_percent_tone = d_tone_count / (len(dis_type_text))*100
            d_percent_combined = d_combined_count / (len(dis_type_text))*100 
        else:
            d_percent_white_noise = 0 
            d_percent_tone = 0
            d_percent_combined = 0
            
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
                
        #print(nd_whitenoise_count, nd_tone_count, nd_combined_count)
        nd_percent_white_noise = nd_whitenoise_count / (len(ndis_type_text))*100
        nd_percent_tone = nd_tone_count / (len(ndis_type_text))*100
        nd_percent_combined =  nd_combined_count / (len(ndis_type_text))*100
        
        print(d_percent_white_noise, nd_percent_white_noise)
        
        percent_distracted_whitenoise = d_whitenoise_count / (d_whitenoise_count + nd_whitenoise_count) *100
        percent_distracted_tone = d_tone_count / (d_tone_count + nd_tone_count) *100
        percent_distracted_combined = d_combined_count / (d_combined_count + nd_combined_count) *100  
        
    # ?    percent_distracted_non_whitenoise = percent_distracted_tone  + percent_distracted_combined
 
        
        percent_dis_whitenoise_group.append(percent_distracted_whitenoise)
        percent_dis_tone_group.append(percent_distracted_tone)
        percent_dis_combined_group.append(percent_distracted_combined)
        
    mean_percent_WHITENOISE = np.mean(percent_dis_whitenoise_group) # the average percentage of JUST whitenoise trials that rats are distracted on 
    mean_percent_TONE = np.mean(percent_dis_tone_group)
    mean_percent_COMBINED = np.mean(percent_dis_combined_group)
    
    
    return discalcgroup, percent_dis_whitenoise_group, percent_dis_tone_group, \
            percent_dis_combined_group, mean_percent_WHITENOISE, mean_percent_TONE, \
            mean_percent_COMBINED

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
# LICK ANALYSIS DPCP (1,2,3)
################################################################################
# (1) Find and extract info from the metafile(s) 
#
metafile_males = '/Volumes/KP_HARD_DRI/kp259/DPCP_ALL/DPCP12Masterfile.csv'
metafile_females = '/Volumes/KP_HARD_DRI/kp259/DPCP_ALL/DPCP3_Metafile.csv'
extract_males = MetaExtractor(metafile_males)
extract_females = MetaExtractor(metafile_females)
# Folder with all medfiles (DPCP1, DPCP2, DPCP3)
medfolder = '/Volumes/KP_HARD_DRI/kp259/DPCP_ALL/' # sometimes need space1 after DRI

# (2) - Get all lick onset/offset data for all rats/sessions (make this a function later)
''' 
  Info DPCP1(16) and DPCP2(16) males, DPCP3(24) females :
      
  DPCP1 |DPCP2 |DPCP3 |CONDITION
  ------|------|------|-----------
  170417|171006|171124|last lick
  170418|171007|171125|dis
  170419|171008|171126|hab 1
  170420|171009|171127|hab 2
  170423|171012|171128|amphetamine
'''
# Subsetting all data by day / drug 
# MALES ***********************************************************************
# SALINE 

last_lick_sal_M = subsetter(extract_males, ['170417','171006'], 'SAL')
lick_minus1_sal_M = subsetter(extract_males, ['170416', '171004'], 'SAL')
lick_minus2_sal_M = subsetter(extract_males, ['170415', '171003'], 'SAL')
distraction_sal_M = subsetter(extract_males, ['170418','171007'], 'SAL', dis=True)
hab1_sal_M = subsetter(extract_males, ['170419','171008'], 'SAL', dis=True)
hab2_sal_M = subsetter(extract_males, ['170420','171009'], 'SAL', dis=True)
amph_sal_M = subsetter(extract_males, ['170423','171012'], 'SAL', dis=True)
# PCP
last_lick_pcp_M = subsetter(extract_males, ['170417','171006'], 'PCP')
lick_minus1_pcp_M = subsetter(extract_males, ['170416', '171004'], 'PCP')
lick_minus2_pcp_M = subsetter(extract_males, ['170415', '171003'], 'PCP')
distraction_pcp_M = subsetter(extract_males, ['170418','171007'], 'PCP', dis=True)
hab1_pcp_M = subsetter(extract_males, ['170419','171008'], 'PCP', dis=True)
hab2_pcp_M = subsetter(extract_males, ['170420','171009'], 'PCP', dis=True)
amph_pcp_M = subsetter(extract_males, ['170423','171012'], 'PCP', dis=True)
# FEMALES **********************************************************************
# SALINE
last_lick_sal_F = subsetter(extract_females, ['171124'], 'SAL')
lick_minus1_sal_F = subsetter(extract_females, ['171122'], 'SAL')
lick_minus2_sal_F = subsetter(extract_females, ['171121'], 'SAL')
lick_minus3_sal_F = subsetter(extract_females, ['171120'], 'SAL')
distraction_sal_F = subsetter(extract_females, ['171125'], 'SAL',dis=True)
hab1_sal_F = subsetter(extract_females, ['171126'], 'SAL',dis=True)
hab2_sal_F = subsetter(extract_females, ['171127'], 'SAL',dis=True)
amph_sal_F = subsetter(extract_females, ['171128'], 'SAL',dis=True)
# PCP
last_lick_pcp_F = subsetter(extract_females, ['171124'], 'PCP')
lick_minus1_pcp_F = subsetter(extract_females, ['171122'], 'PCP')
lick_minus2_pcp_F = subsetter(extract_females, ['171121'], 'PCP')
lick_minus3_pcp_F = subsetter(extract_females, ['171120'], 'PCP')
distraction_pcp_F = subsetter(extract_females, ['171125'], 'PCP',dis=True)
hab1_pcp_F = subsetter(extract_females, ['171126'], 'PCP',dis=True)
hab2_pcp_F = subsetter(extract_females, ['171127'], 'PCP',dis=True)
amph_pcp_F = subsetter(extract_females, ['171128'], 'PCP',dis=True)
 
# (3) Lick calc for last lick day (by group) for male PCP and SAL, for female PCP and SAL
# Licking analysis section, just last lick day 
# assign empty variables to store outputs from lick calc (to find means/groups stats)
# lists where each item is a dictionary (25) derived from lickCalc for each rat / day
lick_analysis_sal_M = lickanalysis(last_lick_sal_M)
lick_analysis_pcp_M = lickanalysis(last_lick_pcp_M)
lick_analysis_sal_F = lickanalysis(last_lick_sal_F)
lick_analysis_pcp_F = lickanalysis(last_lick_pcp_F)
# ***********************************************************************************!!!!!    
## LICK DAY ANALYSIS - BY GROUP 
# Produce medians/means for individual rats and group means 
# Males
# Saline
sal_M_mean_n_bursts, sal_M_mean_n_runs, sal_M_mean_mean_IBI, sal_M_mean_mean_IRI,\
all_n_bursts_sal_M, all_n_runs_sal_M, all_mean_IBI_sal_M, all_mean_IRI_sal_M, \
all_mean_burst_length_sal_M, all_mean_run_length_sal_M = grouped_lickanalysis(lick_analysis_sal_M)
# PCP
pcp_M_mean_n_bursts, pcp_M_mean_n_runs, pcp_M_mean_mean_IBI, pcp_M_mean_mean_IRI,\
all_n_bursts_pcp_M, all_n_runs_pcp_M, all_mean_IBI_pcp_M, all_mean_IRI_pcp_M, \
all_mean_burst_length_pcp_M, all_mean_run_length_pcp_M = grouped_lickanalysis(lick_analysis_pcp_M)
# Females
# Saline
sal_F_mean_n_bursts, sal_F_mean_n_runs, sal_F_mean_mean_IBI, sal_F_mean_mean_IRI,\
all_n_bursts_sal_F, all_n_runs_sal_F, all_mean_IBI_sal_F, all_mean_IRI_sal_F, \
all_mean_burst_length_sal_F, all_mean_run_length_sal_F = grouped_lickanalysis(lick_analysis_sal_F)
# PCP
pcp_F_mean_n_bursts, pcp_F_mean_n_runs, pcp_F_mean_mean_IBI, pcp_F_mean_mean_IRI,\
all_n_bursts_pcp_F, all_n_runs_pcp_F, all_mean_IBI_pcp_F, all_mean_IRI_pcp_F, \
all_mean_burst_length_pcp_F, all_mean_run_length_pcp_F = grouped_lickanalysis(lick_analysis_pcp_F)
################################################################################
# DISTRACTION ANALYSIS DPCP (1,2,3)
################################################################################
# Distraction day analysis (including modalities)
modalitykey = {'whitenoise':[1,4], 'tone':[2,5], 'combined3':[3,6]}
# MALES
# Saline
discalc_sal_M, percent_dis_whitenoise_sal_M, percent_dis_tone_sal_M,\
percent_dis_combined_sal_M, mean_percent_WHITENOISE_sal_M, mean_percent_TONE_sal_M,\
mean_percent_COMBINED_sal_M = discalc_modalities(distraction_sal_M, modalitykey)
# PCP
discalc_pcp_M, percent_dis_whitenoise_pcp_M, percent_dis_tone_pcp_M,\
percent_dis_combined_pcp_M, mean_percent_WHITENOISE_pcp_M, mean_percent_TONE_pcp_M,\
mean_percent_COMBINED_pcp_M = discalc_modalities(distraction_pcp_M, modalitykey)
# FEMALES
# Saline 
discalc_sal_F, percent_dis_whitenoise_sal_F, percent_dis_tone_sal_F,\
percent_dis_combined_sal_F, mean_percent_WHITENOISE_sal_F, mean_percent_TONE_sal_F,\
mean_percent_COMBINED_sal_F = discalc_modalities(distraction_sal_F, modalitykey)
# PCP
discalc_pcp_F, percent_dis_whitenoise_pcp_F, percent_dis_tone_pcp_F,\
percent_dis_combined_pcp_F, mean_percent_WHITENOISE_pcp_F, mean_percent_TONE_pcp_F,\
mean_percent_COMBINED_pcp_F = discalc_modalities(distraction_pcp_F, modalitykey)

# Habituation day 
# SALINE
discalc_hab1_sal_M, percent_dis_whitenoise_hab1_sal_M, percent_dis_tone_hab1_sal_M,\
percent_dis_combined_hab1_sal_M, mean_percent_WHITENOISE_hab1_sal_M, mean_percent_TONE_hab1_sal_M,\
mean_percent_COMBINED_hab1_sal_M = discalc_modalities(hab1_sal_M, modalitykey)

## ••• Need to calculate the mean of ALL others (combined plus tone)
# remember these are means of all trials for each rat then the main value is 
# the mean of all rats 
'''percent_dis_non_whitenoise_hab1_sal_M = percent_dis_tone_hab1_sal_M + percent_dis_combined_hab1_sal_M 
not correct because not the mean !! --> need to add to function 
'''


#PCP
discalc_hab1_pcp_M, percent_dis_whitenoise_hab1_pcp_M, percent_dis_tone_hab1_pcp_M,\
percent_dis_combined_hab1_pcp_M, mean_percent_WHITENOISE_hab1_pcp_M, mean_percent_TONE_hab1_pcp_M,\
mean_percent_COMBINED_hab1_pcp_M = discalc_modalities(hab1_pcp_M, modalitykey)



aa,bb,cc,dd,ee,ff,gg = discalc_modalities(hab1_pcp_M, modalitykey)


aa,bb,cc,dd,ee,ff,gg = discalc_modalities(hab2_sal_M, modalitykey)
aa,bb,cc,dd,ee,ff,gg = discalc_modalities(amph_sal_M, modalitykey)

# Habituation day 2

aa,bb,cc,dd,ee,ff,gg = discalc_modalities(hab2_pcp_M, modalitykey)
aa,bb,cc,dd,ee,ff,gg = discalc_modalities(amph_pcp_M, modalitykey)
# Amphetamine day 
aa,bb,cc,dd,ee,ff,gg = discalc_modalities(hab1_sal_M, modalitykey)
aa,bb,cc,dd,ee,ff,gg = discalc_modalities(hab2_sal_M, modalitykey)
aa,bb,cc,dd,ee,ff,gg = discalc_modalities(amph_sal_M, modalitykey)




# Modelled distractors ˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚
# Not including modalities just where distractors occur and dis vs notdis
# Modelled distractors by group (last lick day)    
mod_dis_sal_M = disbygroup(last_lick_sal_M)
mod_dis_pcp_M = disbygroup(last_lick_pcp_M)
mod_dis_sal_F = disbygroup(last_lick_sal_F)
mod_dis_pcp_F = disbygroup(last_lick_pcp_F)
# Habituation days by group 
hab1_dis_sal_M = disbygroup(hab1_sal_M)
hab2_dis_sal_M = disbygroup(hab2_sal_M)
hab1_dis_pcp_M = disbygroup(hab1_pcp_M)
hab2_dis_pcp_M = disbygroup(hab2_pcp_M)
hab1_dis_sal_F = disbygroup(hab1_sal_F)
hab2_dis_sal_F = disbygroup(hab2_sal_F)
hab1_dis_pcp_F = disbygroup(hab1_pcp_F)
hab2_dis_pcp_F = disbygroup(hab2_pcp_F)
# Amphetamine days by group 
amph_dis_sal_M = disbygroup(amph_sal_M)
amph_dis_pcp_M = disbygroup(amph_pcp_M)
amph_dis_sal_F = disbygroup(amph_sal_F)
amph_dis_pcp_F = disbygroup(amph_pcp_F)

# POST DISTRACTION PAUSES 
# For both distracted and non-distracted trials 

###################################################################################
'''
Info: Structure of the calculated distractors lists 
discalc_sal_M[0][0][0] # [rat][list][licktimestamp]
'''
# CALCULATE THE PDP FOR SPECIFIC GROUPS - finds where in the licks the distractor 
# occurred and then finds the pause before the next lick (ignoring distractors occurring
# on the final lick in a session)

# SALINE MALES - distraction day 
pdps_dis_sal_M, med_pdps_dis_sal_M, preDPs_dis_sal_M,\
pdps_notdis_sal_M, med_pdps_notdis_sal_M, preDPs_notdis_sal_M,\
 = pdpbygroup(discalc_sal_M, distraction_sal_M) 
# ISSUES HERE!!! NAN BECAUSE PDPS_MOD_DIS_SAL_M HAS 2 RATS WITH ZERO so
# med_pdps_mod_dis_sal_M is nan,nan, then fine  
#pdps_mod_sal_M - PDPs for modelled distraction day, modelled distractors
pdps_mod_dis_sal_M, med_pdps_mod_dis_sal_M, preDPs_mod_dis_sal_M,\
pdps_mod_notdis_sal_M, med_pdps_mod_notdis_sal_M, preDPs_mod_notdis_sal_M,\
= pdpbygroup(mod_dis_sal_M, last_lick_sal_M)
#pdps_hab1_sal_M
pdps_hab1_dis_sal_M, med_pdps_hab1_dis_sal_M, preDPs_hab1_dis_sal_M,\
pdps_hab1_notdis_sal_M, med_pdps_hab1_notdis_sal_M, preDPs_hab1_notdis_sal_M,\
= pdpbygroup(hab1_dis_sal_M, hab1_sal_M) 
#pdps_hab2_sal_M
pdps_hab2_dis_sal_M, med_pdps_hab2_dis_sal_M, preDPs_hab2_dis_sal_M,\
pdps_hab2_notdis_sal_M, med_pdps_hab2_notdis_sal_M, preDPs_hab2_notdis_sal_M,\
= pdpbygroup(hab2_dis_sal_M, hab2_sal_M)  
#pdps_amph_sal_M
pdps_amph_dis_sal_M, med_pdps_amph_dis_sal_M, preDPs_amph_dis_sal_M,\
pdps_amph_notdis_sal_M, med_pdps_amph_notdis_sal_M, preDPs_amph_notdis_sal_M,\
= pdpbygroup(amph_dis_sal_M, amph_sal_M)  

 
# PCP MALES - distraction day
pdps_dis_pcp_M, med_pdps_dis_pcp_M, preDPs_dis_pcp_M,\
pdps_notdis_pcp_M, med_pdps_notdis_pcp_M, preDPs_notdis_pcp_M,\
 = pdpbygroup(discalc_pcp_M, distraction_pcp_M)
#pdps_mod_pcp_M
pdps_mod_dis_pcp_M, med_pdps_mod_dis_pcp_M, preDPs_mod_dis_pcp_M,\
pdps_mod_notdis_pcp_M, med_pdps_mod_notdis_pcp_M, preDPs_mod_notdis_pcp_M,\
= pdpbygroup(mod_dis_pcp_M, last_lick_pcp_M)
#pdps_hab1_pcp_M
pdps_hab1_dis_pcp_M, med_pdps_hab1_dis_pcp_M, preDPs_hab1_dis_pcp_M,\
pdps_hab1_notdis_pcp_M, med_pdps_hab1_notdis_pcp_M, preDPs_hab1_notdis_pcp_M,\
= pdpbygroup(hab1_dis_pcp_M, hab1_pcp_M) 
#pdps_hab2_pcp_M
pdps_hab2_dis_pcp_M, med_pdps_hab2_dis_pcp_M, preDPs_hab2_dis_pcp_M,\
pdps_hab2_notdis_pcp_M, med_pdps_hab2_notdis_pcp_M, preDPs_hab2_notdis_pcp_M,\
= pdpbygroup(hab2_dis_pcp_M, hab2_pcp_M)  
#pdps_amph_pcp_M
pdps_amph_dis_pcp_M, med_pdps_amph_dis_pcp_M, preDPs_amph_dis_pcp_M,\
pdps_amph_notdis_pcp_M, med_pdps_amph_notdis_pcp_M, preDPs_amph_notdis_pcp_M,\
= pdpbygroup(amph_dis_pcp_M, amph_pcp_M) 
 
 
# SALINE FEMALES - distraction day
pdps_dis_sal_F, med_pdps_dis_sal_F, preDPs_dis_sal_F,\
pdps_notdis_sal_F, med_pdps_notdis_sal_F, preDPs_notdis_sal_F,\
 = pdpbygroup(discalc_sal_F, distraction_sal_F) 
#pdp_mod_sal_F 
pdps_mod_dis_sal_F, med_pdps_mod_dis_sal_F, preDPs_mod_dis_sal_F,\
pdps_mod_notdis_sal_F, med_pdps_mod_notdis_sal_F, preDPs_mod_notdis_sal_F,\
= pdpbygroup(mod_dis_sal_F, last_lick_sal_F)
#pdps_hab1_sal_F
pdps_hab1_dis_sal_F, med_pdps_hab1_dis_sal_F, preDPs_hab1_dis_sal_F,\
pdps_hab1_notdis_sal_F, med_pdps_hab1_notdis_sal_F, preDPs_hab1_notdis_sal_F,\
= pdpbygroup(hab1_dis_sal_F, hab1_sal_F) 
#pdps_hab2_sal_F
pdps_hab2_dis_sal_F, med_pdps_hab2_dis_sal_F, preDPs_hab2_dis_sal_F,\
pdps_hab2_notdis_sal_F, med_pdps_hab2_notdis_sal_F, preDPs_hab2_notdis_sal_F,\
= pdpbygroup(hab2_dis_sal_F, hab2_sal_F)  
#pdps_amph_sal_F
pdps_amph_dis_sal_F, med_pdps_amph_dis_sal_F, preDPs_amph_dis_sal_F,\
pdps_amph_notdis_sal_F, med_pdps_amph_notdis_sal_F, preDPs_amph_notdis_sal_F,\
= pdpbygroup(amph_dis_sal_F, amph_sal_F)

 
# PCP FEMALES - distraction day 
pdps_dis_pcp_F, med_pdps_dis_pcp_F, preDPs_dis_pcp_F,\
pdps_notdis_pcp_F, med_pdps_notdis_pcp_F, preDPs_notdis_pcp_F,\
 = pdpbygroup(discalc_pcp_F, distraction_pcp_F) 
 
#pdps_mod_pcp_F
pdps_mod_dis_pcp_F, med_pdps_mod_dis_pcp_F, preDPs_mod_dis_pcp_F,\
pdps_mod_notdis_pcp_F, med_pdps_mod_notdis_pcp_F, preDPs_mod_notdis_pcp_F,\
= pdpbygroup(mod_dis_pcp_F, last_lick_pcp_F)
#pdps_hab1_pcp_F
pdps_hab1_dis_pcp_F, med_pdps_hab1_dis_pcp_F, preDPs_hab1_dis_pcp_F,\
pdps_hab1_notdis_pcp_F, med_pdps_hab1_notdis_pcp_F, preDPs_hab1_notdis_pcp_F,\
= pdpbygroup(hab1_dis_pcp_F, hab1_pcp_F) 
#pdps_hab2_pcp_F
pdps_hab2_dis_pcp_F, med_pdps_hab2_dis_pcp_F, preDPs_hab2_dis_pcp_F,\
pdps_hab2_notdis_pcp_F, med_pdps_hab2_notdis_pcp_F, preDPs_hab2_notdis_pcp_F,\
= pdpbygroup(hab2_dis_pcp_F, hab2_pcp_F)  
#pdps_amph_pcp_F
pdps_amph_dis_pcp_F, med_pdps_amph_dis_pcp_F, preDPs_amph_dis_pcp_F,\
pdps_amph_notdis_pcp_F, med_pdps_amph_notdis_pcp_F, preDPs_amph_notdis_pcp_F,\
= pdpbygroup(amph_dis_pcp_F, amph_pcp_F) 
  
'''
# Corelations 
Could have a plot with blue when distracted and grey when not, then correlate with 
the pre dp and the pdp? 
'''

##################################################################
# Percentage distracted    ##################################################################
# Saline Males
percent_dis_dis_sal_M = percentdisgroup(discalc_sal_M)
percent_dis_modelled_sal_M = percentdisgroup(mod_dis_sal_M)
percent_dis_hab1_sal_M = percentdisgroup(hab1_dis_sal_M)
percent_dis_hab2_sal_M = percentdisgroup(hab2_dis_sal_M)
percent_dis_amph_sal_M = percentdisgroup(amph_dis_sal_M)
## PCP Males
percent_dis_dis_pcp_M = percentdisgroup(discalc_pcp_M)
percent_dis_modelled_pcp_M = percentdisgroup(mod_dis_pcp_M)
percent_dis_hab1_pcp_M = percentdisgroup(hab1_dis_pcp_M)
percent_dis_hab2_pcp_M = percentdisgroup(hab2_dis_pcp_M)
percent_dis_amph_pcp_M = percentdisgroup(amph_dis_pcp_M)

############# FEMALES - percent distracted 
## Remember might have an issue with the last PCP rat on amphetamine day
# division by zero potential problem/ Might remove this rat from ALL days for the plots?? 
# SALINE
percent_dis_dis_sal_F = percentdisgroup(discalc_sal_F)
percent_dis_modelled_sal_F = percentdisgroup(mod_dis_sal_F)
percent_dis_hab1_sal_F = percentdisgroup(hab1_dis_sal_F)
percent_dis_hab2_sal_F = percentdisgroup(hab2_dis_sal_F)
percent_dis_amph_sal_F = percentdisgroup(amph_dis_sal_F)
# PCP

### EDITED - 18/06/18 removed [0:-1] for all but amphetamine and added a 0 for amph
percent_dis_dis_pcp_F = percentdisgroup(discalc_pcp_F)
percent_dis_modelled_pcp_F = percentdisgroup(mod_dis_pcp_F)
percent_dis_hab1_pcp_F = percentdisgroup(hab1_dis_pcp_F)
percent_dis_hab2_pcp_F = percentdisgroup(hab2_dis_pcp_F)
percent_dis_amph_pcp_F = percentdisgroup(amph_dis_pcp_F[0:-1])

percent_dis_amph_pcp_F.extend([0]) 


# Using subset lists, finds N licks for the last 3 saccharin training days 
# Stores as list variables for later plotting     
# MALES
nlicks_sal_M, nlicks_minus1_sal_M, nlicks_minus2_sal_M = nlicksgrouped(last_lick_sal_M, lick_minus1_sal_M, lick_minus2_sal_M)
nlicks_pcp_M, nlicks_minus1_pcp_M, nlicks_minus2_pcp_M = nlicksgrouped(last_lick_pcp_M, lick_minus1_pcp_M, lick_minus2_pcp_M)
# FEMALES
nlicks_sal_F, nlicks_minus1_sal_F, nlicks_minus2_sal_F = nlicksgrouped(last_lick_sal_F, lick_minus1_sal_F, lick_minus2_sal_F)
nlicks_pcp_F, nlicks_minus1_pcp_F, nlicks_minus2_pcp_F = nlicksgrouped(last_lick_pcp_F, lick_minus1_pcp_F, lick_minus2_pcp_F)

######################## INDIVIDUAL DIFFERENCES #######################
#sb.jointplot(x=df['nRuns'], y=df['nBursts'], kind='hex')) #or type 'reg' for kernel estimation and regression

################################################################################
# NOR score data, Males and Females saline and PCP

# Added the data in manually and created variables  
# MALES SALINE
# Acquisition
exploration_left_sal_M = [27.18,26.69,19.61,18.17,16.09,30.06,36.54,38.87,24.76,26.86,34.15,0,30.57,38.28]
exploration_right_sal_M = [52.44,21.21,33,35.8,37.41,25.45,31.2,15.71,20.07,28.24,27.82,6.97,40.93,14.52]
# Retention
exploration_fam_sal_M = [28.8,13.56,27.5,11.81,9.06,14.66,21.8,20.62,5.76,17.1,15.16,0,6.18,17.46]
exploration_nov_sal_M = [29.5,17.99,28.96,26.71,16.97,31.82,22.29,20.03,29.02,19.79,18.7,25.74,57.62,29.96]
DI_sal_M = [0.01,0.14,0.03,0.39,0.3,0.37,0.01,-0.01,0.67,0.07,0.1,1,0.81,0.26]
# MALES PCP
# Acquisition
exploration_left_pcp_M = [19.85,29.53,25.52,36.45,19.56,10.28,31.46,12.52,16.28,19.52,27.91,25.59,23.98,29.35,19.89,23.88]
exploration_right_pcp_M = [30.4,22.96,13.76,18.12,20.22,21.21,24.68,21.6,49.08,25.87,9.76,39.46,11.43,20.13,16.41,27.11]
# Retention
exploration_fam_pcp_M = [6.49,7.89,13.81,26.66,29.22,9.28,17.39,18.82,22.95,20.86,11.55,35.84,0.9,13.38,5.68,7.84]
exploration_nov_pcp_M = [13.32,23.52,35.97,8.81,18.53,9.59,38.83,18.02,22.75,44.65,29.64,24.91,2,27.41,15.21,27.35]
DI_pcp_M = [0.34,0.5,0.45,-0.5,-0.22,0.02,0.38,-0.02,0,0.36,0.44,-0.18,0.38,0.34,0.46,0.55]

# FEMALES SALINE
# Acquisition
exploration_left_sal_F = [15.71,16.69,20.77,19.44,30.25,28.73,20.86,14.02,31.13,17.37,29.88,15.51]
exploration_right_sal_F = [27.49,18.7,22.91,18.7,41.41,23.47,16.43,15.47,36.04,23.84,42.94,32.35]
# Retention
exploration_fam_sal_F = [13.47,26.76,9.17,13.63,15.12,7.13,9.29,8.98,14.19,12.2,14.92,14.53]
exploration_nov_sal_F = [13.88,27.04,19.13,16.82,22.45,23.01,12.09,28.08,16.4,17.04,26.5,24.03]
DI_sal_F = [0.02,0.01,0.35,0.1,0.19,0.53,0.13,0.52,0.07,0.17,0.28,0.25]

# FEMALES PCP
# Acquisition
exploration_left_pcp_F = [13.01,27.87,15.54,22.77,16.69,24.87,15.75,20.4,14.23,14.28,11.7,17.37]
exploration_right_pcp_F = [13.83,31.05,40.14,19.66,16.74,15.4,2.02,12.96,23.94,21.65,12.81,11.65]
# Retention
exploration_fam_pcp_F = [13.37,27.41,42.4,10.4,14.05,7.33,14.97,1.94,21.93,3.52,6.69,7.63]
exploration_nov_pcp_F = [25.95,37.09,17.11,11.08,14.53,17.44,8.93,14.44,16.25,10.63,11.04,13.44]
DI_pcp_F = [0.32,0.15,-0.42,0.03,0.02,0.41,-0.25,0.76,-0.15,0.5,0.25,0.28]
   
    
 



   
## Correlations / regressions 
# Males SALINE
# (1) DI score  - nb for males which ones have scores and which do not ? 
# (2) % distracted on distraction day  
# (3) Palatability (licks per burst) 
# (4) % distracted on distraction day 

from scipy import stats
from scipy.stats.stats import pearsonr

import matplotlib.pyplot as plt
pearsonr(np.asarray(DI_sal_M), np.asarray(percent_dis_dis_sal_M[2:]))
pearsonr(all_mean_burst_length_sal_M, percent_dis_dis_sal_M)





'''

import matplotlib as mpl 

import seaborn as sn

sn.set_style("darkgrid")
mpl.rcParams['font.size'] = 14
# MALES DI VS % DISTRACTED
slope, intercept, r_value, p_value, std_err = stats.linregress(DI_sal_M, percent_dis_dis_sal_M[2:])
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(DI_pcp_M, percent_dis_dis_pcp_M)
## Plot the scatter of Male Saline data 
plt.plot(DI_sal_M, percent_dis_dis_sal_M[2:],'o', color='darkgrey', label='saline')
## Add line of best fit for Male Saline data
plt.plot(np.asarray(DI_sal_M), intercept+slope*np.asarray(DI_sal_M), 'darkgrey', label='saline fitted')
## Plot scatters for Male PCP
plt.plot(DI_pcp_M, percent_dis_dis_pcp_M,'o', color='#FFBA08', label='pcp')
## Plot line of best fit for Male PCP
plt.plot(np.asarray(DI_pcp_M), intercept2+slope2*np.asarray(DI_pcp_M), '#FFBA08', label='pcp fitted')

plt.legend()
sn.despine(offset=10, trim=True); 
plt.xlabel('Discrimination index - NOR', fontsize=14)
plt.ylabel('Percentage distracted', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/Corr_DIvs%_M.pdf", bbox_inches='tight')
plt.show()
print('Linear regression, DI vs Percent distracted - Males')
print('SALINE')
print('R squared = ',r_value**2, ', p value = ', p_value)
print('PCP')
print('R squared = ',r_value2**2, ', p value = ', p_value2)



## Palatability MALES vs dis
slope, intercept, r_value, p_value, std_err = stats.linregress(all_mean_burst_length_sal_M, percent_dis_dis_sal_M)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(all_mean_burst_length_pcp_M, percent_dis_dis_pcp_M)
## Plot the scatter of Male Saline data 
plt.plot(all_mean_burst_length_sal_M, percent_dis_dis_sal_M,'o', color='darkgrey', label='saline')
## Add line of best fit for Male Saline data
plt.plot(np.asarray(all_mean_burst_length_sal_M), intercept+slope*np.asarray(all_mean_burst_length_sal_M), 'darkgrey', label='saline fitted')
## Plot scatters for Male PCP
plt.plot(all_mean_burst_length_pcp_M, percent_dis_dis_pcp_M,'o', color= '#FFBA08', label='pcp')
## Plot line of best fit for Male PCP
plt.plot(np.asarray(all_mean_burst_length_pcp_M), intercept2+slope2*np.asarray(all_mean_burst_length_pcp_M), '#FFBA08', label='pcp fitted')

plt.legend()
sn.despine(offset=10, trim=True); 

plt.xlabel('Mean burst length', fontsize=14)
plt.ylabel('Percentage distracted', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/Corr_Burstvs%_M.pdf", bbox_inches='tight')
plt.show()
print('Linear regression, Palatability vs Percent distracted - Males')
print('SALINE')
print('R squared = ',r_value**2, ', p value = ', p_value)
print('PCP')
print('R squared = ',r_value2**2, ', p value = ', p_value2)

# Females   
# Females DI VS % DISTRACTED
slope, intercept, r_value, p_value, std_err = stats.linregress(DI_sal_F, percent_dis_dis_sal_F)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(DI_pcp_F, percent_dis_dis_pcp_F)
## Plot the scatter of Male Saline data 
plt.plot(DI_sal_F, percent_dis_dis_sal_F,'o', color='darkgrey', label='saline')
## Add line of best fit for Male Saline data
plt.plot(np.asarray(DI_sal_F), intercept+slope*np.asarray(DI_sal_F), 'darkgrey', label='saline fitted')
## Plot scatters for Male PCP
plt.plot(DI_pcp_F, percent_dis_dis_pcp_F,'o', color='#249E8D', label='pcp')
## Plot line of best fit for Male PCP
plt.plot(np.asarray(DI_pcp_F), intercept2+slope2*np.asarray(DI_pcp_F), '#249E8D', label='pcp fitted')

plt.legend()
sn.despine(offset=10, trim=True); 
plt.xlabel('Discrimination index - NOR', fontsize=14)
plt.ylabel('Percentage distracted', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/Corr_DIvs%_F.pdf", bbox_inches='tight')
plt.show()
print('Linear regression, DI vs Percent distracted - Females')
print('SALINE')
print('R squared = ',r_value**2, ', p value = ', p_value)
print('PCP')
print('R squared = ',r_value2**2, ', p value = ', p_value2)

## Palatability FEMALES vs dis
slope, intercept, r_value, p_value, std_err = stats.linregress(all_mean_burst_length_sal_F, percent_dis_dis_sal_F)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(all_mean_burst_length_pcp_F, percent_dis_dis_pcp_F)
## Plot the scatter of Male Saline data 
plt.plot(all_mean_burst_length_sal_F, percent_dis_dis_sal_F,'o', color='darkgrey',label='saline')
## Add line of best fit for Male Saline data
plt.plot(np.asarray(all_mean_burst_length_sal_F), intercept+slope*np.asarray(all_mean_burst_length_sal_F), 'darkgrey', label='saline fitted')
## Plot scatters for Male PCP
plt.plot(all_mean_burst_length_pcp_M, percent_dis_dis_pcp_M,'o', color='#249E8D', label='pcp')
## Plot line of best fit for Male PCP
plt.plot(np.asarray(all_mean_burst_length_pcp_F), intercept2+slope2*np.asarray(all_mean_burst_length_pcp_F), '#249E8D', label='pcp fitted')

plt.legend()
#for spine in plt.gca().spines.values():
#   spine.set_visible(False)
#   
sn.despine(offset=10, trim=True);
plt.xlabel('Mean burst length', fontsize=14)
plt.ylabel('Percentage distracted', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/Corr_Burstvs%_F.pdf", bbox_inches='tight')
plt.show()   
print('Linear regression, Palatability vs Percent distracted - Females')
print('SALINE')
print('R squared = ',r_value**2, ', p value = ', p_value)
print('PCP')
print('R squared = ',r_value2**2, ', p value = ', p_value2)   
    
#FFBA08 male pcp
#FFE5A5 male saline 

#249E8D female pcp
#AFDBD5 female saline 
    


    
'''    
    
    
    
    
    
    