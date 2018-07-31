#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 11:25:34 2018

@author: u1490431
"""


## Currently calculates the PERCENTAGE OF DISTRACTED TRIALS that were white noise
## Needs to calculate percentage distracted for JUST the whitenoise trials 


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
      #  print(len(distracted), len(notdistracted))
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
