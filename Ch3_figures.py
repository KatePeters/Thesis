#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 09:12:11 2018

@author: u1490431

#############################################################################

CHAPTER 3 - DISTRACTION FROM ONGOING SACCHARIN CONSUMPTION: EXPERIMENTS
            IN SALINE AND PHENCYCLIDINE TREATED RATS

#############################################################################

This file contains the code for producing FIGURES 
All functions and imports are at the beginning and no others are required

Imports:
     
Functions: 
 
"""
# Imports for generating figures / functions that require them 


# Functions 






# FIGURE 1 


# saline males
# already have these somewhere but easier to index this way! 
 
# Want to plot licking across days 1,2,3 
# Also compare means on each day (mainly the last)
   
nlicks_sal_M, nlicks_minus1_sal_M, nlicks_minus2_sal_M = nlicksgrouped(last_lick_sal_M, lick_minus1_sal_M, lick_minus2_sal_M)
nlicks_pcp_M, nlicks_minus1_pcp_M, nlicks_minus2_pcp_M = nlicksgrouped(last_lick_pcp_M, lick_minus1_pcp_M, lick_minus2_pcp_M)
nlicks_sal_F, nlicks_minus1_sal_F, nlicks_minus2_sal_F = nlicksgrouped(last_lick_sal_F, lick_minus1_sal_F, lick_minus2_sal_F)
nlicks_pcp_F, nlicks_minus1_pcp_F, nlicks_minus2_pcp_F = nlicksgrouped(last_lick_pcp_F, lick_minus1_pcp_F, lick_minus2_pcp_F)

dataM = [[nlicks_minus2_sal_M,nlicks_minus1_sal_M, nlicks_sal_M], [nlicks_minus2_pcp_M,nlicks_minus1_pcp_M, nlicks_pcp_M]]
colors1 = ['#5FADEA','#FFD700','#5FADEA','#FFD700','#5FADEA','#FFD700']
colors2 = ['#0054A5','#0054A5','#0054A5', '#FFBC00','#FFBC00','#FFBC00']
dataF = [[nlicks_minus2_sal_F,nlicks_minus1_sal_F, nlicks_sal_F], [nlicks_minus2_pcp_F,nlicks_minus1_pcp_F, nlicks_pcp_F]]

# Males licking on 3 lick days 
barscatter(dataM, transpose=True, paired=False, barfacecolor=colors1, barfacecoloroption='individual',  ylabel='Licks', grouplabel=['day 1', 'day 2', 'day 3'])
# Females licking on 3 lick days 
barscatter(dataF, transpose=True, paired=False, barfacecolor=colors1, barfacecoloroption='individual',  ylabel='Licks', grouplabel=['day 1', 'day 2', 'day 3'])
# put into data frame --> format so that rat1day1, rat1day2, rat1day3
colors = ['lightblue', 'blue', 'darkblue', 'pink', 'hotpink','red']


# Figure 2 VIOLIN PLOTS FOR LICK PARAMTERS 


# Figure 1 

# Bar scatter over 3 days with saline and PCP (males) = 6 bars
# Bar scatter over 3 days with saline and PCP (females) = 6 bars

# Add significance stars (Python or manualy) to show if there are difs
    # last time I ran analysis NOT significant, but not compared all cohorts
    

# Figure 2

# Lick analysis violin plots for saline vs pcp 
# Violin plots to show no difference 
    # Back this up with stats
    # What stats? MANOVA with just sal and pcp as the IV? 
    # Look up in Field how to do this in SPSS


# Figure 3

# Barscatter percentage distracted
# Males, pcp/sal 
# Females pcp/sal 

