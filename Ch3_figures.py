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

# Do I need to re-format the data so that it is grouped by rat over the 3
# days and not ALL the rats over 3 days? How does it know which rat is which
# to join? Within-subjects format?
# pcp males

# saline females

# pcp females

# put into data frame --> format so that rat1day1, rat1day2, rat1day3
import seaborn as sb
import pandas as pd

sb.swarmplot(nlicks_sal_M, nlicks_minus1_sal_M, nlicks_minus2_sal_M)

ax = sb.swarmplot(x="day", y="total_bill", hue="smoker",\
                   data=tips, palette="Set2", dodge=True)


df = pd.DataFrame({'day1':nlicks_minus2_sal_M, 'day2':nlicks_minus1_sal_M, \
                   'day3':nlicks_sal_M})

rat = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
sb.swarmplot(data=df)

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



