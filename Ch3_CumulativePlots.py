#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 16:32:14 2018

@author: u1490431
"""

# CUMULATIVE PDPS plots 

## SALINE, FEMALES -  LAST LICK DAY AND DISTRACTION DAY PDPS ALL 
fig = plt.figure()
plt.title('Lickday SAL_F', **Calibri, **Size)
ax = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Plots all for lick day with average
###### Issue that on modelled day, no distracted trials so len(pdps) was 0
## For loop to combine distracted and non distracted (as separated these in func)
## MODELLED SAL MALES 
all_pdps_mod_sal_F = []
for index, pdplists in enumerate(pdps_mod_dis_sal_F):    
    C = pdps_mod_dis_sal_F[index] + pdps_mod_notdis_sal_F[index]
    all_pdps_mod_sal_F.append(C)
# Plots all for MODELLED day with average 
for index, licklist in enumerate(all_pdps_mod_sal_F):
    plot = cumulativelickFig(ax, all_pdps_mod_sal_F[index], normed=True, color='lightgrey', log=True)
avg = [item for rat in all_pdps_mod_sal_F for item in rat] 
cumulativelickFig(ax, avg, normed=True, color='dimgrey', log=True)

## Distraction day SAL MALES  
# CUMULATIVE PDPS
fig = plt.figure()
plt.title('Distraction SAL_F', **Calibri, **Size)
ax = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
  
all_pdps_sal_F = []
for index, pdplists in enumerate(pdps_dis_sal_F):    
    C = pdps_dis_sal_F[index] + pdps_notdis_sal_F[index]
    all_pdps_sal_F.append(C)
# Plots all for MODELLED day with average 
for index, licklist in enumerate(all_pdps_sal_F):
    plot = cumulativelickFig(ax, all_pdps_sal_F[index], normed=True, color='lightgrey', log=True)
avg = [item for rat in all_pdps_sal_F for item in rat] 
cumulativelickFig(ax, avg, normed=True, color='darkturquoise', log=True)
#avg2 = [item for rat in pdps_dis_sal_M for item in rat] 
#cumulativelickFig(ax, avg2, normed=True, color='green', log=True)
#avg3 = [item for rat in pdps_notdis_sal_M for item in rat] 
#cumulativelickFig(ax, avg3, normed=True, color='blue', log=True)

fig = plt.figure()
plt.title('Lickday PCP_F', **Calibri, **Size)
ax = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
  
## MODELLED PCP MALES 
all_pdps_mod_pcp_F = []
for index, pdplists in enumerate(pdps_mod_dis_pcp_F):    
    C = pdps_mod_dis_pcp_F[index] + pdps_mod_notdis_pcp_F[index]
    all_pdps_mod_pcp_F.append(C)
# Plots all for MODELLED day with average 
for index, licklist in enumerate(all_pdps_mod_pcp_F):
    plot = cumulativelickFig(ax, all_pdps_mod_pcp_F[index], normed=True, color='lightgrey', log=True)
avg = [item for rat in all_pdps_mod_pcp_F for item in rat] 
cumulativelickFig(ax, avg, normed=True, color='black', log=True)

## Distraction day PCP MALES  
# CUMULATIVE PDPS
fig = plt.figure()
plt.title('Distraction PCP_F', **Calibri, **Size)
ax = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
  
all_pdps_pcp_F = []
for index, pdplists in enumerate(pdps_dis_pcp_F):    
    C = pdps_dis_pcp_F[index] + pdps_notdis_pcp_F[index]
    all_pdps_pcp_F.append(C)
# Plots all for MODELLED day with average 
for index, licklist in enumerate(all_pdps_pcp_F):
    plot = cumulativelickFig(ax, all_pdps_pcp_F[index], normed=True, color='lightgrey', log=True)
avg = [item for rat in all_pdps_pcp_F for item in rat] 
cumulativelickFig(ax, avg, normed=True, color='dodgerblue', log=True)



### All four lines on one graph no individual rats at all 

fig = plt.figure()
plt.title('PDPs by day and group F', **Calibri, **Size)
ax = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
  

avg = [item for rat in all_pdps_mod_sal_F for item in rat] 
avg1 = [item for rat in all_pdps_sal_F for item in rat] 
avg2 = [item for rat in all_pdps_mod_pcp_F for item in rat] 
avg3 = [item for rat in all_pdps_pcp_F for item in rat] 


cumulativelickFig(ax, avg, normed=True, color='dimgrey', log=True)
cumulativelickFig(ax, avg1, normed=True, color='darkturquoise', log=True)
cumulativelickFig(ax, avg2, normed=True, color='black', log=True) 
cumulativelickFig(ax, avg3, normed=True, color='dodgerblue', log=True)