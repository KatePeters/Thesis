#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 17:42:32 2018

@author: u1490431
"""

# Within session habituation (or lack of it)

import itertools  
# Merge the distracted and not distracted presentation times into a single list 
# for each rat
# then, order that list (smallest to largest) so that distracted and not are inter
# mixed
# then, flatten this list (8 rats into one lone list) and compare the PDPs (by flattening this list of lists too)
all_distractors_sorted = []  
for index, value in enumerate(discalc): 
    print(len(value[0]), len(value[1]))
    merged_distractors = list(itertools.chain.from_iterable(value))  
#    sorted_merge_distractors = sorted(merged_distractors)  
    sorted_merge_distractors = merged_distractors  

    all_distractors_sorted.append(sorted_merge_distractors)
    print(len(merged_distractors))
    

all_distractors_sorted_flat = list(itertools.chain.from_iterable(all_distractors_sorted))   
all_pdps_dis_flat = list(itertools.chain.from_iterable(all_pdps_dis))

'''
import seaborn as sn

sn.set_style("darkgrid")
mpl.rcParams['font.size'] = 14
'''
slope, intercept, r_value, p_value, std_err = stats.linregress(all_distractors_sorted_flat, all_pdps_dis_flat)
## Plot the scatter of Male Saline data 

for index, value in enumerate(all_pdps_dis_flat):
    if value > 1:
        plt.plot(all_distractors_sorted_flat[index], all_pdps_dis_flat[index],'o', color='blue', label='saline')
    else:
        plt.plot(all_distractors_sorted_flat[index], all_pdps_dis_flat[index],'o', color='black', label='saline')
        

## Add line of best fit for Male Saline data
plt.plot(np.asarray(all_distractors_sorted_flat), intercept+slope*np.asarray(all_distractors_sorted_flat), 'darkgrey', label='saline fitted')
#plt.legend()
#sn.despine(offset=10, trim=True); 
plt.xlabel('Time in session (sec)', fontsize=14)
plt.ylabel('PDP (sec)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim([-1,50])
#plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/LinReg_pdpsVStime.pdf', bbox_inches='tight')


plt.show()
print('')
print('')
print('R squared = ',r_value**2, ', p value = ', p_value)

########
# GOT RID OF SORTING CODE FOR THE DISTRACTORS, THE LISTS DIS AND NOT DIS 
# ARE CORRECT BECAUSE THE ALL PDPS LIST IS NOT ORDERED BY TIME 
# IT WAS MADE BY MERGING DIS AND NOT DIS IN THAT ORDER 

## could go through and limit to only PDPs less than 20 seconds (as very long???)
# but a bit arbitrary
########
#######################################################################################

### Now - convert the time in session to a number 


all_distractor_number = []
for value in all_distractors_sorted: 
    print(len(np.argsort(value)))
    
    distractor_number = np.argsort(value)
    all_distractor_number.append(distractor_number)

all_dis_num_flat = list(itertools.chain.from_iterable(all_distractor_number))


slope, intercept, r_value, p_value, std_err = stats.linregress(all_dis_num_flat, all_pdps_dis_flat)
## Plot the scatter of Male Saline data 

for index, value in enumerate(all_pdps_dis_flat):
    if value > 1:
        plt.plot(all_dis_num_flat[index], all_pdps_dis_flat[index],'o', color='blue', label='saline')
    else:
        plt.plot(all_dis_num_flat[index], all_pdps_dis_flat[index],'o', color='black', label='saline')
        

## Add line of best fit for Male Saline data
plt.plot(np.asarray(all_dis_num_flat), intercept+slope*np.asarray(all_dis_num_flat), 'darkgrey', label='saline fitted')
#plt.legend()
#sn.despine(offset=10, trim=True); 
plt.xlabel('Distractor number', fontsize=14)
plt.ylabel('PDP (sec)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim([-1,50])
#plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/LinReg_pdpsVSnumber.pdf', bbox_inches='tight')

plt.show()
print('')
print('')
print('R squared = ',r_value**2, ', p value = ', p_value)

################################################################################################
################################################################################################
################################################################################################



# Habituation with percentage distracted? Compare percent distracted in first half and second (or in quarters)
# 900 second bins  

# COuld have smaller bins 

   
def binnedpercentdis(distractiondict):

    ''' Takes a dictionary of distracted and not distracted trials'''  
    
    final1 = [] ## final 1 to 4 = distracted trials in all rats (8 lists) for each quarter
    final2 = []
    final3 = []
    final4 = []
    final5 = [] ## final 5-8 = not dis trials 
    final6 = []
    final7 = []
    final8 = []
    
    for rat in distractiondict:
        dis1st = []      
        dis2nd = []
        dis3rd = []
        dis4th = []
        
        for item in rat[0]:
                if (item > 0) and (item < 900):
                    dis1st.append(item)    
                else:      
                    if (item > 900) and (item < 1800):
                        dis2nd.append(item)
                    else:
                        if (item > 1800) and (item < 2700):
                            dis3rd.append(item)
                        else:
                            if (item > 2700) and (item < 3600):
                                dis4th.append(item)
        
        final1.append(dis1st) # makes 8 lists of the dis trials from fisrt bin
        final2.append(dis2nd)
        final3.append(dis3rd)
        final4.append(dis4th)                           
                                
    for rat in distractiondict:
        notdis1st = []
        notdis2nd = []
        notdis3rd = []
        notdis4th = []
        for item in rat[1]:
                if (item > 0) and (item < 900):
                    notdis1st.append(item)    
                else:      
                    if (item > 900) and (item < 1800):
                        notdis2nd.append(item)
                    else:
                        if (item > 1800) and (item < 2700):
                            notdis3rd.append(item)
                        else:
                            if (item > 2700) and (item < 3600):
                                notdis4th.append(item)
        final5.append(notdis1st)
        final6.append(notdis2nd)
        final7.append(notdis3rd)
        final8.append(notdis4th)
        
    # same indices for all lists here 
    # 8 times, all lists have 8 lists in them     
    percent1st = []
    percent2nd = []
    percent3rd = []
    percent4th = []
    
    
    for index, value in enumerate(final1):
        percent1 = 0
        percent2 = 0
        percent3 = 0 
        percent4 = 0 
        
#        try:
            
        percent1 = len(final1[index]) / (len(final1[index]) + len(final5[index])) * 100 
        percent2 = len(final2[index]) / (len(final2[index]) + len(final6[index])) * 100
#        
        try:
            percent3 = len(final3[index]) / (len(final3[index]) + len(final7[index])) * 100
        except ZeroDivisionError:
            percent3 = np.nan
            pass
        
        percent4 = len(final4[index]) / (len(final4[index]) + len(final8[index])) * 100

            
        percent1st.append(percent1)
        percent2nd.append(percent2)
        percent3rd.append(percent3) 
        percent4th.append(percent4)

               
    return(percent1st, percent2nd, percent3rd, percent4th)
#    return(final1,final2,final3,final4, final5, final6, final7, final8)


quart1, quart2, quart3, quart4 = binnedpercentdis(discalc)
#
#np.nanmean(quart1)
#
##exclude rats 1 and 3 (0 and 2) from the plotting as these have no distractors in 

#del quart1[0], quart1[2], quart2[0], quart2[2], quart3[0], quart3[2],quart4[0], quart4[2]

# Percentage distracted 4 quarts of session: 
quart1, quart2, quart3, quart4 = binnedpercentdis(discalc)    
binned_percent = [[quart1,quart2,quart3,quart4]]
col = ['lightgrey','powderblue','blue','darkblue']
labels = ['15 min','30 min','45 min','60 min']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(binned_percent, transpose=False, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Percent distracted (%)', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40               
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/PercentDisBinned.pdf', bbox_inches='tight')


