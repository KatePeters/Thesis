#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 15:49:43 2018

@author: u1490431
"""
## Correlation between PERCENT DISTRACTED and DISTRCTOR PEAK 
import matplotlib as mpl 
from scipy import stats
import seaborn as sn

#####################

### CORRELATION BETWEEN RELATIVE DISTRACTED AND NOT DISTRACTED PEAK 
relative_dis_peak = []

for index, value in enumerate(peak_distracted):
    relative_dis_peak.append(peak_distracted[index] - peak_notdistracted[index])

## Correlation between PERCENT DISTRACTED and NOT DISTRACTED PEAK 
slope, intercept, r_value, p_value, std_err = stats.linregress(percentdisDis[2:], MultBy100(relative_dis_peak))
## Plot the scatter
plt.plot(percentdisDis[2:], MultBy100(relative_dis_peak),'o', color='#632de9')
## Add line of best fit
plt.plot(np.asarray(percentdisDis[2:]), intercept+slope*np.asarray(percentdisDis[2:]) ,'#632de9')

plt.legend()
sn.despine(offset=10, trim=True); 

plt.xlabel('Percent distracted', fontsize=14)
plt.ylabel('Realtive peak response (% ΔF)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/RelativeResponseCorr.pdf', bbox_inches="tight")


plt.show()
print('Linear regression, Percent distracted VS relative peak')
print('R squared = ',r_value**2, ', p value = ', p_value)



#### 


## Manually input here 12 data points for the 12 rats of interest (minus first2)
## data from GFP number of cells / TH number of cells

percent_dis_11 = percentdisDis[2:-1]
GFP_positive = [69,78,112,54,129,154,152,172,110,93,64]
TH_positive = [75,151,70,130,172,178,133,152,115,127,121]


## Correlation between PERCENT DISTRACTED and GFP expression total
slope, intercept, r_value, p_value, std_err = stats.linregress(percentdisDis[2:-1], GFP_positive)
## Plot the scatter
plt.plot(percentdisDis[2:-1], GFP_positive,'o', color='green')
## Add line of best fit
plt.plot(np.asarray(percentdisDis[2:-1]), intercept+slope*np.asarray(percentdisDis[2:-1]) ,'green')

## Correlation between PERCENT DISTRACTED and TH expression total

slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(percentdisDis[2:-1], TH_positive)
## Plot the scatter
plt.plot(percentdisDis[2:-1], TH_positive,'o', color='red')
## Add line of best fit
plt.plot(np.asarray(percentdisDis[2:-1]), intercept2+slope2*np.asarray(percentdisDis[2:-1]) ,'red')




plt.legend()
sn.despine(offset=10, trim=True); 

plt.xlabel('Percent distracted', fontsize=14)
plt.ylabel('Total number of stained cells in VTA', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Corr_Percent_DistractedPeak.pdf', bbox_inches="tight")
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/GFP_TH_Corr.pdf', bbox_inches="tight")

plt.show()
print('Linear regression, Percent distracted VS GFP')
print('R squared = ',r_value**2, ', p value = ', p_value)
print('Linear regression, Percent distracted VS TH')
print('R squared = ',r_value2**2, ', p value = ', p_value2)


