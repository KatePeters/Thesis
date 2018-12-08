#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 11:12:32 2018

@author: u1490431
"""

nlicksDataM = [[nlicks_sal_M], [nlicks_pcp_M]]
nlicksDataF = [[nlicks_sal_F], [nlicks_pcp_F]]

percent_dis_dis_sal_F - 

percent_dis_dis_sal_M - 

dataMdi = [[DI_sal_M],[DI_pcp_M]]



linear regression / correltation between 


#MALES
# First 2 are excluded as corrupted videos in these rats (DPCP1 rats 1 and 2 sal)
DI_sal_M
nlicks_sal_M[2:]
#FEMALES
DI_sal_F
nlicks_sal_F




# MALES DI VS % LICKS (last lick day) 
import seaborn as sn

slope, intercept, r_value, p_value, std_err = stats.linregress(DI_sal_M, nlicks_sal_M[2:])
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(DI_pcp_M, nlicks_pcp_M)
## Plot the scatter of Male Saline data 
plt.plot(DI_sal_M, nlicks_sal_M[2:],'o', color='darkgrey', label='saline')
## Add line of best fit for Male Saline data
plt.plot(np.asarray(DI_sal_M), intercept+slope*np.asarray(DI_sal_M), 'darkgrey', label='saline fitted')
## Plot scatters for Male PCP
plt.plot(DI_pcp_M, nlicks_pcp_M,'o', color='#FFBA08', label='pcp')
## Plot line of best fit for Male PCP
plt.plot(np.asarray(DI_pcp_M), intercept2+slope2*np.asarray(DI_pcp_M), '#FFBA08', label='pcp fitted')
#plt.legend()
sn.despine(offset=10, trim=True); 
plt.xlabel('Discrimination index - NOR', fontsize=14)
plt.ylabel('Licks on last lick day', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig("/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/Corr_DIvsLicks_M.pdf", bbox_inches='tight')

plt.show()
print('Linear regression, DI vs Licks - Males')
print('SALINE')
print('R squared = ',r_value**2, ', p value = ', p_value)
print('PCP')
print('R squared = ',r_value2**2, ', p value = ', p_value2)



# FEMALES DI vs LICKS (last lick day)
slope, intercept, r_value, p_value, std_err = stats.linregress(DI_sal_F, nlicks_sal_F)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(DI_pcp_F, nlicks_pcp_F)
## Plot the scatter of Male Saline data 
plt.plot(DI_sal_F, nlicks_sal_F,'o', color='darkgrey', label='saline')
## Add line of best fit for Male Saline data
plt.plot(np.asarray(DI_sal_F), intercept+slope*np.asarray(DI_sal_F), 'darkgrey', label='saline fitted')
## Plot scatters for Male PCP
plt.plot(DI_pcp_F, nlicks_pcp_F,'o', color='#249E8D', label='pcp')
## Plot line of best fit for Male PCP
plt.plot(np.asarray(DI_pcp_F), intercept2+slope2*np.asarray(DI_pcp_F), '#249E8D', label='pcp fitted')
#plt.legend()
sn.despine(offset=10, trim=True); 
plt.xlabel('Discrimination index - NOR', fontsize=14)
plt.ylabel('Licks on last lick day', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig("/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/Corr_DIvsLicks_F.pdf", bbox_inches='tight')

plt.show()
print('Linear regression, DI vs Licks - Females')
print('SALINE')
print('R squared = ',r_value**2, ', p value = ', p_value)
print('PCP')
print('R squared = ',r_value2**2, ', p value = ', p_value2)

