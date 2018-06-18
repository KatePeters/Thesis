#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 10:12:21 2018

@author: u1490431
"""

# Function for offsetting violin plot halves
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 16:55:02 2018

@author: u1490431
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
import seaborn as sns
import pandas as pd

def offset_violinplot_halves(ax, delta, width, inner, direction):
    """
    This function offsets the halves of a violinplot to compare tails
    or to plot something else in between them. This is specifically designed
    for violinplots by Seaborn that use the option `split=True`.

    For lines, this works on the assumption that Seaborn plots everything with
     integers as the center.

    Args:
     <ax>    The axis that contains the violinplots.
     <delta> The amount of space to put between the two halves of the violinplot
     <width> The total width of the violinplot, as passed to sns.violinplot()
     <inner> The type of inner in the seaborn
     <direction> Orientation of violinplot. 'hotizontal' or 'vertical'.

    Returns:
     - NA, modifies the <ax> directly
    """
    # offset stuff
    if inner == 'sticks':
        lines = ax.get_lines()
        for line in lines:
            if direction == 'horizontal':
                data = line.get_ydata()
                print(data)
                if int(data[0] + 1)/int(data[1] + 1) < 1:
                    # type is top, move neg, direction backwards for horizontal
                    data -= delta
                else:
                    # type is bottom, move pos, direction backward for hori
                    data += delta
                line.set_ydata(data)
            elif direction == 'vertical':
                data = line.get_xdata()
                print(data)
                if int(data[0] + 1)/int(data[1] + 1) < 1:
                    # type is left, move neg
                    data -= delta
                else:
                    # type is left, move pos
                    data += delta
                line.set_xdata(data)


    for ii, item in enumerate(ax.collections):
        # axis contains PolyCollections and PathCollections
        if isinstance(item, matplotlib.collections.PolyCollection):
            # get path
            path, = item.get_paths()
            vertices = path.vertices
            half_type = _wedge_dir(vertices, direction)
            # shift x-coordinates of path
            if half_type in ['top','bottom']:
               if inner in ["sticks", None]:
                    if half_type == 'top': # -> up
                        vertices[:,1] -= delta
                    elif half_type == 'bottom': # -> down
                        vertices[:,1] += delta
            elif half_type in ['left', 'right']:
                if inner in ["sticks", None]:
                    if half_type == 'left': # -> left
                        vertices[:,0] -= delta
                    elif half_type == 'right': # -> down
                        vertices[:,0] += delta

def _wedge_dir(vertices, direction):
    """
    Args:
      <vertices>  The vertices from matplotlib.collections.PolyCollection
      <direction> Direction must be 'horizontal' or 'vertical' according to how
                   your plot is laid out.
    Returns:
      - a string in ['top', 'bottom', 'left', 'right'] that determines where the
         half of the violinplot is relative to the center.
    """
    if direction == 'horizontal':
        result = (direction, len(set(vertices[1:5,1])) == 1)
    elif direction == 'vertical':
        result = (direction, len(set(vertices[-3:-1,0])) == 1)
    outcome_key = {('horizontal', True): 'bottom',
                   ('horizontal', False): 'top',
                   ('vertical', True): 'left',
                   ('vertical', False): 'right'}
    # if the first couple x/y values after the start are the same, it
    #  is the input direction. If not, it is the opposite
    return outcome_key[result]



# (1) Bar scatter plots for licking 

'''
MALES

which groups? 
what variables?
individual data points with means 


barscatter(data)


  (2)  data - lick day not paired. Bursts, lengths 
  Bursts IBIs
  Runs 
  
  BY GROUP - saline and pcp (male and female )
    
    AND 
    
  (1)  percent distracted lick day, distraction, habituation, hab2, apmetamine IP 
    
    
    
DATA FORMAT :
    
'''


'''
data = [all_mean_IBI_sal_M, all_mean_IBI_sal_M]
distractionData = np.empty((2,), dtype=np.object)
distractionData[0] = np.array(all_mean_IBI_sal_M)
distractionData[1] = np.array(all_mean_IBI_sal_M)
'''
# PUT ALL OF THE BURST LENGTHS IN HERE, THEN SEE AND ADD IN MEAN MEAN BURST LENGTH
meanburstlength = np.empty((2,), dtype=np.object)
meanburstlength[0] = np.array(all_mean_burst_length_sal_M)
meanburstlength[1] = np.array(all_mean_burst_length_pcp_M)
#meanburstlength[2] = np.array(all_mean_burst_length_sal_F)
#meanburstlength[3] = np.array(all_mean_burst_length_pcp_F)

meanrunlength = np.empty((2,), dtype=np.object)
meanrunlength[0] = np.array(all_mean_run_length_sal_M)
meanrunlength[1] = np.array(all_mean_run_length_pcp_M)
#meanrunlength[2] = np.array(all_mean_run_length_sal_F)
#meanrunlength[3] = np.array(all_mean_run_length_pcp_F)

nbursts = np.empty((2,), dtype=np.object)
nbursts[0] = np.array(all_n_bursts_sal_M)
nbursts[1] = np.array(all_n_bursts_pcp_M)
#nbursts[2] = np.array(all_n_bursts_sal_F)
#nbursts[3] = np.array(all_n_bursts_pcp_F)

nruns = np.empty((4,), dtype=np.object)
nruns[0] = np.array(all_n_runs_sal_M)
nruns[1] = np.array(all_n_runs_pcp_M)
#nruns[2] = np.array(all_n_runs_sal_F)
#nruns[3] = np.array(all_n_runs_pcp_F)

meanIBI = np.empty((2,), dtype=np.object)
meanIBI[0] = np.array(all_mean_IBI_sal_M)
meanIBI[1] = np.array(all_mean_IBI_pcp_M)
#meanIBI[2] = np.array(all_mean_IBI_sal_F)
#meanIBI[3] = np.array(all_mean_IBI_pcp_F)

meanIRI = np.empty((2,), dtype=np.object)
meanIRI[0] = np.array(all_mean_IRI_sal_M)
meanIRI[1] = np.array(all_mean_IRI_pcp_M)
#meanIRI[2] = np.array(all_mean_IRI_sal_F)
#meanIRI[3] = np.array(all_mean_IRI_pcp_F)



#### nRuns
a = []
a.extend(nruns[0])
a.extend(nruns[1])
b = ['sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal',\
     'sal','sal','sal','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp',\
     'pcp','pcp','pcp','pcp','pcp','pcp']
# Arbitrary value (not numeric is horizontal) with n values for data points
c = ['a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a']
# c = ['1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1']

df = pd.DataFrame()
df['nRuns'] = a
df['Group'] = b
df['Drug treatment'] = c

median_nRuns_sal_M = np.median(nruns[0])
median_nRuns_pcp_M = np.median(nruns[1])

sb.set_style("white")
fig, ax = plt.subplots(1,1)

plt.yticks()
ax = sb.violinplot(hue=df['Group'], x=df['nRuns'], y=df['Drug treatment'], bw = 0.4, palette=['dodgerblue', 'hotpink'], split=True, saturation=1, scale="width", inner=None)
ax.legend().set_visible(False)
ax.set_xlabel("doesnt change!", fontsize=14)
ax.tick_params(labelsize=14)
ax = sb.swarmplot(x=df["nRuns"], y=df["Drug treatment"], hue=df["Group"], palette=['dodgerblue', 'hotpink'], size=6)
ax.bar(median_nRuns_pcp_M,0.3, color='white') # USE THIS TO ADD IN MEAN/MEDIAN  - MAKE BARS NARROW AND CLEANER
ax.bar(median_nRuns_sal_M,-0.3, color='white')
sb.despine(offset=10, trim=True)
ax.set_ylabel("")
ax.set(yticks=[]) 
# Choose the legend you want from the 4 options plotted
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[2:4], labels[2:4], fontsize=14)
# Make a gap between the distributions for easier comparison, offset determined by delta
inner=None
delta =0.05
final_width = 0.6
inner=None
offset_violinplot_halves(ax, delta, final_width, inner, 'horizontal') ## Add this function to all funcs. be careful with import names sns vs sb

#fig.savefig('/Volumes/KPMSB352/Distraction photometry paper/BehaviourFigs/nRun Violin.pdf', bbox_inches="tight") 

#### nBursts
a = []
a.extend(nbursts[0])
a.extend(nbursts[1])
b = ['sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal',\
     'sal','sal','sal','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp',\
     'pcp','pcp','pcp','pcp','pcp','pcp']
# Arbitrary value (not numeric is horizontal) with n values for data points
c = ['a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a',]

df = pd.DataFrame()
df['nBursts'] = a
df['Group'] = b
df['Drug treatment'] = c

median_nBursts_sal_M = np.median(nbursts[0])
median_nBursts_pcp_M = np.median(nbursts[1])

sb.set_style("white")
fig, ax2 = plt.subplots(1,1)

plt.yticks()
ax2 = sb.violinplot(hue=df['Group'], x=df['nBursts'], y=df['Drug treatment'], bw = 0.4, palette=['dodgerblue', 'hotpink'], split=True, saturation=1, scale="width", inner=None)
ax2.legend().set_visible(False)
ax2.set_xlabel("doesnt change!", fontsize=14)
ax2.tick_params(labelsize=14)
ax2 = sb.swarmplot(x=df["nBursts"], y=df["Drug treatment"], hue=df["Group"], palette=['dodgerblue', 'hotpink'], size=6)
ax2.bar(median_nBursts_pcp_M,0.3, color='white') # USE THIS TO ADD IN MEAN/MEDIAN  - MAKE BARS NARROW AND CLEANER
ax2.bar(median_nBursts_sal_M,-0.3, color='white') 

sb.despine(offset=10, trim=True)
ax2.set_ylabel("")
ax2.set(yticks=[]) 
# Choose the legend you want from the 4 options plotted
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[2:4], labels[2:4], fontsize=14)
# Make a gap between the distributions for easier comparison, offset determined by delta
inner=None
delta =0.05
final_width = 0.6
inner=None
offset_violinplot_halves(ax2, delta, final_width, inner, 'horizontal') ## Add this function to all funcs. be careful with import names sns vs sb

#fig.savefig('/Volumes/KPMSB352/Distraction photometry paper/BehaviourFigs/nBurst Violin.pdf', bbox_inches="tight") 


## Mean burst length 
a = []
a.extend(meanburstlength[0])
a.extend(meanburstlength[1])
b = ['sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal',\
     'sal','sal','sal','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp',\
     'pcp','pcp','pcp','pcp','pcp','pcp']
# Arbitrary value (not numeric is horizontal) with n values for data points
c = ['a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a']
# c = ['1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1']

df = pd.DataFrame()
df['meanBurstlen'] = a
df['Group'] = b
df['Drug treatment'] = c

median_meanBurstlen_sal_M = np.median(meanburstlength[0])
median_meanBurstlen_pcp_M = np.median(meanburstlength[1])

sb.set_style("white")
fig, ax3 = plt.subplots(1,1)

plt.yticks()
ax3 = sb.violinplot(hue=df['Group'], x=df['meanBurstlen'], y=df['Drug treatment'], bw = 0.4, palette=['dodgerblue', 'hotpink'], split=True, saturation=1, scale="width", inner=None)
ax3.legend().set_visible(False)
ax3.set_xlabel("doesnt change!", fontsize=14)
ax3.tick_params(labelsize=14)
ax3 = sb.swarmplot(x=df["meanBurstlen"], y=df["Drug treatment"], hue=df["Group"], palette=['dodgerblue', 'hotpink'], size=6)
ax3.bar(median_meanBurstlen_pcp_M,0.3, color='white') # USE THIS TO ADD IN MEAN/MEDIAN  - MAKE BARS NARROW AND CLEANER
ax3.bar(median_meanBurstlen_sal_M,-0.3, color='white')
sb.despine(offset=10, trim=True)
ax3.set_ylabel("")
ax3.set(yticks=[]) 
# Choose the legend you want from the 4 options plotted
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[2:4], labels[2:4], fontsize=14)
# Make a gap between the distributions for easier comparison, offset determined by delta
inner=None
delta =0.05
final_width = 0.6
inner=None
offset_violinplot_halves(ax3, delta, final_width, inner, 'horizontal') ## Add this function to all funcs. be careful with import names sns vs sb

#fig.savefig('/Volumes/KPMSB352/Distraction photometry paper/BehaviourFigs/MeanBurstLength Violin.pdf', bbox_inches="tight") 


# Mean run length
a = []
a.extend(meanrunlength[0])
a.extend(meanrunlength[1])
b = ['sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal',\
     'sal','sal','sal','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp',\
     'pcp','pcp','pcp','pcp','pcp','pcp']
# Arbitrary value (not numeric is horizontal) with n values for data points
c = ['a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a']
# c = ['1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1']

df = pd.DataFrame()
df['meanRunlen'] = a
df['Group'] = b
df['Drug treatment'] = c

median_meanRunlen_sal_M = np.median(meanrunlength[0])
median_meanRunlen_pcp_M = np.median(meanrunlength[1])

sb.set_style("white")
fig, ax4 = plt.subplots(1,1)

plt.yticks()
ax4 = sb.violinplot(hue=df['Group'], x=df['meanRunlen'], y=df['Drug treatment'], bw = 0.4, palette=['dodgerblue', 'hotpink'], split=True, saturation=1, scale="width", inner=None)
ax4.legend().set_visible(False)
ax4.set_xlabel("doesnt change!", fontsize=14)
ax4.tick_params(labelsize=14)
ax4 = sb.swarmplot(x=df["meanRunlen"], y=df["Drug treatment"], hue=df["Group"], palette=['dodgerblue', 'hotpink'], size=6)
ax4.bar(median_meanRunlen_pcp_M,0.3, color='white') # USE THIS TO ADD IN MEAN/MEDIAN  - MAKE BARS NARROW AND CLEANER
ax4.bar(median_meanRunlen_sal_M,-0.3, color='white')
sb.despine(offset=10, trim=True)
ax4.set_ylabel("")
ax4.set(yticks=[]) 
# Choose the legend you want from the 4 options plotted
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[2:4], labels[2:4], fontsize=14)
# Make a gap between the distributions for easier comparison, offset determined by delta
inner=None
delta =0.05
final_width = 0.6
inner=None
offset_violinplot_halves(ax4, delta, final_width, inner, 'horizontal') ## Add this function to all funcs. be careful with import names sns vs sb
#fig.savefig('/Volumes/KPMSB352/Distraction photometry paper/BehaviourFigs/MeanRunLength Violin.pdf', bbox_inches="tight") 


# Mean IBI 
a = []
a.extend(meanIBI[0])
a.extend(meanIBI[1])
b = ['sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal',\
     'sal','sal','sal','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp',\
     'pcp','pcp','pcp','pcp','pcp','pcp']
# Arbitrary value (not numeric is horizontal) with n values for data points
c = ['a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a']
# c = ['1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1']

df = pd.DataFrame()
df['meanIBI'] = a
df['Group'] = b
df['Drug treatment'] = c

median_IBI_sal_M = np.median(meanIBI[0])
median_IBI_pcp_M = np.median(meanIBI[1])

sb.set_style("white")
fig, ax5= plt.subplots(1,1)

plt.yticks()
ax5 = sb.violinplot(hue=df['Group'], x=df['meanIBI'], y=df['Drug treatment'], bw = 0.4, palette=['dodgerblue', 'hotpink'], split=True, saturation=1, scale="width", inner=None)
ax5.legend().set_visible(False)
ax5.set_xlabel("doesnt change!", fontsize=14)
ax5.tick_params(labelsize=14)

ax5 = sb.swarmplot(x=df["meanIBI"], y=df["Drug treatment"], hue=df["Group"], palette=['dodgerblue', 'hotpink'], size=6)
ax5.bar(median_IBI_pcp_M,0.3, color='white') # USE THIS TO ADD IN MEAN/MEDIAN  - MAKE BARS NARROW AND CLEANER
ax5.bar(median_IBI_sal_M,-0.3, color='white')
sb.despine(offset=10, trim=True)
ax5.set_ylabel("")
ax5.set(yticks=[]) 
# Choose the legend you want from the 4 options plotted
handles, labels = ax5.get_legend_handles_labels()
l = plt.legend(handles[2:4], labels[2:4], fontsize=14)
# Make a gap between the distributions for easier comparison, offset determined by delta
inner=None
delta =0.05
final_width = 0.6
inner=None
offset_violinplot_halves(ax5, delta, final_width, inner, 'horizontal') ## Add this function to all funcs. be careful with import names sns vs sb
#fig.savefig('/Volumes/KPMSB352/Distraction photometry paper/BehaviourFigs/MeanIBI Violin.pdf', bbox_inches="tight") 


# Mean IRI 
a = []
a.extend(meanIRI[0])
a.extend(meanIRI[1])
b = ['sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal','sal',\
     'sal','sal','sal','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp','pcp',\
     'pcp','pcp','pcp','pcp','pcp','pcp']
# Arbitrary value (not numeric is horizontal) with n values for data points
c = ['a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a']
# c = ['1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1']

df = pd.DataFrame()
df['meanIRI'] = a
df['Group'] = b
df['Drug treatment'] = c

median_IRI_sal_M = np.median(meanIRI[0])
median_IRI_pcp_M = np.median(meanIRI[1])

sb.set_style("white")
fig, ax6= plt.subplots(1,1)

plt.yticks()
ax6 = sb.violinplot(hue=df['Group'], x=df['meanIRI'], y=df['Drug treatment'], bw = 0.4, palette=['dodgerblue', 'hotpink'], split=True, saturation=1, scale="width", inner=None)
ax6.legend().set_visible(False)
ax6.set_xlabel("doesnt change!", fontsize=14)
ax6.tick_params(labelsize=14)
ax6 = sb.swarmplot(x=df["meanIRI"], y=df["Drug treatment"], hue=df["Group"], palette=['dodgerblue', 'hotpink'], size=6)
ax6.bar(median_IRI_pcp_M,0.3, color='white') # USE THIS TO ADD IN MEAN/MEDIAN  - MAKE BARS NARROW AND CLEANER
ax6.bar(median_IRI_sal_M,-0.3, color='white')
sb.despine(offset=10, trim=True)
ax6.set_ylabel("")
ax6.set(yticks=[]) 
# Choose the legend you want from the 4 options plotted
handles, labels = ax6.get_legend_handles_labels()
l = plt.legend(handles[2:4], labels[2:4], fontsize=14)
# Make a gap between the distributions for easier comparison, offset determined by delta
inner=None
delta =0.05
final_width = 0.6
inner=None
offset_violinplot_halves(ax6, delta, final_width, inner, 'horizontal') ## Add this function to all funcs. be careful with import names sns vs sb

#fig.savefig('/Volumes/KPMSB352/Distraction photometry paper/BehaviourFigs/MeanIRI Violin.pdf', bbox_inches="tight") 


### CODE OVERWRITES DATAFRAME EACH TIME, CANNOT ACCESS PREVIOUS VARIABLES. RRE NEEDED


##################################################################################

# Barscatter plots for grouped data. PCP vs SAL



# PLOT PERCENTAGE DISTRACTED, by modality - run ANOVA (saline vs pcp)
# PLOT percentage distracted saline vs pcp (males) (female) just distraction day 

# Large plot of all days for just saline MALES
# Large plot of all days for just saline MALES 


# Individual differences, linear regression for PREDP and POSTDP 
    # Try this plot with 16 averages one from each rat. Average pause and average pre??
    # Really want individual data here not averages
    # Very variable though , issues with wild outliers messing up the regression

# Maybe limit to the PDPs under 2 seconds (then same window of 1 second either way and reduced variablity)
'''

compare distracted and not distracted pdps (all together) with all pre-dps too 

plot these 

flatten the arrays 

can go through and find if distracted or not based on the pdp

## Read in ALL of the licks?? Or bursts 



### All bar scatters for percent distracted - or even violin plots 

### Correlation using seaborn model - regression 

## VERY few are less than 2 seconds, when they are distracted they are REALLY distracted 
## Look at the pdps compare the mean or even median for distracted and not 
count = 0
for trial in pdps_dis_sal_M:
    for tril in trial[0]:
        if tril < 10:
            count += 1
            
            # cuts out the extremes. 114 are under 10 seconds. 66 are over 10 seconds 

'''