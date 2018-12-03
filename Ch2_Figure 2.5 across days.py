#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 14:51:10 2018

@author: u1490431
"""
# Run after Ch2_analysis and Ch2_figures

lick_analysis7 = lickanalysis2(last_lick)

lick_analysis6 = lickanalysis2(lick_minus1)

lick_analysis5 = lickanalysis2(lick_minus2)

lick_analysis4 = lickanalysis2(lick_minus3)

lick_analysis3 = lickanalysis2(lick_minus4)

lick_analysis2 = lickanalysis2(lick_minus5)

lick_analysis1 = lickanalysis2(lick_minus6)

# Make burst / cluster / interval measures across days and make bar scatters for these too 

## LICK DAY ANALYSIS  
# Produce medians/means for individual rats and group means (just one group DIS1) 
# Day 1
mean_n_bursts1, mean_n_runs1, mean_mean_IBI1, mean_mean_IRI1,\
all_n_bursts1, all_n_runs1, all_mean_IBI1, all_mean_IRI1, \
all_mean_burst_length1, all_mean_run_length1 = grouped_lickanalysis(lick_analysis1)
# Day 2
mean_n_bursts2, mean_n_runs2, mean_mean_IBI2, mean_mean_IRI2,\
all_n_bursts2, all_n_runs2, all_mean_IBI2, all_mean_IRI2, \
all_mean_burst_length2, all_mean_run_length2 = grouped_lickanalysis(lick_analysis2)
# Day 3
mean_n_bursts3, mean_n_runs3, mean_mean_IBI3, mean_mean_IRI3,\
all_n_bursts3, all_n_runs3, all_mean_IBI3, all_mean_IRI3, \
all_mean_burst_length3, all_mean_run_length3 = grouped_lickanalysis(lick_analysis3)
# Day 4
mean_n_bursts4, mean_n_runs4, mean_mean_IBI4, mean_mean_IRI4,\
all_n_bursts4, all_n_runs4, all_mean_IBI4, all_mean_IRI4, \
all_mean_burst_length4, all_mean_run_length4 = grouped_lickanalysis(lick_analysis4)
# Day 5
mean_n_bursts5, mean_n_runs5, mean_mean_IBI5, mean_mean_IRI5,\
all_n_bursts5, all_n_runs5, all_mean_IBI5, all_mean_IRI5, \
all_mean_burst_length5, all_mean_run_length5 = grouped_lickanalysis(lick_analysis5)
# Day 6
mean_n_bursts6, mean_n_runs6, mean_mean_IBI6, mean_mean_IRI6,\
all_n_bursts6, all_n_runs6, all_mean_IBI6, all_mean_IRI6, \
all_mean_burst_length6, all_mean_run_length6 = grouped_lickanalysis(lick_analysis6)
# Day 7
mean_n_bursts7, mean_n_runs7, mean_mean_IBI7, mean_mean_IRI7,\
all_n_bursts7, all_n_runs7, all_mean_IBI7, all_mean_IRI7, \
all_mean_burst_length7, all_mean_run_length7 = grouped_lickanalysis(lick_analysis7)



##Licking parameters - N BURSTS

nbursts_data = [[all_n_bursts1],[all_n_bursts2],[all_n_bursts3],[all_n_bursts4],[all_n_bursts5],[all_n_bursts6],[all_n_bursts7]]
col = ['powderblue','powderblue','powderblue','powderblue','powderblue','powderblue','powderblue']
labels = ['day 1', 'day 2', 'day 3', 'day 4', 'day 5', 'day 6', 'day 7']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(nbursts_data, transpose=True, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Number of bursts', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45, unequal=False) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/NBurstsBarscatter_7days.pdf',  bbox_inches='tight')

#Licking paramters - N CLUSTERS
nruns_data = [[all_n_runs1],[all_n_runs2],[all_n_runs3],[all_n_runs4],[all_n_runs5],[all_n_runs6],[all_n_runs7]]
col = ['lightpink','lightpink','lightpink','lightpink','lightpink','lightpink','lightpink']
labels = ['day 1', 'day 2', 'day 3', 'day 4', 'day 5', 'day 6', 'day 7']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(nruns_data, transpose=True, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Number of clusters', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45, unequal=False) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/NClustersBarscatter_7days.pdf',  bbox_inches='tight')

#Licking parameters - LICKS PER BURST
nburst_len_data = [[all_mean_burst_length1],[all_mean_burst_length2],[all_mean_burst_length3],[all_mean_burst_length4],[all_mean_burst_length5],[all_mean_burst_length6],[all_mean_burst_length7]]
col = ['powderblue','powderblue','powderblue','powderblue','powderblue','powderblue','powderblue']
labels = ['day 1', 'day 2', 'day 3', 'day 4', 'day 5', 'day 6', 'day 7']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(nburst_len_data, transpose=True, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Mean licks per burst', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40               
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/MeanBurstLenBarscatter_7days.pdf',  bbox_inches='tight')

##Licking parameters - LICKS PER CLUSTER
nrun_len_data = [[all_mean_run_length1],[all_mean_run_length2],[all_mean_run_length3],[all_mean_run_length4],[all_mean_run_length5],[all_mean_run_length6],[all_mean_run_length7]]
col = ['lightpink','lightpink','lightpink','lightpink','lightpink','lightpink','lightpink']
labels = ['day 1', 'day 2', 'day 3', 'day 4', 'day 5', 'day 6', 'day 7']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(nrun_len_data, transpose=True, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Mean licks per cluster', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40               
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/MeanClusterLenBarscatter_7days.pdf',  bbox_inches='tight')

#### Inter-burst (IBI) - could have all the IBIs on the plot if wanted 
IBI_data = [[all_mean_IBI1],[all_mean_IBI2],[all_mean_IBI3],[all_mean_IBI4],[all_mean_IBI5],[all_mean_IBI6],[all_mean_IBI7]]
col = ['powderblue','powderblue','powderblue','powderblue','powderblue','powderblue','powderblue']
labels = ['day 1', 'day 2', 'day 3', 'day 4', 'day 5', 'day 6', 'day 7']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(IBI_data, transpose=True, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Mean inter burst interval (s)', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40               
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/MeanIBIBarscatter_7days.pdf',  bbox_inches='tight')

###### Inter-cluster (IRI) ICI - could have all the IRIs on the plot if wanted (need to extract though)
IRI_data = [[all_mean_IRI1],[all_mean_IRI2],[all_mean_IRI3],[all_mean_IRI4],[all_mean_IRI5],[all_mean_IRI6],[all_mean_IRI7]]
col = ['lightpink','lightpink','lightpink','lightpink','lightpink','lightpink','lightpink']
labels = ['day 1', 'day 2', 'day 3', 'day 4', 'day 5', 'day 6', 'day 7']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(IRI_data, transpose=True, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Mean inter cluster interval (s)', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40               
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/MeanICIBarscatter_7days.pdf',  bbox_inches='tight')


