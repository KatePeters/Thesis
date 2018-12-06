#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 15:47:14 2018

@author: u1490431
"""

# Data for numbers of licks on distraction days (not licking days)

licks_lickday = []
licks_dis = []
licks_hab1 = []
licks_hab2 = []
licks_sal = []
licks_amph = [] 

for index, rat in enumerate(last_lick):
    licks_lickday.append(len(rat[0]))
print(np.mean(licks_lickday))
for index, rat in enumerate(distraction):
    licks_dis.append(len(rat[0]))
print(np.mean(licks_dis))
for index, rat in enumerate(hab1):
    licks_hab1.append(len(rat[0]))
print(np.mean(licks_hab1))
for index, rat in enumerate(hab2):
    licks_hab2.append(len(rat[0]))
print(np.mean(licks_hab2))
for index, rat in enumerate(salineIP):
    licks_sal.append(len(rat[0]))
print(np.mean(licks_sal))
for index, rat in enumerate(amph):
    licks_amph.append(len(rat[0]))
print(np.mean(licks_amph))


## Barscatter of the licking days / distraction / saline / amphetamine   
total_licks = [[licks_lickday, licks_dis, licks_hab1, licks_hab2, licks_sal, licks_amph]]
col = ['gray','#03A9F4','#03A9F4','#03A9F4','#E91E63','#880E4F']
labels = ['lick day','dis','hab1','hab2', 'sal', 'amph']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(total_licks, transpose=False, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Total licks', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40               
#plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/TotalLicks_amph.pdf', bbox_inches='tight')

# assign empty variables to store outputs from lick calc (to find means/groups stats)
# lists where each item is a dictionary (25) derived from lickCalc for each rat / day
lick_analysis = lickanalysis2(last_lick)
micro_lick = lickanalysis2(last_lick)
micro_dis = lickanalysis2(distraction)
micro_hab1 = lickanalysis2(hab1)
micro_hab2 = lickanalysis2(hab2)
micro_sal = lickanalysis2(salineIP)
micro_amph = lickanalysis2(amph) 

## REPEAT THIS FOR ALL DAYS!!! --> then make barscatters for the days 
# LICK DAY
mean_n_bursts1, mean_n_runs1, mean_mean_IBI1, mean_mean_IRI1,\
all_n_bursts1, all_n_runs1, all_mean_IBI1, all_mean_IRI1, \
all_mean_burst_length1, all_mean_run_length1 = grouped_lickanalysis(micro_lick)
# DISTRACTION 
mean_n_bursts2, mean_n_runs2, mean_mean_IBI2, mean_mean_IRI2,\
all_n_bursts2, all_n_runs2, all_mean_IBI2, all_mean_IRI2, \
all_mean_burst_length2, all_mean_run_length2 = grouped_lickanalysis(micro_dis)
# HABITUATION 1
mean_n_bursts3, mean_n_runs3, mean_mean_IBI3, mean_mean_IRI3,\
all_n_bursts3, all_n_runs3, all_mean_IBI3, all_mean_IRI3, \
all_mean_burst_length3, all_mean_run_length3 = grouped_lickanalysis(micro_hab1)
# HABITUATION 2
mean_n_bursts4, mean_n_runs4, mean_mean_IBI4, mean_mean_IRI4,\
all_n_bursts4, all_n_runs4, all_mean_IBI4, all_mean_IRI4, \
all_mean_burst_length4, all_mean_run_length4 = grouped_lickanalysis(micro_hab2)
# SALINE IP 
mean_n_bursts5, mean_n_runs5, mean_mean_IBI5, mean_mean_IRI5,\
all_n_bursts5, all_n_runs5, all_mean_IBI5, all_mean_IRI5, \
all_mean_burst_length5, all_mean_run_length5 = grouped_lickanalysis(micro_sal)
# AMPHETAMINE 
mean_n_bursts6, mean_n_runs6, mean_mean_IBI6, mean_mean_IRI6,\
all_n_bursts6, all_n_runs6, all_mean_IBI6, all_mean_IRI6, \
all_mean_burst_length6, all_mean_run_length6 = grouped_lickanalysis(micro_amph)


#Licking paramters - N BURSTS
nbursts_data = [[all_n_bursts1],[all_n_bursts2],[all_n_bursts3],[all_n_bursts4],[all_n_bursts5],[all_n_bursts6]]
col = ['gray','#03A9F4','#03A9F4','#03A9F4','#E91E63','#880E4F']
labels = ['lick day','dis','hab1','hab2', 'sal', 'amph']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(nbursts_data, transpose=True, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Number of bursts', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40 
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/N_bursts_amph_comparison.pdf', bbox_inches='tight')

#Licking parameters - LICKS PER BURST
nburst_len_data = [[all_mean_burst_length1],[all_mean_burst_length2],[all_mean_burst_length3],[all_mean_burst_length4],[all_mean_burst_length5],[all_mean_burst_length6]]
col = ['gray','#03A9F4','#03A9F4','#03A9F4','#E91E63','#880E4F']
labels = ['lick day','dis','hab1','hab2', 'sal', 'amph']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(nburst_len_data, transpose=True, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Mean burst length', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40 
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/Mean_burst_len_amph_comparison.pdf', bbox_inches='tight')

#Licking paramters - N CLUSTERS
nruns_data = [[all_n_runs1],[all_n_runs2],[all_n_runs3],[all_n_runs4],[all_n_runs5],[all_n_runs6]]
col = ['gray','#03A9F4','#03A9F4','#03A9F4','#E91E63','#880E4F']
labels = ['lick day','dis','hab1','hab2', 'sal', 'amph']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(nruns_data, transpose=True, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Number of clusters', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40 
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/N_cluster_amph_comparison.pdf', bbox_inches='tight')

#Licking parameters - LICKS PER CLUSTER
nrun_len_data = [[all_mean_run_length1],[all_mean_run_length2],[all_mean_run_length3],[all_mean_run_length4],[all_mean_run_length5],[all_mean_run_length6]]
col = ['gray','#03A9F4','#03A9F4','#03A9F4','#E91E63','#880E4F']
labels = ['lick day','dis','hab1','hab2', 'sal', 'amph']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(nrun_len_data, transpose=True, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Mean cluster length', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40 
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/Mean_cluster_len_amph_comparison.pdf', bbox_inches='tight')

##### Inter-cluster (IRI) ICI - could have all the IRIs on the plot if wanted (need to extract though)
IRI_data = [[all_mean_IRI1],[all_mean_IRI2],[all_mean_IRI3],[all_mean_IRI4],[all_mean_IRI5],[all_mean_IRI6]]
col = ['gray','#03A9F4','#03A9F4','#03A9F4','#E91E63','#880E4F']
labels = ['lick day','dis','hab1','hab2', 'sal', 'amph']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(IRI_data, transpose=True, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Mean inter-cluster interval', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40 
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/ICI_amph_comparison.pdf', bbox_inches='tight')

##### Inter-cluster (IBI)
IBI_data = [[all_mean_IBI1],[all_mean_IBI2],[all_mean_IBI3],[all_mean_IBI4],[all_mean_IBI5],[all_mean_IBI6]]
col = ['gray','#03A9F4','#03A9F4','#03A9F4','#E91E63','#880E4F']
labels = ['lick day','dis','hab1','hab2', 'sal', 'amph']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(IBI_data, transpose=True, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Mean inter-burst interval', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40 
plt.savefig('/Volumes/KPMSB352/Thesis/FINAL THESIS/CORRECTIONS/Figures/IBI_amph_comparison.pdf', bbox_inches='tight')
