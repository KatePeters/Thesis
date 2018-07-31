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

Imports: numpy=np, matplotlib.pyplot=plt, chain
     
Functions: barscatter(), setcolors(), data2obj1D(), data2obj2D()
 
"""
# Imports for generating figures / functions that require them 
import numpy as np
import matplotlib.pyplot as plt
from itertools import chain

# Functions 

"""
This function will create bar+scatter plots when passed a 1 or 2 dimensional
array. Data needs to be passed in as a numpy object array, e.g.
data = np.empty((2), dtype=np.object)
data[0] = np.array(allData['nCasLicks'][index])
data[1] = np.array(allData['nMaltLicks'][index])
Various options allow specification of colors and paired/unpaired plotting.
It can return the figures, axes, bars, and scatters for further modification.
e.g.
fig1, ax1, barlist1, sc1 = jmf.barscatter(data)
for i in barlist1[1].get_children():
    i.set_color('g')
"""
def barscatter(data, transpose = False, unequal=False,
                groupwidth = .75,
                barwidth = .9,
                paired = False,
                spaced = False,
                barfacecoloroption = 'same', # other options 'between' or 'individual'
                barfacecolor = ['white'],
                baredgecoloroption = 'same',
                baredgecolor = [''],
                baralpha = 1,
                scatterfacecoloroption = 'same',
                scatterfacecolor = ['white'],
                scatteredgecoloroption = 'same',
                scatteredgecolor = ['grey'],
                scatterlinecolor = 'grey', # Don't put this value in a list
                scattersize = 80,
                scatteralpha = 1,
                linewidth=1,
                ylabel = 'none',
                xlabel = 'none',
                grouplabel = 'auto',
                itemlabel = 'none',
                barlabels = [],
                barlabeloffset=0.1,
                grouplabeloffset=0.2,
                yaxisparams = 'auto',
                show_legend = 'none',
                xrotation=0,
                legendloc='upper right',
                ax=[]):

    if unequal == True:
        dims = np.ndim(data)
        data_obj = np.ndarray((np.shape(data)), dtype=np.object)
        for i1, dim1 in enumerate(data):
            for i2, dim2 in enumerate(dim1):
                data_obj[i1][i2] = np.array(dim2, dtype=np.object)
        data = data_obj
    
    if type(data) != np.ndarray or data.dtype != np.object:
        dims = np.shape(data)
        if len(dims) == 2 or len(dims) == 1:
            data = data2obj1D(data)

        elif len(dims) == 3:
            data = data2obj2D(data)
              
        else:
            print('Cannot interpret data shape. Should be 2 or 3 dimensional array. Exiting function.')
            return

    # Check if transpose = True
    if transpose == True:
        data = np.transpose(data)
        
    # Initialize arrays and calculate number of groups, bars, items, and means
    
    barMeans = np.zeros((np.shape(data)))
    items = np.zeros((np.shape(data)))
    
    nGroups = np.shape(data)[0]
    groupx = np.arange(1,nGroups+1)

    if len(np.shape(data)) > 1:
        grouped = True
        barspergroup = np.shape(data)[1]
        barwidth = (barwidth * groupwidth) / barspergroup
        
        for i in range(np.shape(data)[0]):
            for j in range(np.shape(data)[1]):
                barMeans[i][j] = np.mean(data[i][j])
                items[i][j] = len(data[i][j])
        
    else:
        grouped = False
        barspergroup = 1
        
        for i in range(np.shape(data)[0]):
            barMeans[i] = np.mean(data[i])
            items[i] = len(data[i])
    
    # Calculate x values for bars and scatters
    
    xvals = np.zeros((np.shape(data)))
    barallocation = groupwidth / barspergroup
    k = (groupwidth/2) - (barallocation/2)
    
    if grouped == True:
        
        for i in range(np.shape(data)[0]):
            xrange = np.linspace(i+1-k, i+1+k, barspergroup)
            for j in range(barspergroup):
                xvals[i][j] = xrange[j]
    else:
        xvals = groupx
    
    # Set colors for bars and scatters
     
    barfacecolorArray = setcolors(barfacecoloroption, barfacecolor, barspergroup, nGroups, data)
    baredgecolorArray = setcolors(baredgecoloroption, baredgecolor, barspergroup, nGroups, data)
     
    scfacecolorArray = setcolors(scatterfacecoloroption, scatterfacecolor, barspergroup, nGroups, data, paired_scatter = paired)
    scedgecolorArray = setcolors(scatteredgecoloroption, scatteredgecolor, barspergroup, nGroups, data, paired_scatter = paired)
    
    # Initialize figure
    if ax == []:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    # Make bars
    barlist = []
    barx = []
    for x, y, bfc, bec in zip(xvals.flatten(), barMeans.flatten(),
                              barfacecolorArray, baredgecolorArray):
        barx.append(x)
        barlist.append(ax.bar(x, y, barwidth,
                         facecolor = bfc, edgecolor = bec,
                         zorder=-1))
    
    # Uncomment these lines to show method for changing bar colors outside of
    # function using barlist properties
    #for i in barlist[2].get_children():
    #    i.set_color('r')
    
    # Make scatters
    sclist = []
    if paired == False:
        for x, Yarray, scf, sce  in zip(xvals.flatten(), data.flatten(),
                                        scfacecolorArray, scedgecolorArray):
            for y in Yarray:
                if spaced == True:
                    sclist.append(ax.scatter(x+np.random.random(size=1)*barallocation, y, s = scattersize,
                             c = scf,
                             edgecolors = sce,
                             zorder=20))
                else:
                     sclist.append(ax.scatter(x, y, s = scattersize,
                                     c = scf,
                                     edgecolors = sce,
                                     zorder=20))

    elif grouped == True:
        for x, Yarray, scf, sce in zip(xvals, data, scfacecolorArray, scedgecolorArray):
            for y in np.transpose(Yarray.tolist()):
                sclist.append(ax.plot(x, y, '-o', markersize = scattersize/10,
                         color = scatterlinecolor,
                         linewidth=linewidth,
                         markerfacecolor = scf,
                         markeredgecolor = sce,
                         zorder=20))
    elif grouped == False:
        for n,_ in enumerate(data[0]):
            y = [y[n-1] for y in data]
            sclist.append(ax.plot(xvals, y, '-o', markersize = scattersize/10,
                         color = 'grey',
                         linewidth=linewidth,
                         markerfacecolor = scfacecolorArray[0],
                         markeredgecolor = scedgecolorArray[0],
                         zorder=20))
    
    # Label axes
    if ylabel != 'none':
        plt.ylabel(ylabel)
    
    if xlabel != 'none':
        plt.xlabel(xlabel)
    
    # Set range and tick values for Y axis
    if yaxisparams != 'auto':
        ax.set_ylim(yaxisparams[0])
        plt.yticks(yaxisparams[1])
       
    # X ticks
    ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off') # labels along the bottom edge are off

    if grouplabel == 'auto':
        plt.tick_params(labelbottom='off')
    else:
        if len(barlabels) > 0:
            plt.tick_params(labelbottom='off')
            yrange = ax.get_ylim()[1] - ax.get_ylim()[0]
            offset = ax.get_ylim()[0] - yrange*grouplabeloffset
            for idx, label in enumerate(grouplabel):
                ax.text(idx+1, offset, label, va='top', ha='center')
        else:
            plt.xticks(range(1,nGroups+1), grouplabel)
        
    if len(barlabels) > 0:
        if len(barlabels) != len(barx):
            print('Wrong number of bar labels for number of bars!')
        else:
            yrange = ax.get_ylim()[1] - ax.get_ylim()[0]
            offset = ax.get_ylim()[0] - yrange*barlabeloffset
            for x, label in zip(barx, barlabels):
                ax.text(x, offset, label, va='top', ha='center',rotation=xrotation)
    
    # Hide the right and top spines and set bottom to zero
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_position('zero')
    
    if show_legend == 'within':
        if len(itemlabel) != barspergroup:
            print('Not enough item labels for legend!')
        else:
            legendbar = []
            legendtext = []
            for i in range(barspergroup):
                legendbar.append(barlist[i])
                legendtext.append(itemlabel[i])
            plt.legend(legendbar, legendtext, loc=legendloc)
    
    return ax, barx, barlist, sclist

#plt.savefig('foo.png')
        
# To do
# check if n's are the same for paired and if not default to unpaired
# add color options for scatters
# add alpha options etc
# add axis options
# remove white background
# work out how to export or save as pdf, tiff, eps etc
# work out how to return handles to scatters so can be altered outside of function
# make help doc
# make html file to show usage using ijupyt
      
def setcolors(coloroption, colors, barspergroup, nGroups, data, paired_scatter = False):
            
    nColors = len(colors)
    
    if (paired_scatter == True) & (coloroption == 'within'):
        print('Not possible to make a Paired scatter plot with Within setting.')
        coloroption = 'same'
        
    if coloroption == 'within':
        if nColors < barspergroup:
            print('Not enough colors for this option! Reverting to one color.')
            coloroption = 'same'
        elif nColors > barspergroup:
            colors = colors[:barspergroup]
        coloroutput = [colors for i in data]
        coloroutput = list(chain(*coloroutput))
        
    if coloroption == 'between':
        if nColors < nGroups:
            print('Not enough colors for this option! Reverting to one color.')
            coloroption = 'same'
        elif nColors > nGroups:
            colors = colors[:nGroups]
        if paired_scatter == False:
            coloroutput = [[c]*barspergroup for c in colors]
            coloroutput = list(chain(*coloroutput))
        else:
            coloroutput = colors
            
    if coloroption == 'individual':
        if nColors < nGroups*barspergroup:
            print('Not enough colors for this color option')
            coloroption = 'same'
        elif nColors > nGroups*barspergroup:
            coloroutput = colors[:nGroups*barspergroup]
        else: 
            coloroutput = colors
    
    if coloroption == 'same':
        coloroutput = [colors[0] for x in range(len(data.flatten()))]

    return coloroutput

def data2obj1D(data):
    obj = np.empty(len(data), dtype=np.object)
    for i,x in enumerate(data):
        obj[i] = np.array(x)  
    return obj

def data2obj2D(data):
    obj = np.empty((np.shape(data)[0], np.shape(data)[1]), dtype=np.object)
    for i,x in enumerate(data):
        for j,y in enumerate(x):
            obj[i][j] = np.array(y)
    return obj



import matplotlib as mpl  # import to modify graphs 
# FIGURE 1a and 1b 
# Saline and pcp treated rats across the three lick days before distraction 
# a - Males (n=16 per group)
# b - Females (n=16 per group)

# saline males
# already have these somewhere but easier to index this way! 
 
# Want to plot licking across days 1,2,3 
# Also compare means on each day (mainly the last)
   
nlicks_sal_M, nlicks_minus1_sal_M, nlicks_minus2_sal_M = nlicksgrouped(last_lick_sal_M, lick_minus1_sal_M, lick_minus2_sal_M)
nlicks_pcp_M, nlicks_minus1_pcp_M, nlicks_minus2_pcp_M = nlicksgrouped(last_lick_pcp_M, lick_minus1_pcp_M, lick_minus2_pcp_M)
nlicks_sal_F, nlicks_minus1_sal_F, nlicks_minus2_sal_F = nlicksgrouped(last_lick_sal_F, lick_minus1_sal_F, lick_minus2_sal_F)
nlicks_pcp_F, nlicks_minus1_pcp_F, nlicks_minus2_pcp_F = nlicksgrouped(last_lick_pcp_F, lick_minus1_pcp_F, lick_minus2_pcp_F)

dataM = [[nlicks_minus2_sal_M,nlicks_minus1_sal_M, nlicks_sal_M], [nlicks_minus2_pcp_M,nlicks_minus1_pcp_M, nlicks_pcp_M]]
col3 = ['#FFE5A5','#FFE5A5','#FFE5A5','#FFBA08','#FFBA08','#FFBA08']
dataF = [[nlicks_minus2_sal_F,nlicks_minus1_sal_F, nlicks_sal_F], [nlicks_minus2_pcp_F,nlicks_minus1_pcp_F, nlicks_pcp_F]]
labels = ['-3','-2','-1','-3','-2','-1']
# Males licking on 3 lick days 
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(dataM, transpose=False, ax=ax,paired=True, barfacecolor=col3, barfacecoloroption='individual',  ylabel='Licks', xlabel='Lick days before distraction', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40
plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/LickDays_M.pdf", bbox_inches='tight')



col3 = ['#AFDBD5','#AFDBD5','#AFDBD5','#249E8D','#249E8D','#249E8D']
# Females licking on 3 lick days 
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) 
ax, barx, barlist, sclist = barscatter(dataF, transpose=False, ax=ax,paired=True, barfacecolor=col3, barfacecoloroption='individual',  ylabel='Licks', xlabel='Lick days before distraction', barlabels=labels, barlabeloffset=0.05)#, grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40
plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/LickDays_F.pdf", bbox_inches='tight')

# This section has been replaced with panel figures of all paramters
## Barscatters to replace violin plots here 
#nlicksData =  [[nlicks_sal_M], [nlicks_pcp_M]]
#meanburstlenData = [[all_mean_burst_length_sal_M], [all_mean_burst_length_pcp_M]]
#nburstsData = [[all_n_bursts_sal_M], [all_n_bursts_pcp_M]]
#ax, barx, barlist, sclist = barscatter(nlicksData, transpose=False, paired=False, barfacecolor=['#FFE5A5','#FFBA08'], barfacecoloroption='individual',  ylabel='nLicks', barlabels=labels, itemlabel=['1','2']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
#plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/nlicks2barMale.pdf", bbox_inches='tight')
#ax, barx, barlist, sclist = barscatter(meanburstlenData, transpose=False, paired=False, barfacecolor=['#FFE5A5','#FFBA08'], barfacecoloroption='individual',  ylabel='Licks per burst', barlabels=labels, itemlabel=['1','2']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
#plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/meanburstlen2barMale.pdf", bbox_inches='tight')
#ax, barx, barlist, sclist = barscatter(nburstsData, transpose=False, paired=False, barfacecolor=['#FFE5A5','#FFBA08'], barfacecoloroption='individual',  ylabel='nBursts', barlabels=labels, itemlabel=['1','2']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
#plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/nburst2barMale.pdf", bbox_inches='tight')
#
#
#


## Licking paramters
nlicksDataM = [[nlicks_sal_M], [nlicks_pcp_M]]
meanburstlenDataM = [[all_mean_burst_length_sal_M], [all_mean_burst_length_pcp_M]]
nburstsDataM = [[all_n_bursts_sal_M], [all_n_bursts_pcp_M]]

nlicksDataF = [[nlicks_sal_F], [nlicks_pcp_F]]
meanburstlenDataF = [[all_mean_burst_length_sal_F], [all_mean_burst_length_pcp_F]]
nburstsDataF = [[all_n_bursts_sal_F], [all_n_bursts_pcp_F]]
 
## Initialise figures 

#Ω# Supervisory meeting 
 ## make sure mpl.  
#Males
mpl.rcParams['figure.subplot.wspace'] = 0.6
mpl.rcParams['figure.subplot.right'] = 1
mpl.rcParams['font.size'] = 14
figureA, ax = plt.subplots(nrows=1, ncols=3, figsize=(8,4)) ### x,y 

figureA.tight_layout(pad=3, w_pad=3, h_pad=1.0)

labels = []
ax[0], barx, barlist, sclist = barscatter(nlicksDataM, ax=ax[0],transpose=False, paired=False, barfacecolor=['#FFE5A5','#FFBA08'], barfacecoloroption='individual',  ylabel='nLicks', barlabels=labels, baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[1], barx, barlist, sclist = barscatter(meanburstlenDataM,ax=ax[1], transpose=False, paired=False, barfacecolor=['#FFE5A5','#FFBA08'], barfacecoloroption='individual',  ylabel='Licks per burst', barlabels=labels, baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[2], barx, barlist, sclist = barscatter(nburstsDataM, ax=ax[2],transpose=False, paired=False, barfacecolor=['#FFE5A5','#FFBA08'], barfacecoloroption='individual',  ylabel='nBursts', barlabels=labels, show_legend='within', itemlabel=['SAL', 'PCP'], baredgecolor=[''] )#,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[0].set_ylabel('Licks')
ax[1].set_ylabel('Mean burst length')
ax[2].set_ylabel('Mean number of bursts')

ax[0].set_xticks([])
ax[0].set_ylim([0,4000])
ax[1].set_xticks([])
ax[1].set_ylim([0,25])
ax[2].set_xticks([])
ax[2].set_ylim(0,1200)

ax[0].spines['bottom'].set_visible(False)



plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/LickParams_M.pdf", bbox_inches='tight')



#Females    
mpl.rcParams['figure.subplot.wspace'] = 0.6
figureB, ax = plt.subplots(nrows=1, ncols=3, figsize=(8,4)) ### x,y 
figureB.tight_layout(pad=3, w_pad=3, h_pad=1.0)
ax[0], barx, barlist, sclist = barscatter(nlicksDataF, ax=ax[0],transpose=False, paired=False, barfacecolor=["#AFDBD5", "#249E8D"], barfacecoloroption='individual',  ylabel='nLicks', barlabels=labels, itemlabel=['1','2'],baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[1], barx, barlist, sclist = barscatter(meanburstlenDataF,ax=ax[1], transpose=False, paired=False, barfacecolor=["#AFDBD5", "#249E8D"], barfacecoloroption='individual',  ylabel='Licks per burst', barlabels=labels, itemlabel=['1','2'],baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[2], barx, barlist, sclist = barscatter(nburstsDataF, ax=ax[2],transpose=False, paired=False, barfacecolor=["#AFDBD5", "#249E8D"], barfacecoloroption='individual',  ylabel='nBursts', barlabels=labels, itemlabel=['1','2'],baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])

ax[0].set_ylabel('Licks')
ax[1].set_ylabel('Mean burst length')
ax[2].set_ylabel('Mean number of bursts')

#Change these limits and add in the barcolour removal parameter
ax[0].set_xticks([])
ax[0].set_ylim([0,6000])
ax[1].set_xticks([])
ax[1].set_ylim([0,25])
ax[2].set_xticks([])
ax[2].set_ylim(0,500)

ax[0].spines['bottom'].set_visible(False)


plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/LickParams_F.pdf", bbox_inches='tight')


#### ax[1].set_ylabel # etc



## Colours specified by HEX code

#Selective yellow --> #FFBA08
#Navajo white --> #FFE5A5
#Jungle green --> #249E8D
#Light blue --> #AFDBD5
#Lapis lazuli --> #346699



# Males, pcp/sal 
# Females pcp/sal 
# SALINE

dataMdis = [[percent_dis_modelled_sal_M,percent_dis_dis_sal_M,\
         percent_dis_hab1_sal_M,percent_dis_hab2_sal_M,\
         percent_dis_amph_sal_M], [percent_dis_modelled_pcp_M,\
         percent_dis_dis_pcp_M,percent_dis_hab1_pcp_M,\
         percent_dis_hab2_pcp_M,percent_dis_amph_pcp_M]]

tenbarcolors = ['#FFE5A5','#FFE5A5','#FFE5A5','#FFE5A5','#FFE5A5','#FFBA08','#FFBA08','#FFBA08','#FFBA08','#FFBA08']
labels = ['mod','dis','hab1','hab2','amph','mod','dis','hab1','hab2','amph']
ax, barx, barlist, sclist = barscatter(dataMdis, transpose=False, paired=True, barfacecolor=tenbarcolors, barfacecoloroption='individual',scatterlinecolor = 'lightgrey',  ylabel='Mean percent distracted',  barlabels=labels, itemlabel=labels,xrotation=45, barlabeloffset=0.05) 

ax.spines['bottom'].set_visible(False)

plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/PercentDisM.pdf", bbox_inches='tight')
           
#
tenbarcolors = ['#AFDBD5','#AFDBD5','#AFDBD5','#AFDBD5','#AFDBD5','#249E8D','#249E8D','#249E8D','#249E8D','#249E8D']
dataFdis = [[percent_dis_modelled_sal_F,percent_dis_dis_sal_F,\
         percent_dis_hab1_sal_F,percent_dis_hab2_sal_F,\
         percent_dis_amph_sal_F], [percent_dis_modelled_pcp_F,\
         percent_dis_dis_pcp_F,percent_dis_hab1_pcp_F,\
         percent_dis_hab2_pcp_F,percent_dis_amph_pcp_F]]

labels = ['mod','dis','hab1','hab2','amph','mod','dis','hab1','hab2','amph']
ax, barx, barlist, sclist = barscatter(dataFdis, transpose=False,paired=True, barfacecolor=tenbarcolors, barfacecoloroption='individual', scatterlinecolor = 'lightgrey', ylabel='Mean percent distracted', barlabels=labels, xrotation=45, barlabeloffset=0.05) 
ax.spines['bottom'].set_visible(False)

plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/PercentDisF.pdf", bbox_inches='tight')

## Colours specified by HEX code

#Selective yellow --> #FFBA08
#Navajo white --> #FFE5A5
#Jungle green --> #249E8D
#Light blue --> #AFDBD5
#Lapis lazuli --> #346699


## TESTING PDP PLOTS??????
tenbarcolors = ['#AFDBD5','#AFDBD5','#AFDBD5','#AFDBD5','#AFDBD5','#249E8D','#249E8D','#249E8D','#249E8D','#249E8D']
dataPDPsM = [[med_pdps_mod_dis_sal_M, med_pdps_dis_sal_M, med_pdps_hab1_dis_sal_M,\
              med_pdps_hab2_dis_sal_M, med_pdps_amph_dis_sal_M], [med_pdps_mod_dis_pcp_M, med_pdps_dis_pcp_M, med_pdps_hab1_dis_pcp_M,\
              med_pdps_hab2_dis_pcp_M, med_pdps_amph_dis_pcp_M]]


labels = ['mod','dis','hab1','hab2','amph','mod','dis','hab1','hab2','amph']
ax, barx, barlist, sclist = barscatter(dataPDPsM, transpose=False,paired=True, barfacecolor=tenbarcolors, barfacecoloroption='individual',  ylabel='Mean percent distracted', barlabels=labels)#, xrotation=45) 
#plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/figure2b.pdf", bbox_inches='tight')

## TESTING PDP PLOTS?????? female 
tenbarcolors = ['#AFDBD5','#AFDBD5','#AFDBD5','#AFDBD5','#AFDBD5','#249E8D','#249E8D','#249E8D','#249E8D','#249E8D']
dataPDPsF = [[med_pdps_mod_dis_sal_F, med_pdps_dis_sal_F, med_pdps_hab1_dis_sal_F,\
              med_pdps_hab2_dis_sal_F, med_pdps_amph_dis_sal_F], [med_pdps_mod_dis_pcp_F, med_pdps_dis_pcp_F, med_pdps_hab1_dis_pcp_F,\
              med_pdps_hab2_dis_pcp_F, med_pdps_amph_dis_pcp_F]]

labels = ['mod','dis','hab1','hab2','amph','mod','dis','hab1','hab2','amph']
ax, barx, barlist, sclist = barscatter(dataPDPsF, transpose=False,paired=True, barfacecolor=tenbarcolors, barfacecoloroption='individual',  ylabel='Mean percent distracted', barlabels=labels)#,xrotation=45) 
#plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/figure2b.pdf", bbox_inches='tight')



## NOR plots - Normal bar plots with mean and SEM 
## For DI show all of the data points to show variability 


### MAKE THESE INTO A MULTIPANEL FIGURE WITH ALL 3 ON ONE LINE
## MALES
# Acquisition
barcolors = ['#FFE5A5','#FFE5A5','#FFBA08','#FFBA08']
dataMnorAcq = [[exploration_left_sal_M, exploration_right_sal_M], [exploration_left_pcp_M, exploration_right_pcp_M]]
labels = ['salL','salR', 'pcpL', 'pcpR']
ax, barx, barlist, sclist = barscatter(dataMnorAcq, transpose=False, paired=True, barfacecolor=barcolors, barfacecoloroption='individual',  ylabel='Time (s)', barlabels=labels, unequal=True) 
# Retention 
barcolors = ['#FFE5A5','#FFE5A5','#FFBA08','#FFBA08']
dataMnorRet = [[exploration_fam_sal_M, exploration_nov_sal_M],[exploration_fam_pcp_M, exploration_nov_pcp_M]]
labels = ['salF','salN', 'pcpF', 'pcpN']
ax, barx, barlist, sclist = barscatter(dataMnorRet, transpose=False, paired=True, barfacecolor=barcolors, barfacecoloroption='individual',  ylabel='Time (s)', barlabels=labels, unequal=True) 

# Discrimination index 
barcolors = ['#FFE5A5','#FFBA08']
dataMdi = [[DI_sal_M],[DI_pcp_M]]
labels = ['sal','pcp']
ax, barx, barlist, sclist = barscatter(dataMdi, transpose=False, paired=True, barfacecolor=barcolors, barfacecoloroption='individual',  ylabel='DI', barlabels=labels, unequal=True) 


#*******************************

## FEMALES
# ACQUISITION
barcolors = ['#AFDBD5','#AFDBD5','#249E8D','#249E8D']
dataFnorAcq = [[exploration_left_sal_F, exploration_right_sal_F], [exploration_left_pcp_F, exploration_right_pcp_F]]
labels = ['salL','salR', 'pcpL', 'pcpR']
ax, barx, barlist, sclist = barscatter(dataFnorAcq, transpose=False, paired=True, barfacecolor=barcolors, barfacecoloroption='individual',  ylabel='Time (s)', barlabels=labels) 
# RETENTION
barcolors = ['#AFDBD5','#AFDBD5','#249E8D','#249E8D']
dataFnorRet = [[exploration_fam_sal_F, exploration_nov_sal_F],[exploration_fam_pcp_F, exploration_nov_pcp_F]]
labels = ['salF','salN', 'pcpF', 'pcpN']
ax, barx, barlist, sclist = barscatter(dataFnorRet, transpose=False, paired=True, barfacecolor=barcolors, barfacecoloroption='individual',  ylabel='Time (s)', barlabels=labels) 
# DISCRIMINATION INDEX 
barcolors = ['#AFDBD5','#249E8D']
dataFdi = [[DI_sal_F],[DI_pcp_F]]
labels = ['sal','pcp']
ax, barx, barlist, sclist = barscatter(dataFdi, transpose=False, paired=True, barfacecolor=barcolors, barfacecoloroption='individual',  ylabel='DI', barlabels=labels, itemlabel=['1','2']) 


#plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/FILENAME.pdf", bbox_inches='tight')



## FIX THESE PROBLEMS HERE 

## Multipanel Males NOR
mpl.rcParams['figure.subplot.wspace'] = 0.3
figureB, ax = plt.subplots(nrows=1, ncols=3, figsize=(10,4)) ### x,y 
figureB.tight_layout(pad=3, w_pad=3, h_pad=1.0)

# Set out the different colours here and add manually in each call of bar scatter
barcolors = ['#FFE5A5','#FFE5A5','#FFBA08','#FFBA08']
labels = ['','','','']
ax[0], barx, barlist, sclist = barscatter(dataMnorAcq, ax=ax[0], transpose=False, paired=False, barfacecolor=barcolors, barfacecoloroption='individual',  ylabel='Time (s)', barlabels=labels, unequal=True) 
ax[1], barx, barlist, sclist = barscatter(dataMnorRet, ax=ax[1], transpose=False, paired=False, barfacecolor=barcolors, barfacecoloroption='individual',  ylabel='Time (s)', barlabels=labels, unequal=True) 
barcolors = ['#FFE5A5','#FFBA08']
ax[2], barx, barlist, sclist = barscatter(dataMdi, ax=ax[2],transpose=False, paired=False, barfacecolor=barcolors, barfacecoloroption='individual',  ylabel='Discrimination Index', barlabels=labels, unequal=True) 

# Multipanel Females NOR
barcolors = ['#AFDBD5','#AFDBD5','#249E8D','#249E8D']
mpl.rcParams['figure.subplot.wspace'] = 0.3
figureC, ax = plt.subplots(nrows=1, ncols=3, figsize=(10,4)) ### x,y 
figureC.tight_layout(pad=3, w_pad=3, h_pad=1.0)
ax[0], barx, barlist, sclist = barscatter(dataFnorAcq, ax=ax[0], transpose=False, paired=False, barfacecolor=barcolors, barfacecoloroption='individual',  ylabel='Time (s)', barlabels=labels) 
ax[1], barx, barlist, sclist = barscatter(dataFnorRet, ax=ax[1], transpose=False, paired=False, barfacecolor=barcolors, barfacecoloroption='individual',  ylabel='Time (s)', barlabels=labels) 
barcolors =  ['#AFDBD5','#249E8D']
ax[2], barx, barlist, sclist = barscatter(dataFdi, ax=ax[2], transpose=False, paired=False, barfacecolor=barcolors, barfacecoloroption='individual',  ylabel='Discrimination Index', barlabels=labels) 

## POST editing 
ax[0].set_ylabel('Time (s)')
ax[1].set_ylabel('Time (s)')
ax[2].set_ylabel('Discrimination index')
#Change these limits and add in the barcolour removal parameter
ax[0].set_xticks([])
ax[1].set_xticks([])
ax[2].set_xticks([])

ax[0].spines['bottom'].set_visible(False)


# Add this arg to all plots 
baredgecolor = ['']