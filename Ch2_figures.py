#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 09:11:26 2018

@author: u1490431
"""

"""
Chapter 2 Plotting figures

Barscatters with individual data 


"""

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
                barMeans[i][j] = np.nanmean(data[i][j])
                items[i][j] = len(data[i][j])
        
    else:
        grouped = False
        barspergroup = 1
        
        for i in range(np.shape(data)[0]):
            barMeans[i] = np.nanmean(data[i])
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
                         color = scatterlinecolor,
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



import matplotlib as mpl 

mpl.rcParams['font.size'] = 14
label_size = 14
plt.rcParams['xtick.labelsize'] = label_size 
plt.rcParams['ytick.labelsize'] = label_size 

# Barscatter, licking across days 1,2,3 (all rats DIS1)
lickdata = [nlicks_minus6, nlicks_minus5, nlicks_minus4, nlicks_minus3, nlicks_minus2, nlicks_minus1, nlicks]
col = ['#81D4FA','#03A9F4','#0277BD']
labels = ['1','2','3','4','5','6','7']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(lickdata, transpose=False, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Licks', xlabel='Lick training days', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/LickDaysBarscatter.pdf', bbox_inches='tight')

#  !!!!!!    REMEMBER TO STOP THE IMPORT OF SEABORN BEFORE MAKING BARSCATTERS 
# Light blue - sky 
#81D4FA
#4FC3F7  
#03A9F4
#0277BD 

# Pinks 

#F06292  
#E91E63
#880E4F


# Barscatter, percent distracted across days (all rats, DIS1)
percentdis_data = [[percent_dis_modelled,percent_dis_dis,percent_dis_hab1,percent_dis_hab2,percent_dis_sal,percent_dis_amph ]]
col = ['gray','#03A9F4','#03A9F4','#03A9F4','#E91E63','#880E4F']
labels = ['mod','dis','hab1','hab2','sal','amph']
# Males licking on 3 lick days 
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(percentdis_data, transpose=False, ax=ax,paired=True, barfacecolor=col, barfacecoloroption='individual',  ylabel='Percent distracted', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/PercentDistractedBarscatter.pdf', bbox_inches='tight')


# Bar scatter, distraction day % modalitites (white noise, tone, combined)
# Or white noise and non white noise 

# Barscatter different licking parameters (actually appears first in the chapter)

# Here:
    
# Licks across days already done so not needed 

# Mean licks per burst correlation against mean licks per clusters


## (1) Number of bursts by number of clusters 
#import matplotlib as mpl 
#import seaborn as sn
#sn.set_style("darkgrid")
#mpl.rcParams['font.size'] = 14
## Number of bursts by number of clusters
#slope, intercept, r_value, p_value, std_err = stats.linregress(all_n_bursts, all_n_runs)
### Plot the scatter 
#plt.plot(all_n_bursts, all_n_runs,'o', color='darkgrey')
### Add line of best fit
#plt.plot(np.asarray(all_n_bursts), intercept+slope*np.asarray(all_n_bursts), 'darkgrey')
#
#plt.legend()
#sn.despine(offset=10, trim=True); 
#plt.xlabel('Number of bursts', fontsize=14)
#plt.ylabel('Number of clusters', fontsize=14)
#plt.xticks(fontsize=14)
#plt.yticks(fontsize=14)
##plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/Corr_DIvs%_M.pdf", bbox_inches='tight')
#plt.show()
#print('Linear regression, bursts vs clusters')
#print('SALINE')
#print('R squared = ',r_value**2, ', p value = ', p_value)
#
#
## (2) Licks per burst, licks per cluster (means) - have the means bar at the end too
#
## Licks per burst by licks per cluster
#slope, intercept, r_value, p_value, std_err = stats.linregress(all_mean_burst_length, all_mean_run_length)
### Plot the scatter 
#plt.plot(all_mean_burst_length, all_mean_run_length,'o', color='darkgrey')
### Add line of best fit
#plt.plot(np.asarray(all_mean_burst_length), intercept+slope*np.asarray(all_mean_burst_length), 'darkgrey')
#
#plt.legend()
#sn.despine(offset=10, trim=True); 
#plt.xlabel('Mean licks per burst', fontsize=14)
#plt.ylabel('Mean licks per cluster', fontsize=14)
#plt.xticks(fontsize=14)
#plt.yticks(fontsize=14)
##plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/Corr_DIvs%_M.pdf", bbox_inches='tight')
#plt.show()
#print('Linear regression, licks per burst vs cluster')
#print('SALINE')
#print('R squared = ',r_value**2, ', p value = ', p_value)
#
## (3) Inter burst interval vs intercluster interval 
#
## Number of bursts by mber of nuclusters
#slope, intercept, r_value, p_value, std_err = stats.linregress(all_mean_IBI, all_mean_IRI)
### Plot the scatter 
#plt.plot(all_mean_IBI, all_mean_IRI,'o', color='darkgrey')
### Add line of best fit
#plt.plot(np.asarray(all_mean_IBI), intercept+slope*np.asarray(all_mean_IBI), 'darkgrey')
#
#plt.legend()
#sn.despine(offset=10, trim=True); 
#plt.xlabel('Mean inter-burst interval', fontsize=14)
#plt.ylabel('Mean inter-cluster interval', fontsize=14)
#plt.xticks(fontsize=14)
#plt.yticks(fontsize=14)
##plt.savefig("/Volumes/KPMSB352/Thesis/Chapter 3 - Distraction pcp model/Figures/Corr_DIvs%_M.pdf", bbox_inches='tight')
#plt.show()
#print('Linear regression, bursts vs clusters')
#print('SALINE')
#print('R squared = ',r_value**2, ', p value = ', p_value)
#




#Licking paramters - N LICKS
mean_licks = np.mean(nlicks)
nlick_data = [[nlicks[0]],[nlicks[1]],[nlicks[2]],[nlicks[3]],[nlicks[4]],[nlicks[5]],[nlicks[6]],[nlicks[7]], nlicks]
col = ['lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','powderblue']
labels = ['rat 1', 'rat 2', 'rat 3', 'rat 4', 'rat 5', 'rat 6', 'rat 7', 'rat 8', 'mean']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(nlick_data, transpose=False, ax=ax,paired=False, barfacecolor=col, barfacecoloroption='individual',  ylabel='Licks', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45, unequal=False) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/NLicksBarscatter.pdf', bbox_inches='tight')

#Licking paramters - N BURSTS
mean_n_bursts = np.mean(all_n_bursts)
nbursts_data = [[all_n_bursts[0]],[all_n_bursts[1]],[all_n_bursts[2]],[all_n_bursts[3]],[all_n_bursts[4]],[all_n_bursts[5]],[all_n_bursts[6]],[all_n_bursts[7]], all_n_bursts]
col = ['lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','powderblue']
labels = ['rat 1', 'rat 2', 'rat 3', 'rat 4', 'rat 5', 'rat 6', 'rat 7', 'rat 8', 'mean']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(nbursts_data, transpose=False, ax=ax,paired=False, barfacecolor=col, barfacecoloroption='individual',  ylabel='Number of bursts', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45, unequal=False) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/NBurstsBarscatter.pdf', bbox_inches='tight')

#Licking parameters - LICKS PER BURST
# Maybe here you should have all of the burst lengths for each rat? 
mean_burstlen = np.mean(all_mean_burst_length)
nburst_len_data = [[all_mean_burst_length[0]],[all_mean_burst_length[1]],[all_mean_burst_length[2]],[all_mean_burst_length[3]],[all_mean_burst_length[4]],[all_mean_burst_length[5]],[all_mean_burst_length[6]],[all_mean_burst_length[7]], all_mean_burst_length]
col = ['lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','powderblue']
labels = ['rat 1', 'rat 2', 'rat 3', 'rat 4', 'rat 5', 'rat 6', 'rat 7', 'rat 8', 'mean']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(nburst_len_data, transpose=False, ax=ax,paired=False, barfacecolor=col, barfacecoloroption='individual',  ylabel='Mean licks per burst', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40               
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/MeanBurstLenBarscatter.pdf', bbox_inches='tight')



# Percentage distracted modalities 
modality_data = [[percent_dis_whitenoise,percent_dis_tone, percent_dis_combined]]
col = ['powderblue','blue','darkblue']
labels = ['white noise', 'tone', 'combined']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(modality_data, transpose=False, ax=ax,paired=False, barfacecolor=col, barfacecoloroption='individual',  ylabel='Mean licks per burst', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40               
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/ModalityBarscatter.pdf', bbox_inches='tight')


#Licking paramters - N CLUSTERS
mean_n_runs = np.mean(all_n_runs)
nruns_data = [[all_n_runs[0]],[all_n_runs[1]],[all_n_runs[2]],[all_n_runs[3]],[all_n_runs[4]],[all_n_runs[5]],[all_n_runs[6]],[all_n_runs[7]], all_n_runs]
col = ['lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightpink']
labels = ['rat 1', 'rat 2', 'rat 3', 'rat 4', 'rat 5', 'rat 6', 'rat 7', 'rat 8', 'mean']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(nruns_data, transpose=False, ax=ax,paired=False, barfacecolor=col, barfacecoloroption='individual',  ylabel='Number of clusters', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45, unequal=False) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/NRunsBarscatter.pdf', bbox_inches='tight')


#Licking parameters - LICKS PER CLUSTER
## Maybe here you should have all of the burst lengths for each rat? 
mean_runlen = np.mean(all_mean_run_length)
nrun_len_data = [[all_mean_run_length[0]],[all_mean_run_length[1]],[all_mean_run_length[2]],[all_mean_run_length[3]],[all_mean_run_length[4]],[all_mean_run_length[5]],[all_mean_run_length[6]],[all_mean_run_length[7]], all_mean_run_length]
col = ['lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightpink']
labels = ['rat 1', 'rat 2', 'rat 3', 'rat 4', 'rat 5', 'rat 6', 'rat 7', 'rat 8', 'mean']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(nrun_len_data, transpose=False, ax=ax,paired=False, barfacecolor=col, barfacecoloroption='individual',  ylabel='Mean licks per cluster', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40               
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/MeanRunLenBarscatter.pdf', bbox_inches='tight')


##### Inter-cluster (IRI) ICI - could have all the IRIs on the plot if wanted (need to extract though)
mean_IRI = np.mean(all_mean_IRI)
IRI_data = [[all_mean_IRI[0]],[all_mean_IRI[1]],[all_mean_IRI[2]],[all_mean_IRI[3]],[all_mean_IRI[4]],[all_mean_IRI[5]],[all_mean_IRI[6]],[all_mean_IRI[7]], all_mean_IRI]
col = ['lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightpink']
labels = ['rat 1', 'rat 2', 'rat 3', 'rat 4', 'rat 5', 'rat 6', 'rat 7', 'rat 8', 'mean']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(IRI_data, transpose=False, ax=ax,paired=False, barfacecolor=col, barfacecoloroption='individual',  ylabel='Mean inter cluster interval (s)', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40               
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/MeanIRIBarscatter.pdf', bbox_inches='tight')

### Inter-burst (IBI) - could have all the IBIs on the plot if wanted 
##### Inter-cluster (IRI) ICI - could have all the IRIs on the plot if wanted (need to extract though)
mean_IBI = np.mean(all_mean_IBI)
IBI_data = [[all_mean_IBI[0]],[all_mean_IBI[1]],[all_mean_IBI[2]],[all_mean_IBI[3]],[all_mean_IBI[4]],[all_mean_IBI[5]],[all_mean_IBI[6]],[all_mean_IBI[7]], all_mean_IBI]
col = ['lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','powderblue']
labels = ['rat 1', 'rat 2', 'rat 3', 'rat 4', 'rat 5', 'rat 6', 'rat 7', 'rat 8', 'mean']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(IBI_data, transpose=False, ax=ax,paired=False, barfacecolor=col, barfacecoloroption='individual',  ylabel='Mean inter burst interval (s)', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40               
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/MeanIBIBarscatter.pdf', bbox_inches='tight')

## Licking frequency (check that it is 6-7Hz)
freq_data = [[lick_analysis[0]['freq']], [lick_analysis[1]['freq']],[lick_analysis[2]['freq']],[lick_analysis[3]['freq']],[lick_analysis[4]['freq']],[lick_analysis[5]['freq']],[lick_analysis[6]['freq']],[lick_analysis[7]['freq']]]
freq_data2 = [[lick_analysis[0]['freq']], [lick_analysis[1]['freq']],[lick_analysis[2]['freq']],[lick_analysis[3]['freq']],[lick_analysis[4]['freq']],[lick_analysis[5]['freq']],[lick_analysis[6]['freq']],[lick_analysis[7]['freq']], freq_data]
col = ['lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','k']
labels = ['rat 1', 'rat 2', 'rat 3', 'rat 4', 'rat 5', 'rat 6', 'rat 7', 'rat 8', 'mean']
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4)) ### x,y
ax, barx, barlist, sclist = barscatter(freq_data2, transpose=False, ax=ax,paired=False, barfacecolor=col, barfacecoloroption='individual',  ylabel='Licking frequency (Hz)', barlabels=labels, itemlabel=['1','2'], barlabeloffset=0.05, xrotation=45) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.xaxis.labelpad = 40               
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 2 - Distraction pilot/Figures/FreqBarscatter.pdf', bbox_inches='tight')








