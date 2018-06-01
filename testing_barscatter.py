#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 16:02:33 2018

@author: u1490431
"""

def barscatter(data, transpose = False,
                groupwidth = .75,
                barwidth = .9,
                paired = False,
                barfacecoloroption = 'same', # other options 'between' or 'individual'
                barfacecolor = ['white'],
                baredgecoloroption = 'same',
                baredgecolor = ['black'],
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
                title = 'none',
                grouplabel = 'auto',
                itemlabel = 'none',
                yaxisparams = 'auto',
                show_legend = 'none',
                legendloc='upper right',
                ax=[]):
#
#    if type(data) == float
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
        paired = True
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
    #colors = ['#1abc9c', '#f1c40f', '#d35400', '#3498db', '#8e44ad']
    colors = ['dodgerblue', 'hotpink', 'dodgerblue', 'hotpink']
    colors2 = ['k','k','k', 'k', 'k']
    colors3 = ['white', 'white','white']
    
    barfacecolorArray = setcolors("individual", colors, 1, 2, data, paired_scatter = False)
    baredgecolorArray = setcolors("individual", colors, 1, 2, data, paired_scatter = False)
     
    scfacecolorArray = setcolors("individual", colors3, 1, 2, data, paired_scatter = False)
    scedgecolorArray = setcolors("individual", colors2, 1, 2, data, paired_scatter = False)
 #   scfacecolorArray = setcolors("between", colors3, nGroups=nGroups, barspergroup=barspergroup, data=dataX, paired_scatter = True)
    
# Initialize figure
    if ax == []:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.tick_params(axis='both', which='major', labelsize=14)
        fig.tight_layout()
    
    # Make bars
    barlist = []
    barx = []
    for x, y, bfc, bec in zip(xvals.flatten(), barMeans.flatten(),
                              barfacecolorArray, baredgecolorArray):
        barx.append(x)
        barlist.append(ax.bar(x, y, barwidth,
                         facecolor = bfc, edgecolor = bec,
                         zorder=-1))
    
    # Make scatters
    sclist = []
    if paired == False:
        for x, Yarray, scf, sce  in zip(xvals.flatten(), data.flatten(),
                                        scfacecolorArray, scedgecolorArray):
            for y in Yarray:
                sclist.append(ax.scatter(x, y, s = scattersize,
                         c = scf,
                         edgecolors = sce,
                         zorder=1))

    else:
        try:
            np.shape(data)[1]
            for x, Yarray, scf, sce in zip(xvals, data, scfacecolorArray, scedgecolorArray):
                for y in np.transpose(Yarray.tolist()):
                    sclist.append(ax.plot(x, y, '-o', markersize = scattersize/10,
                             color = scatterlinecolor,
                             linewidth=linewidth,
                             markerfacecolor = scf,
                             markeredgecolor = sce))

# Explicitly added color here, issue with assignment of scf and sce 
        except IndexError:                    
            print(len(data[0]))
            for n,_ in enumerate(data[0]):
                y = [y[n-1] for y in data]
                sclist.append(ax.plot(xvals, y, '-o', markersize = scattersize/10,
                             color = 'grey',
                             linewidth=linewidth,
                             markerfacecolor = 'white',
                             markeredgecolor = 'k'))
                
    # Label axes
    if ylabel != 'none':
        plt.ylabel(ylabel, fontsize=14)
    
    if xlabel != 'none':
        plt.xlabel(xlabel)
        
    if title != 'none':
        plt.title(title, fontsize=14)
    
    # Set range and tick values for Y axis
    if yaxisparams != 'auto':
        ax.set_ylim(yaxisparams[0])
        plt.yticks(yaxisparams[1])
       
    # X ticks
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off') # labels along the bottom edge are off
    
    if grouplabel == 'auto':
        plt.tick_params(labelbottom='off')
    else:
        plt.xticks(range(1,nGroups+1), grouplabel)
        plt.tick_params(top='off')
    
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

    ax.set(ylabel='Mean pdp - notdistracted trials')
    #ax.set(ylabel='Percent distracted trials')
    ax.yaxis.label.set_size(14)      
#    fig.savefig('/Volumes/KPMSB352/Distraction photometry paper/BehaviourFigs/PDP_notdis_salvpcp_M.pdf', bbox_inches="tight") 

    return ax, barx, barlist, sclist


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

