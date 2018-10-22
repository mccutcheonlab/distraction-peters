#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 20:17:44 2018

@author: u1490431
"""

''' Makes all of the bar scatters, not just for distraction, for licking 
    Need to add labels and titles so we know which plts are which 
    Check which scripts need to be run first and make sure to comment out 
    the seaborn import that causes the plotting problems (no black outline)
    and strange grid 
    
'''   
    


## run this script AFTER the ch4 analysis script (photometry licking days and dis, not barscatter)
def MetaExtractorTHPH (metafile):
    f = open(metafile, 'r')
    f.seek(0)
    Metafilerows = f.readlines()[1:]
    tablerows = []

    for row in Metafilerows: 
        items = row.split(',')
        tablerows.append(items)

    TDTFile, MedFilenames, RatID, Date, Session, Include, Licks, Totdistractors, Distracted, \
    Percentdistracted, Note, Endcol = [],[],[],[],[],[],[],[],[],[],[],[]

    for i, lst in enumerate(tablerows):
        
        
        
       TDTFile = TDTFile + [lst[0]]
       MedFilenames = MedFilenames + [lst[1]]
       RatID = RatID + [lst[2]]
       Date = Date + [lst[3]]
       Session = Session + [lst[4]]
       Include = Include + [lst[5]]
       Licks = Licks + [lst[6]]
       Totdistractors = Totdistractors + [lst[7]]  
       Distracted = Distracted + [lst[8]]
       Percentdistracted = Percentdistracted + [lst[9]]
       Note = Note + [lst[10]]
       Endcol = Endcol + [lst[11]]

 
    return ({'MedFilenames':MedFilenames, 'RatID':RatID, 'Date':Date, 'Session':Session, \
             'Include':Include, 'Licks':Licks, 'PercentDistracted':Percentdistracted})



def subsetter2(dictionary, dates, dis=False, verbose=False):
    '''
    SUBSETTER KP
    # Subsets data according to date, reads in dictionnary produced from metafile
    # and subsets into variable based on date(s) and drug condition 
    # if distraction day argument is given as True adds the distractor type 
    # to the output lists for later processing 
    
    Adapted to include & include == 1 (because not all had the same date for last lick day)
    In metafile do not include 2.7 and 2.8 on the lick 3 day but do include them on lick 6
     
    
    '''
    subset = []
    
    for ind, filename in enumerate(dictionary['MedFilenames']):
        path = medfolder + filename
        onsets, offsets, med_dis_times, dis_type = medfilereader(path, ['b', 'c', 'i', 'j'], remove_var_header = True)  # e onset, f offset

        if dis == True:
            if dictionary['Date'][ind] in dates and dictionary['Include'][ind] == '1':
                subset.append([onsets, offsets, dis_type, dictionary['RatID'][ind]])
                
        elif dis==False:   
            if dictionary['Date'][ind] in dates and dictionary['Include'][ind] == '1':
                subset.append([onsets, offsets, dictionary['RatID'][ind]])
            
        if verbose: #assumes true
            print('filename, or comment ...') 
    return subset



def percentdisgroup(distractiondict):
    ''' Discalc_sal_M == distractiondict '''
    
    percent_dis_group = []
    for rat in distractiondict: 
        percentage = len(rat[0]) / (len(rat[0])+len(rat[1])) * 100
        percent_dis_group.append(percentage)
    return percent_dis_group



def disbygroup(dictionary):
    ''' Prodcues times of distracted and not distracted as 2 lists
    takes a dictionary of grouped rat data 
    '''

    dis = []
    for rat in dictionary:
        
        discalc = distractionCalc2(rat[0])         
        distracted, notdistracted = distractedOrNot(discalc, rat[0])
        dis.append([distracted, notdistracted])
        
    return dis

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




################################################################################################
################################################################################################

## BEHAVIOUR SUBSETTING, PERCENT DISTRACTED AND PLOTS [THPH1, THPH2]

metafile = '/Volumes/KP_HARD_DRI/kp259/THPH1AND2/THPH1&2Metafile.csv'
extract_data = MetaExtractorTHPH(metafile)
# Folder with all medfiles (THPH1 and THPH2)
medfolder = '/Volumes/KP_HARD_DRI/kp259/THPH1AND2/med/'



''' 
  Info DPCP1(16) and DPCP2(16) males, DPCP3(24) females :
      
  THPH1      |THPH2.1 - 2.6|THPH2.7 - 2.8|CONDITION
  -----------|-------------|-------------|-----------
  170620     |170809       |170814       |last lick
  170621     |170810       |170815       |dis
  170622     |170811       |170816       |hab 1
  
  All included, can later only subset the first 12 rats and ignore the 
  unequal numbers from no habituation on rat 2.8

'''


last_lick = subsetter2(extract_data, ['170620', '170809', '170814'])
distraction= subsetter2(extract_data, ['170621', '170810', '170815'], dis=True)
habituation = subsetter2(extract_data, ['170622', '170811', '170806'], dis=True)


discalcLick = []
discalcDis = []
discalcHab = []

for rat in last_lick:
    discalc = distractionCalc2(rat[0])
    distracted, notdistracted = distractedOrNot(discalc, rat[0])
    discalcLick.append([distracted, notdistracted])

for rat in distraction:
    discalc = distractionCalc2(rat[0])
    distracted, notdistracted = distractedOrNot(discalc, rat[0])
    discalcDis.append([distracted, notdistracted])
    
for rat in habituation:
    discalc = distractionCalc2(rat[0])
    distracted, notdistracted = distractedOrNot(discalc, rat[0])
    discalcHab.append([distracted, notdistracted])    

percentdisLick = []
percentdisDis = []
percentdisHab = [] 

for rat in discalcLick:
    percent = len(rat[0]) / (len(rat[0])+len(rat[1])) * 100
    percentdisLick.append(percent)

for rat in discalcDis:
    percent = len(rat[0]) / (len(rat[0])+len(rat[1])) * 100
    percentdisDis.append(percent)
    
for rat in discalcHab:
    percent = len(rat[0]) / (len(rat[0])+len(rat[1])) * 100
    percentdisHab.append(percent)    
    

##########################################################################################

# PERCENT DISTRACTED - THPH1 AND THPH2
    
# For barscatter plots excluded rats 2.7 and 2.8 (as not all point on all days)
data = [[percentdisLick[0:12], percentdisDis[0:12], percentdisHab]]

col3 = ['#FFE5A5','#FFE5A5','#FFE5A5']
col3 = ['#a2e283','#4cbb17','#257200']
labels = ['mod', 'dis', 'hab']
mpl.rcParams['font.size'] = 14
figureA, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,5)) ### x,y
ax, barx, barlist, sclist = barscatter(data, transpose=False, ax=ax,paired=True, barfacecolor=col3, barlabels=labels,barfacecoloroption='individual',  ylabel='Percent distracted (%)', itemlabel=['1','2'], barlabeloffset=0.05) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax.spines['bottom'].set_visible(False)
#figureA.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/PercentDisBarScatter.pdf', bbox_inches="tight")

   
 # Barscatters to make: 

## PERHAPS only care about distractORS as not many distracted trials on this habituation day? 

## UNEQUAL NUMBERS HERE, ISSUE FOR BARSCATTER - Maybe just index missing the LAST rat for ALL comparisons 

################################################################################################
################################################################################################

## Short versus long runs of licks ORANGES X 2 (ORANGE AND BRICK OR GREY AND BRICK)
def MultBy100(list):
    output = [x*100 for x in list]
    
    return output
shortVlongPeak = [MultBy100(peak_short_runs), MultBy100(peak_long_runs)]
shortVlongt = [t_short_runs, t_long_runs]
shortVlongPre = [MultBy100(pre_short_runs), MultBy100(pre_long_runs)]
shortVlongPost = [MultBy100(post_short_runs), MultBy100(post_long_runs)]

mpl.rcParams['figure.subplot.wspace'] = 0.6
mpl.rcParams['figure.subplot.right'] = 1
mpl.rcParams['font.size'] = 14
figureA, ax = plt.subplots(nrows=1, ncols=4, figsize=(10,4)) ### x,y 
figureA.tight_layout(pad=3, w_pad=3, h_pad=1.0)

labels = []
ax[0], barx, barlist, sclist = barscatter(shortVlongPeak, ax=ax[0],transpose=False, paired=True, barfacecolor=['#ffb349','#ef7700'], barfacecoloroption='individual',  ylabel='Peak (%)', barlabels=labels, baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[1], barx, barlist, sclist = barscatter(shortVlongt, ax=ax[1], transpose=False, paired=True, barfacecolor=['#ffb349','#ef7700'], barfacecoloroption='individual',  ylabel='t (s)', barlabels=labels, baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[2], barx, barlist, sclist = barscatter(shortVlongPre, ax=ax[2],transpose=False, paired=True, barfacecolor=['#ffb349','#ef7700'], barfacecoloroption='individual',  ylabel='Pre-event period (mean %)', barlabels=labels,  baredgecolor=[''] )#,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[3], barx, barlist, sclist = barscatter(shortVlongPost, ax=ax[3],transpose=False, paired=True, barfacecolor=['#ffb349','#ef7700'], barfacecoloroption='individual',  ylabel='Post-event period (mean %)', barlabels=labels, baredgecolor=[''] )#,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])

ax[0].set_ylabel('Peak (% ΔF)')
ax[1].set_ylabel('t (s)')
ax[2].set_ylabel('Pre-event period (mean % ΔF)')
ax[3].set_ylabel('Post-event period (mean % ΔF)')

ax[0].set_xticks([])
#ax[0].set_ylim([0,4000])
ax[1].set_xticks([])
#ax[1].set_ylim([0,25])
ax[2].set_xticks([])
#ax[2].set_ylim(0,1200)
ax[3].set_xticks([])
#ax[3].set_ylim(0,1200)

ax[0].spines['bottom'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[2].spines['bottom'].set_visible(False)
ax[3].spines['bottom'].set_visible(False)
#figureA.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/ShortVsLongPeaksBarScatter.pdf', bbox_inches="tight")

################################################################################################
################################################################################################

## All runs versus all distractors (BRICK AND GREEN or ORANGE AND GREEN)
runVdisPeak = [MultBy100(peak_runs), MultBy100(peak_distractor)]
runVdist = [t_runs, t_distractor]
runVdisPre = [MultBy100(pre_runs), MultBy100(pre_distractor)]
runVdisPost = [MultBy100(post_runs), MultBy100(post_distractor)]

figureA, ax = plt.subplots(nrows=1, ncols=4, figsize=(10,4)) ### x,y 
figureA.tight_layout(pad=3, w_pad=3, h_pad=1.0)

labels = []
ax[0], barx, barlist, sclist = barscatter(runVdisPeak, ax=ax[0],transpose=False, paired=True, barfacecolor=['#ef7700','#4cbb17'], barfacecoloroption='individual',  ylabel='Peak (%)', barlabels=labels, baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[1], barx, barlist, sclist = barscatter(runVdist, ax=ax[1], transpose=False, paired=True, barfacecolor=['#ef7700','#4cbb17'], barfacecoloroption='individual',  ylabel='t (s)', barlabels=labels, baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[2], barx, barlist, sclist = barscatter(runVdisPre, ax=ax[2],transpose=False, paired=True, barfacecolor=['#ef7700','#4cbb17'], barfacecoloroption='individual',  ylabel='Pre-event period (mean %)', barlabels=labels,  baredgecolor=[''] )#,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[3], barx, barlist, sclist = barscatter(runVdisPost, ax=ax[3],transpose=False, paired=True, barfacecolor=['#ef7700','#4cbb17'], barfacecoloroption='individual',  ylabel='Post-event period (mean %)', barlabels=labels, baredgecolor=[''] )#,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])

ax[0].set_ylabel('Peak (% ΔF)')
ax[1].set_ylabel('t (s)')
ax[2].set_ylabel('Pre-event period (mean % ΔF)')
ax[3].set_ylabel('Post-event period (mean % ΔF)')

ax[0].set_xticks([])
#ax[0].set_ylim([0,4000])
ax[1].set_xticks([])
#ax[1].set_ylim([0,25])
ax[2].set_xticks([])
#ax[2].set_ylim(0,1200)
ax[3].set_xticks([])
#ax[3].set_ylim(0,1200)

ax[0].spines['bottom'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[2].spines['bottom'].set_visible(False)
ax[3].spines['bottom'].set_visible(False)
#figureA.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/RunVsDisPeaksBarScatter.pdf', bbox_inches="tight")


################################################################################################
################################################################################################

## Distracted vs not distracted (GREEN X 2)
disVnotPeak = [MultBy100(peak_notdistracted),MultBy100(peak_distracted)] 
disVnott = [t_notdistracted,t_distracted]
disVnotPre = [MultBy100(pre_notdistracted),MultBy100(pre_distracted)]
disVnotPost = [MultBy100(post_notdistracted),MultBy100(post_distracted)]

figureA, ax = plt.subplots(nrows=1, ncols=4, figsize=(10,4)) ### x,y 
figureA.tight_layout(pad=3, w_pad=3, h_pad=1.0)

labels = []
ax[0], barx, barlist, sclist = barscatter(disVnotPeak, ax=ax[0],transpose=False, paired=True, barfacecolor=['#257200','#4cbb17'], barfacecoloroption='individual',  ylabel='Peak (%)', barlabels=labels, baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[1], barx, barlist, sclist = barscatter(disVnott, ax=ax[1], transpose=False, paired=True, barfacecolor=['#257200','#4cbb17'], barfacecoloroption='individual',  ylabel='t (s)', barlabels=labels, baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[2], barx, barlist, sclist = barscatter(disVnotPre, ax=ax[2],transpose=False, paired=True, barfacecolor=['#257200','#4cbb17'], barfacecoloroption='individual',  ylabel='Pre-event period (mean %)', barlabels=labels,  baredgecolor=[''] )#,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[3], barx, barlist, sclist = barscatter(disVnotPost, ax=ax[3],transpose=False, paired=True, barfacecolor=['#257200','#4cbb17'], barfacecoloroption='individual',  ylabel='Post-event period (mean %)', barlabels=labels, baredgecolor=[''] )#,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])

ax[0].set_ylabel('Peak (% ΔF)')
ax[1].set_ylabel('t (s)')
ax[2].set_ylabel('Pre-event period (mean % ΔF)')
ax[3].set_ylabel('Post-event period (mean % ΔF)')

ax[0].set_xticks([])
#ax[0].set_ylim([0,4000])
ax[1].set_xticks([])
#ax[1].set_ylim([0,25])
ax[2].set_xticks([])
#ax[2].set_ylim(0,1200)
ax[3].set_xticks([])
#ax[3].set_ylim(0,1200)

ax[0].spines['bottom'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[2].spines['bottom'].set_visible(False)
ax[3].spines['bottom'].set_visible(False)
figureA.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/DisVsNotDisPeaksBarScatter.pdf', bbox_inches="tight")


################################################################################################
################################################################################################

## Modelled versus distracion day presented distractors (GREY and GREEN)
modVdisPeak = [MultBy100(peak_distractorMOD), MultBy100(peak_distractor)]
modVdist = [t_distractorMOD, t_distractor]
modVdisPre = [MultBy100(pre_distractorMOD), MultBy100(pre_distractor)]
modVdisPost = [MultBy100(post_distractorMOD), MultBy100(post_distractor)]

figureA, ax = plt.subplots(nrows=1, ncols=4, figsize=(10,4)) ### x,y 
figureA.tight_layout(pad=3, w_pad=3, h_pad=1.0)

labels = []
ax[0], barx, barlist, sclist = barscatter(modVdisPeak, ax=ax[0],transpose=False, paired=True, barfacecolor=['lightgray','#4cbb17'], barfacecoloroption='individual',  ylabel='Peak (%)', barlabels=labels, baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[1], barx, barlist, sclist = barscatter(modVdist, ax=ax[1], transpose=False, paired=True, barfacecolor=['lightgray','#4cbb17'], barfacecoloroption='individual',  ylabel='t (s)', barlabels=labels, baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[2], barx, barlist, sclist = barscatter(modVdisPre, ax=ax[2],transpose=False, paired=True, barfacecolor=['lightgray','#4cbb17'], barfacecoloroption='individual',  ylabel='Pre-event period (mean %)', barlabels=labels,  baredgecolor=[''] )#,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[3], barx, barlist, sclist = barscatter(modVdisPost, ax=ax[3],transpose=False, paired=True, barfacecolor=['lightgray','#4cbb17'], barfacecoloroption='individual',  ylabel='Post-event period (mean %)', barlabels=labels, baredgecolor=[''] )#,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])

ax[0].set_ylabel('Peak (% ΔF)')
ax[1].set_ylabel('t (s)')
ax[2].set_ylabel('Pre-event period (mean % ΔF)')
ax[3].set_ylabel('Post-event period (mean % ΔF)')

ax[0].set_xticks([])
#ax[0].set_ylim([0,4000])
ax[1].set_xticks([])
#ax[1].set_ylim([0,25])
ax[2].set_xticks([])
#ax[2].set_ylim(0,1200)
ax[3].set_xticks([])
#ax[3].set_ylim(0,1200)

ax[0].spines['bottom'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[2].spines['bottom'].set_visible(False)
ax[3].spines['bottom'].set_visible(False)
figureA.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/ModVsDisPeaksBarScatter.pdf', bbox_inches="tight")

################################################################################################
################################################################################################
## Distraction day vs habituation day (GREEN, light and dark)
disVhabPeak = [MultBy100(peak_distractor), MultBy100(peak_distractorHAB)]
disVhabt = [t_distractor, t_distractorHAB]
disVhabPre = [MultBy100(pre_distractor), MultBy100(pre_distractorHAB)]
disVhabPost = [MultBy100(post_distractor), MultBy100(post_distractorHAB)] 

figureA, ax = plt.subplots(nrows=1, ncols=4, figsize=(10,4)) ### x,y 
figureA.tight_layout(pad=3, w_pad=3, h_pad=1.0)

labels = []
ax[0], barx, barlist, sclist = barscatter(disVhabPeak, ax=ax[0],transpose=False, paired=True, barfacecolor=['#a2e283','#257200'], barfacecoloroption='individual',  ylabel='Peak (%)', barlabels=labels, baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[1], barx, barlist, sclist = barscatter(disVhabt, ax=ax[1], transpose=False, paired=True, barfacecolor=['#a2e283','#257200'], barfacecoloroption='individual',  ylabel='t (s)', barlabels=labels, baredgecolor=['']) #,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[2], barx, barlist, sclist = barscatter(disVhabPre, ax=ax[2],transpose=False, paired=True, barfacecolor=['#a2e283','#257200'], barfacecoloroption='individual',  ylabel='Pre-event period (mean %)', barlabels=labels,  baredgecolor=[''] )#,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])
ax[3], barx, barlist, sclist = barscatter(disVhabPost, ax=ax[3],transpose=False, paired=True, barfacecolor=['#a2e283','#257200'], barfacecoloroption='individual',  ylabel='Post-event period (mean %)', barlabels=labels, baredgecolor=[''] )#,grouplabel=['Sal', 'Pcp', 'day -2', 'day -1'])

ax[0].set_ylabel('Peak (% ΔF)')
ax[1].set_ylabel('t (s)')
ax[2].set_ylabel('Pre-event period (mean % ΔF)')
ax[3].set_ylabel('Post-event period (mean % ΔF)')

ax[0].set_xticks([])
#ax[0].set_ylim([0,4000])
ax[1].set_xticks([])
#ax[1].set_ylim([0,25])
ax[2].set_xticks([])
#ax[2].set_ylim(0,1200)
ax[3].set_xticks([])
#ax[3].set_ylim(0,1200)

ax[0].spines['bottom'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[2].spines['bottom'].set_visible(False)
ax[3].spines['bottom'].set_visible(False)
#figureA.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/DisVsHabPeaksBarScatter.pdf', bbox_inches="tight")

#???????

# ANALYSIS FOR WHITE NOISE VS NON-WHITE NOISE HERE 
# (7) 4 bars --> wn vs nwn dis day vs hab day - peak, t, pre, post 



## Correlation between PERCENT DISTRACTED and DISTRCTOR PEAK 
import matplotlib as mpl 
from scipy import stats
import seaborn as sn
slope, intercept, r_value, p_value, std_err = stats.linregress(percentdisDis[2:], MultBy100(peak_distractor))
## Plot the scatter
plt.plot(percentdisDis[2:], MultBy100(peak_distractor),'o', color='limegreen')
## Add line of best fit
plt.plot(np.asarray(percentdisDis[2:]), intercept+slope*np.asarray(percentdisDis[2:]), 'limegreen')

plt.legend()
sn.despine(offset=10, trim=True); 

plt.xlabel('Percent distracted', fontsize=14)
plt.ylabel('Peak response to distractors (% ΔF)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Corr_Percent_DistractorPeak.pdf', bbox_inches="tight")
plt.show()
print('Linear regression, Percent distracted VS distractor peak')
print('R squared = ',r_value**2, ', p value = ', p_value)


## Correlation between PERCENT DISTRACTED and DISTRACTED PEAK 
slope, intercept, r_value, p_value, std_err = stats.linregress(percentdisDis[2:], MultBy100(peak_distracted))
## Plot the scatter
plt.plot(percentdisDis[2:], MultBy100(peak_distracted),'o', color='#257200')
## Add line of best fit
plt.plot(np.asarray(percentdisDis[2:]), intercept+slope*np.asarray(percentdisDis[2:]) ,'#257200')

plt.legend()
sn.despine(offset=10, trim=True); 

plt.xlabel('Percent distracted', fontsize=14)
plt.ylabel('Peak response distracted (% ΔF)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Corr_Percent_DistractedPeak.pdf', bbox_inches="tight")
plt.show()
print('Linear regression, Percent distracted VS distracted peak')
print('R squared = ',r_value**2, ', p value = ', p_value)




## Correlation between PERCENT DISTRACTED and RUNS PEAK 
slope, intercept, r_value, p_value, std_err = stats.linregress(percentdisDis[2:], MultBy100(peak_runs))
## Plot the scatter
plt.plot(percentdisDis[2:], MultBy100(peak_runs),'o', color='darkorange')
## Add line of best fit
plt.plot(np.asarray(percentdisDis[2:]), intercept+slope*np.asarray(percentdisDis[2:]), 'darkorange')

plt.legend()
sn.despine(offset=10, trim=True); 

plt.xlabel('Percent distracted', fontsize=14)
plt.ylabel('Peak response runs (% ΔF)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('/Volumes/KPMSB352/Thesis/Chapter 4 - Photometry VTA/Figures/Corr_Percent_RunsPeak.pdf', bbox_inches="tight")

plt.show()
print('Linear regression, Percent distracted VS all runs peak')
print('R squared = ',r_value**2, ', p value = ', p_value)


################################################################################################
################################################################################################


#  White noise separation - thph1 and thph2

def discalc_modalities(dictionary, modalitykey):
    ''' Calculates distractors, distracted and modalities for dictionary of 
    rats by group for just distraction day only 
    
    Calculates a grouped percentage of each modality and how distracted 
    by that modality rats are on average (by group)
    
    '''
    percent_dis_whitenoise_group = []
    percent_dis_tone_group = []
    percent_dis_combined_group = []
    percent_dis_all_non_WN_group = []
    discalcgroup = []
    ## SAL MALES - DISTRACTION DAY ONLY - DISTRACTOR TYPE ANALYSIS INCLUDED
# Finds distracted or not (corrects for med slipping issue)
    for rat in dictionary:
        print('a')
        discalc = distractionCalc2(rat[0])
        distracted, notdistracted = distractedOrNot(discalc, rat[0])
        print(len(distracted), len(notdistracted))
      #  work out percentage and add this too 
        discalcgroup.append([distracted, notdistracted])
    
        dis_numeric = []
        ndis_numeric = []
    # Modality analysis - calculates which distractors contain different features (whitenoise, tone or combination)
    # Then works out on how many of these trials rats are distracted (individual) before creating a mean 
    # Tr and except because there are some cases where there are ZERO distracted or non distracted trials (or both)
        for d in distracted:
                        try:
                            dis_numeric.append([rat[2][idx] for idx, val in enumerate(discalc) if val == d][0])
                        except IndexError:
                            pass
                
      
        for nd in notdistracted:
                        try:
                             ndis_numeric.append([rat[2][idx] for idx, val in enumerate(discalc) if val == nd][0])   
                        except IndexError:
                            pass
                        
        # Makes the distracted trial types into integers 
        dis_numeric = [int(d) for d in dis_numeric]
        # Counts to work out percentages after finding how many are each modality 
        d_whitenoise_count = 0
        d_tone_count = 0
        d_combined_count = 0 
        
        dis_type_text = [] #labels the distypes with text labels and adds to the counts
        for d in dis_numeric:
            if d in modalitykey['whitenoise']:
                dis_type_text.append('whitenoise')
                d_whitenoise_count += 1
            elif d in modalitykey['tone']:
                dis_type_text.append('tone')
                d_tone_count += 1
            elif d in modalitykey['combined3']:
                dis_type_text.append('combined3')
                d_combined_count += 1 
   # Non-distracted trials by modality 
        ndis_numeric = [int(d) for d in ndis_numeric]
        nd_whitenoise_count = 0
        nd_tone_count = 0
        nd_combined_count = 0 
        
        ndis_type_text = []
        for d in ndis_numeric:
            if d in modalitykey['whitenoise']:
                ndis_type_text.append('whitenoise')
                nd_whitenoise_count += 1
            elif d in modalitykey['tone']:
                ndis_type_text.append('tone')
                nd_tone_count += 1
            elif d in modalitykey['combined3']:
                ndis_type_text.append('combined3')
                nd_combined_count += 1 
 
## DEBUGGING SECTION        

# Not an issue if there are zero distracted trials, because calculates 0%
# without problem, will be an issue if both 0 dis and 0 non-dis (because no 
# trials of that modality, need to exclude)
         
#        if d_whitenoise_count < 1:
#            print('y', d_whitenoise_count, nd_whitenoise_count)
#        if d_tone_count < 1:
#            print('x', d_tone_count, nd_tone_count)
#        if d_combined_count < 1:
#            print('z', d_combined_count, nd_combined_count)
 
# WHITENOISE        
        if nd_whitenoise_count > 0:            
            percent_distracted_whitenoise = d_whitenoise_count / (d_whitenoise_count + nd_whitenoise_count) *100
        elif d_whitenoise_count < 1:
            try:
                percent_distracted_whitenoise = d_whitenoise_count / (d_whitenoise_count + nd_whitenoise_count) *100
            except ZeroDivisionError:
                print('oops')
                percent_distracted_whitenoise = -1
        elif d_whitenoise_count > 1:
            percent_distracted_whitenoise = 100 
# TONE            
        ## what to do if both zero 
        if nd_tone_count > 0:
            percent_distracted_tone = d_tone_count / (d_tone_count + nd_tone_count) *100
        elif d_tone_count < 1:
            try:
                percent_distracted_tone = d_tone_count / (d_tone_count + nd_tone_count) *100
            except ZeroDivisionError:
                print('oopsT')
                percent_distracted_tone = -1
        elif d_tone_count > 1:
            percent_distracted_tone = 100 

# COMBINED        
        if nd_combined_count > 0:
            percent_distracted_combined = d_combined_count / (d_combined_count + nd_combined_count) *100  

        elif d_combined_count < 1:
            try:
                percent_distracted_combined = d_combined_count / (d_combined_count + nd_combined_count) *100  
            except ZeroDivisionError:
                print('oopsC')
                percent_distracted_combined = -1
        elif d_combined_count > 1:
            percent_distracted_combined = 100 

# NON WHITE NOISE

        if nd_combined_count > 0 or nd_tone_count > 0:
            percent_distracted_all_non_WN = (d_tone_count + d_combined_count) / ((d_tone_count + d_combined_count )+ (nd_tone_count + nd_combined_count)) *100  
                
        elif d_combined_count < 1 and d_tone_count < 1:
            print('all less')
            try:
                percent_distracted_all_non_WN = (d_tone_count + d_combined_count) / ((d_tone_count + d_combined_count )+ (nd_tone_count + nd_combined_count)) *100  
            except ZeroDivisionError:
                print('oopsNWN')
                percent_distracted_all_non_WN = -1     
        elif d_combined_count or d_tone_count > 1:
            percent_distracted_all_non_WN = 100
     
      
        percent_dis_whitenoise_group.append(percent_distracted_whitenoise)
        percent_dis_tone_group.append(percent_distracted_tone)
        percent_dis_combined_group.append(percent_distracted_combined)
        percent_dis_all_non_WN_group.append(percent_distracted_all_non_WN)
    
## Removing the occurrences of -1 which were added to account for zero distractors
    
    percent_dis_whitenoise_group = [x for x in percent_dis_whitenoise_group if x != -1]
    percent_dis_tone_group = [x for x in percent_dis_tone_group if x != -1]
    percent_dis_combined_group = [x for x in percent_dis_combined_group if x != -1]
    percent_dis_all_non_WN_group = [x for x in percent_dis_all_non_WN_group if x != -1]       

        
    mean_percent_WHITENOISE = np.mean(percent_dis_whitenoise_group) # the average percentage of JUST whitenoise trials that rats are distracted on 
    mean_percent_TONE = np.mean(percent_dis_tone_group)
    mean_percent_COMBINED = np.mean(percent_dis_combined_group)
    mean_percent_ALL_NON_WN = np.mean(percent_dis_all_non_WN_group)
    
    return discalcgroup, percent_dis_whitenoise_group, percent_dis_tone_group, \
            percent_dis_combined_group, percent_dis_all_non_WN_group, mean_percent_WHITENOISE, mean_percent_TONE, \
            mean_percent_COMBINED, mean_percent_ALL_NON_WN



################################################################################################
################################################################################################

# Distraction day analysis (including modalities)
modalitykey = {'whitenoise':[1,4], 'tone':[2,5], 'combined3':[3,6]}


discalcgroup, percent_dis_whitenoise_group, percent_dis_tone_group, \
percent_dis_combined_group, percent_dis_all_non_WN_group, mean_percent_WHITENOISE, \
mean_percent_TONE, mean_percent_COMBINED, mean_percent_ALL_NON_WN = discalc_modalities(distraction, modalitykey)


## Report that they are no more distracting in these animals so no need to do the photometry!!!
## Maybe report this early with the behavioural stuff? 