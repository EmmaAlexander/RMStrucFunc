# Emma Alexander MPhys Project
# Sept 2016
# RM map -> plot

from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import pylab as pl
import math as mth

datafile = 'randerrorSF2.txt'
datafile2 = 'rightarm.txt'
datafile3 = 'southarm.txt'
#rndfile = 'RTOrandom.txt'
RTfile = 'struct_map2_new.txt'
#halfbm = 'RTOhalfbeam.txt'
figname = 'strucmap2pix.png'

distance = 9300. # distance to Galaxy in kpc

data = np.loadtxt(datafile) # text file of data
# take log of angular size and corresponding data values
logdata = np.log10(data) 

data2 = np.loadtxt(datafile2) # text file of data
# take log of angular size and corresponding data values
logdata2 = np.log10(data2) 

data3 = np.loadtxt(datafile3) # text file of data
# take log of angular size and corresponding data values
logdata3 = np.log10(data3) 

#randdata = np.loadtxt(rndfile)
#randlogdata = np.log10(randdata)

RTdata = np.loadtxt(RTfile)
RTlogdata = np.log10(RTdata)
RTpixels = np.multiply(RTdata[:,0],1200)
RTpixelslog = np.log10(RTpixels)

#hlfbmdat = np.loadtxt(halfbm)
#hlfbmlogdat = np.log10(hlfbmdat)

# start plotting 
fig = plt.figure()
ax1 = fig.add_subplot(111)
#ax1.scatter(logdata[:,0],logdata[:,1], c='k',marker='x')
#ax1.scatter(randlogdata[:,0],randlogdata[:,1],c='g',marker='o')
#ax1.scatter(randlogdata[:,0],randlogdata[:,2],c='k',marker='o')
plt1=ax1.scatter(RTpixelslog,RTlogdata[:,1],c='r',marker='x',s=20,label='Radial')
plt2=ax1.scatter(RTpixelslog,RTlogdata[:,2],c='b',marker='+',s=30,label='Azimuthal')
plt3=ax1.scatter(RTpixelslog,RTlogdata[:,3],c='k',marker='s',s=6,label='Neither/mix')
#ax1.scatter(hlfbmlogdat[:,0],hlfbmlogdat[:,1],c='r',marker='+')
#ax1.scatter(hlfbmlogdat[:,0],hlfbmlogdat[:,2],c='b',marker='+')
#ax1.scatter(hlfbmlogdat[:,0],hlfbmlogdat[:,3],c='k',marker='+')
#ax1.scatter(randlogdata[:,0],randlogdata[:,1],c='r',marker='+')
#ax1.scatter(randlogdata[:,0],randlogdata[:,2],c='b',marker='+')
#ax1.scatter(randlogdata[:,0],randlogdata[:,3],c='k',marker='+')
#plt1=ax1.scatter(logdata[:,0],logdata[:,1],c='r',marker='x',s=20,label='Left arm')
#plt2=ax1.scatter(logdata2[:,0],logdata2[:,1],c='k',marker='s',s=6,label='Right arm')
#plt3=ax1.scatter(logdata3[:,0],logdata3[:,1],c='b',marker='+',s=30,label='South arm')
ax1.set_xlabel(r"$log_{10}(pixels)$")
ax1.set_ylabel(r"$log_{10}(SF value)$")
lims = ax1.get_xlim()
#plt.figtext(0,0,datafile)
distticks = [0.2,0.5,1,2,5,10,20]
minor_ticks = [0.3,0.4,0.6,0.7,0.8,0.9,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.]
scaleddistticks = np.log10(np.rad2deg(np.arctan(np.divide(distticks,distance))))
scaled_minor = np.log10(np.rad2deg(np.arctan(np.divide(minor_ticks,distance))))
#parkerline = np.log10(np.rad2deg(np.arctan(np.divide(3.7,distance))))
#plt.plot((parkerline, parkerline), (ax1.get_ybound()), 'k--',color = '0.75')
#plt1=ax1.scatter(logdata[:,0],logdata[:,1],c='r',marker='x',s=20)
#plt2=ax1.scatter(logdata2[:,0],logdata2[:,1],c='k',marker='s',s=6)
#plt3=ax1.scatter(logdata3[:,0],logdata3[:,1],c='b',marker='+',s=30)

#ax1.set_ybound(2.65,3.0)

#ax2 = ax1.twiny()
#ax2.set_xlabel("Distance scale at galaxy [kpc]")
#ax2.set_xticks(scaleddistticks)                                                     
#ax2.set_xticks(scaled_minor, minor=True) 
#ax2.set_xticklabels(distticks)                                                      
#ax2.set_xbound(ax1.get_xbound())
#red_patch = pat.Patch(color='red', label=r'Radial $\delta \theta$')
#blue_patch = pat.Patch(color='blue', label=r'Azimuthal $\delta \theta$')
#green_patch = pat.Patch(color='green', label = 'Randomised radial')
#black_patch = pat.Patch(color='black', label = 'Neither/ Mix')
#red_patch = pat.Patch(color='red', label=r'Left arm')
#blue_patch = pat.Patch(color='blue', label=r'Right arm')
#black_patch = pat.Patch(color='black', label = r'South arm')
#plt.legend(handles=[red_patch,blue_patch,black_patch],loc='upper left',prop={'size':11})

handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels,prop={'size':11},loc='upper left')
plt.savefig(figname,dpi=500)
plt.show()

