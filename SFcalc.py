# Emma Alexander MPhys Project
# Jan 2017
# Calculates structure function values for a map

from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt
import pylab as pl
import math as mth
import sys

# functions
# function to calculate distance between two points
def dist(x1,x2,y1,y2):
    xdif = mth.pow(x2-x1,2)
    ydif = mth.pow(y2-y1,2)
    distance = mth.pow(xdif+ydif,0.5)
    return float(distance)

# function to return angle between two points and x axis
# polar coordinates convention
def ang(x1,x2,y1,y2):
    xdif = x2-x1
    ydif = y2-y1
    angle = np.arctan2(ydif,xdif) 
    if (angle >= mth.pi/2.) and (angle <= mth.pi):
        correction = -mth.pi/2.
    else:
        correction = 3.*mth.pi/2.
    return angle + correction


# function to return midpoints
def midpnt(x1,x2,y1,y2):
	x = int((x1+x2)/2.)
	y = int((y1+y2)/2.)
	return x,y

def printhelp():
    # print help message with options
    print "\nCode to calculate the values of the structure function"
    print "Usage: python SFcalc.py options infile outfile \n"
    print "Options:"
    print "-h   Display these options"
    print "-b   Specify beam size in arcseconds"
    print "     Usage: -b number, where number is the beam width in arcseconds"
    print "     Default: 18  arcsec"
    print "-s   Specify pixel scale in arcseconds"
    print "     Usage: -b number, where number width of pixels in arcseconds"
    print "     Default: 3  arcsec"
    print "-c   Crop the map into a square"
    print "     Usage: -c number, where number is the size lenth of the square in pixels"
    print "     Default: 220 pixels"
    print "-m   Specify maximum distance scale in beam widths"
    print "     Usage: -m number, where number is the number of beam widths"
    print "     Default: 10"
    print "-o   Split into radial and tangential componets" 
    print "-i   Calculate intermediate values at half beam widths in addition to full beam separations"
    print "-r   Randomise the map for test purposes"   
    print "-p   Show plot of structure function" 
    print "     Usage: -p figname, where figname is where figure wil be saved"

#---------------------get input parameters------------------------#

#parse arguments 
args = sys.argv[1:]
#assert len(args) >= 3
# EDIT JULY 2018: WHAT WAS I THINKING WRITING THE OPTIONS LIKE THIS?????
if (len(args) == 1 and args[0] == '-h'):
    # print help message with options
    printhelp()
    sys.exit(0)

else:
    options = args[0:-2]
    infile = str(args[-2])
    outfile = str(args[-1])

if any('-h' in s for s in options):
    printhelp()
if any('-b' in s for s in options):
    # beam width specified 
    beamsize =  float(options[options.index('-b')+1])
else:
    beamsize = 18.
if any('-s' in s for s in options):
    # cell scale size in arcsecs
    cellsize =  float(options[options.index('-s')+1])
else:
    cellsize = 3.
if any('-c' in s for s in options):
    # want to crop the image into a square
    cropparam =  int(options[options.index('-c')+1])
else:
    cropparam = 220
# redefine to match code below
cropparam = cropparam/2
if any('-m' in s for s in options):
    # max points specified
    maxpoints =  int(options[options.index('-m')+1])
else:
    maxpoints = 10
if any('-r' in s for s in options):
    #randomise
    print "Will Randomise map"
    randomise = True
else: randomise = False
if any('-o' in s for s in options):
    #randomise
    print "Will spilt into radial and tangential components"
    split = True 
else: split = False
if any('-i' in s for s in options):
    # want half beam values
    print "Will calculate values at half the beam width" 
    halfbeam = True
else: halfbeam = False
if any('-p' in s for s in options):
    print "Will plot figure"
    figname = options[options.index('-p')+1]
    plot = True
else: plot = False

# values calculated from parameters
# size of the beam in pixels
beampix = beamsize/cellsize

#---------------------------begin-------------------------------#
print infile
# deal with files differently depending on input type 
fitpos = ['.FITS','.fits']
txtpos = ['.txt','.TXT','.dat','.DAT']

if any(x in infile for x in fitpos):
    # open file of input values
    hdu = fits.open(infile)
    header = hdu[0].header
    data = hdu[0].data
    # close HDU
    hdu.close()

    # get dimensions of input image in pixels
    xlength = data.shape[1]
    ylength = data.shape[0]
    xcentre = round(xlength/2.)
    ycentre = round(ylength/2.) -1

    # crop the data to be centred and square
    datacropped = data[int(ycentre-cropparam):int(ycentre+cropparam),int(xcentre-cropparam):int(xcentre+cropparam)]
    datacropped = np.ma.masked_invalid(datacropped)
    cropxlength = datacropped.shape[1]
    cropylength = datacropped.shape[0]

elif any(x in infile for x in txtpos):
    vals = np.loadtxt(infile)
    datacropped = vals
    cropxlength = vals.shape[0]
    cropylength = vals.shape[1] 
else:
    print 'Could not determine input file type (please use .fits, .dat or .txt)'

plt.imshow(datacropped, origin='lower')
plt.show()

# test thing
if (randomise == True):
    datashape = datacropped.shape
    print datashape
    datacropped = np.reshape(datacropped,-1)
    print datacropped.shape
    np.random.shuffle(datacropped)
    print datacropped.shape
    datacroppedr = np.reshape(datacropped,datashape)
    datacropped = datacroppedr

# array of distances from a central pixel, in pixels
fulldistarray = np.empty([2*cropxlength+1,2*cropylength+1])
#fullangarray = np.empty([2*cropxlength+1,2*cropylength+1])

for x in range(0,2*cropxlength+1):
    for y in range(0,2*cropylength+1):
        fulldistarray[x,y] = dist(cropxlength,x,cropylength,y)
        #fullangarray[x,y] = ang(cropxlength,x,cropylength,y)

# get index values for each point
index = np.empty(maxpoints)
if halfbeam == True:
    for i in range(0,maxpoints):
        index[i] = float((i+1)/2.)
        binwidth = beampix/4.
else:
    for i in range(0,maxpoints):
        index[i] = float((i+1))
        binwidth = beampix/2.

print binwidth

S = np.empty(maxpoints)
SR=np.empty(maxpoints)
ST=np.empty(maxpoints)
Sother=np.empty(maxpoints)
dtheta = np.empty(maxpoints)

# file to write numerical values out to
f = open(outfile, 'w')
header1 = "# Structure Function values for "+str(infile)+"\n"
f.write(header1)

if split == True: 
    # splitting into radial and tangential
    header2 = "# Columns are Theta [degrees], SF vals [rad^2 m^-4] (Radial, Tangential, Mix)\n"
else:
    header2 = "# Columns are Theta [degrees], SF vals [rad^2 m^-4] \n"
f.write(header2)

anglist = []

if split == True:
    for i in range(0,maxpoints):
        print index[i]
        SarrT = np.empty(datacropped.shape)
        SarrR = np.empty(datacropped.shape)
        Sarrother = np.empty(datacropped.shape)

        tanweight = np.empty(datacropped.shape)
        radweight = np.empty(datacropped.shape)
        othweight = np.empty(datacropped.shape)

        beampixlower = index[i]*beampix - binwidth
        beampixhigher = index[i]*beampix + binwidth
        dtheta[i] = (index[i]*beamsize)/3600. #dtheta in degrees
        for xval in range(0,cropxlength):
            for yval in range(0,cropylength):
                if datacropped[xval,yval] != 0:
                    distarray = fulldistarray[cropxlength-xval:2*cropxlength-xval,cropylength-yval:2*cropylength-yval]
                    # find points the requred distance away
                    distlocs1 = np.where(distarray < beampixhigher,0,1)
                    distlocs2 = np.where(distarray > beampixlower,0,2)

                    distflags = np.where(distlocs1 == distlocs2,1,0)
 
                    wanteddata = np.multiply(datacropped,distflags)
                    wanteddata = np.ma.masked_where(wanteddata==0,wanteddata)
                    # split points into radial, tangential or neither
                    if np.ma.MaskedArray.count(wanteddata) != 0:
                        tandata = np.empty(wanteddata.shape)
                        raddata = np.empty(wanteddata.shape)
                        othdata = np.empty(wanteddata.shape)
                        indices = np.ma.MaskedArray.nonzero(wanteddata)
                        for j in range(0,len(indices[0])):
                            x = indices[0][j] 
                            y = indices[1][j]
                            xmid,ymid = midpnt(x,xval,y,yval)
                            midang = ang(xmid,xval,ymid,yval)
                            organg = ang(xmid,cropxlength/2,ymid,cropylength/2)
                            angdif = organg - midang
                            if angdif < 0:
                                angdif = angdif + 2*mth.pi
                            if angdif > 2*mth.pi:
                                angdif = angdif - 2*mth.pi
                            if (angdif < mth.pi/8.) or (angdif > 15*mth.pi/8. and angdif < 17* mth.pi/8.) or (angdif > 
7*mth.pi/8. and angdif < 9*mth.pi/8.):
                                 # radial
                                 raddata[x,y] = wanteddata[x,y]
                                 tandata[x,y] = 0
                                 othdata[x,y] = 0
                            elif (angdif > 3*mth.pi/8. and angdif < 5*mth.pi/8.) or (angdif > 11*mth.pi/8. and angdif <
 13* mth.pi/8.):
                                 #tangential
                                 raddata[x,y] = 0
                                 tandata[x,y] = wanteddata[x,y]
                                 othdata[x,y] = 0
                            else:
                                 #other 
                                 tandata[x,y] = 0
                                 raddata[x,y] = 0
                                 othdata[x,y] = wanteddata[x,y]

                    else: 
                        tandata = np.zeros(wanteddata.shape)
                        raddata = np.zeros(wanteddata.shape)
                        othdata = np.zeros(wanteddata.shape) 

                    # count the number of points, so mean can be weighted
                    tanweight[xval,yval] = np.count_nonzero(tandata)
                    radweight[xval,yval] = np.count_nonzero(raddata)
                    othweight[xval,yval] = np.count_nonzero(othdata)

                    tandata = np.ma.masked_where(tandata==0,tandata)
                    raddata = np.ma.masked_where(raddata==0,raddata)
                    othdata = np.ma.masked_where(othdata==0,othdata)

                    SarrT[xval,yval] = np.nanmean(np.square(tandata - datacropped[xval,yval]))
                    SarrR[xval,yval] = np.nanmean(np.square(raddata - datacropped[xval,yval]))
                    Sarrother[xval,yval] = np.nanmean(np.square(othdata - datacropped[xval,yval]))

        SarrR = np.ma.masked_invalid(SarrR)
        SarrT = np.ma.masked_invalid(SarrT)
        Sarrother = np.ma.masked_invalid(Sarrother)
    
        SR[i] = np.ma.average(SarrR,weights=radweight)
        ST[i] = np.ma.average(SarrT,weights=tanweight)
        Sother[i] = np.ma.average(Sarrother,weights=othweight)
        print("{} {} {}".format(SR[i],ST[i],Sother[i]))

        s = str(dtheta[i]) + ' ' + str(SR[i]) + ' ' + str(ST[i]) + ' ' + str(Sother[i]) + '\n'
        f.write(s)
else:
    for i in range(0,maxpoints):
        print index[i]
        Sarr = np.empty(datacropped.shape)
        weight = np.empty(datacropped.shape)
        beampixlower = index[i]*beampix - binwidth
        beampixhigher = index[i]*beampix + binwidth
        dtheta[i] = (index[i]*beamsize)/3600. #dtheta in degrees

        for xval in range(0,cropxlength):
            for yval in range(0,cropylength):
                distarray = fulldistarray[cropxlength-xval:2*cropxlength-xval,cropylength-yval:2*cropylength-yval]
                datalocs1 = np.where(distarray < beampixhigher,0,1)
                datalocs2 = np.where(distarray > beampixlower,0,2)
                flags = np.where(datalocs1 == datalocs2,1,0)
                wanteddata = np.multiply(datacropped,flags)    
                wanteddata = np.ma.masked_where(wanteddata==0,wanteddata)
                Sarr[xval,yval] = np.nanmean(np.square(wanteddata - datacropped[xval,yval]))
                weight[xval,yval] = np.count_nonzero(Sarr[xval,yval])
        Sarr = np.ma.masked_invalid(Sarr)
        S[i] = np.ma.average(Sarr,weights=weight)
        print("{} {}".format(dtheta[i],S[i]))
    s = str(dtheta[i]) + ' ' + str(S[i]) + '\n'
    f.write(s)

f.close()

# copied code from seperate plot function
if (plot == True):
    distance = 9300. # distance to Galaxy in kpc
    beamsize = 0.005 # beam size in degrees
    data = np.loadtxt(outfile) # text file of data
    rows = data.shape[0]
    columns = data.shape[1]
    # take log of angular size and corresponding data values
    logdata = np.log10(data) 

    m = ['x','+','s'] # marker style
    ms = [20,30,6] # marker size
    mc= ['r','b','k'] #marker colour

    # start plotting 
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for i in range (1,columns):    
        j = i-1
        ax1.scatter(logdata[:,0],logdata[:,j],marker=m[j],c=mc[j],s=ms[j])

    ax1.set_xlabel(r"$log(\delta \theta [deg])$")
    ax1.set_ylabel(r"$log(SF_{RM} [rad^2 m^{-4}])$")
    plt.figtext(0,0,filename) 
    lims = ax1.get_xlim()

    distticks = [0.2,0.5,1,2,5,10,20]
    minor_ticks = [0.3,0.4,0.6,0.7,0.8,0.9,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.]
    scaleddistticks = np.log10(np.rad2deg(np.arctan(np.divide(distticks,distance))))
    scaled_minor = np.log10(np.rad2deg(np.arctan(np.divide(minor_ticks,distance))))

    ax2 = ax1.twiny()
    ax2.set_xlabel("Distance scale at galaxy [kpc]")
    ax2.set_xticks(scaleddistticks)                                                     
    ax2.set_xticks(scaled_minor, minor=True) 
    ax2.set_xticklabels(distticks)                                                      
    ax2.set_xbound(ax1.get_xbound())
    red_patch = pat.Patch(color='red', label=r'Radial $\delta \theta$')
    blue_patch = pat.Patch(color='blue', label=r'Tangential $\delta \theta$')
    green_patch = pat.Patch(color='green', label = 'Randomised radial')
    black_patch = pat.Patch(color='black', label = 'Neither/ Mix')
    plt.legend(handles=[red_patch,blue_patch,black_patch],loc='upper left',prop={'size':11})
    plt.savefig(figname,dpi=500)
