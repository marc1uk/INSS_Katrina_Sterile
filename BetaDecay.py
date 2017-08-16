# vim:set noexpandtab tabstop=4 wrap
import numpy as np
#import pylab as P
import math 
import matplotlib.pyplot
import os, sys, csv
import ROOT
from array import array
# could replace with from ROOT import *
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TF1
from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double

# A short program to estimate the sensitivity of the Katrin experiment to sterile neutrinos. 
# The program is intended to produce contours in the sterile mass - sterile mixing angle plane
# that encompass regions that could be excluded at a 3 sigma level, when taking into account
# the expected errors expected in Katrin.

# The analysis flow is as follows:
# 1. Implement a function representing the Tritium beta decay spectrum, including the effects
#    of a single sterile of a given mass and mixing angle.
#    The no-sterile spectrum - null hypothesis - is this spectrum evaluated with a mixing angle of 0.
#    This will be used as the underlying pdf for generating simulated data, and for calculating
#    expected rates for the calculation of chi^2 of a fit to the null hypothesis.
# 2. Implement a function that will take a given mass and mixing angle and will generate a simulated
#    spectrum, representing the expected counts obtained over a nominal 3-year measurement period.
# 3. Implement a function that will take a simulated spectrum and will evaluate the chi^2 between 
#    the data and the null hypothesis.
# 4. Loop over a vector of prospective neutrino masses. For each mass:
#    +  Loop over a vector of prospective mixing angles. For each mixing angle:
#       - Generate a set of datapoints representing a measured spectrum for that mass and mixing angle.
#       - Calculate the chi^2 of the null hypothesis fit to the data. 
#       - Store the result in a vector for this mass.
#    +  Fit a curve to the span of chi^2 vs mixing angle values, for this mass
#    +  Find the mixing angles at which chi^2 = 3 sigma.
#    +  Add the pair of values to a vector - one entry per mass, each entry the max/min mixing angles.
# 5. Implement a function that takes the results - a vector of pairs - and produces a contour plot
#    Each pair corresponds to the range of mixing angles that are consistent with the null hypothesis
#    at that mass. 

# ===========================================================================================
# Global constants
# ===========================================================================================
# Use eV as as universal unit(?)
Ns = 1.47e-13						# Rate normalisation, value for Katrin. [s^-1 ev^-5]
me = 511000							# Mass of electron [eV]
kb = 8.76e-5						# Boltzmann constant [eV K^-1]
a0 = 2.68e-4						# Bohr radius [fm]
m1 = 0.								# lightest neutrino mass [eV]
deltamsq21 = 7.53e-5				# solar mass splitting 
deltamsq32 = 2.45e-3				# atmospheric mass splitting
m2 = math.sqrt(m1**2 + deltamsq21)
m3 = math.sqrt(m2**2 + deltamsq32)
ActiveMasses = [m1,m2,m3]
theta12 = np.arcsin(np.sqrt(0.297))	# from PDG 2016 mixing review
theta13 = np.arcsin(np.sqrt(0.216))	# 
theta23 = np.arcsin(np.sqrt(0.500))	# 
deltacp = 1.35*np.pi				# 
Q = 18571.8							# Q value for Tritium [eV]
numkebins = 100						# number of bins for simulated data
nummixingangs = 30					# number of steps in mixing angle range scan
nummasses = 20						# number of steps in mass range scan
kemin = 1.							# min of the electron energy spectrum [eV]
kemax = Q							# max of the electron energy spectrum [eV]
kelist = []
mixingangmin = 6					# we'll scan a range of mixings such that sin^2(theta_s)
mixingangmax = 10					# spans the range 10^-(mixingangmax) -> 10^-(mixingangmin)
massmin = 1000.						# minimum sterile mass (eV)
massmax = 20000.					# maximum sterile mass (eV)
#alpha = 2.7x10^-3 = 0.0027 for 3 sigma -> chi^2 ~8
chi2for3sigma = 8					# what chi2 value corresponds to a 3 sigma deviation from the null hypothesis
integrationtime = 3*360*24*3600 	# 3 years of running, assuming spectrum is defined in counts/s FIXME if necessary

excitationfile = 'HeTtable.csv'
with open(excitationfile,'rb') as file:
	reader = csv.reader(file,quoting=csv.QUOTE_NONNUMERIC)
	daughtertable = list(reader)
# ===========================================================================================
# 1. Function for representing the tritium beta decay spectrum, including a potential sterile
# ===========================================================================================
# Support functions
# -------------------------------
def Energy(mass,Ke):			# Energy of electron
	return mass + Ke
# -------------------------------
def Momentum(mass,Ke):			# Momentum of electrom
	return math.sqrt(Energy(mass,Ke)**2 - mass**2)
# -------------------------------
def Fermi(mass,Ke):					# The fermi function
	return (1/(1+np.exp((Ke-Ef(Ke))/((2./3.)*Ke))))
# -------------------------------
def Ef(Ke): 					# The fermi level of H-3
	eta = 1/(a0*Momentum(me, Ke))
	return 4*np.pi*eta/(1-np.exp(-4*np.pi*eta))
# -------------------------------
def GetPMNS():
	Ue1   =  np.cos(theta12)*np.cos(theta13)
	Ue2   =  np.sin(theta12)*np.cos(theta13)
	Ue3   =  np.sin(theta13)*cmath.exp(-1j*deltacp)
	Umu1  = -np.sin(theta12)*np.cos(theta23) - np.cos(theta12)*np.sin(theta23)*np.sin(theta13)*exp(-1j*deltacp)
	Umu2  =  np.cos(theta12)*np.cos(theta23) - np.sin(theta12)*np.sin(theta23)*np.sin(theta13)*exp(-1j*deltacp)
	Umu3  =  np.sin(theta23)*np.cos(theta13)
	Utau1 =  np.sin(theta12)*np.sin(theta23) - np.cos(theta12)*np.cos(theta23)*np.sin(theta13)*exp(-1j*deltacp)
	Utau2 = -np.cos(theta12)*np.sin(theta23) - np.sin(theta12)*np.cos(theta23)*np.sin(theta13)*exp(-1j*deltacp)
	Utau3 =  np.cos(theta23)*np.cos(theta13)
	
	Upmns= [[Ue1   ,  Ue2   ,  Ue3   ],
			[Umu1  ,  Umu2  ,  Umu3  ],
			[Utau1 ,  Utau2 ,  Utau3 ]]
	
	return Upmns

def beta(Ke,Sinsq14=0.,m4=-1.):
	Mixing = [math.cos(theta12)**2 * math.cos(theta13)**2 * (1.-Sinsq14), math.cos(theta13)**2 * math.sin(theta12)**2 * (1.-Sinsq14), math.sin(theta13)**2 * (1.-Sinsq14)]
	Masses = []
	if Sinsq14 >=0.:
		Masses = ActiveMasses + [m4]
		Mixing.append(Sinsq14)
	
	rate = Ns*Fermi(me,Ke)*Energy(me,Ke)*Momentum(me,Ke)
	
	asum = 0
	
	for excitpair in daughtertable:
		
		for (mix,mass) in zip(Mixing,Masses):
			if ((Q-excitpair[0]-Ke) - mass) >= 0. and ((Q-excitpair[0]-Ke)**2 - mass**2) > 0:
				temp = excitpair[1]*(Q-excitpair[0]-Ke) * mix * math.sqrt((Q-excitpair[0]-Ke)**2 - mass**2)
				asum += temp
	
	rate *= asum
	
	return rate

#def BetaHist():  # Generate the beta histogram
#    hist = []    # Declare the content of the histogram
#    binWidth = (kemin-kemax)/numbins
#    for ki in range(kemin, kemax, binWidth):
#        leftBinEdge = beta(ki, Q, 0, 0)  # FIXME assuming Q is a variable and null scenario.
#        rightBinEdge = beta(ki+binWidth, Q, 0, 0)
#        binContent = binWidth*(leftBinEdge+rightBinEdge)/2
#        binContent *= integrationtime    # Unnecesarry if bin content is event rate.
#        hist.append(binContent);
#    return hist
          
# -------------------------------
# To create a TF1 object with a custom function, the function must take it's arguments via arrays
# We could do away with this function by suitably defining beta - this is just a wrapper.
def BetaSpectrum(variables,parameters):
	Ke = variables[0]
	sin2thetae4 = parameters[0]
	sterilemass = parameters[1]
	return beta(Ke, sin2thetae4, sterilemass)
# -------------------------------
# A TF1 object based on the spectrum is needed to generate the expected counts for simulated data
npars = 2						# Num parameters required for our custom function
betatf1 = TF1('betaspectrum', BetaSpectrum, kemin, kemax, npars) # this is the rate per unit time per unit energy at a given energy

# ===========================================================================================
# 2. Function for generating fake data spectrum for given mass and mixing angle
# ===========================================================================================
def DataSpectrum(mass,mixing):	# simulated count spectrum for a given mass and mixing angle

	spectrum=[]
#	# method 1: generate a suitable number of 'hits' point by point and bin into a histogram
#	spectrumhist = TH1F('spectrum', 'Simulated Data', numkebins, kemin, kemax)
#	betatf1.SetParameters(mass,mixing)
#	for hit in range(0,totalcounts):
#		KeOfThisHit = betatf1.GetRandom()
#		spectrumhist.Fill(KeOfThisHit)
#	spectrumbuffer = spectrumhist.GetArray()				# this returns a pointer to a ROOT buffer
#	for bini in range(0,numkebins):							# convert histogram to a normal list
#		spectrum.append(spectrumbuffer[bini])
	
#	# method 2: generate the number of entries in each bin directly, taking into account statistics
#	kebinwidth = (kemax-kemin)/numkebins
#	for i in range(0,numkebins-1):							# spectrum runs from 0.1 to 18.574keV over 100 bins.
#		Ke = kemin + (i/numkebins)*(kemax-kemin)			# Ke of this bin centre in eV
#		expectedbincount = BetaSpectrum(Ke)*kebinwidth*integrationtime
#		bincontents = np.random.poisson(expectedbincount)	# pull from a poisson with mean @ expected. Is binoimial preferable?
#		spectrum.append(bincontents)
	
	# method 3: statistical fluctuations in simulated data are not incorporated (...)
	# so 'data' is exactly the expected number of counts for a given mass, mixing
	kebinwidth = (kemax-kemin)/numkebins
	for i in range(0,numkebins+1):							# spectrum runs from 0.1 to 18.574keV over 100 bins.
		Ke = kemin + (1.*i)*kebinwidth							# Ke of this bin centre in eV
		if (len(kelist)!=(numkebins+1)):
			kelist.append(Ke)
		parameters = [mixing,mass]
		expectedbincount = BetaSpectrum([Ke],parameters)*kebinwidth*integrationtime
		spectrum.append(expectedbincount)					# same as above just skip the Poisson step

	return spectrum											# return simulated data

# ===========================================================================================
# 3. Function for calculating the chi^2 between a given dataset and the no-sterile hypothesis
# ===========================================================================================
def Chi2Test(observe, expect):								# chi2 test for each bin
     sigma = np.sqrt(expect)
     value=0
     if (sigma):
         value = (observe**2 - expect**2)/sigma**2
     return value

def Chi2TestPoisson(observe, expect):						# chi2 test for each bin using the Poisson formula
	if (observe == 0) or (expect == 0):
		return 0
	chi2forthispoint = 2.* (expect - observe + observe*np.log(observe/expect))
#	if(chi2forthispoint<0):
#		print "\nOH NO, CHI2<0! expect=" + str(expect) + ", observe=" + str(observe),
#		print ", observe/expect=" + str(observe/expect) + "log(o/e)=" + str(np.log(observe/expect)),
#		print ", observe(1+log(o/e))=" + str(observe*(1+np.log(observe/expect))) + ", chi2 =" + str(chi2forthispoint) + "\n"
	return chi2forthispoint

def Chi2FitToNull(nullSpectrum, modelSpectrum):				# Sum the elemental chi2 value.
    assert len(nullSpectrum) == len(modelSpectrum), "The bin numbers of the spectra do not match."    
    summation = 0
    for i in range (0, len(modelSpectrum)):
        chi2 = Chi2Test(nullSpectrum[i], modelSpectrum[i])
        summation += chi2
    return summation

def Chi2FitToNullIntegral(nullSpectrum, modelSpectrum):			# Sum the elemental chi2 value. Counts based on retarding potl.
    assert len(nullSpectrum) == len(modelSpectrum), "The bin numbers of the spectra do not match."    
    summation = 0
    nullintegral=0
    modelintegral=0
    for i in range (0, len(modelSpectrum)):
        nullintegral+=nullSpectrum[numkebins-i]
        modelintegral+=modelSpectrum[numkebins-i]
        chi2 = Chi2Test(nullintegral, modelintegral)
        summation += chi2
    return summation
    
# ===========================================================================================
# 4. Function for calculating the matrix of chi^2 over a range of masses and mixing angles
# ===========================================================================================
def CalculateChi2Matrix():
	chi2array = []		# matrix of chi2. rows are masses, columns are mixing angles
	mixinganglelist = []
	masslist = []
	nullSpectrum = DataSpectrum(0, 0) # Generate a 3-mixing-only scenario
# 4. Loop over a vector of prospective neutrino masses.

#	matplotlib.pyplot.figure()
	for i in range(0,nummasses+1):
		mass = massmin + ((1.*i)/nummasses)*(massmax-massmin)
		chi2forthismass = []
#    +  Loop over a vector of prospective mixing angles.
		for j in range(0,nummixingangs+1):
			power = mixingangmin + ((j*1.)/nummixingangs)*(mixingangmax-mixingangmin)
			sin2mixang = 10 ** -power
			print "mass = " + str(mass) + ", sin2mixang = " + str(sin2mixang)
			if (j==0):
				masslist.append(mass)
			if (i==0):
				mixinganglelist.append(sin2mixang)
#       - Generate a set of datapoints representing a measured spectrum for that mass and mixing angle.
			dataset = DataSpectrum(mass, sin2mixang)
			
			# to take the difference we need to convert to np arrays
#			nullarray = np.array(nullSpectrum)
#			datasetarray = np.array(dataset)
#			matplotlib.pyplot.plot(nullarray,'k-',label='null')
#			matplotlib.pyplot.plot(datasetarray,color='b',linestyle='dashed',label='data')
#			matplotlib.pyplot.show()
#			matplotlib.pyplot.figure()
#			ratioarray = datasetarray/nullarray
#			ratioplot = matplotlib.pyplot.plot(kelist,ratioarray,label=str(sin2mixang))
#			title = "Sterile/Null for mass " + str(mass)# + " and sin2thetae4 " + str(sin2mixang)
#			matplotlib.pyplot.title(title)
#			matplotlib.pyplot.ylim([0.999999,1.0000001])
#			matplotlib.pyplot.show()
			
#       - Calculate the chi^2 of the null hypothesis fit to the data
			chi2 = Chi2FitToNullIntegral(nullSpectrum, dataset)
#       - Store the result in a vector for this mass.
			chi2forthismass.append(chi2)
		#matplotlib.pyplot.show() # will halt until closed
		# plot the chi2 curve, for checking
#		thismassplot = matplotlib.pyplot.plot(mixinganglelist,chi2forthismass,label=str(mass))
#		matplotlib.pyplot.title('Chi2 vs Mixing Angle for Mass ' + str(mass))
#		matplotlib.pyplot.xlabel('Sin^2 Mixing angle')
#		matplotlib.pyplot.ylabel('Chi2')
#		matplotlib.pyplot.show()
		# turns out to draw contours we can just use the curves directly, we don't need to find the values at 3 sigma
		chi2array.append(chi2forthismass)
#	matplotlib.pyplot.legend()
#	matplotlib.pyplot.show()
	return chi2array, masslist, mixinganglelist
# Obselete code, originally to fit a polynomial to the range of chi2's and find the mixing angles where chi2 corresponds to 3 sigma
#		# TGraph takes array.array datatypes, so we need to convert from python lists
#		mixinganglelistar = array('d',mixinganglelist)
#		chi2forthismassar = array('d',chi2forthismass)
#		chi2vsmasstgraph = TGraph(nummixingangs,mixinganglelistar,chi2forthismassar)
##    +  Fit a curve so we can lookup the mixing angle for a given chi2
#		fitfunc = TF1( 'fitfunc', 'pol2',  mixinganglelist[0],  mixinganglelist[-1]) # PRESUME quadratic
#		thefit = chi2vsmasstgraph.Fit(fitfunc)		# should check the status of the fit
##    +  Find the mixing angles at which chi^2 = 3 sigma.
#		# To do this, we can shift the parabola down such that the 3-sigma line is at 0, then find the roots
#		fitparameters = thefit.GetParameters()
#		fitparametersaslist = [fitparameters[2],fitparameters[1],fitparameters[0]-chi2for3sigma] # reverse order!
#		theroots = np.roots(fitparametersaslist)
##    +  Store the pair of values to a vector - one entry per mass, each entry the max/min mixing angles.
#		chi2array.append([theroots[0], theroots[1]])
# End obselete code -------

# ===========================================================================================
# 5. Function for producing the contour plots
# ===========================================================================================
#http://matplotlib.org/examples/pylab_examples/contour_demo.html
#https://stackoverflow.com/questions/25206580/in-python-can-matplotlibs-contour-function-return-the-points-of-a-particular-c?rq=1
# 5. Implement a function that takes the results and produces a contour plot. Results are stored in
#    a matrix, where each row stores the range of mixing angles that are consistent with
#    the null hypothesis at that mass
def MakeContourPlots(mixings, masses, chi2vals, contourlimits):
	matplotlib.pyplot.figure()
	contourplot = matplotlib.pyplot.contourf(mixings, masses, chi2vals, contourlimits) # contour gives lines, contourf fills
	matplotlib.pyplot.title('3 Sigma Exclusion Limits over 3 Years, Statistical Only')
	matplotlib.pyplot.xlabel('Sterile-Electron Mixing Angle (rads)')
	matplotlib.pyplot.ylabel('Sterile Mass (eV)')
	matplotlib.pyplot.xscale('log')
	matplotlib.pyplot.clabel(contourplot, inline=1, fontsize=10)
	matplotlib.pyplot.show()

# ===========================================================================================
# 6. The Main Function
# ===========================================================================================
def main():
	print "doing main"
#	nullSpectrum = DataSpectrum(0, 0); # Generate a 3-mixing-only scenario.
#	thismassplot = matplotlib.pyplot.plot(nullSpectrum) # visual check of spectrum
#	matplotlib.pyplot.show()
	chi2array, masslist, mixinganglelist = CalculateChi2Matrix()
	chi2array2 = np.array(chi2array)
	chi2arraytranspose = chi2array2.transpose()
	chi2limits = [2.30,4.61,5.99,6.18,9.21]
	# 1 - 5 sigma for 2 DOF
	MakeContourPlots(mixinganglelist, masslist, chi2array, 5 )
#	# wait for input to keep the plot alive - not necessary, it waits for you anyway.
#	if __name__ == '__main__':
#		rep = ''
#		while not rep in [ 'q', 'Q' ]:
#			rep = raw_input( 'enter "q" to quit: ' )
#			if 1 < len(rep):
#				rep = rep[0]

main()

# a function to generate a plot of the effect on the beta spectrum, using artificially large mixing to enhance visibility
def effectPlot():
	matplotlib.pyplot.figure()
	nullSpectrum = DataSpectrum(0, 0); # Generate a 3-mixing-only scenario.
	anexamplespectrum = DataSpectrum(10000, 0.2)
	nullplot = matplotlib.pyplot.plot(kelist,nullSpectrum, label='no sterile')
	sterileplot = matplotlib.pyplot.plot(kelist,anexamplespectrum, label='with sterile')
	matplotlib.pyplot.title('Effect of a 10keV Sterile with Artificially Large Mixing')
	matplotlib.pyplot.xlabel('Electron Energy (eV)')
	matplotlib.pyplot.ylabel('Rate per unit Energy ($s^-1 eV^-1$)')
	matplotlib.pyplot.legend()
	matplotlib.pyplot.show()
	
# a function to plot the ratio for a realistic pair of values to show the shape
def ratioPlot():
	matplotlib.pyplot.figure()
	nullSpectrum = DataSpectrum(0, 0); # Generate a 3-mixing-only scenario.
	nullSpectrumArray = np.array(nullSpectrum)
	exampleSpectrum = DataSpectrum(10000, 10**-8)
	exampleSpectrumArray = np.array(exampleSpectrum)
	ratiospectrum = exampleSpectrumArray/nullSpectrumArray
	ratioplot = matplotlib.pyplot.plot(kelist,ratiospectrum)
	matplotlib.pyplot.title('Effect of a 2keV Sterile with')
	matplotlib.pyplot.xlabel('Electron Energy (eV)')
	matplotlib.pyplot.ylabel('Ratio of Rates with and Without Sterile ($s^-1 eV^-1$)')
	matplotlib.pyplot.show()
	
#effectPlot()
#ratioPlot()
