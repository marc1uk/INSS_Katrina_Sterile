# vim:set noexpandtab tabstop=4 wrap
import numpy as np
import pylab as P
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
Ns = 1.47e-13					# Total number of emitted electron for KATRIN
me = 511						# Mass of electron
kb = 8.76e-5					# Boltzmann radius
a0 = 2.68e-4					# Bohr radius
deltamsq21 = 7.53e-5                            # solar mass splitting 
deltamsq32 = 2.45e-3                            # atmospheric mass splitting
numkebins = 100					# number of bins for simulated data
nummixingangs = 20				# number of steps in mixing angle range scan
nummasses = 20					# number of steps in mass range scan
kemin = 0						# min of the electron energy spectrum (eV)
kemax = 18575					# max of the electron energy spectrum (eV)
mixingangmin = 6				# we'll scan a range of mixings such that sin^2(theta_s)
mixingangmax = 10				# spans the range 10^-(mixingangmax) -> 10^-(mixingangmin)
massmin = 1000					# minimum sterile mass (eV)
massmax = 20000					# maximum sterile mass (eV)
#alpha = 2.7x10^-3 = 0.0027 for 3 sigma -> chi^2 ~8
chi2for3sigma = 8				# what chi2 value corresponds to a 3 sigma deviation from the null hypothesis
integrationtime = 3*360*24*3600 # 3 years of running, assuming spectrum is defined in counts/s FIXME if necessary

totalrate = 1.47*(10.**-13)		# total rate in counts second^-1 eV^-5???
totalcounts = totalrate * integrationtime * (kemax*kemin) # FIXME

excitationfile = 'HeTtable.csv'
with open(excitationfile,'rb') as file:
	reader = csv.reader(f,quoting=csv.QUOTE_NONNUMERIC)
	daughthertable = list(reader)
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
def Fermi(Ke):					# The fermi function
	return 1/(1+np.exp(Ke+mass-Ef(Ke))/Ke*kb))
# -------------------------------
def Ef(Ke): 					# The fermi level of H-3
	eta = 1/(a0*Momentrum(me, Ke))
	return 4*np.pi*eta/(1-np.exp(-4*np.pi*eta))
# -------------------------------
def beta(Ke,Q,Mixing,masses): # (NOT DONE)
	### FIXME
	### isn't Q a constant? Why is it needed as an argument?
	### Why is 'Masses' plural? doesn't the spectrum only depend on the (one) sterile mass?
	
	m1 = masses[0]
	m2 = math.sqrt(m1**2 + deltamsq21)
	m3 = math.sqrt(m2**2 + deltamsq32)
	Masses = [m1,m2,m3]
	if len(masses)==2:
		m4 = masses[1]
		Masses.append(m4)


	rate = Ns*Fermi(Ke)*Energy(me,Ke)*Momentum(me,Ke)
	
	sum = 0
	
	for excitpair in daughtertable:
		
		for (mix,mass) in zip(Mixing,Masses):
			temp = excitpair[1]*(Q-excitpair[0]-Ke) * mix**2 * math.sqrt((Q-excitpair[0]-Ke)**2 - mass**2)
	
			if ((Q-excitpair[0]-Ke) - mass) >= 0.:
				sum += temp
	
	rate *= sum
	
	return rate


def betaBACKUP(Ke,Q,Mixing,ms): # (NOT DONE) A backup option of implementing the 3+1 model with keV sterile neutrion.
    # Mixing element in this function is only U_e4, instead of all other U_ei.
    
    Uei = [0.82, 0.54, -0.15] # Declare three U_ei elements.
	m1 = 0
	m2 = math.sqrt(m1**2 + deltamsq21)
	m3 = math.sqrt(m2**2 + deltamsq32)
	Masses = [m1,m2,m3] # 3-gen masses
    
    sin2e4 = Mixing
    cos2e4 = 1-sin2e4
   
	rate = Ns*Fermi(Ke)*Energy(me,Ke)*Momentum(me,Ke)
	
	sum = 0
	
	for pair in len(daughtertable):
        
		for (mix,mass) in zip(Uei,Masses):
			lightPart = daughtertable[pair][1]*(Q-daughtertable[pair][0]-Ke) * mix**2 * math.sqrt((Q-daughtertable[pair][0]-Ke)**2 - mass**2)
        ## The sterile component 
        heavyPart = daughtertable[pair][1]*(Q-daughtertable[pair][0]-Ke)*math.sqrt((Q-daughtertable[pair][0]-Ke)**2 - ms**2)

			if ((Q-daughtertable[pair][0]-Ke) - mass) >= 0.:
				lightSum += lightPart
                heavySum += heavyPart
                
	rate *= cos2e4*lightSum + sin2e4*heavySum
    

def BetaHist():  # Generate the beta histogram
    hist = []    # Declare the content of the histogram
    binWidth = (kemin-kemax)/numbins
    for ki in range(kemin, kemax, binWidth):
        leftBinEdge = beta(ki, Q, 0, 0)  # FIXME assuming Q is a variable and null scenario.
        rightBinEdge = beta(ki+binWidth, Q, 0, 0)
        binContent = binWidth*(leftBinEdge+rightBinEdge)/2
        binContent *= integrationtime    # Unnecesarry if bin content is event rate.
        hist.append(binContent);
    return hist
          
# -------------------------------
# To create a TF1 object with a custom function, the function must take it's arguments via arrays
# We could do away with this function by suitably defining beta - this is just a wrapper.
def BetaSpectrum(variables,parameters):
	Ke = variables[0]
	Mixing = parameters[0]
	Mass = parameters[1]
	return beta(Ke, 0, Mixing, Mass) ### FIXME remove 0 if not passing Q
# -------------------------------
# A TF1 object based on the spectrum is needed to generate the expected counts for simulated data
npars = 2						# Num parameters required for our custom function
betatf1 = TF1('betaspectrum', BetaSpectrum, kemin, kemax, npars) # this is the rate per unit time per unit energy at a given energy

# ===========================================================================================
# 2. Function for generating fake data spectrum for given mass and mixing angle
# ===========================================================================================
def DataSpectrum(mass,Mixing):	# simulated count spectrum for a given mass and mixing angle
	### FIXME it would be easier/nicer to step over Sin^2(theta_s) directly than theta_s
	### so ensure DataSpectrum takes this as an argument

#	# method 1: generate a suitable number of 'hits' point by point and bin into a histogram
#	spectrumhist = TH1F('spectrum', 'Simulated Data', numkebins, kemin, kemax)
#	betatf1.SetParameters(mass,Mixing)
#	for hit in range(0,totalcounts):
#		KeOfThisHit = betatf1.GetRandom()
#		spectrumhist.Fill(KeOfThisHit)
#	spectrumbuffer = spectrumhist.GetArray()				# this returns a pointer to a ROOT buffer
#	spectrum=[]
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
	for i in range(0,numkebins-1):							# spectrum runs from 0.1 to 18.574keV over 100 bins.
		Ke = kemin + i*kebinwidth							# Ke of this bin centre in eV
		expectedbincount = BetaSpectrum(Ke)*kebinwidth*integrationtime
		spectrum.append(expectedbincount)					# same as above just skip the Poisson step

	return spectrum											# return simulated data

# ===========================================================================================
# 3. Function for calculating the chi^2 between a given dataset and the no-sterile hypothesis
# ===========================================================================================
def Chi2Test(observe, expect):	# chi2 test for each bin
    sigma = np.sqrt(expect)
    return value = (observe**2 - expect**2)/sigma**2

def Chi2FitToNull(dataSpectrum):	# Sum the elemental chi2 value.
    nullSpectrum = DataSpectrum(0, 0)
    assert len(nullSpectrum) == len(dataSpectrum), "The bin numbers of the spectra do not match."    
    summation = 0
    for i in range (0, len(dataSpectrum)):
        chi2 = Chi2Test(nullSpectrum[i], dataspectrum[i]) ### XXX fixed to pass as 2 args rather than passing diff
        summation += chi2
    return summation
    
# ===========================================================================================
# 4. Function for calculating the matrix of chi^2 over a range of masses and mixing angles
# ===========================================================================================
def CalculateChi2Matrix:
	chi2array = []				# make 
	mixinganglelist = []
# 4. Loop over a vector of prospective neutrino masses.
	for i in range(0,nummasses-1):
		mass = massmin + (i/nummasses)*(massmax-massmin)
		chi2forthismass = []
#    +  Loop over a vector of prospective mixing angles.
		for j in range(0,nummixingangs-1):
			power = mixingangmin + (j/nummixingangs)*(mixingangmax-mixingangmin)
			sin2mixang = 10 ** power
			if (i==0):
				mixinganglelist.append(sin2mixang)
#       - Generate a set of datapoints representing a measured spectrum for that mass and mixing angle.
			dataset = DataSpectrum(mass, sin2mixang)
#       - Calculate the chi^2 of the null hypothesis fit to the data
			chi2 = Chi2FitToNull(dataset)
#       - Store the result in a vector for this mass.
			chi2forthismass.append(chi2)
		# plot the chi2 curve, for checking
		thismassplot = matplotlib.pyplot.plot(chi2forthismass)
		matplotlib.pyplot.show()
		# turns out to draw contours we can just use the curves directly, we don't need to find the values at 3 sigma
		chi2array.append(chi2forthismass)
	return chi2array
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
def MakeContourPlots(masses, mixings, chi2vals):
	matplotlib.pyplot.figure()
	contourplot = matplotlib.pyplot.contourf(masses, mixings, chi2vals, chi2for3sigma) # contour gives lines, contourf fills
	matplotlib.pyplot.title('3 Sigma Exclusion Limits over 3 Years, Statistical Only')
	matplotlib.pyplot.xlabel('Sterile Mass (eV)')
	matplotlib.pyplot.ylabel('Sterile-Electron Mixing Angle (rads)')
	matplotlib.pyplot.show()

# ===========================================================================================
# 6. The Main Function
# ===========================================================================================
def main():
	thechi2matrix = CalculateChi2Matrix()
	MakeContourPlots()
	## wait for input to keep the plot alive. 
	if __name__ == '__main__':
		rep = ''
		while not rep in [ 'q', 'Q' ]:
			rep = raw_input( 'enter "q" to quit: ' )
			if 1 < len(rep):
				rep = rep[0]
