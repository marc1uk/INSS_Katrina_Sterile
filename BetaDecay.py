# vim:set noexpandtab tabstop=4 wrap
import numpy as np
import pylab as P
import math 
import matplotlib.pyplot
import os, sys
import ROOT
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
# 1. Function for representing the tritium beta decay spectrum, including a potential sterile
# ===========================================================================================
# constants and support functions
# -------------------------------
# Use eV as as universal unit(?)
Ns = 1.47e-13					# Total number of emitted electron for KATRIN
me = 511						# Mass of electron
kb = 8.76e-5					# Boltzmann radius
a0 = 2.68e-4;					# Bohr radius
numbins = 100					# number of bins for simulated data
kemin = 1000					# min of the energy range - 1keV
kemax = 20000					# max of the energy range - 20keV
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
def beta(Ke,Q,Mixing,Masses): # (NOT DONE)
	### FIXME
	### isn't Q a constant? Why is it needed as an argument?
	### Why is 'Masses' plural? doesn't the spectrum only depend on the (one) sterile mass?
	
	p = []						# excitation probabilities
	Exe = []					# excitation energies
	
	rate = Ns*Fermi(Ke)*Energy(me,Ke)*Momentum(me,Ke)
	
	sum = 0
	
	for (prob,excit) in zip(p,Exe):
		
		for (mix,mass) in zip(Mixing,Masses):
			temp = prob*(Q-excit-Ke) * mix**2 * math.sqrt((Q-excit-Ke)**2 - mass**2)
	
			if ((Q-excit-Ke) - mass) >= 0.:
				sum += temp
	
	rate *= sum
	
	return rate
# -------------------------------
# To create a TF1 object with a custom function, the function must take it's arguments via arrays
# We could do away with this function by suitably defining beta - this is just a wrapper.
def BetaSpectrum(variables,parameters):
	Ke = variables[0]
	Mixing = parameters[0]
	Mass = parameters[1]
	return beta(Ke, 0, Mixing, Mass) ### FIXME remove 0 if not passing Q
# -------------------------------
# A TF1 object based on the spectrum is needed to do fits with PyROOT, find chi2 etc.
npars = 2						# Num parameters required for our custom function
betatf1 = TF1('betaspectrum', BetaSpectrum, kemin, kemax, npars)

# ===========================================================================================
# 2. Function for generating fake data spectrum for given mass and mixing angle
# ===========================================================================================
def DataSpectrum(mass,Mixing):	# simulated count spectrum for a given mass and mixing angle
	### FIXME it would be easier/nicer to step over Sin^2(theta_s) directly than theta_s
	### so ensure DataSpectrum takes this as an argument
	betatf1.SetParameters(mass,Mixing)

	spectrumhist = TH1F('spectrum', 'Simulated Data', numbins, kemin, kemax)
	# not sure how to do this - do we want to generate a suitable number of 'hits'
	for hit in range(0,numhits):
		KeOfThisHit = betatf1.GetRandom()
		spectrumhist.Fill(KeOfThisHit)
	# convert histogram to array...?
	spectrum = spectrumhist.GetArray()			# this returns a pointer in c++...??? what about in python??
	return spectrum								# return simulated data
	
	# or can we somehow just generate the number of entries in each bin directly?
	for i in range(0,numbins-1):				# spectrum runs from 0.1 to 18.574keV over 100 bins.
		Ke = kemin + (i/numbins)*(kemax-kemin)	# Ke of this bin centre in eV
		bincount = ??? ### FIXME
		# normalise to the number of counts over a 3 year period?
		datacount = totalcounts*unitcount 		### FIXME necesary? need to define totalcounts.
		spectrum.append(datacount)
	
	return spectrum								# return simulated data

# ===========================================================================================
# 3. Function for calculating the chi^2 between a given dataset and the no-sterile hypothesis
# ===========================================================================================
def Chi2Test(observe, expect):	# chi2 test for each bin
    sigma = np.sqrt(expect)
    return value = (observe**2 - expect**2)/sigma**2

def Chi2FitToNull(dataSpectrum):	# Sum the elemental chi2 value.
    nullSpectrum = DataSpectrum(0, 0);
    assert len(nullSpectrum) == len(dataSpectrum), "The bin numbers of the spectra do not match."    
    summation = 0
	for i in range (0, len(dataSpectrum)):
        chi2 = Chi2Test(nullSpectrum[i], dataspectrum[i]); ### XXX pass as 2 args rather than passing diff
        summation += chi2
    return summation
    
# ===========================================================================================
# 4. Function for calculating the matrix of chi^2 over a range of masses and mixing angles
# ===========================================================================================
def CalculateChi2Matrix:
	nummixingangs = 20			# number of steps in mixing angle range scan
	mixingangmin = 6			# we'll scan a range of mixings such that
	mixingangmax = 10			# sin^2(theta_s) = 10^-(mixingangmax) -> 10^-(mixingangmin)
	nummasses = 20				# number of steps in mass range scan
	massmin = 1					# minimum mass in keV
	massmax = 20				# maximum mass in keV
	#alpha = 2.7x10^-3 = 0.0027 for 3 sigma -> chi^2 ~8
	chi2for3sigma = 8			# what chi2 value corresponds to a 3 sigma deviation from the null hypothesis
	chi2matrix = 
	sin2mixangmin = 0			# range of sin^2(theta_s) values from mixingangmax and mixingangmin
	sin2mixangmax = 0			# this range is needed for fitting a curve to the chi2 values
# 4. Loop over a vector of prospective neutrino masses.
	for i in range(0,nummasses-1):
		mass = massmin*1000 + (i/nummasses)*((massmax-massmin)*1000)
		chi2forthismass = range(1,nummixingangs)
#    +  Loop over a vector of prospective mixing angles.
		for j in range(0,nummixingangs-1):
			power = mixingangmin + (j/nummixingangs)*(mixingangmax-mixingangmin)
			sin2mixang = 10 ** power
#       - Note the max and min sin^2(mixingang) values; we'll use them for fitting 
#         our range of chi2 vs sin2mixang
			if (j==0)
				sin2mixangmin = sin2mixang
			else if (j==nummixingangs-1)
				sin2mixangmax = sin2mixang
#       - Generate a set of datapoints representing a measured spectrum for that mass and mixing angle.
			dataset = DataSpectrum(mass, sin2mixang)
#       - Calculate the chi^2 of the null hypothesis fit to the data
			chi2 = Chi2FitToNull(dataset)
#       - Store the result in a vector for this mass.
			chi2forthismass[j] = chi2
		# plot the chi2 curve, for checking
		thismassplot = matplotlib.pyplot.plot(chi2forthismass)
		matplotlib.pyplot.show()
#    +  Fit a curve to the span of chi^2 vs mixing angle values, for this mass
		fitfunc = TF1( 'fitfunc', 'pol2',  sin2mixangmin,  sin2mixangmax) # PRESUME quadratic is suitable
		thefit = massgraph.Fit(fitfunc)		### FIXME Fit is a method of TGraph or TH1
											# - need chi2 vs mixingang as a TGraph... 
#    +  Find the mixing angles at which chi^2 = 3 sigma.
		exclusionlimits = thefit.eval(chi2for3sigma)			### FIXME How do i lookup curve values?
#    +  Add the pair of values to a vector - one entry per mass, each entry the max/min mixing angles.
		chi2matrix.append( [exclusionlimits.first, exclusionlimits.second] )

# ===========================================================================================
# 5. Function for producing the contour plots
# ===========================================================================================
# 5. Implement a function that takes the results and produces a contour plot. Results are stored in
#    a matrix, where each row stores the range of mixing angles that are consistent with
#    the null hypothesis at that mass


