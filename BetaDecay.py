# vim:set noexpandtab tabstop=4 wrap
import numpy as np
import pylab as P
import math 

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
def beta(Ke,Q,Mixing,Masses):	# The beta spectrum (NOT DONE)
	### FIXME
	### isn't Q a constant? Why is it needed as an argument?
	### is this returning the data count spectrum returning (a vector),
	### or the function specifying the underlying distribution (a TF1)?
	
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

# ===========================================================================================
# 2. Function for generating fake data spectrum for given mass and mixing angle
# ===========================================================================================
def DataSpectrum(mass,Mixing):	# simulated count spectrum for a given mass and mixing angle
	### FIXME it would be easier/nicer to step over Sin^2(theta_s) directly than theta_s
	### so ensure DataSpectrum takes this as an argument
	spectrum = []
	numbins = 100
	for i in range(0,numbins-1):		# spectrum runs from 0.1 to 18.574keV over 100 bins.
		Ke = 0.1 + 186.6*i		# gives Ke in eV.
		spectrum.append(beta(Ke,0,Mixing,mass))  # FIXME Q=???
	return spectrum				# return the spectrum binned into 100 bins
# ===========================================================================================
# 3. Function for calculating the chi^2 between a given dataset and the no-sterile hypothesis
# ===========================================================================================
def Chi2FitToNull(dataspectrum):
### TODO to be filled in

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
# 4. Loop over a vector of prospective neutrino masses.
	for i in range(0,nummasses-1):
		mass = massmin*1000 + (i/nummasses)*((massmax-massmin)*1000)
		chi2forthismass = range(1,nummixingangs)
#    +  Loop over a vector of prospective mixing angles.
		for j in range(0,nummixingangs-1):
			power = mixingangmin + (j/nummixingangs)*(mixingangmax-mixingangmin)
			sin2mixang = 10 ** power
#       - Generate a set of datapoints representing a measured spectrum for that mass and mixing angle.
			dataset = DataSpectrum(mass, sin2mixang)
#       - Calculate the chi^2 of the null hypothesis fit to the data
			chi2 = Chi2FitToNull(dataset)
#       - Store the result in a vector for this mass.
			chi2forthismass[j]=chi2
#    +  Fit a curve to the span of chi^2 vs mixing angle values, for this mass
		thefit = SomeFit(chi2forthismass)	### FIXME How do i fit curve?
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


