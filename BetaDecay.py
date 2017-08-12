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
#    The no-sterile spectrum - null hypothesis - is just this spectrum evaluated with a mixing angle of 0.
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


# -------------------------------------------------------------------------
# 1. Function for representing the tritium beta decay spectrum, including a potential sterile
# -------------------------------------------------------------------------
def BetaSpectrum(mass,Mixing):

Ns = 1.47e-13

def Energy(mass,Ke):
    return mass + Ke

def Momentum(mass,Ke):
    return math.sqrt(Energy(mass,Ke)**2 - mass**2)

def beta(Ke,Q,Mixing,Masses):
    Ns*Fermi(2,Ke)*Ee*pe*

  return spectrum          # this should contain the binned spectrum, with 125 bins
                           #  (based on Katrin's energy resolution)
                           
                           
