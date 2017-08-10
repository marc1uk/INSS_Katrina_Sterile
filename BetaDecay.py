import numpy as np
import pylab as P
import math 

# Use eV as as universal unit(?)
Ns = 1.47e-13	#Total number of emitted electron for KATRIN
me = 511			#Mass of electron
kb = 8.76e-5	#Boltzmann radius
a0 = 2.68e-4;

def Energy(mass,Ke):	#Energy of electron
    return mass + Ke

def Momentum(mass,Ke):	#Momentum of electrom
    return math.sqrt(Energy(mass,Ke)**2 - mass**2)

def beta(Ke,Q,Mixing,Masses):	#The beta spectrum (NOT DONE)
    Ns*Fermi(2,Ke)*Ee*pe*

#testing
def Fermi(Ke):	# The fermi function
	return 1/(1+np.exp(Ke+mass-Ef(Ke))/Ke*kb))

def Ef(Ke): 		# The fermi level of H-3
	eta = 1/(a0*Momentrum(me, Ke))
	return 4*np.pi*eta/(1-np.exp(-4*np.pi*eta))
