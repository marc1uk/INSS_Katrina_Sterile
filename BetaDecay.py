import numpy as np
import pylab as P
import math 

Ns = 1.47e-13
mass = 0.511
kb = 8.76e-8

def Energy(mass,Ke):
    return mass + Ke

def Momentum(mass,Ke):
    return math.sqrt(Energy(mass,Ke)**2 - mass**2)

def beta(Ke,Q,Mixing,Masses):
    Ns*Fermi(2,Ke)*Ee*pe*

#testing
def Fermi(Z,Ke):
	1/(1+np.exp(Ke+mass-Ef(Z))/Ke*kb))
