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
<<<<<<< HEAD
def Fermi(Z,Ke)
		1/(1+np.exp(Ke+mass-Ef(Z))/Ke*kb))
=======
>>>>>>> 785383b416f64c7323e51807bf105f35c6eae15c
