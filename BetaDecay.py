import numpy as np
import pylab as P
import math 

Ns = 1.47e-13

def Energy(mass,Ke):
    return mass + Ke

def Momentum(mass,Ke):
    return math.sqrt(Energy(mass,Ke)**2 - mass**2)

def beta(Ke,Q,Mixing,Masses):
    Ns*Fermi(2,Ke)*Energy(me,Ke)*Momentum(me,Ke)*


