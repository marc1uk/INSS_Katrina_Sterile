import os, sys
import ROOT
from ROOT import TCanvas, TF1, TGraph
from ROOT import gROOT, gSystem, Double

x = ( 1.913521, 1.953769, 2.347435, 2.883654, 3.493567,
      4.047560, 4.337210, 4.364347, 4.563004, 5.054247,
      5.194183, 5.380521, 5.303213, 5.384578, 5.563983,
      5.728500, 5.685752, 5.080029, 4.251809, 3.372246,
      2.207432, 1.227541, 0.8597788,0.8220503,0.8046592,
      0.7684097,0.7469761,0.8019787,0.8362375,0.8744895,
      0.9143721,0.9462768,0.9285364,0.8954604,0.8410891,
      0.7853871,0.7100883,0.6938808,0.7363682,0.7032954,
      0.6029015,0.5600163,0.7477068,1.188785, 1.938228,
      2.602717, 3.472962, 4.465014, 5.177035 )

def myFunc(x,par):
	# x is the variable, par is an array of parameters, whose size we need to specify when creating the TF1
	a=par[0]
	b=par[1]
	c=par[2]
	xval = x[0]
	returnval = a*(xval**2) + b*xval + c
	return returnval

fitrangelow  = -10.
fitrangehigh =  10.
fa3 = TF1("afit",myFunc,fitrangelow,fitrangehigh,3)
fa3.SetParameters(3.,0.,-5.)
fa3.Draw()
rep = ''
while not rep in [ 'q', 'Q' ]:
	rep = raw_input( 'enter "q" to quit: ' )
	if 1 < len(rep):
		rep = rep[0]
