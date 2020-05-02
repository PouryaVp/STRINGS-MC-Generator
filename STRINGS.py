#!/usr/bin/python

###################################################################################################################################
# STRINGS version 1-00 - Monte Carlo Event Generator for the Prtoduction and Decay of String Resonances in Proton-Proton Collisions
# Authors:  P. Vakilipourtakalou, D. M. Gingrich. October 2018
###################################################################################################################################

from __future__ import division    
import lhapdf                     
import math                       
import random                     
import argparse                   
from subprocess import call       

############################################################################################################
# Event Generator Main Code
# All of the input parameters and their default values are defined after the definition of the main function
############################################################################################################

def main(RandGenSEED, COME, Ms, ymax, Counter, MinM, MaxM, stringCoeff, SecondStringCoeff, QCDCoeff, PDFSet, PDFScale, Coupling, CouplingScale, dmass, umass, smass, cmass, bmass, tmass, gg2gg, gg2qqbar, gq2gq, gqbar2gqbar, qqbar2gg, gg2gGamma, gq2qGamma):

	# A seed for the random number generator
	# By changing this seed, different sequences of random numbers are generated
	random.seed(RandGenSEED)

	# Printing the input parameters with their values on the screen
	print("####                                                                                                                              ")
	print("#### Inputs:                                                                                                                      ")
	print("####                                                                                                                              ")
	print("#### {0:<10}{1:<22}{2:<30}".format(RandGenSEED, 'RandGenSEED', 'Seed for the the Random Number Generator'))
	print("#### {0:<10}{1:<22}{2:<30}".format(COME, 'COME', 'Centre of Mass Energy (GeV)'))
	print("#### {0:<10}{1:<22}{2:<30}".format(Counter, 'Number', 'Number of Events'))
	print("#### {0:<10}{1:<22}{2:<30}".format(Ms, 'Ms', 'String Scale (GeV)'))
	print("#### {0:<10}{1:<22}{2:<30}".format(MinM, 'MinMass', 'Minimum Invariant Mass (GeV)'))
	print("#### {0:<10}{1:<22}{2:<30}".format(MaxM, 'MaxMass', 'Maximum Invariant Mass (GeV)'))
	print("#### {0:<10}{1:<22}{2:<30}".format(ymax, 'ymax', 'Upper Bound for the Rapidity of the Outgoing Partons'))
	print("#### {0:<10}{1:<22}{2:<30}".format(PDFSet, 'PDFSet', 'PDF Set of the LHAPDF'))
	if(PDFScale == -22):
		PDFScale = Ms
	print("#### {0:<10}{1:<22}{2:<30}".format(PDFScale, 'PDFScale', 'Scale at Which the PDF Set is Evaluated (GeV)'))
	if(Coupling == -1):
		print("#### {0:<10}{1:<22}{2:<30}".format(Coupling, 'Coupling', 'Running Coupling Constant (alpha_s without 4*pi Factor)'))
		if(CouplingScale == -22):
			CouplingScale = Ms
		print("#### {0:<10}{1:<22}{2:<30}".format(CouplingScale, 'CouplingScale', 'Scale at Which the Running Coupling is Calculated (GeV)'))
	else:
		print("#### {0:<10}{1:<22}{2:<30}".format(Coupling, 'Coupling', 'Coupling Constant (alpha_s without 4*pi Factor)'))
	
	print("#### {0:<10}{1:<22}{2:<30}".format(dmass, 'dMass', 'Mass of the Down Quark (GeV)'))
	print("#### {0:<10}{1:<22}{2:<30}".format(umass, 'dMass', 'Mass of the Up Quark (GeV)'))
	print("#### {0:<10}{1:<22}{2:<30}".format(smass, 'dMass', 'Mass of the Strange Quark (GeV)'))
	print("#### {0:<10}{1:<22}{2:<30}".format(cmass, 'dMass', 'Mass of the Charm Quark (GeV)'))
	print("#### {0:<10}{1:<22}{2:<30}".format(bmass, 'dMass', 'Mass of the Bottom Quark (GeV)'))
	print("#### {0:<10}{1:<22}{2:<30}".format(tmass, 'dMass', 'Mass of the Top Quark (GeV)'))


	if(QCDCoeff == "true" or QCDCoeff == "True"):
		print("#### {0:<10}{1:<22}{2:<30}".format(QCDCoeff, 'QCDCoeff', '(Enabled)  Production of QCD tree-level diparton'))
		QCDCoeff11 = 1
	else:
		print("#### {0:<10}{1:<22}{2:<30}".format(QCDCoeff, 'QCDCoeff', '(Disabled) Production of QCD tree-level diparton'))
		QCDCoeff11 = 0

	if(stringCoeff == "true" or stringCoeff == "True"):
		print("#### {0:<10}{1:<22}{2:<30}".format(stringCoeff, 'FirstStringCoeff', '(Enabled)  Production of First  String Resonance ( 2 --> 2 Partonic Scattering and 2-Parton --> Parton-Gamma Scattering )'))
		stringCoeff11 = 1
	else:
		print("#### {0:<10}{1:<22}{2:<30}".format(stringCoeff, 'FirstStringCoeff', '(Disabled) Production of First  String Resonance  ( 2 --> 2 Partonic Scattering and 2-Parton --> Parton-Gamma Scattering )'))
		stringCoeff11 = 0

	if(SecondStringCoeff == "true" or SecondStringCoeff == "True"):
		print("#### {0:<10}{1:<22}{2:<30}".format(SecondStringCoeff, 'SecondStringCoeff', '(Enabled)  Production of Second String Resonance ( 2 --> 2 Partonic Scattering )'))
		SecondStringCoeff11 = 1
	else:
		print("#### {0:<10}{1:<22}{2:<30}".format(SecondStringCoeff, 'SecondStringCoeff', '(Disabled) Production of Second String Resonance ( 2 --> 2 Partonic Scattering )'))
		SecondStringCoeff11 = 0

	NPRUP  = 0

	if(gg2gg == "True" or gg2gg == "true"):
		print("#### {0:<10}{1:<22}{2:<30}".format(gg2gg, 'gg2gg', '(Enabled)  gg --> gg Subprocess       (ID = 1)'))
		gg2gg11 = 1
		NPRUP = NPRUP + 1
	else:
		print("#### {0:<10}{1:<22}{2:<30}".format(gg2gg, 'gg2gg', '(Disabled) gg --> gg Subprocess       (ID = 1)'))
		gg2gg11 = 0

	if(gg2qqbar == "True" or gg2qqbar == "true"):
		print("#### {0:<10}{1:<22}{2:<30}".format(gg2qqbar, 'gg2qqbar', '(Enabled)  gg --> qqbar Subprocess    (ID = 2)'))
		gg2qqbar11 = 1
		NPRUP = NPRUP + 1
	else:
		print("#### {0:<10}{1:<22}{2:<30}".format(gg2qqbar, 'gg2qqbar', '(Disabled) gg --> qqbar Subprocess    (ID = 2)'))
		gg2qqbar11 = 0

	if(gq2gq == "True" or gq2gq == "true"):
		print("#### {0:<10}{1:<22}{2:<30}".format(gq2gq, 'gq2gq', '(Enabled)  gq --> gq Subprocess       (ID = 3)'))
		gq2gq11 = 1
		NPRUP = NPRUP + 1
	else:
		print("#### {0:<10}{1:<22}{2:<30}".format(gq2gq, 'gq2gq', '(Disabled) gq --> gq Subprocess       (ID = 3)'))
		gq2gq11 = 0

	if(gqbar2gqbar == "True" or gqbar2gqbar == "true"):
		print("#### {0:<10}{1:<22}{2:<30}".format(gqbar2gqbar, 'gqbar2gqbar', '(Enabled)  gqbar --> gqbar Subprocess (ID = 4)'))
		gqbar2gqbar11 = 1
		NPRUP = NPRUP + 1
	else:
		print("#### {0:<10}{1:<22}{2:<30}".format(gqbar2gqbar, 'gqbar2gqbar', '(Disabled) gqbar --> gqbar Subprocess (ID = 4)'))
		gqbar2gqbar11 = 0

	if(qqbar2gg == "True" or qqbar2gg == "true"):
		print("#### {0:<10}{1:<22}{2:<30}".format(qqbar2gg, 'qqbar2gg', '(Enabled)  qqbar --> gg Subprocess    (ID = 5)'))
		qqbar2gg11 = 1
		NPRUP = NPRUP + 1
	else:
		print("#### {0:<10}{1:<22}{2:<30}".format(qqbar2gg, 'qqbar2gg', '(Disabled) qqbar --> gg Subprocess    (ID = 5)'))
		qqbar2gg11 = 0

	if(gg2gGamma == "True" or gg2gGamma == "true"):
		print("#### {0:<10}{1:<22}{2:<30}".format(gg2gGamma, 'gg2gGamma', '(Enabled)  gg --> gGamma Subprocess   (ID = 6)'))
		gg2gGamma11 = 1
		NPRUP = NPRUP + 1
	else:
		print("#### {0:<10}{1:<22}{2:<30}".format(gg2gGamma, 'gg2gGamma', '(Disabled) gg --> gGamma Subprocess   (ID = 6)'))
		gg2gGamma11 = 0

	if(gq2qGamma == "True" or gq2qGamma == "true"):
		print("#### {0:<10}{1:<22}{2:<30}".format(gq2qGamma, 'gq2qGamma', '(Enabled)  gq --> qGamma Subprocess   (ID = 7)'))
		gq2qGamma11 = 1
		NPRUP = NPRUP + 1
	else:
		print("#### {0:<10}{1:<22}{2:<30}".format(gq2qGamma, 'gq2qGamma', '(Disabled) gq --> qGamma Subprocess   (ID = 7)'))
		gq2qGamma11 = 0

	print("####                                                ")
	print("########################################################")
	print("########################################################")


	# Specifying the Parton Distribution Functions

	DistFunc = lhapdf.mkPDF(PDFSet, 0)

	print("\nInitializing...\n")

	print("\n(Each dot represents generation of 50 events)\n")

	# Centre of mass energy

	s = COME**2

	N = 3
	Nf = 6

	
	if(Coupling == -1):
		# Runnung Coupling Constant. It is calculated at the CouplingScale
		alpha3 = 1/(1/0.118 + 7/2/math.pi*math.log(CouplingScale/91.2))
	else:
		# Fixed coupling constant
		alpha3 = Coupling

	g3 = math.sqrt(4*math.pi*alpha3)

	# Decay widths that will be used in the scattering amplitudes

	# Decay Widths of First Resonances
	GJ0gS   = g3**2/4/math.pi*Ms*N/4
	GJ0cS   = g3**2/4/math.pi*Ms*N/2
	GJhfqS  = g3**2/4/math.pi*Ms*N/8
	GJ2gS   = g3**2/4/math.pi*Ms*(N/10 + Nf/40)
	GJ2cS   = g3**2/4/math.pi*Ms*(N/5  + Nf/40)
	GJ3hfqS = g3**2/4/math.pi*Ms*N/16

	# Decay Widths of Second Resonances
	GssJhfqS  = g3**2/4/math.pi*math.sqrt(2)*Ms*N/24
	GssJ3hfqS = g3**2/4/math.pi*math.sqrt(2)*Ms*19*N/240
	GssJ5hfqS = g3**2/4/math.pi*math.sqrt(2)*Ms*N/60

	def uniformlist(xmi, xma, n):
		a = []
		for i in range(n+1):
			a.append(xmi + (xma-xmi)*i/n)
		return(a)

	# Definition of the seven functions corresponding to seven subprocesses
	# Each event is specified by (M, Y, y). The ID of the first and second incoming partons are determined by FstpartonID and ScndpartonID, respectively.

	# gg->gg 
	def MonteCarlogggg(parameters, *args):

	    Y, y                           = parameters
	    M, FstpartonID, ScndpartonID   = args
	    
	    tau = M**2/s

	    xa  = math.sqrt(tau)*math.exp( Y)
	    xb  = math.sqrt(tau)*math.exp(-Y)

	    shat = M**2
	    that = -0.5*M**2*math.exp(-y)/math.cosh(y)
	    uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

	    # Production of the first string resonance
	    MSquareString = stringCoeff11*(8/N**2*g3**4/Ms**4*( (N**2-4)**2/4/(N**2-1)*( Ms**8/((shat-Ms**2)**2+(Ms*GJ0gS)**2) + (uhat**4+that**4)/((shat-Ms**2)**2 + (Ms*GJ2gS)**2)) + Ms**8/((shat-Ms**2)**2+(Ms*GJ0cS)**2) + (uhat**4+that**4)/((shat-Ms**2)**2 + (Ms*GJ2cS)**2)))

	    # Production of the QCD diparton
	    MSquareQCD    = QCDCoeff11*(g3**4*(1/shat**2+1/that**2+1/uhat**2)*(2*N**2/(N**2-1)*(shat**2+that**2+uhat**2)+4*(3-N**2)/N**2/(N**2-1)*(shat+that+uhat)**2))

	    # Total scattering amplitude
	    MSquareTotal  = MSquareString + MSquareQCD

	    
	    #Convolution of the scattering amplitudes with parton distribution function
	    return(gg2gg11*DistFunc.xfxQ(FstpartonID, xa, PDFScale)*DistFunc.xfxQ(ScndpartonID, xb, PDFScale)*(MSquareTotal)/(math.cosh(y)**2*16*math.pi*M*shat))

	#gg->qqbar
	def MonteCarloggqqbar(parameters, *args):
	    
	    Y, y                           = parameters
	    M, FstpartonID, ScndpartonID   = args
	    
	    tau = M**2/s

	    xa  = math.sqrt(tau)*math.exp( Y)
	    xb  = math.sqrt(tau)*math.exp(-Y)

	    shat = M**2
	    that = -0.5*M**2*math.exp(-y)/math.cosh(y)
	    uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

	    # Production of the first string resonance
	    MSquareString = stringCoeff11*(2*2/N/(N**2-1)*Nf*g3**4/Ms**4*(0.5*(N**2-4)*uhat*that*(uhat**2+that**2)/((shat-Ms**2)**2+(Ms*GJ2gS)**2) + uhat*that*(uhat**2+that**2)/((shat-Ms**2)**2+(Ms*GJ2cS)**2)))

	    # Production of the QCD diparton
	    MSquareQCD    = QCDCoeff11*(2*g3**4*Nf*(that**2+uhat**2)/shat**2*(1/2/N/uhat/that*(that+uhat)**2-N/(N**2-1)))

	    # Total scattering amplitude
	    MSquareTotal  = MSquareString + MSquareQCD

	    
	    #Convolution of the scattering amplitudes with parton distribution function
	    return(gg2qqbar11*2*DistFunc.xfxQ(FstpartonID, xa, PDFScale)*DistFunc.xfxQ(ScndpartonID, xb, PDFScale)*(MSquareTotal)/(math.cosh(y)**2*16*math.pi*M*shat))

	#gq->gq
	def MonteCarlogqgq(parameters, *args):
	    
	    Y, y                           = parameters
	    M, FstpartonID, ScndpartonID   = args
	    
	    tau = M**2/s

	    xa  = math.sqrt(tau)*math.exp( Y)
	    xb  = math.sqrt(tau)*math.exp(-Y)

	    shat = M**2
	    that = -0.5*M**2*math.exp(-y)/math.cosh(y)
	    uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

	    # Production of the first string resonance
	    MSquareString       = stringCoeff11*((N**2-1)/2/N**2*g3**4/Ms**2*(-Ms**4*that/((shat-Ms**2)**2+(Ms*GJhfqS)**2) - that**3/((shat-Ms**2)**2+(Ms*GJ3hfqS)**2)) + (N**2-1)/2/N**2*g3**4/Ms**2*(-Ms**4*uhat/((shat-Ms**2)**2+(Ms*GJhfqS)**2) - uhat**3/((shat-Ms**2)**2+(Ms*GJ3hfqS)**2)))

	    # Production of the QCD diparton
	    MSquareQCD          = QCDCoeff11*(g3**4*(shat**2+uhat**2)/that**2*(1-(N**2-1)/2/N**2/uhat/shat*(shat+uhat)**2) + g3**4*(shat**2+that**2)/uhat**2*(1-(N**2-1)/2/N**2/that/shat*(shat+that)**2))

	    # Production of the second string resonance
	    MSquareSecondString = SecondStringCoeff11*(2*(N**2-1)/N**2*(g3**4/2/Ms**2*(1/9*Ms**4*(-uhat)/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJhfqS)**2) + (-1/9)*uhat*(3*that+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2)) + g3**4/8/Ms**6*(9/25*Ms**4*(-uhat)**3/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2) + (-1/25)*uhat**3*(5*that+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ5hfqS)**2))) + 2*(N**2-1)/N**2*(g3**4/2/Ms**2*(1/9*Ms**4*(-that)/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJhfqS)**2) + (-1/9)*that*(3*uhat+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2)) + g3**4/8/Ms**6*(9/25*Ms**4*(-that)**3/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2) + (-1/25)*that**3*(5*uhat+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ5hfqS)**2))))

	    #Total scattering amplitude
	    MSquareTotal        = MSquareString + MSquareQCD + MSquareSecondString

	    
	    #Convolution of the scattering amplitudes with parton distribution function
	    return(gq2gq11*DistFunc.xfxQ(FstpartonID, xa, PDFScale)*DistFunc.xfxQ(ScndpartonID, xb, PDFScale)*(MSquareTotal)/(math.cosh(y)**2*16*math.pi*M*shat) + gq2gq11*DistFunc.xfxQ(FstpartonID, xb, PDFScale)*DistFunc.xfxQ(ScndpartonID, xa, PDFScale)*(MSquareTotal)/(math.cosh(y)**2*16*math.pi*M*shat))

	#gqbar->gqbar
	def MonteCarlogqbargqbar(parameters, *args):
	    
	    Y, y                           = parameters
	    M, FstpartonID, ScndpartonID   = args
	    
	    tau = M**2/s

	    xa  = math.sqrt(tau)*math.exp( Y)
	    xb  = math.sqrt(tau)*math.exp(-Y)

	    shat = M**2
	    that = -0.5*M**2*math.exp(-y)/math.cosh(y)
	    uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

	    # Production of the first string resonance
	    MSquareString       = stringCoeff11*((N**2-1)/2/N**2*g3**4/Ms**2*(-Ms**4*that/((shat-Ms**2)**2+(Ms*GJhfqS)**2) - that**3/((shat-Ms**2)**2+(Ms*GJ3hfqS)**2)) + (N**2-1)/2/N**2*g3**4/Ms**2*(-Ms**4*uhat/((shat-Ms**2)**2+(Ms*GJhfqS)**2) - uhat**3/((shat-Ms**2)**2+(Ms*GJ3hfqS)**2)))

	    # Production of the QCD diparton
	    MSquareQCD          = QCDCoeff11*(g3**4*(shat**2+uhat**2)/that**2*(1-(N**2-1)/2/N**2/uhat/shat*(shat+uhat)**2) + g3**4*(shat**2+that**2)/uhat**2*(1-(N**2-1)/2/N**2/that/shat*(shat+that)**2))

	    # Production of the second string resonance
	    MSquareSecondString = SecondStringCoeff11*(2*(N**2-1)/N**2*(g3**4/2/Ms**2*(1/9*Ms**4*(-uhat)/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJhfqS)**2) + (-1/9)*uhat*(3*that+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2)) + g3**4/8/Ms**6*(9/25*Ms**4*(-uhat)**3/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2) + (-1/25)*uhat**3*(5*that+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ5hfqS)**2))) + 2*(N**2-1)/N**2*(g3**4/2/Ms**2*(1/9*Ms**4*(-that)/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJhfqS)**2) + (-1/9)*that*(3*uhat+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2)) + g3**4/8/Ms**6*(9/25*Ms**4*(-that)**3/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2) + (-1/25)*that**3*(5*uhat+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ5hfqS)**2))))

	    # Total scattering amplitude
	    MSquareTotal        = MSquareString + MSquareQCD + MSquareSecondString

	    
	    #Convolution of the scattering amplitudes with parton distribution function
	    return(gqbar2gqbar11*DistFunc.xfxQ(FstpartonID, xa, PDFScale)*DistFunc.xfxQ(ScndpartonID, xb, PDFScale)*(MSquareTotal)/(math.cosh(y)**2*16*math.pi*M*shat) + gqbar2gqbar11*DistFunc.xfxQ(FstpartonID, xb, PDFScale)*DistFunc.xfxQ(ScndpartonID, xa, PDFScale)*(MSquareTotal)/(math.cosh(y)**2*16*math.pi*M*shat))

	#ggbar->gg
	def MonteCarloqqbargg(parameters, *args):

	    Y, y                           = parameters
	    M, FstpartonID, ScndpartonID   = args
	    
	    tau = M**2/s
	    
	    xa  = math.sqrt(tau)*math.exp( Y)
	    xb  = math.sqrt(tau)*math.exp(-Y)

	    shat = M**2
	    that = -0.5*M**2*math.exp(-y)/math.cosh(y)
	    uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

	    # Production of the first string resonance
	    MSquareString = stringCoeff11*(2*2*(N**2-1)/N**3*g3**4/Ms**4*((N**2-4)/2*uhat*that*(uhat**2+that**2)/((shat-Ms**2)**2+(Ms*GJ2gS)**2) + uhat*that*(uhat**2+that**2)/((shat-Ms**2)**2+(Ms*GJ2cS)**2)))

	    # Production of the QCD diparton
	    MSquareQCD    = QCDCoeff11*(2*g3**4*(that**2+uhat**2)/shat**2*((N**2-1)**2/2/N**3/uhat/that*(that+uhat)**2-(N**2-1)/N))

	    # Total scattering amplitude
	    MSquareTotal  = MSquareString + MSquareQCD

	    #Convolution of the scattering amplitudes with parton distribution function
	    return(qqbar2gg11*DistFunc.xfxQ(FstpartonID, xa, PDFScale)*DistFunc.xfxQ(ScndpartonID, xb, PDFScale)*(MSquareTotal)/(math.cosh(y)**2*16*math.pi*M*shat) + qqbar2gg11*DistFunc.xfxQ(FstpartonID, xb, PDFScale)*DistFunc.xfxQ(ScndpartonID, xa, PDFScale)*(MSquareTotal)/(math.cosh(y)**2*16*math.pi*M*shat))

	#gg->gGamma
	def MonteCarloggGamma(parameters, *args):

	    Y, y                           = parameters
	    M, FstpartonID, ScndpartonID   = args
	    
	    tau = M**2/s
	    
	    xa  = math.sqrt(tau)*math.exp( Y)
	    xb  = math.sqrt(tau)*math.exp(-Y)

	    shat = M**2
	    that = -0.5*M**2*math.exp(-y)/math.cosh(y)
	    uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

	    # Production of the first string resonance
	    MSquareTotal  = stringCoeff11*2*5*g3**4/math.sqrt(6)**2/3/Ms**4*(Ms**8/((shat*Ms**2)**2 + (GJ0gS*Ms)**2) + (that**4 + uhat**4)/((shat-Ms**2)**2 + (GJ2gS*Ms)**2))

	    #Convolution of the scattering amplitudes with parton distribution function
	    return(2*gg2gGamma11*DistFunc.xfxQ(FstpartonID, xa, PDFScale)*DistFunc.xfxQ(ScndpartonID, xb, PDFScale)*(MSquareTotal)/(math.cosh(y)**2*16*math.pi*M*shat))

	#gq->qGamma
	def MonteCarlogqGamma(parameters, *args):

	    Y, y                           = parameters
	    M, FstpartonID, ScndpartonID   = args
	    
	    tau = M**2/s
	    
	    xa  = math.sqrt(tau)*math.exp( Y)
	    xb  = math.sqrt(tau)*math.exp(-Y)

	    shat = M**2
	    that = -0.5*M**2*math.exp(-y)/math.cosh(y)
	    uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

	    # Production of the first string resonance
	    MSquareTotal  = -stringCoeff11*g3**4/math.sqrt(6)**2/3/Ms**4*(Ms**4*uhat/((shat*Ms**2)**2 + (GJhfqS*Ms)**2) + (uhat**3)/((shat-Ms**2)**2 + (GJ3hfqS*Ms)**2)) - g3**4/math.sqrt(6)**2/3/Ms**4*(Ms**4*that/((shat*Ms**2)**2 + (GJhfqS*Ms)**2) + (that**3)/((shat-Ms**2)**2 + (GJ3hfqS*Ms)**2))

	    #Convolution of the scattering amplitudes with parton distribution function
	    return(gq2qGamma11*DistFunc.xfxQ(FstpartonID, xa, PDFScale)*DistFunc.xfxQ(ScndpartonID, xb, PDFScale)*(MSquareTotal)/(math.cosh(y)**2*16*math.pi*M*shat) + gq2qGamma11*DistFunc.xfxQ(FstpartonID, xb, PDFScale)*DistFunc.xfxQ(ScndpartonID, xa, PDFScale)*(MSquareTotal)/(math.cosh(y)**2*16*math.pi*M*shat))

	MC = 4
	MaxPar = [0, 0]

	# Quark Masses
	QuarkMasses = [dmass, umass, smass, cmass, bmass, tmass]

	# Assigning the parameters that are used in the LHE file
	IDBMUP = 2212              # Proton beams
	EBMUP  = COME/2            # Energy of each proton beam
	PDFGUP = 0
	PDFSUP = DistFunc.lhapdfID # ID of the PDF set
	IDWTUP = 3                 # Unweighted events

	QED = 1/128  # QED coupling constant


	# Cross-sections of the individual subprocesses
	CrossSectiongggg      = 0
	CrossSectionErrorgggg = 0

	CrossSectionggqqbar      = 0
	CrossSectionErrorggqqbar = 0

	CrossSectiongqgq      = 0
	CrossSectionErrorgqgq = 0

	CrossSectiongqbargqbar      = 0
	CrossSectionErrorgqbargqbar = 0

	CrossSectionqqbargg      = 0
	CrossSectionErrorqqbargg = 0

	CrossSectionggGamma      = 0
	CrossSectionErrorggGamma = 0

	CrossSectiongqGamma      = 0
	CrossSectionErrorgqGamma = 0



	# Finding the minimum and maximum of the functions
	sumggggMin      = 0
	sumggqqbarMin   = 0
	sumggGammaMin   = 0

	sumGammadMin  = 0
	sumGammauMin  = 0
	sumGammasMin  = 0
	sumGammacMin  = 0
	sumGammabMin  = 0
	sumGammatMin  = 0

	sumgdMin = 0
	sumguMin = 0
	sumgsMin = 0
	sumgcMin = 0
	sumgbMin = 0
	sumgtMin = 0

	sumgdbarMin = 0
	sumgubarMin = 0
	sumgsbarMin = 0
	sumgcbarMin = 0
	sumgbbarMin = 0
	sumgtbarMin = 0

	sumddbarMin = 0
	sumuubarMin = 0
	sumssbarMin = 0
	sumccbarMin = 0
	sumbbbarMin = 0
	sumttbarMin = 0

	argsggMin = (MinM, 21, 21)

	argsgdMin = (MinM, 21, 1)
	argsguMin = (MinM, 21, 2)
	argsgsMin = (MinM, 21, 3)
	argsgcMin = (MinM, 21, 4)
	argsgbMin = (MinM, 21, 5)
	argsgtMin = (MinM, 21, 6)

	argsgdbarMin = (MinM, 21, -1)
	argsgubarMin = (MinM, 21, -2)
	argsgsbarMin = (MinM, 21, -3)
	argsgcbarMin = (MinM, 21, -4)
	argsgbbarMin = (MinM, 21, -5)
	argsgtbarMin = (MinM, 21, -6)

	argsddbarMin = (MinM, 1, -1)
	argsuubarMin = (MinM, 2, -2)
	argsssbarMin = (MinM, 3, -3)
	argsccbarMin = (MinM, 4, -4)
	argsbbbarMin = (MinM, 5, -5)
	argsttbarMin = (MinM, 6, -6)


	sumggggMax      = 0
	sumggqqbarMax   = 0
	sumggGammaMax   = 0

	sumGammadMax  = 0
	sumGammauMax  = 0
	sumGammasMax  = 0
	sumGammacMax  = 0
	sumGammabMax  = 0
	sumGammatMax  = 0

	sumgdMax = 0
	sumguMax = 0
	sumgsMax = 0
	sumgcMax = 0
	sumgbMax = 0
	sumgtMax = 0

	sumgdbarMax = 0
	sumgubarMax = 0
	sumgsbarMax = 0
	sumgcbarMax = 0
	sumgbbarMax = 0
	sumgtbarMax = 0

	sumddbarMax = 0
	sumuubarMax = 0
	sumssbarMax = 0
	sumccbarMax = 0
	sumbbbarMax = 0
	sumttbarMax = 0

	if(Ms < MinM):
		Msss = MinM
	elif(Ms > MaxM):
		Msss = MinM
	else:
		Msss = Ms

	argsggMax = (Msss, 21, 21)

	argsgdMax = (Msss, 21, 1)
	argsguMax = (Msss, 21, 2)
	argsgsMax = (Msss, 21, 3)
	argsgcMax = (Msss, 21, 4)
	argsgbMax = (Msss, 21, 5)
	argsgtMax = (Msss, 21, 6)

	argsgdbarMax = (Msss, 21, -1)
	argsgubarMax = (Msss, 21, -2)
	argsgsbarMax = (Msss, 21, -3)
	argsgcbarMax = (Msss, 21, -4)
	argsgbbarMax = (Msss, 21, -5)
	argsgtbarMax = (Msss, 21, -6)

	argsddbarMax = (Msss, 1, -1)
	argsuubarMax = (Msss, 2, -2)
	argsssbarMax = (Msss, 3, -3)
	argsccbarMax = (Msss, 4, -4)
	argsbbbarMax = (Msss, 5, -5)
	argsttbarMax = (Msss, 6, -6)

	tauMin = MinM**2/s
	YmaxMin = min([math.log(1/math.sqrt(tauMin)), ymax])

	tauMax = Msss**2/s
	YmaxMax = min([math.log(1/math.sqrt(tauMax)), ymax])

	# Defining the range of the integral over y and Y
	# Finding the maximum and minimum of the functions
	for Y in uniformlist(-YmaxMin+10**-8, YmaxMin-10**-8, MC):
		for y in uniformlist(-(ymax-math.fabs(Y)), ymax-math.fabs(Y), MC):

			ActualParameters = [Y, y]

			sumggggMin     += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogggg(ActualParameters, *argsggMin)
			sumggqqbarMin  += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarloggqqbar(ActualParameters, *argsggMin)

			sumggGammaMin  += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarloggGamma(ActualParameters, *argsggMin)

			sumGammadMin   += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgdMin)
			sumGammauMin   += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsguMin)
			sumGammasMin   += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgsMin)
			sumGammacMin   += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgcMin)
			sumGammabMin   += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgbMin)
			sumGammatMin   += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgtMin)

			sumgdMin       += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgdMin)
			sumguMin       += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsguMin)
			sumgsMin       += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgsMin)
			sumgcMin       += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgcMin)
			sumgbMin       += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgbMin)
			sumgtMin       += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgtMin)

			sumgdbarMin    += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgdbarMin)
			sumgubarMin    += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgubarMin)
			sumgsbarMin    += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgsbarMin)
			sumgcbarMin    += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgcbarMin)
			sumgbbarMin    += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgbbarMin)
			sumgtbarMin    += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgtbarMin)

			sumddbarMin    += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsddbarMin)
			sumuubarMin    += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsuubarMin)
			sumssbarMin    += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsssbarMin)
			sumccbarMin    += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsccbarMin)
			sumbbbarMin    += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsbbbarMin)
			sumttbarMin    += 4*YmaxMin*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsttbarMin)

	for Y in uniformlist(-YmaxMax+10**-8, YmaxMax-10**-8, MC):
		for y in uniformlist(-(ymax-math.fabs(Y)), ymax-math.fabs(Y), MC):

			ActualParameters = [Y, y]

			sumggggMax     += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogggg(ActualParameters, *argsggMax)
			sumggqqbarMax  += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarloggqqbar(ActualParameters, *argsggMax)

			sumggGammaMax  += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarloggGamma(ActualParameters, *argsggMax)

			sumGammadMax   += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgdMax)
			sumGammauMax   += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsguMax)
			sumGammasMax   += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgsMax)
			sumGammacMax   += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgcMax)
			sumGammabMax   += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgbMax)
			sumGammatMax   += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgtMax)

			sumgdMax       += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgdMax)
			sumguMax       += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsguMax)
			sumgsMax       += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgsMax)
			sumgcMax       += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgcMax)
			sumgbMax       += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgbMax)
			sumgtMax       += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgtMax)

			sumgdbarMax    += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgdbarMax)
			sumgubarMax    += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgubarMax)
			sumgsbarMax    += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgsbarMax)
			sumgcbarMax    += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgcbarMax)
			sumgbbarMax    += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgbbarMax)
			sumgtbarMax    += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgtbarMax)

			sumddbarMax    += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsddbarMax)
			sumuubarMax    += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsuubarMax)
			sumssbarMax    += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsssbarMax)
			sumccbarMax    += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsccbarMax)
			sumbbbarMax    += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsbbbarMax)
			sumttbarMax    += 4*YmaxMax*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsttbarMax)

	ArrayOfMaxes = []

	ArrayOfMaxes.append(max(sumgdMin, sumgdMax))
	ArrayOfMaxes.append(max(sumguMin, sumguMax))
	ArrayOfMaxes.append(max(sumgsMin, sumgsMax))
	ArrayOfMaxes.append(max(sumgcMin, sumgcMax))
	ArrayOfMaxes.append(max(sumgbMin, sumgbMax))
	ArrayOfMaxes.append(max(sumgtMin, sumgtMax))

	ArrayOfMaxes.append(max(sumddbarMin, sumddbarMax))
	ArrayOfMaxes.append(max(sumuubarMin, sumuubarMax))
	ArrayOfMaxes.append(max(sumssbarMin, sumssbarMax))
	ArrayOfMaxes.append(max(sumccbarMin, sumccbarMax))
	ArrayOfMaxes.append(max(sumbbbarMin, sumbbbarMax))
	ArrayOfMaxes.append(max(sumttbarMin, sumttbarMax))

	ArrayOfMaxes.append(max(sumggggMin, sumggggMax))

	ArrayOfMaxes.append(max(sumgdbarMin, sumgdbarMax))
	ArrayOfMaxes.append(max(sumgubarMin, sumgubarMax))
	ArrayOfMaxes.append(max(sumgsbarMin, sumgsbarMax))
	ArrayOfMaxes.append(max(sumgcbarMin, sumgcbarMax))
	ArrayOfMaxes.append(max(sumgbbarMin, sumgbbarMax))
	ArrayOfMaxes.append(max(sumgtbarMin, sumgtbarMax))

	ArrayOfMaxes.append(max(sumggqqbarMin, sumggqqbarMax))

	ArrayOfMaxes.append(max(sumggGammaMin, sumggGammaMax))

	ArrayOfMaxes.append(max(sumGammadMin, sumGammadMax))
	ArrayOfMaxes.append(max(sumGammauMin, sumGammauMax))
	ArrayOfMaxes.append(max(sumGammasMin, sumGammasMax))
	ArrayOfMaxes.append(max(sumGammacMin, sumGammacMax))
	ArrayOfMaxes.append(max(sumGammabMin, sumGammabMax))
	ArrayOfMaxes.append(max(sumGammatMin, sumGammatMax))

	# Creating the LHE file for saving the events
	fh = open("events.lhe", "w")

	DeltaMass = (MaxM - MinM)/Counter

	numberofevents = 0

	# Event generation loop
	while(numberofevents < Counter):

		# A list that saves the weight of each subprocess at a random invariant mass. This is used to specify the type of subprocess
		FuncAtMass = []

		bol = 1

		while(bol == 1):

			# Random invariant mass in the specified interval			
			M = random.uniform(MinM,MaxM)

			if(M == COME):
				break
			
			tau = M**2/s
			Ymax = min([math.log(1/math.sqrt(tau)), ymax])

			# Initial values for the calculation of the cross-sections of the subprocesses
			sumgggg      = 0
			sumggqqbar   = 0
			sumggGamma   = 0

			sumGammad  = 0
			sumGammau  = 0
			sumGammas  = 0
			sumGammac  = 0
			sumGammab  = 0
			sumGammat  = 0

			sumgd = 0
			sumgu = 0
			sumgs = 0
			sumgc = 0
			sumgb = 0
			sumgt = 0

			sumgdbar = 0
			sumgubar = 0
			sumgsbar = 0
			sumgcbar = 0
			sumgbbar = 0
			sumgtbar = 0

			sumddbar = 0
			sumuubar = 0
			sumssbar = 0
			sumccbar = 0
			sumbbbar = 0
			sumttbar = 0

			argsgg = (M, 21, 21)

			argsgd = (M, 21, 1)
			argsgu = (M, 21, 2)
			argsgs = (M, 21, 3)
			argsgc = (M, 21, 4)
			argsgb = (M, 21, 5)
			argsgt = (M, 21, 6)

			argsgdbar = (M, 21, -1)
			argsgubar = (M, 21, -2)
			argsgsbar = (M, 21, -3)
			argsgcbar = (M, 21, -4)
			argsgbbar = (M, 21, -5)
			argsgtbar = (M, 21, -6)

			argsddbar = (M, 1, -1)
			argsuubar = (M, 2, -2)
			argsssbar = (M, 3, -3)
			argsccbar = (M, 4, -4)
			argsbbbar = (M, 5, -5)
			argsttbar = (M, 6, -6)

			for Y in uniformlist(-Ymax+10**-8, Ymax-10**-8, MC):
				for y in uniformlist(-(ymax-math.fabs(Y)), ymax-math.fabs(Y), MC):

					ActualParameters = [Y, y]

					sumgggg     += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogggg(ActualParameters, *argsgg)
					sumggqqbar  += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarloggqqbar(ActualParameters, *argsgg)

					sumggGamma  += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarloggGamma(ActualParameters, *argsgg)

					sumGammad   += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgd)
					sumGammau   += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgu)
					sumGammas   += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgs)
					sumGammac   += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgc)
					sumGammab   += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgb)
					sumGammat   += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqGamma(ActualParameters, *argsgt)

					sumgd       += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgd)
					sumgu       += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgu)
					sumgs       += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgs)
					sumgc       += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgc)
					sumgb       += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgb)
					sumgt       += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqgq(ActualParameters, *argsgt)

					sumgdbar    += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgdbar)
					sumgubar    += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgubar)
					sumgsbar    += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgsbar)
					sumgcbar    += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgcbar)
					sumgbbar    += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgbbar)
					sumgtbar    += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarlogqbargqbar(ActualParameters, *argsgtbar)

					sumddbar    += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsddbar)
					sumuubar    += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsuubar)
					sumssbar    += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsssbar)
					sumccbar    += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsccbar)
					sumbbbar    += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsbbbar)
					sumttbar    += 4*Ymax*(ymax-math.fabs(Y))/MC**2*MonteCarloqqbargg(ActualParameters, *argsttbar)


			FuncAtMass.append(sumgd)
			FuncAtMass.append(sumgu)
			FuncAtMass.append(sumgs)
			FuncAtMass.append(sumgc)
			FuncAtMass.append(sumgb)
			FuncAtMass.append(sumgt)

			FuncAtMass.append(sumddbar)
			FuncAtMass.append(sumuubar)
			FuncAtMass.append(sumssbar)
			FuncAtMass.append(sumccbar)
			FuncAtMass.append(sumbbbar)
			FuncAtMass.append(sumttbar)

			FuncAtMass.append(sumgggg)

			FuncAtMass.append(sumgdbar)
			FuncAtMass.append(sumgubar)
			FuncAtMass.append(sumgsbar)
			FuncAtMass.append(sumgcbar)
			FuncAtMass.append(sumgbbar)
			FuncAtMass.append(sumgtbar)

			FuncAtMass.append(sumggqqbar)

			FuncAtMass.append(sumggGamma)

			FuncAtMass.append(sumGammad)
			FuncAtMass.append(sumGammau)
			FuncAtMass.append(sumGammas)
			FuncAtMass.append(sumGammac)
			FuncAtMass.append(sumGammab)
			FuncAtMass.append(sumGammat)

			FuncMax = FuncAtMass[0]

			InteractionIndex = 0

			# Determination of the subprocess type

			for k in range(26):
				if (FuncAtMass[k+1] > FuncMax):
					FuncMax = FuncAtMass[k+1]
					InteractionIndex = k+1

			RandomSig = random.uniform(0, FuncMax)

			for k in range(27):
				if(FuncAtMass[k] > RandomSig):
					if(FuncAtMass[k] < FuncAtMass[InteractionIndex]):
						InteractionIndex = k

			RandomY = random.uniform(-Ymax+10**-9, Ymax-10**-9)
			Randomy = random.uniform(-(ymax-math.fabs(RandomY)), (ymax-math.fabs(RandomY)))

			# Seven if statements for calculating the kinematic variables and filling the LHE files, corresponding to seven different subprocesses
			if(InteractionIndex < 6):

				IDQ = InteractionIndex + 1
				args = (M, 21, IDQ)

				RInvMass = random.uniform(0, ArrayOfMaxes[InteractionIndex])
				if(RInvMass > FuncAtMass[InteractionIndex]):
					break

				MaximumEvaluation = MonteCarlogqgq(MaxPar, *args)
				RandomEval = random.uniform(0, MaximumEvaluation)
				ActualParameters = [RandomY, Randomy]
				ActualEval = MonteCarlogqgq(ActualParameters, *args)

				if(RandomEval > ActualEval):
					break

				numberofevents += 1

				# Calculation of the kinematic variables

				Qmass = QuarkMasses[InteractionIndex]

				y1 = RandomY + Randomy
				y2 = RandomY - Randomy

				xa = math.sqrt(tau)*math.exp( RandomY)
				xb = math.sqrt(tau)*math.exp(-RandomY)

				Ea = xa*EBMUP
				Eb = xb*EBMUP

				phi = random.uniform(0, 2*math.pi)

				pt = 0.5*M/math.cosh(Randomy)

				px1 = pt*math.cos(phi)
				py1 = pt*math.sin(phi)

				px2 = -px1
				py2 = -py1

				E1 = pt*math.cosh(y1)
				E2 = pt*math.cosh(y2)

				pz1 = E1*math.tanh(y1)
				pz2 = E2*math.tanh(y2)

				weight = FuncAtMass[InteractionIndex]*(10**9)/2.56819*DeltaMass

				CrossSectiongqgq += weight
				CrossSectionErrorgqgq += weight**2

				# Filling the LHE file

				fh.write("<event>\n")
				fh.write("{}  {}  {}  {}  {} {}\n".format(4, 3, weight, PDFScale, QED, alpha3))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,-1,0,0,101,102,0,0,Ea,Ea,0,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(IDQ,-1,0,0,102,0,0,0,-Eb,Eb,Qmass,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,1,1,2,101,103,px1,py1,pz1,E1,0,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(IDQ,1,1,2,103,0,px2,py2,pz2,E2,Qmass,0,9))
				fh.write("#pdf  {}  {}  {}  {}  {}  {}  {}\n".format(21,IDQ,xa,xb,PDFScale,DistFunc.xfxQ(21, xa, PDFScale),DistFunc.xfxQ(IDQ, xb, PDFScale)))
				fh.write("</event>\n")
				if(numberofevents == 500): print("\n{} Events are Generated\n".format(numberofevents))
				if(numberofevents%50   == 0): print(".")
				if(numberofevents%1000 == 0): print("\n{} Events are Generated\n".format(numberofevents))
				break

			elif(InteractionIndex < 12):

				IDQ = InteractionIndex - 5
				args = (M, IDQ, -IDQ)

				RInvMass = random.uniform(0, ArrayOfMaxes[InteractionIndex])
				if(RInvMass > FuncAtMass[InteractionIndex]):
					break

				MaximumEvaluation = MonteCarloqqbargg(MaxPar, *args)
				RandomEval = random.uniform(0, MaximumEvaluation)
				ActualParameters = [RandomY, Randomy]
				ActualEval = MonteCarloqqbargg(ActualParameters, *args)
				Qmass = QuarkMasses[InteractionIndex - 6]

				if(RandomEval > ActualEval):
					break

				numberofevents += 1

				# Calculation of the kinematic variables

				y1 = RandomY + Randomy
				y2 = RandomY - Randomy

				xa = math.sqrt(tau)*math.exp( RandomY)
				xb = math.sqrt(tau)*math.exp(-RandomY)

				Ea = xa*EBMUP
				Eb = xb*EBMUP

				phi = random.uniform(0, 2*math.pi)

				pt = 0.5*M/math.cosh(Randomy)

				px1 = pt*math.cos(phi)
				py1 = pt*math.sin(phi)

				px2 = -px1
				py2 = -py1

				E1 = pt*math.cosh(y1)
				E2 = pt*math.cosh(y2)

				pz1 = E1*math.tanh(y1)
				pz2 = E2*math.tanh(y2)

				weight = FuncAtMass[InteractionIndex]*(10**9)/2.56819*DeltaMass

				CrossSectionqqbargg += weight
				CrossSectionErrorqqbargg += weight**2

				# Filling the LHE file

				fh.write("<event>\n")
				fh.write("{}  {}  {}  {}  {} {}\n".format(4, 5, weight, PDFScale, QED, alpha3))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(IDQ,-1,0,0,101,0,0,0,Ea,Ea,Qmass,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(-IDQ,-1,0,0,0,103,0,0,-Eb,Eb,Qmass,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,1,1,2,101,102,px1,py1,pz1,E1,0,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,1,1,2,102,103,px2,py2,pz2,E2,0,0,9))
				fh.write("#pdf  {}  {}  {}  {}  {}  {}  {}\n".format(IDQ,-IDQ,xa,xb,PDFScale,DistFunc.xfxQ(IDQ, xa, PDFScale),DistFunc.xfxQ(-IDQ, xb, PDFScale)))
				fh.write("</event>\n")
				if(numberofevents == 500): print("\n{} Events are Generated\n".format(numberofevents))
				if(numberofevents%50   == 0): print(".")
				if(numberofevents%1000 == 0): print("\n{} Events are Generated\n".format(numberofevents))
				break

			elif(InteractionIndex == 12):

				args = (M, 21, 21)

				RInvMass = random.uniform(0, ArrayOfMaxes[InteractionIndex])
				if(RInvMass > FuncAtMass[InteractionIndex]):
					break

				MaximumEvaluation = MonteCarlogggg(MaxPar, *args)
				RandomEval = random.uniform(0, MaximumEvaluation)
				ActualParameters = [RandomY, Randomy]
				ActualEval = MonteCarlogggg(ActualParameters, *args)

				if(RandomEval > ActualEval):
					break

				numberofevents += 1

				# Calculation of the kinematic variables

				y1 = RandomY + Randomy
				y2 = RandomY - Randomy

				xa = math.sqrt(tau)*math.exp( RandomY)
				xb = math.sqrt(tau)*math.exp(-RandomY)

				Ea = xa*EBMUP
				Eb = xb*EBMUP

				phi = random.uniform(0, 2*math.pi)

				pt = 0.5*M/math.cosh(Randomy)

				px1 = pt*math.cos(phi)
				py1 = pt*math.sin(phi)

				px2 = -px1
				py2 = -py1

				E1 = pt*math.cosh(y1)
				E2 = pt*math.cosh(y2)

				pz1 = E1*math.tanh(y1)
				pz2 = E2*math.tanh(y2)

				weight = FuncAtMass[InteractionIndex]*(10**9)/2.56819*DeltaMass

				CrossSectiongggg += weight
				CrossSectionErrorgggg += weight**2

				# Filling the LHE file

				fh.write("<event>\n")
				fh.write("{}  {}  {}  {}  {} {}\n".format(4, 1, weight, PDFScale, QED, alpha3))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,-1,0,0,101,102,0,0,Ea,Ea,0,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,-1,0,0,102,103,0,0,-Eb,Eb,0,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,1,1,2,101,104,px1,py1,pz1,E1,0,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,1,1,2,104,103,px2,py2,pz2,E2,0,0,9))
				fh.write("#pdf  {}  {}  {}  {}  {}  {}  {}\n".format(21,21,xa,xb,PDFScale,DistFunc.xfxQ(21, xa, PDFScale),DistFunc.xfxQ(21, xb, PDFScale)))
				fh.write("</event>\n")
				if(numberofevents == 500): print("\n{} Events are Generated\n".format(numberofevents))
				if(numberofevents%50   == 0): print(".")
				if(numberofevents%1000 == 0): print("\n{} Events are Generated\n".format(numberofevents))
				break

			elif(InteractionIndex < 19):

				IDQ = -(InteractionIndex - 12)
				args = (M, 21, IDQ)

				RInvMass = random.uniform(0, ArrayOfMaxes[InteractionIndex])
				if(RInvMass > FuncAtMass[InteractionIndex]):
					break

				MaximumEvaluation = MonteCarlogqbargqbar(MaxPar, *args)
				RandomEval = random.uniform(0, MaximumEvaluation)
				ActualParameters = [RandomY, Randomy]
				ActualEval = MonteCarlogqbargqbar(ActualParameters, *args)

				Qmass = QuarkMasses[InteractionIndex - 13]

				if(RandomEval > ActualEval):
					break

				numberofevents += 1

				# Calculation of the kinematic variables

				y1 = RandomY + Randomy
				y2 = RandomY - Randomy

				xa = math.sqrt(tau)*math.exp( RandomY)
				xb = math.sqrt(tau)*math.exp(-RandomY)

				Ea = xa*EBMUP
				Eb = xb*EBMUP

				phi = random.uniform(0, 2*math.pi)

				pt = 0.5*M/math.cosh(Randomy)

				px1 = pt*math.cos(phi)
				py1 = pt*math.sin(phi)

				px2 = -px1
				py2 = -py1

				E1 = pt*math.cosh(y1)
				E2 = pt*math.cosh(y2)

				pz1 = E1*math.tanh(y1)
				pz2 = E2*math.tanh(y2)

				weight = FuncAtMass[InteractionIndex]*(10**9)/2.56819*DeltaMass

				CrossSectiongqbargqbar += weight
				CrossSectionErrorgqbargqbar += weight**2

				# Filling the LHE file

				fh.write("<event>\n")
				fh.write("{}  {}  {}  {}  {} {}\n".format(4, 4, weight, PDFScale, QED, alpha3))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,-1,0,0,101,102,0,0,Ea,Ea,0,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(IDQ,-1,0,0,0,101,0,0,-Eb,Eb,Qmass,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,1,1,2,103,102,px1,py1,pz1,E1,0,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(IDQ,1,1,2,0,103,px2,py2,pz2,E2,Qmass,0,9))
				fh.write("#pdf  {}  {}  {}  {}  {}  {}  {}\n".format(21,IDQ,xa,xb,PDFScale,DistFunc.xfxQ(21, xa, PDFScale),DistFunc.xfxQ(IDQ, xb, PDFScale)))
				fh.write("</event>\n")
				if(numberofevents == 500): print("\n{} Events are Generated\n".format(numberofevents))
				if(numberofevents%50   == 0): print(".")
				if(numberofevents%1000 == 0): print("\n{} Events are Generated\n".format(numberofevents))
				break

			elif(InteractionIndex == 19):

				IDQ = random.randint(1,6)
				args = (M, 21, 21)

				RInvMass = random.uniform(0, ArrayOfMaxes[InteractionIndex])
				if(RInvMass > FuncAtMass[InteractionIndex]):
					break

				MaximumEvaluation = MonteCarloggqqbar(MaxPar, *args)
				RandomEval = random.uniform(0, MaximumEvaluation)
				ActualParameters = [RandomY, Randomy]
				ActualEval = MonteCarloggqqbar(ActualParameters, *args)
				Qmass = QuarkMasses[IDQ-1]

				if(RandomEval > ActualEval):
					break

				numberofevents += 1

				# Calculation of the kinematic variables

				y1 = RandomY + Randomy
				y2 = RandomY - Randomy

				xa = math.sqrt(tau)*math.exp( RandomY)
				xb = math.sqrt(tau)*math.exp(-RandomY)

				Ea = xa*EBMUP
				Eb = xb*EBMUP

				phi = random.uniform(0, 2*math.pi)

				pt = 0.5*M/math.cosh(Randomy)

				px1 = pt*math.cos(phi)
				py1 = pt*math.sin(phi)

				px2 = -px1
				py2 = -py1

				E1 = pt*math.cosh(y1)
				E2 = pt*math.cosh(y2)

				pz1 = E1*math.tanh(y1)
				pz2 = E2*math.tanh(y2)

				weight = FuncAtMass[InteractionIndex]*(10**9)/2.56819*DeltaMass

				CrossSectionggqqbar += weight
				CrossSectionErrorggqqbar += weight**2

				# Filling the LHE file

				fh.write("<event>\n")
				fh.write("{}  {}  {}  {}  {} {}\n".format(4, 2, weight, PDFScale, QED, alpha3))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,-1,0,0,101,102,0,0,Ea,Ea,0,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,-1,0,0,103,101,0,0,-Eb,Eb,0,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(-IDQ,1,1,2,0,102,px1,py1,pz1,E1,Qmass,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(IDQ,1,1,2,103,0,px2,py2,pz2,E2,Qmass,0,9))
				fh.write("#pdf  {}  {}  {}  {}  {}  {}  {}\n".format(21,21,xa,xb,PDFScale,DistFunc.xfxQ(21, xa, PDFScale),DistFunc.xfxQ(21, xb, PDFScale)))
				fh.write("</event>\n")
				if(numberofevents == 500): print("\n{} Events are Generated\n".format(numberofevents))
				if(numberofevents%50   == 0): print(".")
				if(numberofevents%1000 == 0): print("\n{} Events are Generated\n".format(numberofevents))
				break

			elif(InteractionIndex == 20):

				args = (M, 21, 21)

				RInvMass = random.uniform(0, ArrayOfMaxes[InteractionIndex])
				if(RInvMass > FuncAtMass[InteractionIndex]):
					break

				MaximumEvaluation = MonteCarloggGamma(MaxPar, *args)
				RandomEval = random.uniform(0, MaximumEvaluation)
				ActualParameters = [RandomY, Randomy]
				ActualEval = MonteCarloggGamma(ActualParameters, *args)

				if(RandomEval > ActualEval):
					break

				numberofevents += 1

				# Calculation of the kinematic variables

				y1 = RandomY + Randomy
				y2 = RandomY - Randomy

				xa = math.sqrt(tau)*math.exp( RandomY)
				xb = math.sqrt(tau)*math.exp(-RandomY)

				Ea = xa*EBMUP
				Eb = xb*EBMUP

				phi = random.uniform(0, 2*math.pi)

				pt = 0.5*M/math.cosh(Randomy)

				px1 = pt*math.cos(phi)
				py1 = pt*math.sin(phi)

				px2 = -px1
				py2 = -py1

				E1 = pt*math.cosh(y1)
				E2 = pt*math.cosh(y2)

				pz1 = E1*math.tanh(y1)
				pz2 = E2*math.tanh(y2)

				weight = FuncAtMass[InteractionIndex]*(10**9)/2.56819*DeltaMass

				CrossSectionggGamma += weight
				CrossSectionErrorggGamma += weight**2

				# Filling the LHE file

				fh.write("<event>\n")
				fh.write("{}  {}  {}  {}  {} {}\n".format(4, 6, weight, PDFScale, QED, alpha3))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,-1,0,0,101,102,0,0,Ea,Ea,0,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,-1,0,0,103,101,0,0,-Eb,Eb,0,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,1,1,2,103,102,px1,py1,pz1,E1,0,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(22,1,1,2,0,0,px2,py2,pz2,E2,0,0,9))
				fh.write("#pdf  {}  {}  {}  {}  {}  {}  {}\n".format(21,21,xa,xb,PDFScale,DistFunc.xfxQ(21, xa, PDFScale),DistFunc.xfxQ(21, xb, PDFScale)))
				fh.write("</event>\n")
				if(numberofevents == 500): print("\n{} Events are Generated\n".format(numberofevents))
				if(numberofevents%50   == 0): print(".")
				if(numberofevents%1000 == 0): print("\n{} Events are Generated\n".format(numberofevents))
				break

			elif(InteractionIndex < 27):

				IDQ = InteractionIndex - 20
				args = (M, 21, IDQ)

				RInvMass = random.uniform(0, ArrayOfMaxes[InteractionIndex])
				if(RInvMass > FuncAtMass[InteractionIndex]):
					break

				MaximumEvaluation = MonteCarlogqGamma(MaxPar, *args)
				RandomEval = random.uniform(0, MaximumEvaluation)
				ActualParameters = [RandomY, Randomy]
				ActualEval = MonteCarlogqGamma(ActualParameters, *args)
				Qmass = QuarkMasses[IDQ-1]

				if(RandomEval > ActualEval):
					break

				numberofevents += 1

				# Calculation of the kinematic variables

				y1 = RandomY + Randomy
				y2 = RandomY - Randomy

				xa = math.sqrt(tau)*math.exp( RandomY)
				xb = math.sqrt(tau)*math.exp(-RandomY)

				Ea = xa*EBMUP
				Eb = xb*EBMUP

				phi = random.uniform(0, 2*math.pi)

				pt = 0.5*M/math.cosh(Randomy)

				px1 = pt*math.cos(phi)
				py1 = pt*math.sin(phi)

				px2 = -px1
				py2 = -py1

				E1 = pt*math.cosh(y1)
				E2 = pt*math.cosh(y2)

				pz1 = E1*math.tanh(y1)
				pz2 = E2*math.tanh(y2)

				weight = FuncAtMass[InteractionIndex]*(10**9)/2.56819*DeltaMass

				CrossSectiongqGamma += weight
				CrossSectionErrorgqGamma += weight**2

				# Filling the LHE file

				fh.write("<event>\n")
				fh.write("{}  {}  {}  {}  {} {}\n".format(4, 7, weight, PDFScale, QED, alpha3))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(21,-1,0,0,102,101,0,0,Ea,Ea,0,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(IDQ,-1,0,0,101,0,0,0,-Eb,Eb,Qmass,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(IDQ,1,1,2,102,0,px1,py1,pz1,E1,Qmass,0,9))
				fh.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(22,1,1,2,103,0,px2,py2,pz2,E2,0,0,9))
				fh.write("#pdf  {}  {}  {}  {}  {}  {}  {}\n".format(21,21,xa,xb,PDFScale,DistFunc.xfxQ(21, xa, PDFScale),DistFunc.xfxQ(21, xb, PDFScale)))
				fh.write("</event>\n")
				if(numberofevents == 500): print("\n{} Events are Generated\n".format(numberofevents))
				if(numberofevents%50   == 0): print(".")
				if(numberofevents%1000 == 0): print("\n{} Events are Generated\n".format(numberofevents))
				break

	fh.write("</LesHouchesEvents>")
	fh.close()

	print("\nFinalizing...\n")

	# Calculating the error in the cross-sections
	ggggStandardError           = math.sqrt(math.fabs(Counter*CrossSectionErrorgggg - CrossSectiongggg**2))/math.sqrt(Counter)
	gqgqStandardError           = math.sqrt(math.fabs(Counter*CrossSectionErrorgqgq - CrossSectiongqgq**2))/math.sqrt(Counter)
	ggqqbarStandardError        = math.sqrt(math.fabs(Counter*CrossSectionErrorggqqbar - CrossSectionggqqbar**2))/math.sqrt(Counter)
	gqbargqbarStandardError     = math.sqrt(math.fabs(Counter*CrossSectionErrorgqbargqbar - CrossSectiongqbargqbar**2))/math.sqrt(Counter)
	qqbarggStandardError        = math.sqrt(math.fabs(Counter*CrossSectionErrorqqbargg - CrossSectionqqbargg**2))/math.sqrt(Counter)
	ggGammaStandardError        = math.sqrt(math.fabs(Counter*CrossSectionErrorggGamma - CrossSectionggGamma**2))/math.sqrt(Counter)
	gqGammaStandardError        = math.sqrt(math.fabs(Counter*CrossSectionErrorgqGamma - CrossSectiongqGamma**2))/math.sqrt(Counter)

	TotalCrossSection = CrossSectiongggg + CrossSectionggqqbar + CrossSectiongqgq + CrossSectiongqbargqbar + CrossSectionqqbargg + CrossSectionggGamma + CrossSectiongqGamma

	# Printing the Calculated cross-sections on the screen
	print("{0:<20}{1:<25}{2:<30}".format('SubProcess', 'ID', 'Cross-Section (pb)'))

	if(gg2gg11 == 1):
		print("{0:<20}{1:<25}{2:<30}".format('gg2gg', '1', CrossSectiongggg))
	if(gg2qqbar11 == 1):
		print("{0:<20}{1:<25}{2:<30}".format('gg2qqbar', '2', CrossSectionggqqbar))
	if(gq2gq11 == 1):
		print("{0:<20}{1:<25}{2:<30}".format('gq2gq', '3', CrossSectiongqgq))
	if(gqbar2gqbar11 == 1):
		print("{0:<20}{1:<25}{2:<30}".format('gqbar2gqbar', '4', CrossSectiongqbargqbar))
	if(qqbar2gg11 == 1):
		print("{0:<20}{1:<25}{2:<30}".format('qqbar2gg', '5', CrossSectionqqbargg))
	if(gg2gGamma11 == 1):
		print("{0:<20}{1:<25}{2:<30}".format('gg2gGamma', '6', CrossSectionggGamma))
	if(gq2qGamma11 == 1):
		print("{0:<20}{1:<25}{2:<30}".format('gq2qGamma', '7', CrossSectiongqGamma))


	print("\n")
	print("{0:<45}{1:<65}".format('Total Cross-Section (pb)', TotalCrossSection))

	# Writing the input parameters to the LHE file

	fh = open("begin.lhe", "w")

	fh.write("<LesHouchesEvents version=\"%.1f\">\n" % 1.0)
	fh.write("<!--\n")
	fh.write("  File Written by STRINGS - Version 1.00, October 2018\n")
	fh.write("  Generator by Pourya Vakilipourtakalou and Douglas M. Gingrich\n")
	fh.write("  Parameters:\n")
	fh.write("  {0:<10}{1:<22}{2:<30}".format(RandGenSEED, 'RandGenSEED', 'Seed for the the Random Number Generator'))
	fh.write("\n")
	fh.write("  {0:<10}{1:<22}{2:<30}".format(COME, 'COME', 'Centre of Mass Energy (GeV)'))
	fh.write("\n")
	fh.write("  {0:<10}{1:<22}{2:<30}".format(Counter, 'Number', 'Number of Events'))
	fh.write("\n")
	fh.write("  {0:<10}{1:<22}{2:<30}".format(Ms, 'Ms', 'String Scale (GeV)'))
	fh.write("\n")
	fh.write("  {0:<10}{1:<22}{2:<30}".format(MinM, 'MinMass', 'Minimum Invariant Mass (GeV)'))
	fh.write("\n")
	fh.write("  {0:<10}{1:<22}{2:<30}".format(MaxM, 'MaxMass', 'Maximum Invariant Mass (GeV)'))
	fh.write("\n")
	fh.write("  {0:<10}{1:<22}{2:<30}".format(ymax, 'ymax', 'Upper Bound for the Rapidity of the Outgoing Partons'))
	fh.write("\n")
	fh.write("  {0:<10}{1:<22}{2:<30}".format(PDFSet, 'PDFSet', 'PDF Set of the LHAPDF'))
	fh.write("\n")
	fh.write("  {0:<10}{1:<22}{2:<30}".format(PDFScale, 'PDFScale', 'Scale at Which the PDF Set is Evaluated (GeV)'))
	fh.write("\n")
	if(Coupling == -1):
		fh.write("  {0:<10}{1:<22}{2:<30}".format(Coupling, 'Coupling', 'Running Coupling Constant (alpha_s without 4*pi Factor)'))
		fh.write("\n")
		fh.write("  {0:<10}{1:<22}{2:<30}".format(CouplingScale, 'CouplingScale', 'Scale at Which the Running Coupling is Calculated (GeV)'))
		fh.write("\n")
	else:
		fh.write("  {0:<10}{1:<22}{2:<30}".format(Coupling, 'Coupling', 'Coupling Constant (alpha_s without 4*pi Factor)'))
		fh.write("\n")
	
	fh.write("  {0:<10}{1:<22}{2:<30}".format(dmass, 'dMass', 'Mass of the Down Quark (GeV)'))
	fh.write("\n")
	fh.write("  {0:<10}{1:<22}{2:<30}".format(umass, 'dMass', 'Mass of the Up Quark (GeV)'))
	fh.write("\n")
	fh.write("  {0:<10}{1:<22}{2:<30}".format(smass, 'dMass', 'Mass of the Strange Quark (GeV)'))
	fh.write("\n")
	fh.write("  {0:<10}{1:<22}{2:<30}".format(cmass, 'dMass', 'Mass of the Charm Quark (GeV)'))
	fh.write("\n")
	fh.write("  {0:<10}{1:<22}{2:<30}".format(bmass, 'dMass', 'Mass of the Bottom Quark (GeV)'))
	fh.write("\n")
	fh.write("  {0:<10}{1:<22}{2:<30}".format(tmass, 'dMass', 'Mass of the Top Quark (GeV)'))
	fh.write("\n")

	if(QCDCoeff11 == 1):
		fh.write("  {0:<10}{1:<22}{2:<30}".format(QCDCoeff, 'QCDCoeff', '(Enabled)  Production of QCD tree-level diparton'))
		fh.write("\n")
	else:
		fh.write("  {0:<10}{1:<22}{2:<30}".format(QCDCoeff, 'QCDCoeff', '(Disabled) Production of QCD tree-level diparton'))
		fh.write("\n")

	if(stringCoeff11 == 1):
		fh.write("  {0:<10}{1:<22}{2:<30}".format(stringCoeff, 'FirstStringCoeff', '(Enabled)  Production of First  String Resonance ( 2 --> 2 Partonic Scattering and 2-Parton --> Parton-Gamma Scattering )'))
		fh.write("\n")
	else:
		fh.write("  {0:<10}{1:<22}{2:<30}".format(stringCoeff, 'FirstStringCoeff', '(Disabled) Production of First  String Resonance  ( 2 --> 2 Partonic Scattering and 2-Parton --> Parton-Gamma Scattering )'))
		fh.write("\n")

	if(SecondStringCoeff11 == 1):
		fh.write("  {0:<10}{1:<22}{2:<30}".format(SecondStringCoeff, 'SecondStringCoeff', '(Enabled)  Production of Second String Resonance ( 2 --> 2 Partonic Scattering )'))
		fh.write("\n")
	else:
		fh.write("  {0:<10}{1:<22}{2:<30}".format(SecondStringCoeff, 'SecondStringCoeff', '(Disabled) Production of Second String Resonance ( 2 --> 2 Partonic Scattering )'))
		fh.write("\n")

	if(gg2gg11 == 1):
		fh.write("  {0:<10}{1:<22}{2:<30}".format(gg2gg, 'gg2gg', '(Enabled)  gg --> gg Subprocess       (ID = 1)'))
		fh.write("\n")
	else:
		fh.write("  {0:<10}{1:<22}{2:<30}".format(gg2gg, 'gg2gg', '(Disabled) gg --> gg Subprocess       (ID = 1)'))
		fh.write("\n")

	if(gg2qqbar11 == 1):
		fh.write("  {0:<10}{1:<22}{2:<30}".format(gg2qqbar, 'gg2qqbar', '(Enabled)  gg --> qqbar Subprocess    (ID = 2)'))
		fh.write("\n")
	else:
		fh.write("  {0:<10}{1:<22}{2:<30}".format(gg2qqbar, 'gg2qqbar', '(Disabled) gg --> qqbar Subprocess    (ID = 2)'))
		fh.write("\n")

	if(gq2gq11 == 1):
		fh.write("  {0:<10}{1:<22}{2:<30}".format(gq2gq, 'gq2gq', '(Enabled)  gq --> gq Subprocess       (ID = 3)'))
		fh.write("\n")
	else:
		fh.write("  {0:<10}{1:<22}{2:<30}".format(gq2gq, 'gq2gq', '(Disabled) gq --> gq Subprocess       (ID = 3)'))
		fh.write("\n")

	if(gqbar2gqbar11 == 1):
		fh.write("  {0:<10}{1:<22}{2:<30}".format(gqbar2gqbar, 'gqbar2gqbar', '(Enabled)  gqbar --> gqbar Subprocess (ID = 4)'))
		fh.write("\n")
	else:
		fh.write("  {0:<10}{1:<22}{2:<30}".format(gqbar2gqbar, 'gqbar2gqbar', '(Disabled) gqbar --> gqbar Subprocess (ID = 4)'))
		fh.write("\n")

	if(qqbar2gg11 == 1):
		fh.write("  {0:<10}{1:<22}{2:<30}".format(qqbar2gg, 'qqbar2gg', '(Enabled)  qqbar --> gg Subprocess    (ID = 5)'))
		fh.write("\n")
	else:
		fh.write("  {0:<10}{1:<22}{2:<30}".format(qqbar2gg, 'qqbar2gg', '(Disabled) qqbar --> gg Subprocess    (ID = 5)'))
		fh.write("\n")

	if(gg2gGamma11 == 1):
		fh.write("  {0:<10}{1:<22}{2:<30}".format(gg2gGamma, 'gg2gGamma', '(Enabled)  gg --> gGamma Subprocess   (ID = 6)'))
		fh.write("\n")
	else:
		fh.write("  {0:<10}{1:<22}{2:<30}".format(gg2gGamma, 'gg2gGamma', '(Disabled) gg --> gGamma Subprocess   (ID = 6)'))
		fh.write("\n")

	if(gq2qGamma11 == 1):
		fh.write("  {0:<10}{1:<22}{2:<30}".format(gq2qGamma, 'gq2qGamma', '(Enabled)  gq --> qGamma Subprocess   (ID = 7)'))
		fh.write("\n")
	else:
		fh.write("  {0:<10}{1:<22}{2:<30}".format(gq2qGamma, 'gq2qGamma', '(Disabled) gq --> qGamma Subprocess   (ID = 7)'))
		fh.write("\n")

	print("\nSuccessfully Generated {} Events".format(numberofevents))
	print("Events are saved in STRINGSfile.lhe\n")
	print("Thanks for using STRINGS-1.00\n")

	# Updating the LHE file with the calculated cross-sections


	fh.write("-->\n")
	fh.write("<init>\n")
	fh.write("  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(IDBMUP, IDBMUP, EBMUP, EBMUP, PDFGUP, PDFGUP, PDFSUP, PDFSUP, IDWTUP, NPRUP))
	if(gg2gg11 == 1):
		fh.write("  {}  {}  {}  {}\n".format(CrossSectiongggg, ggggStandardError, 1, 1))
	if(gg2qqbar11 == 1):
		fh.write("  {}  {}  {}  {}\n".format(CrossSectionggqqbar, ggqqbarStandardError, 1, 2))
	if(gq2gq11 == 1):
		fh.write("  {}  {}  {}  {}\n".format(CrossSectiongqgq, gqgqStandardError, 1, 3))
	if(gqbar2gqbar11 == 1):
		fh.write("  {}  {}  {}  {}\n".format(CrossSectiongqbargqbar, gqbargqbarStandardError, 1, 4))
	if(qqbar2gg11 == 1):
		fh.write("  {}  {}  {}  {}\n".format(CrossSectionqqbargg, qqbarggStandardError, 1, 5))
	if(gg2gGamma11 == 1):
		fh.write("  {}  {}  {}  {}\n".format(CrossSectionggGamma, ggGammaStandardError, 1, 6))
	if(gq2qGamma11 == 1):
		fh.write("  {}  {}  {}  {}\n".format(CrossSectiongqGamma, gqGammaStandardError, 1, 7))
	fh.write("</init>\n")

	fh.close()

	call("cat events.lhe >> begin.lhe",shell=True)

	call("rm events.lhe",shell=True)
	call("mv begin.lhe STRINGSfile.lhe",shell=True)



# Definition of the input parameters with their default values. 
# Two versions of the arguments are assigned for each input that can be used to change the parameter.


parser = argparse.ArgumentParser(description='Parameters Passed to the STRINGS Generator')

parser.add_argument('-v', '--version', action='version', version='STRINGS-1.00')

parser.add_argument('-RandGenSEED', '--RandGenSEEDvalue', default=123456, help='Seed for Random Number Generator',type=float)

parser.add_argument('-COME', '--COMEvalue', default=13000, help='Centre of Mass Energy (GeV)',type=float)

parser.add_argument('-Ms', '--Msvalue', default=7000, help='String Scale of the String Theory (GeV)',type=float)

parser.add_argument('-ymax', '--ymaxvalue', default=2.5, help='Upper Limit for the Rapidity of the Outgoing Partons',type=float)

parser.add_argument('-Number', '--Numbervalue', default=10000, help='Number of Events',type=int)

parser.add_argument('-MinMass', '--MinMassvalue', default=6000, help='Lower Bound of the Invariant Mass (GeV)',type=float)

parser.add_argument('-MaxMass', '--MaxMassvalue', default=8000, help='Upper Bound of the Invariant Mass (GeV)',type=float)

parser.add_argument('-FirstStringCoeff', '--FirstStringCoeffvalue', default="true", help='First String Resonance',type=str)

parser.add_argument('-SecondStringCoeff', '--SecondStringCoeffvalue', default="false", help='Second String Resonance',type=str)

parser.add_argument('-QCDCoeff', '--QCDCoeffvalue', default="false", help='QCD Tree-Level',type=str)

parser.add_argument('-PDFSet', '--PDFSetvalue', default="cteq6l1", help='Parton Distribution Function Set of the LHAPDF',type=str)

parser.add_argument('-Coupling', '--Couplingvalue', default=-1, help='Coupling Constant',type=float)

parser.add_argument('-CouplingScale', '--CouplingScalevalue', default=-22, help='Scale at Which the Running Coupling Constant is Calculated(GeV)',type=float)

parser.add_argument('-PDFScale', '--PDFScalevalue', default=-22, help='Scale at which the Parton Distribution Function is Evaluated (GeV)',type=float)

parser.add_argument('-dmass', '--dmassvalue', default=5e-3, help='Mass of the Down Quark (GeV)',type=float)

parser.add_argument('-umass', '--umassvalue', default=2e-3, help='Mass of the Up Quark (GeV)',type=float)

parser.add_argument('-smass', '--smassvalue', default=1e-3, help='Mass of the Strange Quark (GeV)',type=float)

parser.add_argument('-cmass', '--cmassvalue', default=1.27, help='Mass of the Charm Quark (GeV)',type=float)

parser.add_argument('-bmass', '--bmassvalue', default=4.4, help='Mass of the Bottom Quark (GeV)',type=float)

parser.add_argument('-tmass', '--tmassvalue', default=172, help='Mass of the Top Quark (GeV)',type=float)

parser.add_argument('-gg2gg', '--gg2ggvalue', default="true", help='gg -> gg Subprocess',type=str)

parser.add_argument('-gg2qqbar', '--gg2qqbarvalue', default="true", help='gg -> qqbar Subprocess',type=str)

parser.add_argument('-gq2gq', '--gq2gqvalue', default="true", help='gq -> gq Subprocess',type=str)

parser.add_argument('-gqbar2gqbar', '--gqbar2gqbarvalue', default="true", help='gqbar -> gqbar Subprocess',type=str)

parser.add_argument('-qqbar2gg', '--qqbar2ggvalue', default="true", help='qqbar -> gg Subprocess',type=str)

parser.add_argument('-gg2gGamma', '--gg2gGammavalue', default="false", help='gg -> gGamma Subprocess',type=str)

parser.add_argument('-gq2qGamma', '--gq2qGammavalue', default="false", help='gq -> qGamma Subprocess',type=str)

args = parser.parse_args()

print("########################################################")
print("################# STRINGS-VERSION 1.00 #################")
print("########################################################")
print("####                                                ####")
print("#### October 2018                                   ####")
print("####                                                ####")
print("#### Authors:                                       ####")
print("####                                                ####")
print("####\t Pourya Vakilipourtakalou                   ####")
print("####\t Email:                                     ####") 
print("####\t \t vakilipo@ualberta.ca               ####")
print("####\t \t pourya.vakilipourtakalou@cern.ch   ####")
print("####                                                ####")
print("####\t Douglas M. Gingrich                        ####")
print("####\t Email:                                     ####")
print("####\t \t gingrich@ualberta.ca               ####")
print("####                                                ####")
print("########################################################")
print("########################################################")


# A checker is set to make sure all of the input variables are correctly assigned
checker = 0

if(args.QCDCoeffvalue == 'false' or args.QCDCoeffvalue == 'False'):
	if(args.FirstStringCoeffvalue == 'false' or args.FirstStringCoeffvalue == 'False'):
		if(args.SecondStringCoeffvalue == 'false' or args.SecondStringCoeffvalue == 'False'):
			checker = 1
			print("Error: Type of Event Generation is not Specified\n")


if(args.gg2ggvalue == 'false' or args.gg2ggvalue == 'False'):
	if(args.gg2qqbarvalue == 'false' or args.gg2qqbarvalue == 'False'):
		if(args.gq2gqvalue == 'false' or args.gq2gqvalue == 'False'):
			if(args.gqbar2gqbarvalue == 'false' or args.gqbar2gqbarvalue == 'False'):
				if(args.qqbar2ggvalue == 'false' or args.qqbar2ggvalue == 'False'):
					if(args.gg2gGammavalue == 'false' or args.gg2gGammavalue == 'False'):
						if(args.gq2qGammavalue == 'false' or args.gq2qGammavalue == 'False'):
							checker = 1
							print("Error: No Subprocess is Turned On\n")
if(args.gg2ggvalue == 'false' or args.gg2ggvalue == 'False'):
	if(args.gg2qqbarvalue == 'false' or args.gg2qqbarvalue == 'False'):
		if(args.gq2gqvalue == 'false' or args.gq2gqvalue == 'False'):
			if(args.gqbar2gqbarvalue == 'false' or args.gqbar2gqbarvalue == 'False'):
				if(args.qqbar2ggvalue == 'false' or args.qqbar2ggvalue == 'False'):
					if(args.gg2gGammavalue == 'true' or args.gg2gGammavalue == 'True'):
						if(args.QCDCoeffvalue == 'true' or args.QCDCoeffvalue == 'True'):
							checker = 1
							print("Error: QCD Diparton Cannot be Produced Through gg->gGamma Scattering\nPlease Turn On Either of the 2->2 Partonic Scatterings")
						if(args.SecondStringCoeffvalue == 'true' or args.SecondStringCoeffvalue == 'True'):
							checker = 1
							print("Error: Second string cannot be produced through gg->gGamma scattering\nPlease Turn On Either of the 2->2 Partonic Scatterings")
					if(args.gq2qGammavalue == 'true' or args.gq2qGammavalue == 'True'):
						if(args.QCDCoeffvalue == 'true' or args.QCDCoeffvalue == 'True'):
							checker = 1
							print("Error: QCD Diparton Cannot be Produced Through gq->qGamma Scattering\nPlease Turn On Either of the 2->2 Partonic Scatterings")
						if(args.SecondStringCoeffvalue == 'true' or args.SecondStringCoeffvalue == 'True'):
							checker = 1
							print("Error: Second string cannot be produced through gq->qGamma scattering\nPlease Turn On Either of the 2->2 Partonic Scatterings")
if(args.COMEvalue <= 0):
	print("Error: Invalid Value for Centre of Mass Energy\nIt Should be Positive\n")
	checker = 1

if(args.Msvalue <= 0):
	print("Error: Invalid Value for the String Scale\nIt Should be Positive\n")
	checker = 1

if(args.ymaxvalue == -1):
	args.ymaxvalue = float("inf")
elif(args.ymaxvalue < 0):
	print("Error: Invalid Value for the Upper Bound of the Rapidity\nIt Should be Positive or for the Upper Bound to be Infinity, Use -1\n")
	checker = 1

if(args.Numbervalue <= 0):
	print("Error: Invalid Value for the Number of Generated Events\nIt Should be Positive\n")
	checker = 1

if(args.MinMassvalue <= 0):
	print("Error: Invalid Value for the Lower Bound of the Invariant Mass\nIt Should be Positive\n")
	checker = 1

if(args.MaxMassvalue <= 0):
	print("Error: Invalid Value for the Upper Bound of the Invariant Mass\nIt Should be Positive\n")
	checker = 1
if(args.MaxMassvalue > args.COMEvalue):
	print("Error: Invalid Value for the Upper Bound of the Invariant Mass\nIt Should Not be larger than the Centre of Mass Energy\n")
	checker=1
elif(args.MaxMassvalue <= args.MinMassvalue):
	print("Error: Invalid Value for the Upper Bound of the Invariant Mass\nIt Should be greater than the Lower Bound of the Invariant Mass\n")
	checker = 1
elif(args.MaxMassvalue > args.COMEvalue):
	print("Error: Invalid Value for the Upper Bound of the Invariant Mass\nIt Should not be greater than the Centre of Mass Energy\n")
	checker = 1

if(args.FirstStringCoeffvalue != "True" and args.FirstStringCoeffvalue != "true" and args.FirstStringCoeffvalue != "false" and args.FirstStringCoeffvalue != "False"):
	print("Error: Invalid Value for the FirstStringCoeff\nIt is a boolean and should be either of the followings: True, true, False, false\n")
	checker = 1

if(args.SecondStringCoeffvalue != "True" and args.SecondStringCoeffvalue != "true" and args.SecondStringCoeffvalue != "false" and args.SecondStringCoeffvalue != "False"):
	print("Error: Invalid Value for the SecondStringCoeff\nIt is a boolean and should be either of the followings: True, true, False, false\n")
	checker = 1

if(args.QCDCoeffvalue != "True" and args.QCDCoeffvalue != "true" and args.QCDCoeffvalue != "False" and args.QCDCoeffvalue != "false"):
	print("Error: Invalid Value for the QCDCoeff\nIt is a boolean and should be either of the followings: True, true, False, false\n")
	checker = 1

if(args.PDFScalevalue <= 0 and args.PDFScalevalue != -22):
	print("Error: Invalid Value for the Scale at Which PDFs are Determined\nIt Should be Positive\n")
	checker = 1

if(args.Couplingvalue == -1):
	args.Couplingvalue = -1
elif(args.Couplingvalue < 0):
	print("Error: Invalid Value for the Coupling Constant\nIt Should be positive or for running coupling use -1\n")
	checker = 1

if(args.CouplingScalevalue <= 0 and args.CouplingScalevalue != -22):
	print("Error: Invalid Value for the Scale at Which Running Coupling is Measured\nIt Should be Positive\n")
	checker = 1

if(args.dmassvalue < 0):
	print("Error: Invalid Value for the Mass of the Down Quark\nIt Should be Positive\n")
	checker = 1

if(args.umassvalue < 0):
	print("Error: Invalid Value for the Mass of the UP Quark\nIt Should be Positive\n")
	checker = 1

if(args.smassvalue < 0):
	print("Error: Invalid Value for the Mass of the Strange Quark\nIt Should be Positive\n")
	checker = 1

if(args.cmassvalue < 0):
	print("Error: Invalid Value for the Mass of the Charm Quark\nIt Should be Positive\n")
	checker = 1

if(args.bmassvalue < 0):
	print("Error: Invalid Value for the Mass of the Bottom Quark\nIt Should be Positive\n")
	checker = 1

if(args.tmassvalue < 0):
	print("Error: Invalid Value for the Mass of the Top Quark\nIt Should be Positive\n")
	checker = 1

if(args.gg2ggvalue != "True" and args.gg2ggvalue != "true" and args.gg2ggvalue != "False" and args.gg2ggvalue != "false"):
	print("Error: Invalid Value; In order to enable or disable gg --> gg subprocess, set gg2gg to true or false, respectively.\n")
	checker = 1

if(args.gg2qqbarvalue != "True" and args.gg2qqbarvalue != "true" and args.gg2qqbarvalue != "False" and args.gg2qqbarvalue != "false"):
	print("Error: Invalid Value; In order to enable or disable gg --> qqbar subprocess, set gg2qqbar to true or false, respectively.\n")
	checker = 1

if(args.gq2gqvalue != "True" and args.gq2gqvalue != "true" and args.gq2gqvalue != "False" and args.gq2gqvalue != "false"):
	print("Error: Invalid Value; In order to enable or disable gq --> gq subprocess, set gq2gq to true or false, respectively.\n")
	checker = 1

if(args.gqbar2gqbarvalue != "True" and args.gqbar2gqbarvalue != "true" and args.gqbar2gqbarvalue != "False" and args.gqbar2gqbarvalue != "false"):
	print("Error: Invalid Value; In order to enable or disable gqbar --> gqbar subprocess, set gqbar2gqbar to true or false, respectively.\n")
	checker = 1

if(args.qqbar2ggvalue != "True" and args.qqbar2ggvalue != "true" and args.qqbar2ggvalue != "false" and args.qqbar2ggvalue != "False"):
	print("Error: Invalid Value; In order to enable or disable qqbar --> gg subprocess, set qqbar2gg to true or false, respectively.\n")
	checker = 1

if(args.gg2gGammavalue != "True" and args.gg2gGammavalue != "true" and args.gg2gGammavalue != "False" and args.gg2gGammavalue != "false"):
	print("Error: Invalid Value; In order to enable or disable gg --> gGamma subprocess, set gg2gGamma to true or false, respectively.\n")
	checker = 1

if(args.gq2qGammavalue != "True" and args.gq2qGammavalue != "true" and args.gq2qGammavalue != "false" and args.gq2qGammavalue != "False"):
	print("Error: Invalid Value; In order to enable or disable gq --> qGamma subprocess, set gq2qGamma to true or false, respectively.\n")
	checker = 1

# If all of the input parameters make sense, the main function is called
if(checker == 0):
	main(args.RandGenSEEDvalue, args.COMEvalue, args.Msvalue, args.ymaxvalue, args.Numbervalue, args.MinMassvalue, args.MaxMassvalue, args.FirstStringCoeffvalue, args.SecondStringCoeffvalue, args.QCDCoeffvalue, args.PDFSetvalue, args.PDFScalevalue, args.Couplingvalue, args.CouplingScalevalue, args.dmassvalue, args.umassvalue, args.smassvalue, args.cmassvalue, args.bmassvalue, args.tmassvalue, args.gg2ggvalue, args.gg2qqbarvalue, args.gq2gqvalue, args.gqbar2gqbarvalue, args.qqbar2ggvalue, args.gg2gGammavalue, args.gq2qGammavalue)