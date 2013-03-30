#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
from numpy import *
from matplotlib.pyplot import *
from scipy.signal import argrelextrema
from scipy.interpolate import *

def generateSignal(t):
	return sin(2 * pi * t) + .5 * cos(4 * pi * t)

def getEnvelope(x, t, comp):
	num = argrelextrema(x, comp)[0]
	tck = splrep(t[num], x[num])
	return splev(t, tck)[70:-70]
	
def getEnvelopeMean(signal, t):
	envelopeMin = getEnvelope(signal, t, less)
	envelopeMax = getEnvelope(signal, t, greater)
	return .5 * (envelopeMin + envelopeMax)	

def recursieveGetMode(signal, t):
	someMode = signal[70:-70]
	mean = getEnvelopeMean(someMode, t)
	while (mean >= 0.4).any():
		someMode -= mean
		mean = getEnvelopeMean(someMode, t)
	return someMode			
		

def main():
	
	t = arange(.0, 10., .01)
	signal = generateSignal(t)
	fMode = recursieveGetMode(signal, t)
	ff = getEnvelopeMean(signal, t)
	subplot(211)
	plot(t, signal)
	subplot(212)
	#plot(t[70:-70], fMode)
	show()
		
	return 0

if __name__ == '__main__':
	main()

