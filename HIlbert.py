#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from numpy import *
from scipy.integrate import odeint
from matplotlib.pyplot import *

def VdP(y, t0, p):
	yp = empty_like(y)
	yp[0] = y[1]
	yp[1] = (p['r'] - y[0]**2)*y[1]- y[0]*p['w0']**2
	return yp

def hilbert(signal, D):

	weight = ones((signal.shape[0]))
	weight[:D] = linspace(0, 1, D)
	weight[-D:] = linspace(1, 0, D)
	signal = signal * weight
	return (fft.irfft((fft.rfft(signal)*1j)), weight) 	
	
def phaseVupr(signalS, signal):
	phiAll = arctan(signalS/signal)
	phiV = empty((phiAll.shape[0]))
	phiV[0] = phiAll[0]
	alpha = 0
	for i in xrange(1, phiAll.shape[0]):
		if (phiAll[i] < phiAll[i - 1]):
			alpha = alpha + pi
		phiV[i] = phiAll[i] + alpha
	return phiAll, phiV					

def main():
	pp = 50000
	x0 = [.1, -.2]
	t = arange(.1, 3000., .05)
	x = odeint(VdP, x0, t, args=({'r':1.4, 'w0':.2}, ))
	sopr = hilbert(x[pp:, 0], 6)
	gs = GridSpec(3, 2)		
	rcParams['figure.figsize'] = [12., 10.]	
	subplot(gs[0, 0])
	plot(t[pp:], x[pp:, 0], "b-")
	xlabel(r"$t$", fontsize=22)
	ylabel(r'$x_1$', fontsize=22)
	subplot(gs[1, 0])
	plot(t[pp:], sopr[0], "r-")
	subplot(gs[0, 1])
	plot(x[pp:, 0]*sopr[1], sopr[0], "g-")
	phase = phaseVupr(x[pp:, 0]*sopr[1], sopr[0])
	subplot(gs[1, 1])
	plot(t[pp:], phase[0], "r.")
	subplot(gs[2, 1])
	plot(t[pp:], phase[1])
	show()
	x = x[pp:, :]
	return 0

if __name__ == '__main__':
	main()

