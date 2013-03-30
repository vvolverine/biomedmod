#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from numpy import *
from scipy.integrate import odeint
from matplotlib.pyplot import *

def VdP(y, t0, p):
	yp = empty_like(y)
	yp[0] = y[1]
	yp[1] = (p['r1'] - y[0]**2)*y[1]- y[0]*p['w01']**2 + 0.1 * (y[2]-y[0])
	yp[2] = y[3]
	yp[3] = (p['r2'] - y[2]**2)*y[3]- y[2]*p['w02']**2 + 0.1 * (y[2]-y[0])
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
	
def koeffSinhr(phase1, phase2):
	return abs(mean(exp((phase1-phase2)*1j)))			

def main():
	pp = 50000
	x0 = [.1, -.2, -.1, .2]
	t = arange(.1, 3000., .05)
	x = odeint(VdP, x0, t, args=({'r1':1.4, 'w01':.2, 'r2':1.2, 'w02':.3}, ))
	soprX = hilbert(x[pp:, 0], 6)
	soprY = hilbert(x[pp:, 2], 6)
	xwind = x[pp:, 0]*soprX[1]
	ywind = x[pp:, 2]*soprY[1]
	phaseX = phaseVupr(xwind, soprX[0])
	phaseY = phaseVupr(ywind, soprY[0])
	gs = GridSpec(4, 3)		
	rcParams['figure.figsize'] = [12., 10.]	
	subplot(gs[0, 0])
	plot(t[pp:], x[pp:, 0], "b-")
	xlabel(r"$t$", fontsize=22)
	ylabel(r'$x_1$', fontsize=22)
	subplot(gs[1, 0])
	plot(t[pp:], soprX[0], "r-")
	subplot(gs[0, 2])
	plot(xwind, soprX[0], "g-")
	subplot(gs[0, 1])
	plot(t[pp:], phaseX[0], "r.")
	subplot(gs[1, 1])
	plot(t[pp:], phaseX[1])
	# -------- Plot second signal ---------
	subplot(gs[2, 0])
	plot(t[pp:], x[pp:, 2], "b-")
	xlabel(r"$t$", fontsize=22)
	ylabel(r'$x_1$', fontsize=22)
	subplot(gs[3, 0])
	plot(t[pp:], soprY[0], "r-")
	subplot(gs[2, 2])
	plot(ywind, soprY[0], "g-")
	subplot(gs[2, 1])
	plot(t[pp:], phaseY[0], "r.")
	subplot(gs[3, 1])
	plot(t[pp:], phaseY[1])
	print koeffSinhr(phaseX[1], phaseY[1])
	show()
	x = x[pp:, :]
	return 0

if __name__ == '__main__':
	main()

