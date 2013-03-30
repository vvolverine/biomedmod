#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
from numpy import *
from matplotlib.pyplot import *
from matplotlib.gridspec import *

def generateSignal(w0, t, phi0, gamma):
	return cos(w0 * t + gamma * (t ** 2) + phi0)

def furieModule(signal, N, dt):
	furie = fft.rfft(signal)
	return (abs(furie) ** 2) / (.5 * (N))**2

def main():
	N = 10000; dt = .01; L = 1000; S = 1; komp = 200;
	t = arange(0, N*dt, dt)
	mas = empty((L/2 + 1, (N-L)/S+1))
	freq = linspace(0., 1./(2*dt), L/2+1)
	signal = generateSignal(1, t, 0, 0.05)
	gs = GridSpec(2, 1)	
	rcParams['figure.figsize'] = [10., 20.]		
	subplot(gs[0 ,0])
	plot(t, signal)
	
	weight = ones((L))
	weight[:komp] = linspace(0, 1, komp)
	weight[-komp:] = linspace(1, 0, komp)
	for i in range(0, ((N-L)/S + 1), S):
		window = signal[i*S:i*S+L] * weight
		ffx = furieModule(window, L, dt)
		mas[:, i] = ffx;
	
	subplot(gs[1, 0])
	imshow(mas[:mas.shape[0]/5, :], extent=[0, (N-L)*dt, 0, .5/dt/5], aspect=1/((.5/dt/5)/((N-L)*dt)), origin='lower')	
	savefig('soknom.png', dpi=96)		
	return 0

if __name__ == '__main__':
	main()

