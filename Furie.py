#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Без имени.py
from numpy import *
import matplotlib.pyplot as pl
from matplotlib.gridspec import *


def signalGeneration(t):	
	signal = sin(2 * pi * t) + random.uniform(-50, 50, t.shape) 
	return signal
	
def furieModule(signal, N, dt):
	furie = fft.rfft(signal)
	s = abs(furie)/(.5 * (N/dt))
	freq = arange(0., 1./(2*dt), 1./N)
	return (freq, s[:-1])
		

def main():
	N = 100000; dt = .01
	t = arange(0, N, dt)
	x = signalGeneration(t)
	fff = []
	ff = furieModule(x[:10000], 10000*dt, dt)
	for i in range(100):
		fff.append(furieModule(x[i*10000:(i+1)*10000], 10000, dt)[1])
	fff = require(fff)
	meanFurie = mean(fff, 0)	
	gs = GridSpec(3, 1)
	pl.subplot(gs[0, 0])
	pl.plot(t[:10000], x[:10000])
	pl.subplot(gs[1, 0])
	pl.plot(ff[0], meanFurie)
	pl.subplot(gs[2, 0])
	pl.plot(ff[0], ff[1])
	pl.show()
	
	return 0

if __name__ == '__main__':
	main()

