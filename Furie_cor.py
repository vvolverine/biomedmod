#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Без имени.py
from numpy import *
import matplotlib.pyplot as pl
from matplotlib.gridspec import *


def signalGeneration(t, ansambleSize):
	ansambleSignal = []
	for i in range(ansambleSize):
		ansambleSignal.append(sin(2 * pi * t) + random.uniform(-50, 50, t.shape))	
	return require(ansambleSignal)
	
def furieModule(signal, N, dt):
	furie = fft.rfft(signal)
	return (abs(furie) ** 2) / (.5 * (N))**2
		
def furieModuleCor(signalX, signalY, N, dt):
	furieX = fft.rfft(signalX)
	furieY = fft.rfft(signalY)
	return (furieX.conj() * furieY) / (.5 * (N))**2


def main():
	N = 10000; dt = .01; ansamble = 500;
	t = arange(0, N*dt, dt)
	x = signalGeneration(t, ansamble)
	y = signalGeneration(t, ansamble)
	fffx = []; fffy = []; fffxy = []
	freq = linspace(0., 1./(2*dt), len(fft.rfft(x[0])))
	for signal in x:
		fffx.append(furieModule(signal, 10000, dt))
	for signal in y:
		fffy.append(furieModule(signal, 10000, dt))
	for i in range(ansamble):
		fffxy.append(furieModuleCor(x[i], y[i], 10000, dt))
	fffx = require(fffx)
	fffy = require(fffy)
	fffxy = require(fffxy)
	meanFurieX = mean(fffx, 0)
	meanFurieY = mean(fffy, 0)
	meanFurieXY = mean(fffxy, 0)
	cxy = abs(meanFurieXY) / sqrt(meanFurieX * meanFurieY)
	print meanFurieX.shape, freq.shape 		
	gs = GridSpec(2, 2)
	pl.subplot(gs[0, 0])
	pl.plot(t[:10000], x[:10000][0])
	pl.subplot(gs[1, 0])
	pl.plot(freq, meanFurieX)
	pl.subplot(gs[1, 1])
	pl.plot(freq, meanFurieY)
	pl.subplot(gs[0, 1])
	pl.plot(freq, cxy)
	pl.show()
	
	return 0

if __name__ == '__main__':
	main()

