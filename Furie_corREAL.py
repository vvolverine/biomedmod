#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Без имени.py
from numpy import *
import matplotlib.pyplot as pl
from matplotlib.gridspec import *


def signalGeneration(path):
	return loadtxt(path)
	
def furieModule(signal, N, dt):
	furie = fft.rfft(signal)
	return (abs(furie) ** 2) / (.5 * (N))**2
		
def furieModuleCor(signalX, signalY, N, dt):
	furieX = fft.rfft(signalX)
	furieY = fft.rfft(signalY)
	return (furieX.conj() * furieY) / (.5 * (N))**2
	
def getFrequence(dt, N, x):
	return linspace(0., 1./(2*dt), len(fft.rfft(x[:N])))


def main():
	# declare variable
	N = 1250; dt = .004; ansamble = 84;
	t = arange(0, N*dt, dt)
	x = signalGeneration("F:\Biomedmod\Toronto_data\one_subject_network_source_1.txt")
	y = signalGeneration("F:\Biomedmod\Toronto_data\one_subject_network_source_3.txt")
	fffx = []; fffy = []; fffxy = []
	freq = linspace(0., 1./(2*dt), len(fft.rfft(x[:1250])))
	for i in range(ansamble):
		fffxy.append(furieModuleCor(x[i*1250:(i+1)*1250], y[i*1250:(i+1)*1250], 1250, dt))
		fffx.append(furieModule(x[i*1250:(i+1)*1250], 1250, dt))
		fffy.append(furieModule(y[i*1250:(i+1)*1250], 1250, dt))
	# calculate average	
	fffx = require(fffx); fffy = require(fffy); fffxy = require(fffxy)
	meanFurieX = mean(fffx, 0); meanFurieY = mean(fffy, 0); meanFurieXY = mean(fffxy, 0)
	cxy = abs(meanFurieXY) / sqrt(meanFurieX * meanFurieY)
	# show it		
	gs = GridSpec(2, 2)
	pl.subplot(gs[0, 0])
	pl.plot(t, x[:1250])
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

