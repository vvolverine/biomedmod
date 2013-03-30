#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from numpy import *
from scipy.integrate import odeint
from matplotlib.pyplot import *
from matplotlib.gridspec import *

def ML(y, t0, p):
	yp = empty_like(y)
	# as v1
	yp[0] = p['I']-p['gl']*(y[0]-p['vl'])-p['gk']*y[1]*(y[0]-p['vk'])-0.5*p['gca']*(1+tanh((y[0]+.01)/.15))*(y[0]-p['vca'])+p['k']*(y[2]-y[0])
	# as w1
	yp[1] = (cosh((y[0]-.1)/(2*.145))/3)*(.5*(1+tanh((y[0]-.1)/.145))-y[1]) 
	# as v2
	yp[2] = p['I']-p['gl']*(y[2]-p['vl'])-p['gk']*y[3]*(y[2]-p['vk'])-0.5*p['gca']*(1+tanh((y[2]+.01)/.15))*(y[2]-p['vca'])+p['k']*(y[0]-y[2])
	# as w2
	yp[3] = (cosh((y[2]-.1)/(2*.145))/3)*(.5*(1+tanh((y[2]-.1)/.145))-y[3]) 
	return yp
	

def main():
	pp = 0
	x0 = [2.57, 3.02, -.67, 2.02]
	t = arange(.1, 700., .05)
	x = odeint(ML, x0, t, args=({'gca':1.33, 'gk':2., 'gl':0.5, 'vca':1., 'vk':-0.7, 'vl':-0.5, 'eps':0.02, 'I':.75, 'k':5.0}, ))
	rcParams['figure.figsize'] = [12., 10.]	
	gs = GridSpec(2, 2)
	subplot(gs[0, 0])
	plot(t[pp:], x[pp:, 0], "b")
	plot(t[pp:], x[pp:, 2], "g")
	xlabel(r"$t$", fontsize=22)
	ylabel(r'$x_1$', fontsize=22)
	subplot(gs[0, 1])
	plot(t[pp:], x[pp:, 1], "r")
	plot(t[pp:], x[pp:, 3], "orange")
	xlabel(r"$t$", fontsize=22)
	ylabel(r'$x_2$', fontsize=22)
	subplot(gs[1, :])
	plot(x[pp:, 0], x[pp:, 1], 'r')
	plot(x[pp:, 2], x[pp:, 3], 'y')
	xlabel(r'$x_1$', fontsize=22)
	ylabel(r'$x_2$', fontsize=22)
	annotate(r'$Crivulka$', xy=(.0, .5), xytext=(.2, .3),
             arrowprops=dict(facecolor='red', shrink=0.3),
             )
	show()
	x = x[pp:, :]
	return 0

if __name__ == '__main__':
	main()


