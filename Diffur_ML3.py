#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from numpy import *
from scipy.integrate import odeint
from matplotlib.pyplot import *
from matplotlib.gridspec import *

def ML(y, t0, p):
	yp = empty_like(y)
	for i in range(0, 12, 2):
		yp[i] = pr1(y[i], y[i+1], p)+p['k']*y[i-2]
		yp[i+1] = pr2(y[i], y[i+1], p)  
	return yp
	
def pr1(a, b, p):
	# a as y[0] -- v
	# b as y[1] -- w
	return p['I']-p['gl']*(a-p['vl'])-p['gk']*b*(a-p['vk'])-0.5*p['gca']*(1+tanh((a+.01)/.15))*(a-p['vca'])	
	
def pr2(a, b, p):
	# a as y[0] -- v
	# b as y[1] -- w
	return (cosh((a-.1)/(2*.145))/3)*(.5*(1+tanh((a-.1)/.145))-b) 	
	

def main():
	pp = 0
	x0 = [1.6, -4.02, 1.4, 2.02, 1.2, 4.02, 1., -2.02, .8, 4.02, .6, 2.02]
	#~ x0 = [2.57, 3.02, -.67, 2.02]
	t = arange(.1, 70., .05)
	x = odeint(ML, x0, t, args=({'gca':1.33, 'gk':2., 'gl':0.5, 'vca':1., 'vk':-0.7, 'vl':-0.5, 'eps':0.02, 'I':.108, 'k':.01}, ))
	rcParams['figure.figsize'] = [12., 10.]	
	gs = GridSpec(2, 2)
	ax = subplot(gs[0, 0])
	for i in range(0, 12, 2):
		plot(t[pp:], x[pp:, i], "b")
		yticks((-2, -1, 0, 1, 2))
		ax.set_yticklabels(('a', 'd', 'd', 'f', 'f'))
	xlabel(r"$t$", fontsize=22)
	ylabel(r'$x_1$', fontsize=22)
	subplot(gs[0, 1])
	for i in range(1, 13, 2):
		plot(t[pp:], x[pp:, i], "r")
		yticks((-5, -2,  0, 2, 5))
	xlabel(r"$t$", fontsize=22)
	ylabel(r'$x_2$', fontsize=22)
	subplot(gs[1, :])
	for i in range(0, 12, 2):
		plot(x[pp:, i], x[pp:, i+1], color=(i*.1, .5, .1))
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


