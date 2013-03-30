#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from numpy import *
from scipy.integrate import odeint
from matplotlib.pyplot import *

def ML(y, t0, p):
	yp = empty_like(y)
	# as v
	yp[0] = p['I']-p['gl']*(y[0]-p['vl'])-p['gk']*y[1]*(y[0]-p['vk'])-0.5*p['gca']*(1+tanh((y[0]+.01)/.15))*(y[0]-p['vca'])
	# as w
	yp[1] = (cosh((y[0]-.1)/(2*.145))/3)*(.5*(1+tanh((y[0]-.1)/.145))-y[1]) 
	return yp
	
def ML_arg(y, t0, p):
	yp = empty_like(y)
	# as v1
	yp[0] = p['I']-p['gl']*(y[0]-p['vl'])-p['gk']*y[1]*(y[0]-p['vk'])-p['gca']*mb(y[0])*(y[0]-p['vca'])
	# as w1
	yp[1] = lambd(y[0])*y[0]*(wb(y[0])*y[0]-y[1])
	
def mb(x):
	y = .5*(1+tanh((x+.01)/.15))
	return y
	
def wb(x):
	y = .5*(1+tanh((x-.1)/.145))
	return y
	
def lambd(x):
	y = (cosh((x-.1)/(2*.145)))/3
	return y	


def main():
	pp = 700
	x0 = [-.07, 1.02]
	t = arange(.1, 3500., .05)
	x = odeint(ML_arg, x0, t, args=({'gca':1.33, 'gk':2., 'gl':0.5, 'vca':1., 'vk':-0.7, 'vl':-0.5, 'eps':0.02, 'I':.075}, ))
	rcParams['figure.figsize'] = [12., 10.]	
	subplot(221)
	plot(t[pp:], x[pp:, 0], "b")
	xlabel(r"$t$", fontsize=22)
	ylabel(r'$x_1$', fontsize=22)
	subplot(222)
	plot(t[pp:], x[pp:, 1], "r")
	xlabel(r"$t$", fontsize=22)
	ylabel(r'$x_2$', fontsize=22)
	subplot(212)
	plot(x[pp:, 0], x[pp:, 1])
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


