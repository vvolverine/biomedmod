#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from numpy import *
from scipy.integrate import odeint
from matplotlib.pyplot import *

def ML(y, t0, p):
	yp = empty_like(y)
	# as x0
	yp[0] = y[1]
	# as x1
	yp[1] = ((2*p['tau']*p['Ps']*y[1])/abs((y[0]+p['x0']+p['tau']*y[1]))-p['r']*y[1]-p['k']*y[0])/p['m']
	return yp
	
def main():
	pp = 700
	x0 = [.0000002, -0.01]
	t = arange(0, 1.5, .000005)
	dc = {'m':4.76, 'd':.3e-2, 'lg':1.4e-2,
		  'Ps':9.745e2, 'x0':0.5e-3, 'r':2e3, 'k':1e6}
	dc['tau'] = (pi/12) * sqrt(dc['m']/dc['k'])
	print dc['tau']  
	x = odeint(ML, x0, t, args=(dc, ))
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


