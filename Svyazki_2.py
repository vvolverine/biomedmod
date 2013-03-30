#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from numpy import *
from scipy.integrate import odeint
from matplotlib.pyplot import *

def ML(y, t0, p):
	# verify and set som logical parameter
	# for k1
	eps = p['ag01'] * 1e-5
	
	ag1 = ag(p['ag01'], p['lg'], y[0])
	ag2 = ag(p['ag02'], p['lg'], y[2])
	
	lg1 = lgn(p['d1'], p['ro'], ag1)
	lg2 = lgn(p['d2'], p['ro'], ag2)
	
	rk1 = frk1(p['ro'], ag1)
	rk2 = frk2(p['ro'], ag2, p['a1'])
	
	er1 = ern(p['nu'], p['lg'], p['d1'], ag1)
	er2 = ern(p['nu'], p['lg'], p['d2'], ag2)
	
	if ((ag1 <= eps) or (ag2 <= eps)):
		ugp = 0
	else:
		ugp = ((lg1 + lg2)**-1) * (p['ps'] - (rk1 + rk2) * y[4] * abs(y[4]) - (er1 + er2) * y[4])
	
	if (ag1 >= eps):
		k1 = k_n1(y[0], p['k1'], p['etta1'])
		r1 = r_n(p['dz1'], p['m1'], p['k1'])
	else:
		k1 = k_n2(y[0], p['k1'], p['etta1'], p['h1'], p['ag01'], p['lg'], p['ettah1'])
		r1 = r_n(p['dzh1'], p['m1'], p['k1'])	
			 		
	if (ag2 >= eps):
		k2 = k_n1(y[2], p['k2'], p['etta2'])
		r2 = r_n(p['dz2'], p['m2'], p['k2'])
	else:		 
		k2 = k_n2(y[2], p['k2'], p['etta2'], p['h2'], p['ag02'], p['lg'], p['ettah2']) 
		r2 = r_n(p['dzh2'], p['m2'], p['k2'])
		
	# Define force f
	if (ag1 > eps) and (ag2 > eps):
		pm1 = Pm1(p['ps'], y[4], ag1, p['nu'], p['lg'], p['d1'], p['ro'], ugp)
		f1 = pm1 * p['d1'] * p['lg']
		f2 = Pm2(pm1, p['nu'], p['lg'], p['d1'], p['d2'], ag1, ag2, p['ro'], y[4], ugp) * p['d2'] * p['lg']
	elif ((ag1 < eps)):
		f1 = p['ps'] * p['d1'] * p['lg']
		f2 = 0
	elif (ag1 > eps) and (ag2 < eps):
		f1 = p['ps'] * p['d1'] * p['lg']
		f2 = p['ps'] * p['d2'] * p['lg']	
						
	yp = empty_like(y)
	# 0 as x1, 1 as y1, 2 as x2, 3 as y2, 4 as ug
	# as x1
	yp[0] = y[1]
	# as y1
	yp[1] = (f1 - p['kc']*(y[0] - y[2]) - k1 - r1*y[1])/p['m1']
	# as x2
	yp[2] = y[3]
	# as y2
	yp[3] = (f2 - p['kc']*(y[2] - y[0]) - k2 - r2*y[3])/p['m2']
	# as ug	
	yp[4] = ugp
	return yp

def Pm1(ps, ug, ag1, nu, lg, d1, ro, ugp):
	er1 = ern(nu, lg, d1, ag1)
	lg1 = lgn(d1, ro, ag1)
	return ps - 0.5 * (1.37 *ro* ((ug/ag1)**2) - (er1 * ug + lg1 * ugp)) 
	
def Pm2(pm1, nu, lg, d1, d2, ag1, ag2, ro, ug, ugp):
	er1 = ern(nu, lg, d1, ag1)
	er2 = ern(nu, lg, d2, ag2)
	lg1 = lgn(d1, ro, ag1)
	lg2 = lgn(d2, ro, ag2)
	return pm1 - 0.5 * (ug * (er1 + er2) + (lg1 + lg2)*ugp + ro * (ug ** 2) * (ag2 ** -2 - ag1 ** -2))

def ag(ago, lg, x):
	return ago + 2 * lg * x

def k_n1(x, k, etta):
	return k * x * (1 + etta * (x ** 2)) 
	
def k_n2(xn, k, ettan, hn, agon, lg, ettahn):
	return k*xn*(1 + ettan*(xn**2)) + hn*(xn + agon/(2*lg))*(1 + ettahn * ((xn + agon/(2*lg)) ** 2))
	
def r_n(dz, m, k):
	return dz * sqrt(m / k)	
	
def frk1(ro, ag1):
	return (.19 * ro) / (ag1 ** 2)
	
def frk2(ro, ag2, a1):
	return (ro * (.5 - ag2/a1)) / (ag2 ** 2)	
	
def lgn(dn, ro, agn):
	return (dn * ro) / (agn)
	
def ern(nu, lg, dn, agn):
	return (12 * nu * (lg ** 2) * dn) / (agn ** 3) 
	
def main():
	pp = 8000
	x0 = [.01, -0.01, .01, 1.4, 0.05]
	t = arange(0, 0.15, .000005)
	dc = {'m1':0.125, 'm2':0.025, 'lg':1.4, 'd1':0.25, 'd2':0.05, 'ag01':0.056, 
		  'ag02':0.056, 'a1':0.168, 'kc':2.5e4, 'k1':8e4, 'k2':8e3, 'etta1':1e2, 'etta2':1e2, 
		  'ettah1':5e2, 'ettah2':5e2, 'dz1':0.3e4, 'dz2':1.e4, 'dzh1':0.1e5, 'dzh2':0.3e5, 
		  'ro':1.14e-3, 'nu':18.2466e-5, 'ps':9.5e3}
	dc['h1'] = 3 * dc['k1']
	dc['h2'] = 3 * dc['k2'] 
	x = odeint(ML, x0, t, args=(dc, ))
	rcParams['figure.figsize'] = [12., 10.]	
	subplot(221)
	plot(t[pp:], x[pp:, 0], "b")
	plot(t[pp:], x[pp:, 2], "g")
	xlabel(r"$t$", fontsize=22)
	ylabel(r'$x_1$', fontsize=22)
	subplot(222)
	plot(t[pp:], x[pp:, 1], "r")
	plot(t[pp:], x[pp:, 3], "m")
	xlabel(r"$t$", fontsize=22)
	ylabel(r'$x_2$', fontsize=22)
	subplot(212)
	plot(x[pp:, 0], x[pp:, 2])
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


