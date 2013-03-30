#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from numpy import *
from scipy.integrate import odeint
from matplotlib.pyplot import *


def VdP(y, t0, p):
    yp = empty_like(y)
    yp[0] = y[1]
    yp[1] = (p['r'] - y[0] ** 2) * y[1] - y[0] * p['w0'] ** 2
    return yp


def main():
    pp = 700
    x0 = [.1, -.2]
    t = arange(.1, 500., .05)
    x = odeint(VdP, x0, t, args=({'r': .4, 'w0': .2}, ))
    rcParams['figure.figsize'] = [12., 10.]
    subplot(221)
    plot(t[pp:], x[pp:, 0], "bo")
    xlabel(r"$t$", fontsize=22)
    ylabel(r'$x_1$', fontsize=22)
    subplot(222)
    plot(t[pp:], x[pp:, 1], "ro")
    xlabel(r"$t$", fontsize=22)
    ylabel(r'$x_2$', fontsize=22)
    subplot(212)
    plot(x[pp:, 0], x[pp:, 1])
    xlabel(r'$x_1$', fontsize=22)
    ylabel(r'$x_2$', fontsize=22)
    annotate(r'$Crivulka$', xy=(-.5, .2), xytext=(-1, .5),
             arrowprops=dict(facecolor='red', shrink=0.003),
             )
    show()
    x = x[pp:, :]
    return 0

if __name__ == '__main__':
    main()

