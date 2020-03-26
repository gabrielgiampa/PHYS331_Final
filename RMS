#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 15:15:34 2020

@author: gabriel_giampa
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sc

#
n = 1500
dom = np.linspace(0.005,10,n)
#the domain for the graphs, in terms of R/a from the center of the galaxy
#

def integrand(x, b, mu):
    first = x**(2*b - 1)
    first = first / (1 + x)**5
    second = x**(2*b - 3)
    second = mu*second/(1+x)**3
    return first + second

vr = np.zeros(n)
err = np.zeros(n)


def rms(b, mu):
    vrms = np.zeros(n)
    for i in range(n):
        vr[i], err[i] = sc.quad(integrand, dom[i], np.inf, (b,mu))
        if(dom[i] != 0):
            vr[i] = vr[i] * (1 + dom[i])**3
            vr[i] = vr[i] / (dom[i])**(2*b - 1)
        else:
            vr[i] = 0
        vrms[i] = np.sqrt(vr[i] - 2*(b - 1)*vr[i])
    return [vrms,vr]   


"""
Use the code below this point to graph RMS velocities
Simply use plt.semilogx(dom, rms(beta, mu)) subbing in whatever values 
you want
In theory, beta cannot exceed 1 but can stretch down to - infinity
Add a label if you want the graph you made to show up in the legend

"""

plt.figure()
plt.title("$v_{RMS}$ of Stars in Hernquist Galaxy w/ Scale Length $ a$")
plt.ylabel("Relative RMS veloicty $v_{RMS} / \sqrt{GM_g /a} $")
plt.xlabel("Relative Distance from Galactic Center (R/$a$)")
plt.semilogx(dom,rms(0.0,0)[0], label = "beta = 0, mu = 0")
plt.semilogx(dom,rms(0,0.002)[0], label = "beta = 0, mu = 0.002")
plt.semilogx(dom,rms(0,0.004)[0], label = "beta = 0, mu = 0.004")
plt.semilogx(dom,rms(0,0.006)[0], label = "beta = 0, mu = 0.006")
plt.ylim(0,1.5)
plt.xlim(0.006, 10)
plt.legend()
plt.show()

"""
THIS SECTION IS NOT FUNCTIONAL YET
plt.figure()
plt.title("Dispersion of LOSV")
plt.ylabel("$\sigma_{||} / \sqrt{GM_g /a}  $")
plt.xlabel("Relative Distance from Galactic Center (R/$a$)")
plt.semilogx(dom,np.sqrt(rms(0.0,0)[1]), label = "beta = 0, mu = 0")
plt.semilogx(dom,np.sqrt(rms(0,0.002)[1]), label = "beta = 0, mu = 0.002")
plt.semilogx(dom,np.sqrt(rms(0,0.004)[1]), label = "beta = 0, mu = 0.004")
plt.ylim(0, 0.6)
plt.xlim(0.006, 10)
plt.legend()
plt.show()
"""
