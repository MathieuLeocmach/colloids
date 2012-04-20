"""Phase behaviour of Colloid+Polymer mixture according to generalized free volume theory
see Gerard J Fleer and Remco Tuinier, Advances in Colloid and Interface Science 143, 1-47 (2008).
"""

import numpy as np
from scipy.optimize import fsolve, fminbound
from scipy import interpolate

A = lambda q: (1+q)**3 -1
B = lambda q: 3*q**2*(q+1.5)
C = lambda q: q**3
Q = lambda f, q: A(q)*f + B(q)*f**2 + C(q)*f**3
Q1 = lambda f, q: A(q) + 2*B(q)*f + 3*C(q)*f**2
Q2 = lambda f, q: 2*B(q) + 6*C(q)*f
Q3 = lambda f, q: 6*C(q)
beta = lambda f,q: np.exp(-Q(f,q))
beta1 = lambda f,q: -beta(f,q)*Q1(f,q)
beta2 = lambda f,q: -beta(f,q)*Q2(f,q) - beta1(f,q)*Q1(f,q)
beta3 = lambda f,q: -beta(f,q)*Q3(f,q) - 2*beta1(f,q)*Q2(f,q)- beta2(f,q)*Q1(f,q)
pv_0 = lambda f: f/(1+f) + 4*f**2 + 2* f**3
pv_0_1 = lambda f: (1+f)**(-2) + 8*f + 6*f**2
pv_0_2 = lambda f: -2*(1+f)**(-3) + 8 + 12*f
mu_0 =  lambda f: np.log(f/(1+f)) + 8*f + 7*f**2 + 2*f**3
g = lambda f,q: (beta(f,q) - (1+f)*beta1(f,q))/(1+q)**3
h = lambda f,q: beta(f,q) - f*beta1(f,q)
pv = lambda f, piv, q: pv_0(f) + piv * h(f,q)
mu = lambda f, piv, q: mu_0(f) + piv * (1+q)**3 * g(f, q)
alpha = lambda f, q: 1/(beta(f,q) * (1+f))
qR2q = lambda qR: 0.9*qR**0.9
piv2y = lambda piv, qR: qR**3*piv
f2vf = lambda f: f/(1+f)


def critical_point(q):
    """Critical point coordinates in the (PIv, f) plane function of the effective size ratio q=delta/a"""
    fc = fsolve(lambda f: 1/f + beta3(f,q)/beta2(f,q) - pv_0_2(f)/pv_0_1(f), 0.5)[0]
    PIvc = pv_0_1(fc)/fc/beta2(fc, q)
    return (fc, PIvc)
    
def binodalGL(piv, q, guess=None):
    """return (f_gas, f_liquid) on the binodal line at a given insersion work piv"""
    if guess is None:
        fc = critical_point(q)[0]
        fm = fminbound(mu, fc, 2*fc, args=(piv, q))
        fM = fminbound(lambda f: -pv(f, piv,q), 0, fc)
        guess = [0.5*fM, fm+fM]
    return fsolve(lambda fs: [
        pv(fs[0], piv, q) - pv(fs[1], piv, q), 
        mu(fs[0], piv, q) - mu(fs[1], piv, q)
        ], guess)

def spinodalGL(q, bins=None):
    fc, pivc = critical_point(q)
    if bins is None:
        bins = 1/np.linspace(1/pivc, 1/(8*pivc))
    return np.column_stack((bins, np.vstack([(
        fminbound(lambda f: -pv(f, piv, 0.06), 0, fc),
        fminbound(lambda f: pv(f, piv, 0.06), fc, 12)
        ) for piv in bins])))
    
        
def all_GL(q, maxpiv=None):
    fc, pivc = critical_point(q)
    #start sensibly above the critical point
    startp = pivc*1.1
    fm = fminbound(mu, fc, 2*fc, args=(startp, q))
    fM = fminbound(lambda f: -pv(f, startp, q), 0, fc)
    initial_guess = [0.5*fM, fm+fM]
    #construct the top of the GL binodal
    if maxpiv is None:
        maxpiv = startp*2
    topp = np.linspace(startp, maxpiv)
    topGL = [initial_guess]
    for piv in topp:
        topGL.append(binodalGL(piv, q, topGL[-1]))
    #construct the GL binodal between the starting piv and the critical point
    botp = np.linspace(startp, pivc)[:-1]
    botGL = [initial_guess]
    for piv in botp:
        botGL.append(binodalGL(piv, q, botGL[-1]))
    #join the two results
    binodal = np.vstack((
        [[pivc, fc, fc]],
        np.column_stack((botp, botGL[1:]))[::-1],
        np.column_stack((topp, topGL[1:]))[1:]
        ))
    #spinodal
    spL = [fminbound(mu, G, L, args=(piv, q)) for piv, G, L in binodal]
    spG = [fminbound(lambda f: -mu(f, piv, 0.06), G, L) for piv, G, L in binodal]
    #join everything
    return np.column_stack((binodal, spG, spL))
    
    
    
