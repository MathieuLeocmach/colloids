"""Phase behaviour of Colloid+Polymer mixture according to generalized free volume theory
see Gerard J Fleer and Remco Tuinier, Advances in Colloid and Interface Science 143, 1-47 (2008).
"""

import numpy as np
from scipy.optimize import fsolve, fminbound
from scipy import interpolate

def logOnePlusX(x):
    """Calculate log(1 + x), preventing loss of precision for small values of x."""
    #assert x > -1.0, "Invalid input argument (%f); must be greater than -1.0"%x
    if np.fabs(x) > 0.375:
        # x is sufficiently large that the obvious evaluation is OK
        return np.log(1.0 + x);
    # For smaller arguments we use a rational approximation
    # to the function log(1+x) to avoid the loss of precision
    # that would occur if we simply added 1 to x then took the log.
    p1 =  -0.129418923021993e+01
    p2 =   0.405303492862024e+00
    p3 =  -0.178874546012214e-01
    q1 =  -0.162752256355323e+01
    q2 =   0.747811014037616e+00
    q3 =  -0.845104217945565e-01
    t = x/(x + 2.0)
    t2 = t*t
    w = (((p3*t2 + p2)*t2 + p1)*t2 + 1.0)/(((q3*t2 + q2)*t2 + q1)*t2 + 1.0);
    return 2.0*t*w;

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
pvS_0 = lambda f: 3/(1.0/f - 1/2.853)
pv_0_1 = lambda f: (1+f)**(-2) + 8*f + 6*f**2
pv_0_2 = lambda f: -2*(1+f)**(-3) + 8 + 12*f
mu_0 =  lambda f: np.log(f) - logOnePlusX(f) + 8*f + 7*f**2 + 2*f**3
muS_0 = lambda f: 2.1306 + 3.0/((1+f)/f - 1/0.741) - 3*np.log((1+f)/f - 1/0.741)
g = lambda f,q: (beta(f,q) - (1+f)*beta1(f,q))/(1+q)**3
h = lambda f,q: beta(f,q) - f*beta1(f,q)
pv = lambda f, piv, q: pv_0(f) + piv * h(f,q)
pvS = lambda f, piv, q: pvS_0(f) + piv * h(f,q)
mu = lambda f, piv, q: mu_0(f) + piv * (1+q)**3 * g(f, q)
muS = lambda f, piv, q: muS_0(f) + piv * (1+q)**3 * g(f, q)
alpha = lambda f, q: 1/(beta(f,q) * (1+f))
qR2q = lambda qR: 0.9*qR**0.9
piv2y = lambda piv, qR: qR**3*piv
f2vf = lambda f: f/(1+f)

def mu_of_log(F, piv, q):
    f = np.exp(F)
    return F - logOnePlusX(f) + 8*f + 7*f**2 + 2*f**3  + piv * (1+q)**3 * g(f, q)
def muS_of_log (F, piv, q):
    f = np.exp(F)
    u = f*(1-1/0.741)
    return 2.1306 + 3.0*f/(1+u) + 3*F - 3*logOnePlusX(u) + piv * (1+q)**3 * g(f, q)

pv_of_log = lambda F, piv, q: pv(np.exp(F), piv, q)
pvS_of_log = lambda F, piv, q: pvS(np.exp(F), piv, q)


def critical_point(q):
    """Critical point coordinates in the (PIv, f) plane function of the effective size ratio q=delta/a"""
    fc = fsolve(lambda f: 1/f + beta3(f,q)/beta2(f,q) - pv_0_2(f)/pv_0_1(f), 0.5)[0]
    PIvc = pv_0_1(fc)/fc/beta2(fc, q)
    return (fc, PIvc)
    
def binodalGL(piv, q, guess=None):
    """return (log(f_gas), log(f_liquid)) on the binodal line at a given insersion work piv"""
    if guess is None:
        fc = critical_point(q)[0]
        fspG = fminbound(lambda f: -pv(f, piv,q), 0, fc)
        fspL = fminbound(pv, fc, 2*fc, args=(piv, q))
        guess = np.log([0.5*fspG, fspL+fspG])
    return fsolve(lambda Fs: [
        pv_of_log(Fs[0], piv, q) - pv_of_log(Fs[1], piv, q), 
        mu_of_log(Fs[0], piv, q) - mu_of_log(Fs[1], piv, q)
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
    Fc = np.log(fc)
    #start sensibly above the critical point
    startp = pivc*1.1
    fm = fminbound(mu, fc, 2*fc, args=(startp, q))
    fM = fminbound(lambda f: -pv(f, startp, q), 0, fc)
    initial_guess = np.log([0.5*fM, fm+fM])
    #construct the top of the GL binodal
    if maxpiv is None:
        maxpiv = startp*2
    topp = 1./np.linspace(1./startp, 1./maxpiv)
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
        np.column_stack((botp, np.exp(botGL[1:])))[::-1],
        np.column_stack((topp, np.exp(topGL[1:])))[1:]
        ))
    #spinodal
    spinodal = spinodalGL(q, binodal[:,0])
    #join everything
    return np.column_stack((binodal, spinodal[:,1:]))
    
def binodalFS(piv, q, guess=np.log([0.970, 1.185])):
    """return (log(f_Fluid), log(f_solid)) on the binodal line at a given insersion work piv"""
    return fsolve(lambda Fs: [
        pv_of_log(Fs[0], piv, q) - pvS_of_log(Fs[1], piv, q), 
        mu_of_log(Fs[0], piv, q) - muS_of_log(Fs[1], piv, q)
        ], guess)
    
