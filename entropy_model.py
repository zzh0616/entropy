#!/usr/bin/env python3
from scipy.optimize import minimize
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

import sys
import math
import numpy
import scipy
import matplotlib
import re
import types
import zzh_model 

def gpob(x,x0,sigma):
    pos=math.log(math.pow(2*3.1416,-0.5)/sigma)-(x-x0)*(x-x0)/2/sigma/sigma
    return pos

def entropy(r,r500,a0,gamma0,k0,alpha):
    k=a0*math.pow(r,gamma0)*math.exp(alpha*(r/r500)*(r/r500))+k0
    return k

def main(R,R500_C,R200_C,T0,N1,N2,RS,A0,GAMMA0,K0,K_OBS,N3,Z,RHO,A,BETA,NT,NM,X,alpha):
    def nfw(r,rs=RS,rho=RHO):
        return rho/(r/rs)/(1+r/rs)/(1+r/rs)
    def ave_den(r,rs=RS,rho=RHO):
        return 4*3.1416*rho*rs*rs*rs*(math.log((rs+r)/rs)-r/(r+rs))/(4/3*3.1416*r*r*r)
    f500 = lambda x: ave_den(x)-500*cri_den
    f200 = lambda x: ave_den(x)-200*cri_den
    
    R200=R200_C
    R500=R500_C
    
#    R200=fsolve(f200,R200_C)
#    R500=fsolve(f500,R500_C)
#    print('R200=%.4e'%R200)
#    print('R500=%.4e'%R500)
    K1=abs(K1)
    NT=abs(NT)
#    K_MOD=entropy(R,R500,A0,GAMMA0,K0,alpha)
    K_MOD=A0*math.pow(R,GAMMA0)+K0
#    K_OBS=A1*math.pow(R,GAMMA1)+K1
    M200 = abs( 4 * 3.1416*RHO*RS*RS*RS*(math.log((RS+R200)/RS)-R200/(RS+R200)))
    ETA=1-A*math.pow(1+Z,BETA)*math.pow(R/R500,NT)*\
            math.pow(M200/3/math.pow(10,14),NM)
    T=(T0+(1-N3*X)*(N1*RS*RHO*RS*RS*math.log((RS+R)/RS)/R))/\
            (1/ETA-N3-abs(N2)*(K_OBS-K_MOD)/K_OBS)
    return T
def beta_model(r,N0,BETA0,RC0,N1,BETA1,RC1,BKG):
    return abs(N0)*math.pow(1+r*r/RC0/RC0,-1.5*abs(BETA0))+abs(N1)*math.pow(1+r*r/RC1/RC1,-1.5*abs(BETA1))+abs(BKG)

def dydr(y,r0,p):
    G=6.67e-11 #m^3 kg^-1 s^-2
    pi=3.1415926
    mp=1.67e-27  #kg
    kev=1.6e-16 #J
    kpc=3.086e19 #m
    Msun=2e30 #kg
    mu=0.61
    rs=p[2]
    rho=p[7]
    c2=p[1]
    c3=p[6]
    c4=1-c3-c2
    a0=p[3]
    gamma0=p[4]
    k0=p[5]
    c1=p[0]
    k=a0*numpy.power(r0,gamma0)+k0
    dkdr=a0*gamma0*numpy.power(r0,gamma0-1)
    z1=-4*pi*rho*rs*rs*rs*(numpy.log((r0+rs)/rs)-r0/(rs+r0))*G*mu*mp/r0/r0\
            *Msun/kpc/kev   #unit in kev/kpc same as z2
    z2=(1-c3)*c1*rho*rs*rs*rs*(-numpy.log((rs+r0)/rs)/r0/r0+1/(r0*(rs+r0)))
    ne=y[1]
    T=y[0]
    dTdr=z1-(z2-c4*z1-c2*dkdr*numpy.power(ne,2/3))/(2/3*k*c2/T*numpy.power(\
            ne,2/3)-c4)
    dnedr=(z2-c4*z1-c2*dkdr*numpy.power(ne,2/3))/(2*k*c2/3/numpy.power(\
            ne,1/3)-c4*T/ne)
    return [dTdr,dnedr]

