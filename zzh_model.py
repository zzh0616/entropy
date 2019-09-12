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


def gpob(x,x0,sigma):
    pos=numpy.log(numpy.power(2*3.1416,-0.5)/sigma)-(x-x0)*(x-x0)/2/sigma/sigma
    return pos

def calc_cri_den(z,H0=2.3e-18,omega_m=0.27):
    G=6.673*1e-8 #cm^3 g^-1 s^2
    E=math.sqrt(omega_m*(1+z)*(1+z)*(1+z)+1-omega_m)
    H=H0*E
    den=3*H*H/8/3.1416/G # g/cm^3
    msun=2e33 #g
    kpc=3.086e21 #cm
    return den*kpc*kpc*kpc/msun #msun/kpc^3

def nfw(r,rs,rho):
    return rho/(r/rs)/(1+r/rs)/(1+r/rs)
def ave_den(r,rs,rho):
    return 4*3.1416*rho*rs*rs*rs*(math.log((rs+r)/rs)-r/(r+rs))/(4/3*3.1416*r*r*r)
#rs in kpc ; rho in msun/kpc^3 N1 in unit kev*kpc/Msun
def main(R,R500_C,R200_C,T0,N1,N2,RS,A0,GAMMA0,K0,A1,GAMMA1,K1,N3,Z,RHO,A,BETA,NT,NM,X):
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
    K_MOD=A0*math.pow(R,GAMMA0)+K0
    K_OBS=A1*math.pow(R,GAMMA1)+K1
    M200 = abs( 4 * 3.1416*RHO*RS*RS*RS*(math.log((RS+R200)/RS)-R200/(RS+R200)))
    ETA=1-A*math.pow(1+Z,BETA)*math.pow(R/R500,NT)*\
            math.pow(M200/3/math.pow(10,14),NM)
    T=(T0+(1-N3*X)*(N1*RS*RHO*RS*RS*math.log((RS+R)/RS)/R))/\
            (1/ETA-N3-abs(N2)*(K_OBS-K_MOD)/K_OBS)
    return T
def dm(R,R500_C,R200_C,T0,N1,N2,RS,A0,GAMMA0,K0,A1,GAMMA1,K1,N3,Z,RHO,A,BETA,NT,NM,X,EFF):
    msun=2e33 #g
    kpc=3.086e21 #cm
    cri_den=calc_cri_den(Z)*kpc*kpc*kpc/msun #msun/kpc^3
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
    K_MOD=A0*math.pow(R,GAMMA0)+K0
    K_OBS=A1*math.pow(R,GAMMA1)+K1
    M200 = 4 * 3.1416*RHO*RS*RS*RS*(math.log((RS+R200)/RS)-R200/(RS+R200))
    ETA=1-A*math.pow(1+Z,BETA)*math.pow(R/R500,NT)*\
            math.pow(M200/3/math.pow(10,14),NM)
    T=(T0+(1+5*EFF-N3*X)*N1*RS*RHO*RS*RS*math.log((RS+R)/RS)/R)/\
            (1/ETA-N3-abs(N2)*(K_OBS-K_MOD)/K_OBS)
    return T

def hydro(R,R500_C,R200_C,T0,N1,N2,RS,A0,GAMMA0,K0,A1,GAMMA1,K1,N3,Z,RHO,A,BETA,NT,NM,X):
    msun=2e33 #g
    K1=abs(K1)
    kpc=3.086e21 #cm
    cri_den=calc_cri_den(Z)*kpc*kpc*kpc/msun #msun/kpc^3
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
    K_MOD=A0*math.pow(R,GAMMA0)+K0
    K_OBS=A1*math.pow(R,GAMMA1)+K1
    M200=0
    NT=abs(NT)
    ETA=1-A*math.pow(1+Z,BETA)*math.pow(R/R500,NT)*\
            math.pow(M200/3/math.pow(10,14),NM)
    T=(T0+(1-N3*X)*N1*RS*RHO*RS*RS*math.log((RS+R)/RS)/R)/\
            (1/ETA-N3-abs(N2)*(K_OBS-K_MOD)/K_OBS)
    return T
def beta_model(r,N0,BETA0,RC0,N1,BETA1,RC1,BKG):
    return abs(N0)*math.pow(1+r*r/RC0/RC0,-1.5*abs(BETA0))+abs(N1)*math.pow(1+r*r/RC1/RC1,-1.5*abs(BETA1))+abs(BKG)


def extra_energy(ROUT,N0,BETA0,RC0,N1,BETA1,RC1,BKG,R500_C,R200_C,T0,C1,C2,RS,A0,GAMMA0,K0,A1,GAMMA1,K1,C3,Z,RHO,A,BETA,NT,NM,X):
    energy=0.0
    tot_num=0
    grav_energy=0
    res_energy=0
    kpc=3.086e21 #cm
    M200 = 4 * 3.1416*RHO*RS*RS*RS*(math.log((RS+R200_C)/RS)-R200_C/(RS+R200_C))
    for i in range(int(ROUT)):
        if i == 0:
            continue
        den=beta_model(i-0.5,N0,BETA0,RC0,N1,BETA1,RC1,BKG) #cm^-3
        den=den*2.2
        vol=4/3*3.1416*(3*i*i-3*i+1) #kpc^3
        num=den*vol*kpc*kpc*kpc
        tot_num=tot_num+num
        T=main(i-0.5,R500_C,R200_C,T0,C1,C2,RS,A0,GAMMA0,K0,A1,GAMMA1,K1,C3,Z,RHO,A,BETA,NT,NM,X)
        K_MOD=A0*math.pow(i-0.5,GAMMA0)+K0
        K_OBS=A1*math.pow(i-0.5,GAMMA1)+K1
#        print('K_MOD=%.4e'%(K_MOD))
#        print('K_OBS=%.4e'%(K_OBS))
        energy=energy+1.5*C2*T*(K_OBS-K_MOD)/K_OBS*num
        grav_energy=C1*1.5*RHO*RS*RS*RS*math.log((RS+i-0.5)/(RS))/(i-0.5)*num+grav_energy
        R=i-0.5
        ETA=1-A*math.pow(1+Z,BETA)*math.pow(R/R500_C,NT)*\
                math.pow(M200/3/math.pow(10,14),NM)
        res_energy=1.5*T*num/ETA+res_energy
#        print(energy)
    return [energy,energy/tot_num,grav_energy,res_energy]

