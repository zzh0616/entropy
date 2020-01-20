#!/usr/bin/env python3
## $1 temperature data; radius in kpc
## $2 temperature parameter file
## $3 sbp data; 'radius err sbp err' radius in kpc
## $4 cooling function file , radius in kpc, radius should match density
## $5 cm_per_pixel for this source !! no need anymore
## $6 cooling function file for chandra fit
## $7 outfile name for fitted surface brightness !not used
## $8 outfile name for fitted entropy ! not used
####
#tempory var define area
####
from scipy.optimize import minimize
from scipy.optimize import fsolve
from scipy.integrate import odeint
from scipy.integrate import ode
from scipy.integrate import quad
from zzh_model import gpob
from numpy.random import random
import zzh_model
import deproject_model
import modnfw_readarray
import matplotlib.pyplot as plt
import sys
import math
import numpy
import numpy as np
import scipy
import matplotlib
import re
import types
import pymc
import json
import cProfile
import pstats
#global T0_0,T0_MIN,T0_MAX,N1_0,N1_MIN,N1_MAX,N2_0,N2_ERR
#global RS_0,RS_ERR,A0_0,A0_ERR,GAMMA0_0,GAMMA0_ERR,K0_0,K0_MIN,K0_MAX
#global A1_0,A1_ERR,GAMMA1_0,GAMMA1_ERR,K1_0,K1_ERR
#global N3_0,N3_MIN,N3_MAX,RHO_0,RHO_ERR,Z_0,Z_MIN,Z_MAX
#global A_0,A_ERR,BETA_0,BETA_MIN,BETA_MAX,NT_0,NT_ERR
#global NM_0,NM_MIN,NM_MAX,X_0,X_MIN,X_MAX
#global r_array,re_array,t_array,te_array
numpy.set_printoptions(linewidth=100000)
pi=3.1415926
cm_per_pixel=np.float(sys.argv[5])
name=sys.argv[1][0:-4]
print(name)
#print(cm_per_pixel)
G=6.67e-11 #m^3 kg^-1 s^-2
mp=1.67e-27  #kg
kev=1.6e-16 #J
kpc=3.086e19 #m
Msun=2e30 #kg
mu=0.61
script_dir=''
tmp=sys.argv[0].split('/')
for i in range(len(tmp)-1):
    script_dir=script_dir+tmp[i]
    script_dir=script_dir+'/'
ta1=numpy.load(script_dir+'lrs_ori.npy')
ta2=numpy.load(script_dir+'lrs_dvr.npy')
ta3=numpy.load(script_dir+'hrs_ori.npy')
ta4=numpy.load(script_dir+'hrs_dvr.npy')
t_total=[ta1,ta2,ta3,ta4]
##p=[rho,rs,delta,delta2]
def mod_nfw(r,p):
    rho=p[0]
    rs=p[1]
    delta=p[2]
    delta2=p[3]
    den=rho/(np.power(r/rs,delta)*(numpy.power(1+r/rs,delta2-delta)))*4*pi*r*r #Msun/kpc^3
    return den

def Modefied_Mnfw(r,p):
    M=modnfw_readarray.calc('n',r,p,t_total)
    return M

def mod_nfw_divbyr(r,p):
    a=mod_nfw(r,p)/r
    return a
def nfw(r,p):
    a=mod_nfw(r,[p[0],p[1],1.0,3.0])
    return a
def nfw_divbyr(r,p):
    a=nfw(r,p)/r
    return a

def calc_den(r,p_den,p_mnfw):
    rs=p_mnfw[1]
    tau=p_den[0]
    fg=p_den[1]
    s=p_den[2]
    rtmp=s*np.power(r/s,tau)
    tmp_const=Msun/kpc/kpc/kpc/1000000/(0.61*mp)
    den=tau*fg*np.power(r/s,3*tau-3)*mod_nfw(rtmp,p_mnfw)/4/pi/rtmp/rtmp*tmp_const #cm^-3
    return den

### length of x(radius) and k_fit(entropy) must be the same
def calc_T(x,den_fit,m_factor,p):
    r0=x
    G=6.67e-11 #m^3 kg^-1 s^-2
    mp=1.67e-27  #kg
    kev=1.6e-16 #J
    kpc=3.086e19 #m
    Msun=2e30 #kg
    mu=0.61
    rs=p[2]
    rho=p[7]
    c2=p[1]
    c3=p[6]
    a0=p[3]
    gamma0=p[4]
    k0=p[5]
    c1=p[0]
    k=a0*numpy.power(r0,gamma0)+k0
    dkdr=a0*gamma0*numpy.power(r0,gamma0-1)
    r200=R200_0
    r500=R500_0
    nth_a=p[8]
    nth_b=p[9]
    nth_gamma=p[10]
    delta=p[11]
    delta2=p[12]
    c4=p[13]
    p_MMnfw=[rho,rs,delta,delta2]
    eta=nth_a*(1+numpy.exp(-numpy.power(r0/R200M/nth_b,nth_gamma)))
    Tx=(T0_0+(1-c4*c3)*c1*rho*rs*rs*rs*numpy.log((rs+x)/rs)/x*m_factor-k*c2*np.power(den_fit,2/3))/(1/eta-c3-c2)
    return Tx

def entropy_model(T_array,ne_array):
    return T_array*np.power(ne_array,-2/3)
def calc_mass(r_array,t_array,ne_array,p_eta):
    mass_array=[]
    a=p_eta[0]
    b=p_eta[1]
    c=p_eta[2]
    tmp_const=1/G/mu/mp/Msun*kpc*kev
    for i in range(len(r_array)):
        r0=r_array[i]
        eta=0.0+1.0*a*(1+math.exp(-math.pow(r0/R200M/b,c)))
        detadr=a*math.exp(-math.pow(r0/R200M/b,c))*\
                c*-math.pow(r0/R200M/b,c-1)/R200M/b
        detadr=1.0*detadr
#        eta=1
#        detadr=0
        if i==0:
            dTdr=(t_array[i+1]-t_array[i])/(r_array[i+1]-r_array[i])
            dnedr=(ne_array[i+1]-ne_array[i])/(r_array[i+1]-r_array[i])
            z1=(dnedr*t_array[i]/ne_array[i]+dTdr-t_array[i]/eta*detadr)/eta #kev/kpc
            mass=-z1*r0*r0*tmp_const #Msun
        elif i==len(r_array)-1:
            dTdr=(t_array[i]-t_array[i-1])/(r_array[i]-r_array[i-1])
            dnedr=(ne_array[i]-ne_array[i-1])/(r_array[i]-r_array[i-1])
            z1=(dnedr*t_array[i]/ne_array[i]+dTdr-t_array[i]/eta*detadr)/eta #kev/kpc
            mass=-z1*r0*r0*tmp_const #Msun
        else:
            dTdr_r=(t_array[i+1]-t_array[i])/(r_array[i+1]-r_array[i])
            dnedr_r=(ne_array[i+1]-ne_array[i])/(r_array[i+1]-r_array[i])
            dTdr_l=(t_array[i]-t_array[i-1])/(r_array[i]-r_array[i-1])
            dnedr_l=(ne_array[i]-ne_array[i-1])/(r_array[i]-r_array[i-1])
            dTdr=(dTdr_l+dTdr_r)/2
            dnedr=(dnedr_l+dnedr_r)/2
            z1=(dnedr*t_array[i]/ne_array[i]+dTdr-t_array[i]/eta*detadr)/eta #kev/kpc
            mass=-z1*r0*r0*tmp_const #Msun
        if mass<=1e7:
            mass=1e7
        mass_array.append(mass)
    return mass_array
def clumping_model(x,n1,n2,n3,n4,n5):
    tmp=numpy.power(1+x,n1)*numpy.exp(x*n2)+n3*numpy.exp(-(x-n4)*(x-n4)/n5)
    for i in range(len(tmp)):
        if tmp[i] < 1 :
            tmp[i]=1
    return tmp
r_array=[]
re_array=[]
t_array=[]
te_array=[]
rsbp_array=[]
rsbpe_array=[]
sbp_array=[]
sbpe_array=[]
rcsbp_array=[]
rcsbpe_array=[]
csbp_array=[]
csbpe_array=[]
flag_tproj_array=[]
f_sbp_array=[]
for i in open(sys.argv[1]):
    r,rer,t,te,flag_tproj=i.split()
    r=float(r)
    rer=float(rer)
    t=float(t)
    te=float(te)
    r_array.append(r)
    re_array.append(rer)
    t_array.append(t)
    te_array.append(te)
    flag_tproj_array.append(flag_tproj)


#r_array=numpy.array(r_array)
#re_array=numpy.array(re_array)
#t_array=numpy.array(t_array)
#te_array=numpy.array(te_array)

for i in open(sys.argv[3]):
    r,rer,sbp,sbpe,f_sbp=i.split()
    r=float(r)
    rer=float(rer)
    sbp=float(sbp)
    sbpe=float(sbpe)
    if f_sbp=='c':
        rcsbp_array.append(r)
        rcsbpe_array.append(rer)
        csbp_array.append(sbp)
        csbpe_array.append(sbpe)
    else:
        rsbp_array.append(r)
        rsbpe_array.append(rer)
        sbp_array.append(sbp)
        sbpe_array.append(sbpe)
    f_sbp_array.append(f_sbp)
#rsbp_array=numpy.array(rsbp_array)
#rsbpe_array=numpy.array(rsbpe_array)
#sbp_array=numpy.array(sbp_array)
#sbpe_array=numpy.array(sbpe_array)
FLAG_SBPC=0
FLAG_FIT=0
FLAG_ITN=0
FLAG_CLUMP=0
C4_0=0.5
C4_ERR=2.0
C4_MIN=0.0
C4_MAX=1
for i in open(sys.argv[2],'r'):
    if re.match(r'^T0\s',i):
        T0_0=float(i.split()[1])
        T0_MIN=float(i.split()[2])
        T0_MAX=float(i.split()[3])
    elif re.match(r'^N1\s',i):
        N1_0=float(i.split()[1])
        N1_MIN=float(i.split()[2])
        N1_MAX=float(i.split()[3])
    elif re.match(r'^N2\s',i):
        N2_0=float(i.split()[1])
        N2_ERR=float(i.split()[2])
    elif re.match(r'^rs\s',i):
        RS_0=float(i.split()[1])
        RS_ERR=float(i.split()[2])
    elif re.match(r'^a0\s',i):
        A0_0=float(i.split()[1])
        A0_ERR=float(i.split()[2])
    elif re.match(r'^gamma0\s',i):
        GAMMA0_0=float(i.split()[1])
        GAMMA0_ERR=float(i.split()[2])
    elif re.match(r'k0\s',i):
        if len(i.split())==4:
            K0_0=float(i.split()[1])
            K0_MIN=float(i.split()[2])
            K0_MAX=float(i.split()[3])
            K0_ERR=K0_MAX-K0_MIN
        elif len(i.split())==5:
            K0_0=float(i.split()[1])
            K0_ERR=float(i.split()[2])
            K0_MIN=float(i.split()[3])
            K0_MAX=float(i.split()[4])
    elif re.match(r'^kmod_a\s',i):
        KMOD_A0=float(i.split()[1])
        KMOD_AERR=float(i.split()[2])
    elif re.match(r'^kmod_b\s',i):
        KMOD_B0=float(i.split()[1])
        KMOD_BERR=float(i.split()[2])
    elif re.match(r'^kmod_c\s',i):
        KMOD_C0=float(i.split()[1])
        KMOD_CERR=float(i.split()[2])
    elif re.match(r'^kmod_k0\s',i):
        KMOD_K0=float(i.split()[1])
        KMOD_K0ERR=float(i.split()[2])
    elif re.match(r'^N3\s',i):
        if len(i.split())==4:
            N3_0=float(i.split()[1])
            N3_C=1.3
            N3_MIN=float(i.split()[2])
            N3_MAX=float(i.split()[3])
            N3_ERR=1.5
        elif len(i.split())==5:
            N3_C=1.2
            N3_0=float(i.split()[1])
            N3_ERR=float(i.split()[2])
            N3_MIN=float(i.split()[3])
            N3_MAX=float(i.split()[4])
    elif re.match(r'^z\s',i):
        Z_0=float(i.split()[1])
#        Z_MIN=float(i.split()[2])
#        Z_MAX=float(i.split()[3])
    elif re.match(r'^rho\s',i):
        RHO_0=float(i.split()[1])
        RHO_ERR=float(i.split()[2])
    elif re.match(r'^a\s',i):
        A_0=float(i.split()[1])
        A_ERR=float(i.split()[2])
    elif re.match(r'^beta\s',i):
        BETA_0=float(i.split()[1])
        BETA_MIN=float(i.split()[2])
        BETA_MAX=float(i.split()[3])
    elif re.match(r'^nt\s',i):
        NT_0=float(i.split()[1])
        NT_ERR=float(i.split()[2])
    elif re.match(r'^nm\s',i):
        NM_0=float(i.split()[1])
        NM_MIN=float(i.split()[2])
        NM_MAX=float(i.split()[3])
    elif re.match(r'^R200\s',i):
        R200_0=float(i.split()[1])
        R200_MIN=float(i.split()[2])
        R200_MAX=float(i.split()[3])
        R200M=R200_0
    elif re.match(r'^R500\s',i):
        R500_0=float(i.split()[1])
        R500_MIN=float(i.split()[2])
        R500_MAX=float(i.split()[3])
    elif re.match(r'^x\s',i):
        X_0=float(i.split()[1])
        X_MIN=float(i.split()[2])
        X_MAX=float(i.split()[3])
    elif re.match(r'^alpha\s',i):
        alpha_0=float(i.split()[1])
        alpha_err=float(i.split()[2])
    elif re.match(r'^nec\s',i):
        NEC_0=float(i.split()[1])
        NEC_ERR=float(i.split()[2])
    elif re.match(r'^TC\s',i):
        TC_0=float(i.split()[1])
        TC_ERR=float(i.split()[2])
    elif re.match(r'nth_A\s',i):
        NTH_A_0=float(i.split()[1])
        NTH_A_ERR=float(i.split()[2])
    elif re.match(r'nth_B\s',i):
        NTH_B_0=float(i.split()[1])
        NTH_B_ERR=float(i.split()[2])
    elif re.match(r'nth_gamma\s',i):
        NTH_GAMMA_0=float(i.split()[1])
        NTH_GAMMA_ERR=float(i.split()[2])
    elif re.match(r'sbp_c\s',i):
        SBP_C0=float(i.split()[1])
        SBP_CERR=float(i.split()[2])
        SBP_CMIN=float(i.split()[3])
        SBP_CMAX=float(i.split()[4])
        if SBP_CMAX > sbp_array[-1]:
            SBP_CMAX=sbp_array[-1]
        FLAG_SBPC=1
    elif re.match(r'fit_factor\s',i):
        T_FACTOR=float(i.split()[1])
        SBP_FACTOR=float(i.split()[2])
        MASS_FACTOR=float(i.split()[3])
        FLAG_FIT=1
    elif re.match(r'^aa\s',i):
        aa=int(i.split()[1])
        bb=int(i.split()[2])
        cc=int(i.split()[3])
        FLAG_ITN=1
    elif re.match(r'^tau\s',i):
        TAU_0=float(i.split()[1])
        TAU_ERR=float(i.split()[2])
        TAU_MIN=float(i.split()[3])
        TAU_MAX=float(i.split()[4])
    elif re.match(r'^fg\s',i):
        FG_0=float(i.split()[1])
        FG_ERR=float(i.split()[2])
        FG_MIN=float(i.split()[3])
        if FG_MIN <=0.05:
            FG_MIN=0.05
        FG_MAX=float(i.split()[4])
    elif re.match(r'^s\s',i):
        S_0=float(i.split()[1])
        S_ERR=float(i.split()[2])
        S_MIN=float(i.split()[3])
        S_MAX=float(i.split()[4])
    elif re.match(r'^c4\s',i):
        C4_0=float(i.split()[1])
        C4_ERR=float(i.split()[2])
        C4_MIN=float(i.split()[3])
        C4_MAX=float(i.split()[4])
    elif re.match(r'clumping\s',i):
        cp_p0,cp_perr,cp_e0,cp_eerr,cp_g00,cp_g0err,cp_xmin,cp_xmax,cp_sigma0,cp_sigmaerr=np.array(i.split()[1:],dtype=float)
        FLAG_CLUMP=1
if FLAG_SBPC==0:
    SBP_C0=1e-12
    SBP_CERR=1
    SBP_CMIN=0
    SBP_CMAX=sbp_array[-1]*1.0
if FLAG_FIT==0:
    T_FACTOR=1
    SBP_FACTOR=1
    MASS_FACTOR=1
if FLAG_ITN==0:
    aa=30000
    bb=10000
    cc=50
if FLAG_CLUMP==0:
    cp_p0=-3.7
    cp_perr=1.0
    cp_e0=3.7
    cp_eerr=1.0
    cp_g00=0.1
    cp_g0err=3.0
    cp_xmin=0
    cp_xmax=0.1
    cp_sigma0=1e-4
    cp_sigmaerr=0.033
cfunc_ori_array=[]
r_cfunc_array=[]
cfunc_use_array=[]
cfunc_ori_array.append(0)
r_cfunc_array.append(0)
cfunc_use_array.append(0)
cfunc_cori_array=[]
r_cfunc_c_array=[]
cfunc_cuse_array=[]
cfunc_cori_array.append(0)
r_cfunc_c_array.append(0)
cfunc_cuse_array.append(0)
tmp1=list(range(1,11))
tmp2=list(range(11,41,3))
tmp3=list(range(41,3001,5))
tmp1.extend(tmp2)
tmp1.extend(tmp3)
rne_array=numpy.array(tmp1)
for i in open(sys.argv[4]):
    r_cfunc_array.append(float(i.split()[0]))
    cfunc_ori_array.append(float(i.split()[1]))
for i in open(sys.argv[6]):
    r_cfunc_c_array.append(float(i.split()[0]))
    cfunc_cori_array.append(float(i.split()[1]))
for i in rne_array:
    if i==0:
        continue
    for j in range(len(r_cfunc_array)):
        if r_cfunc_array[j]>=i:
            cfunc_use_array.append(cfunc_ori_array[j])
            cfunc_cuse_array.append(cfunc_cori_array[j])
            break
cfunc_use_array=numpy.array(cfunc_use_array)
cfunc_cuse_array=np.array(cfunc_cuse_array)
tmp_array=[]
tmp_array.append(0)
i_before=0
for i in rne_array:
    if i==rne_array[0]:
        continue
    tmp_array.append((i+i_before)/2)
    if i==rne_array[-1]:
        tmp_array.append(1.5*i-0.5*i_before)
    i_before=i

T0=pymc.TruncatedNormal('T0',0,1/0.15**2,T0_MIN,T0_MAX,value=T0_0)
n2=pymc.TruncatedNormal('n2',N2_0,1/np.square(N2_ERR),0.1,2.5,value=N2_0)
rs=pymc.TruncatedNormal('rs',RS_0,1/np.square(RS_ERR),1,RS_0*3,value=RS_0)
a0=pymc.TruncatedNormal('a0',A0_0,1/np.square(A0_ERR),0,5,value=A0_0)
gamma0=pymc.TruncatedNormal('gamma0',GAMMA0_0,1/np.square(GAMMA0_ERR),0,2.0,value=GAMMA0_0)
k0=pymc.TruncatedNormal('k0',K0_0,1/np.square(K0_ERR),K0_MIN,K0_MAX,value=K0_0)
n3=pymc.TruncatedNormal('n3',N3_C,1/np.square(N3_ERR),N3_MIN,N3_MAX,value=N3_0)
rho=pymc.TruncatedNormal('rho',RHO_0,1/np.square(RHO_ERR),0,np.inf,value=RHO_0)
ne0=pymc.TruncatedNormal('ne0',NEC_0,1/np.square(NEC_ERR),1e-5,0.5,value=NEC_0)
kmod_a=pymc.TruncatedNormal('kmod_a',KMOD_A0,1/np.square(KMOD_AERR),0,5,value=KMOD_A0)
kmod_b=pymc.TruncatedNormal('kmod_b',KMOD_B0,1/np.square(KMOD_BERR),0.1,2,value=KMOD_B0)
kmod_c=pymc.TruncatedNormal('kmod_c',KMOD_C0,1/np.square(KMOD_CERR),-1,3,value=KMOD_C0)
kmod_k0=pymc.TruncatedNormal('kmod_k0',KMOD_K0,1/np.square(KMOD_K0ERR),0,np.inf,value=KMOD_K0)
sbp_c=pymc.TruncatedNormal('sbp_c',SBP_C0,1/np.square(SBP_CERR),SBP_CMIN,SBP_CMAX,value=SBP_C0)
delta=pymc.TruncatedNormal('delta',1.0,1/(0.12*0.12),0.5,1.5,value=1.0)
delta2=pymc.TruncatedNormal('delta2',3.0,1/0.09,2,4.5,value=3.0)
nth_a=pymc.TruncatedNormal('nth_a',NTH_A_0,1/np.square(NTH_A_ERR),0.4,0.5,value=NTH_A_0)
nth_b=pymc.TruncatedNormal('nth_b',NTH_B_0,1/np.square(NTH_B_ERR),0,np.inf,value=NTH_B_0)
nth_gamma=pymc.TruncatedNormal('nth_gamma',NTH_GAMMA_0,1/np.square(NTH_GAMMA_ERR),0,np.inf,value=NTH_GAMMA_0)
cp_p=pymc.TruncatedNormal('cp_p',cp_p0,1/np.square(cp_perr),-10,10,value=cp_p0)
cp_e=pymc.TruncatedNormal('cp_e',cp_e0,1/np.square(cp_eerr),-10,10,value=cp_e0)
cp_g0=pymc.TruncatedNormal('cp_g0',cp_g00,1/np.square(cp_g0err),0,10,value=cp_g00)
cp_x0=pymc.Uniform('cp_x0',lower=cp_xmin,upper=cp_xmax,value=(cp_xmin+cp_xmax)/2)
cp_sigma=pymc.TruncatedNormal('cp_sigma',cp_sigma0,1/np.square(cp_sigmaerr),1e-5,0.05,value=cp_sigma0)
c4=pymc.TruncatedNormal('c4',C4_0,1/np.square(C4_ERR),C4_MIN,C4_MAX,value=C4_0)
tau=pymc.TruncatedNormal('tau',TAU_0,1/np.square(TAU_ERR),TAU_MIN,TAU_MAX,value=TAU_0)
s=pymc.TruncatedNormal('s',S_0,1/np.square(S_ERR),S_MIN,S_MAX,value=S_0)
fg=pymc.TruncatedNormal('fg',FG_0,1/np.square(FG_ERR),FG_MIN,FG_MAX,value=FG_0)

##############
###MCMC fitting

####################
t_array=numpy.array(t_array)
te_array=numpy.array(te_array)
@pymc.stochastic
def temp_ne2(n2=n2,rs=rs,a0=a0,gamma0=gamma0,delta=delta,k0=k0,n3=n3,\
        rho=rho,T0=T0,ne0=ne0,sbp_c=sbp_c,kmod_a=kmod_a,kmod_b=kmod_b,\
        kmod_c=kmod_c,nth_a=nth_a,nth_b=nth_b,nth_gamma=nth_gamma,\
        kmod_k0=kmod_k0,delta2=delta2,cp_p=cp_p,cp_e=cp_e,cp_g0=cp_g0,\
        cp_x0=cp_x0,cp_sigma=cp_sigma,c4=c4,s=s,tau=tau,fg=fg,value=2):
# there is a bug in pymc2 that using Slice sampling method will sometimes break the upper or lower limit of that parameter. We check it manually here.
    if a0<=0:
        a0=0.01
        return np.nan
    if rs<=0:
        rs=1
        return np.nan
    if rho<=0:
        rho=1
        return np.nan
    if gamma0<=0:
        gamma0=0.1
        return np.nan
    if kmod_a<=0:
        kmod_a=0.01
        return np.nan
    if kmod_b<=0:
        kmod_b=0.1
        return np.nan
    lhood=0
    lhood=lhood+gpob((cp_e+cp_p)/(cp_e-cp_p),-0.01,0.03)
#    if kmod_a-a0<=0:
#        lhood=lhood+gpob((kmod_a-a0)/kmod_a,0,1)
    if rho<RHO_0:
        lhood=lhood+gpob(np.log(rho),np.log(RHO_0),1)
    else:
        lhood=lhood+gpob(np.log(RHO_0),np.log(RHO_0),1)
    if fg>0.17:
        lhood=lhood+gpob(fg,0.16,0.007)
    else:
        lhood=lhood+gpob(0.17,0.16,0.007)
    if fg<0.10:
        lhood=lhood+gpob(fg,0.10,0.007)
    else:
        lhood=lhood+gpob(0.10,0.10,0.007)
    y0=[ne0]
    p=[N1_0,n2,rs,a0,gamma0,k0,n3,rho,nth_a,nth_b,nth_gamma,delta,delta2,c4,tau,fg,s]
    flag_print=0
    if scipy.random.random()<0.02:
        flag_print=1
        print(p)
        print([kmod_a,kmod_b,kmod_c,kmod_k0])
        print([cp_p,cp_e,cp_g0,cp_x0,cp_sigma])
    T_array=[]
#    K_ARRAY_FIT=entropy_model(rne_array,kmod_a,kmod_b,kmod_c,kmod_k0)
    m_factor=[]
    p_MMnfw=[rho,rs,delta,delta2]
    p_den=[tau,fg,s]
    m_factor_tmp=[]
    den_array_fit=calc_den(rne_array,p_den,p_MMnfw)
    for j in rne_array:
        r=j
        e_ori=4*pi*rho*rs*rs*rs/r*np.log((rs+r)/rs)
        e_mod=modnfw_readarray.calc('n',r,p_MMnfw,t_total)/r+modnfw_readarray.calc('r',r,p_MMnfw,t_total)
        m_tmp=e_mod/e_ori
        m_factor_tmp.append(m_tmp)
    m_factor=numpy.array(m_factor_tmp)
    ne_array=den_array_fit/1.93
    T_array=calc_T(rne_array,ne_array,m_factor,p)
    k_fit_array=T_array*np.power(ne_array,-2/3)

    kmod_test=k_fit_array[312]
    kth_test=a0*1501**gamma0+k0
    if kth_test>kmod_test:
        lhood=lhood+gpob((kth_test-kmod_test)/kmod_test,0,0.2)
    kmod_test=k_fit_array[412]
    kth_test=a0*2001**gamma0+k0
    if kth_test>kmod_test:
        lhood=lhood+gpob((kth_test-kmod_test)/kmod_test,0,0.1)
    kmod_test=k_fit_array[512]
    kth_test=a0*2501**gamma0+k0
    if kth_test>kmod_test:
        lhood=lhood+gpob((kth_test-kmod_test)/kmod_test,0,0.03)
    kmod_test=k_fit_array[212]
    kth_test=a0*1001**gamma0+k0
    if kth_test>kmod_test:
        lhood=lhood+gpob((kth_test-kmod_test)/kmod_test,0,0.2)
    kmod_test=k_fit_array[112]
    kth_test=a0*501**gamma0+k0
    if kth_test>kmod_test:
        lhood=lhood+gpob((kth_test-kmod_test)/kmod_test,0,0.1)
    p_nth=[nth_a,nth_b,nth_gamma]
    m_array=calc_mass(rne_array,T_array,ne_array,p_nth)
    for i in range(len(m_array)):
        if i==0:
            continue
        if m_array[i]< m_array[i-1]:
            return np.nan
    ne_cl_array=clumping_model(rne_array/R200_0,cp_p,cp_e,cp_g0,cp_x0,cp_sigma)*ne_array
    for i in range(len(ne_cl_array)):
        if ne_cl_array[i]<ne_array[i]:
            ne_cl_array[i]=ne_array[i]
    ne_array=numpy.insert(ne_array,0,0.0)
    ne_cl_array=numpy.insert(ne_cl_array,0,0.0)
    T_array=numpy.insert(T_array,0,0.0)
    rtmp=numpy.array(list(range(45)))
    rmass_array=30*numpy.power(1.11,rtmp)
    for i in range(len(rmass_array)):
        rmass_array[i]=numpy.int(rmass_array[i])
    rmass_array=numpy.insert(rmass_array,0,10)
    rmass_array=numpy.insert(rmass_array,1,20)
    for i in range(len(rmass_array)):
        if rmass_array[i]>=10:
            for tmp in range(len(rne_array)):
                if rne_array[tmp]>=rmass_array[i]:
                    r_this=tmp
                    break
            Mnfw_model=Modefied_Mnfw(rmass_array[i],[rho,rs,delta,delta2])
            M_this=m_array[r_this]
            M_this_err=M_this*(te_array.sum()/t_array.sum()+0.15)
            if flag_print==1:
                print(rmass_array[i],Mnfw_model,M_this,M_this_err,Mnfw_model-M_this)
            if rmass_array[i]<=rsbp_array[-1]:
                if rmass_array[i]>=50:
                    lhood=lhood+gpob(Mnfw_model,M_this,M_this_err)*1*MASS_FACTOR
                elif rmass_array[i]>=10:
                    lhood=lhood+gpob(Mnfw_model,M_this,M_this_err)*0.95*MASS_FACTOR
            elif rmass_array[i]>R200_0:
                lhood=lhood+gpob(Mnfw_model,M_this,M_this_err)*0.95*MASS_FACTOR
            elif rmass_array[i]>R200_0/1.5:
                lhood=lhood+gpob(Mnfw_model,M_this,M_this_err)*0.95*MASS_FACTOR
            else:
                lhood=lhood+gpob(Mnfw_model,M_this,M_this_err)*1.0*MASS_FACTOR
    for i in range(len(r_array)):
        for j in range(len(rne_array)):
            if rne_array[j]>r_array[i]:
                r_this=j
                break
        if flag_tproj_array[i]=='1':
            t_this=T_array[r_this]
        else:
            t_this=deproject_model.calc_projT(r_array[i],tmp_array,T_array,ne_array)
        te_tmp=te_array[i]+t_array[i]*0.08
        if flag_print==1:
            print(t_this,t_array[i])
        lhood=lhood+gpob(t_this,t_array[i],te_tmp)*1.3*T_FACTOR
    for i in range(len(rsbp_array)):
        sbp_this=deproject_model.calc_sb(rsbp_array[i],tmp_array,ne_cl_array,cfunc_use_array)
        sbp_this=abs(sbp_this)
        sbp_this=sbp_this+sbp_c
        sbp_this=abs(sbp_this)
        tmp_sbp=abs(sbp_array[i])
#        tmp_sbpe=abs(numpy.log(1+sbpe_array[i]/sbp_array[i]))+0.05
        tmp_sbpe=tmp_sbp*(sbpe_array[i]/sbp_array[i]+0.1)
        if flag_print==1:
            print(rsbp_array[i],sbp_this,tmp_sbp,tmp_sbpe,sbp_this-tmp_sbp)
        if rsbp_array[i]>100:
            lhood=lhood+gpob(sbp_this,tmp_sbp,tmp_sbpe)*1*SBP_FACTOR
        else:
            lhood=lhood+gpob(sbp_this,tmp_sbp,tmp_sbpe)*0.5*SBP_FACTOR
    for i in range(len(rcsbp_array)):
        sbp_this=deproject_model.calc_sb(rcsbp_array[i],tmp_array,ne_cl_array,cfunc_cuse_array)
        sbp_this=abs(sbp_this)
#        sbp_this=sbp_this+sbp_c*0
#        sbp_this=abs(sbp_this)
        tmp_sbp=abs(csbp_array[i])
#        tmp_sbpe=abs(numpy.log(1+sbpe_array[i]/sbp_array[i]))+0.05
        tmp_sbpe=tmp_sbp*(csbpe_array[i]/csbp_array[i]+0.1)
        if flag_print==1:
            print(rcsbp_array[i],sbp_this,tmp_sbp,tmp_sbpe,sbp_this-tmp_sbp)
        if rcsbp_array[i]>100:
            lhood=lhood+gpob(sbp_this,tmp_sbp,tmp_sbpe)*1*SBP_FACTOR
        else:
            lhood=lhood+gpob(sbp_this,tmp_sbp,tmp_sbpe)*0.5*SBP_FACTOR
    return lhood
M2=pymc.MCMC(set([T0,n2,rs,a0,gamma0,k0,n3,rho,ne0,sbp_c,temp_ne2,nth_a,nth_b,nth_gamma,kmod_a,kmod_b,kmod_c,kmod_k0,delta,delta2,cp_p,cp_e,cp_g0,cp_x0,cp_sigma,c4,tau,fg,s]),db='pickle',dbname='sampled.pickle')
M2.db
#M2.use_step_method(pymc.Slicer,gamma0,w=0.1,doubling=True)
#M2.use_step_method(pymc.Slicer,rs,w=10,doubling=True)
#M2.use_step_method(pymc.Slicer,rho,w=50000,doubling=True)
#M2.use_step_method(pymc.Slicer,kmod_b,w=0.1,doubling=True)
#M2.use_step_method(pymc.Slicer,kmod_a)
#M2.use_step_method(pymc.Slicer,a0,w=0.1,doubling=True)
#M2.use_step_method(pymc.Slicer,delta2)
#M2.use_step_method(pymc.Slicer,k0)
pr=cProfile.Profile()
pr.enable()
M2.sample(iter=aa,burn=bb,thin=cc)
pr.disable()
p=pstats.Stats(pr)
p.sort_stats('cumulative').print_stats(100)
median=numpy.int((aa-bb)/cc/2)
NAME=sys.argv[1][:-4]
conchk_raft=pymc.raftery_lewis(M2,q=0.025,r=0.02,s=0.95,epsilon=0.003,verbose=2)
fi=open(NAME+'_raftery.txt','w')
print(conchk_raft,file=fi)
fi.close()
#pymc.Matplot.autocorrelation(M2)
#gel=pymc.gelman_rubin(M2)
#fi=open(NAME+'_gelman.txt','w')
#print(gel,file=fi)
#fi.close()
M2.write_csv("result.csv")
t_array=list(t_array)
te_array=list(te_array)
n1_f=N1_0
n2_f=M2.trace('n2')[:]
rs_f=M2.trace('rs')[:]
a0_f=M2.trace('a0')[:]
gamma0_f=M2.trace('gamma0')[:]
k0_f=M2.trace('k0')[:]
n3_f=M2.trace('n3')[:]
rho_f=M2.trace('rho')[:]
ne0_f=M2.trace('ne0')[:]
T0_f=M2.trace('T0')[:]
sbp_c_f=M2.trace('sbp_c')[:]
delta_f=M2.trace('delta')[:]
delta2_f=M2.trace('delta2')[:]
nth_a_f=M2.trace('nth_a')[:]
nth_b_f=M2.trace('nth_b')[:]
nth_gamma_f=M2.trace('nth_gamma')[:]
kmod_a_f=M2.trace('kmod_a')[:]
kmod_b_f=M2.trace('kmod_b')[:]
kmod_c_f=M2.trace('kmod_c')[:]
kmod_k0_f=M2.trace('kmod_k0')[:]
cp_p_f=M2.trace('cp_p')[:]
cp_e_f=M2.trace('cp_e')[:]
cp_g0_f=M2.trace('cp_g0')[:]
cp_x0_f=M2.trace('cp_x0')[:]
cp_sigma_f=M2.trace('cp_sigma')[:]
c4_f=M2.trace('c4')[:]
tau_f=M2.trace('tau')[:]
s_f=M2.trace('s')[:]
fg_f=M2.trace('fg')[:]
M2.db.close()
ne_array=[]
T_array=[]
p_all_f=numpy.array([n2_f,rs_f,a0_f,gamma0_f,k0_f,n3_f,rho_f,ne0_f,T0_f,sbp_c_f,delta_f,delta2_f,nth_a_f,nth_b_f,nth_gamma_f,kmod_a_f,kmod_b_f,kmod_c_f,kmod_k0_f,cp_p_f,cp_e_f,cp_g0_f,cp_x0_f,cp_sigma_f,c4_f,tau_f,fg_f,s_f])
numpy.save('p_all',p_all_f)
##############################
#print result

SUM_T_array=[]
SUM_Tproj_array=[]
SUM_ne_array=[]
SUM_sbp_fit=[]
SUM_csbp_fit=[]
SUM_k_array=[]
SUM_kfit_array=[]
SUM_nfw_fit_array=[]
SUM_mass_array=[]
SUM_ne_cl_array=[]
for i in range(len(ne0_f)):
    print(i)
    p=[n1_f,n2_f[i],rs_f[i],a0_f[i],gamma0_f[i],k0_f[i],n3_f[i],rho_f[i],\
            nth_a_f[i],nth_b_f[i],nth_gamma_f[i],delta_f[i],delta2_f[i],c4_f[i],tau_f[i],fg_f[i],s_f[i]]
#    K_ARRAY_FIT=entropy_model(rne_array,kmod_a_f[i],kmod_b_f[i],kmod_c_f[i],kmod_k0_f[i])
    p_nth=[nth_a_f[i],nth_b_f[i],nth_gamma_f[i]]
#    m_factor=1
    m_factor=[]
    p_MMnfw=[rho_f[i],rs_f[i],delta_f[i],delta2_f[i]]
    p_den=[tau_f[i],fg_f[i],s_f[i]]
    for j in range(len(rne_array)):
        r=rne_array[j]
#        e_ori=quad(nfw,0,r,p_MMnfw)[0]/r+quad(nfw_divbyr,r,numpy.inf,p_MMnfw)[0]
        e_ori=4*pi*rho_f[i]*np.power(rs_f[i],3)/r*np.log((rs_f[i]+r)/rs_f[i])
        e_mod=modnfw_readarray.calc('n',r,p_MMnfw,t_total)/r+modnfw_readarray.calc('r',r,p_MMnfw,t_total)
#        e_mod=quad(mod_nfw,0,r,p_MMnfw)[0]/r+quad(mod_nfw_divbyr,r,numpy.inf,p_MMnfw)[0]
        m_tmp=e_mod/e_ori
        m_factor.append(m_tmp)
    m_factor=numpy.array(m_factor)
    gden_array=calc_den(rne_array,p_den,p_MMnfw)
    ne_array=gden_array/1.93
    T_array=calc_T(rne_array,ne_array,m_factor,p)
    SUM_T_array.append(T_array)
#    ne_array=numpy.power(K_ARRAY_FIT/T_array,-1.5)
    SUM_ne_array.append(ne_array)
    mass_array=calc_mass(rne_array,T_array,ne_array,p_nth)
    SUM_mass_array.append(mass_array)
    nfw_fitted_array=[]
    k_array_fitted=T_array*np.power(ne_array,-2/3)
    for j in range(len(rne_array)):
        nfw_fitted=Modefied_Mnfw(rne_array[j],[rho_f[i],rs_f[i],delta_f[i],delta2_f[i]])
        nfw_fitted_array.append(nfw_fitted)
    SUM_nfw_fit_array.append(nfw_fitted_array)
    ne_cl_array=clumping_model(rne_array/R200_0,cp_p_f[i],cp_e_f[i],cp_g0_f[i],cp_x0_f[i],cp_sigma_f[i])*ne_array
    for j in range(len(ne_cl_array)):
        if ne_cl_array[j]<ne_array[j]:
            ne_cl_array[j]=ne_array[j]
    ne_cl_array=numpy.insert(ne_cl_array,0,0.0)
    SUM_ne_cl_array.append(ne_cl_array)
    ne_array=numpy.insert(ne_array,0,0.0)
    T_array=numpy.insert(T_array,0,0.0)
    sbp_fit=[]
    t2d_array=[]
    csbp_fit=[]
    for j in rne_array:
        a=deproject_model.calc_sb(j,tmp_array,ne_cl_array,cfunc_use_array)
        a=a+sbp_c_f[i]
        t2d=deproject_model.calc_projT(j,tmp_array,T_array,ne_array)
        sbp_fit.append(a)
        t2d_array.append(t2d)
        b=deproject_model.calc_sb(j,tmp_array,ne_cl_array,cfunc_cuse_array)
        csbp_fit.append(b)
    SUM_sbp_fit.append(sbp_fit)
    SUM_Tproj_array.append(t2d_array)
    SUM_csbp_fit.append(csbp_fit)
#    k_array_fitted=entropy_model(rne_array,kmod_a_f[i],kmod_b_f[i],kmod_c_f[i],kmod_k0_f[i])
    SUM_kfit_array.append(k_array_fitted)
IND_10=numpy.int(numpy.round((aa-bb)/cc*0.1))
IND_16=numpy.int(numpy.round((aa-bb)/cc*0.16))
IND_50=numpy.int(numpy.round((aa-bb)/cc*0.5))
IND_84=numpy.int(numpy.round((aa-bb)/cc*0.84))
IND_90=numpy.int(numpy.round((aa-bb)/cc*0.9))
T_10_ARRAY=[]
T_16_ARRAY=[]
T_50_ARRAY=[]
T_84_ARRAY=[]
T_90_ARRAY=[]
TPROJ_10_ARRAY=[]
TPROJ_16_ARRAY=[]
TPROJ_50_ARRAY=[]
TPROJ_84_ARRAY=[]
TPROJ_90_ARRAY=[]
DEN_10_ARRAY=[]
DEN_16_ARRAY=[]
DEN_50_ARRAY=[]
DEN_84_ARRAY=[]
DEN_90_ARRAY=[]
DEN_CL_10_ARRAY=[]
DEN_CL_16_ARRAY=[]
DEN_CL_50_ARRAY=[]
DEN_CL_84_ARRAY=[]
DEN_CL_90_ARRAY=[]
K_10_ARRAY=[]
K_16_ARRAY=[]
K_50_ARRAY=[]
K_84_ARRAY=[]
K_90_ARRAY=[]
KFIT_10_ARRAY=[]
KFIT_16_ARRAY=[]
KFIT_50_ARRAY=[]
KFIT_84_ARRAY=[]
KFIT_90_ARRAY=[]
M_10_ARRAY=[]
M_16_ARRAY=[]
M_50_ARRAY=[]
M_84_ARRAY=[]
M_90_ARRAY=[]
MFIT_10_ARRAY=[]
MFIT_16_ARRAY=[]
MFIT_50_ARRAY=[]
MFIT_84_ARRAY=[]
MFIT_90_ARRAY=[]
SBP_10_ARRAY=[]
SBP_16_ARRAY=[]
SBP_50_ARRAY=[]
SBP_84_ARRAY=[]
SBP_90_ARRAY=[]
CSBP_10_ARRAY=[]
CSBP_16_ARRAY=[]
CSBP_50_ARRAY=[]
CSBP_84_ARRAY=[]
CSBP_90_ARRAY=[]
for i in range(len(rne_array)):
    tmp_t_array=[]
    tmp_k_array=[]
    tmp_den_array=[]
    tmp_den_cl_array=[]
    tmp_sbp_array=[]
    tmp_csbp_array=[]
    tmp_kfit_array=[]
    tmp_m_array=[]
    tmp_mfit_array=[]
    tmp_t2d_array=[]
    for j in range(len(ne0_f)):
        tmp_t_array.append(SUM_T_array[j][i])
        tmp_t2d_array.append(SUM_Tproj_array[j][i])
        tmp_den_array.append(SUM_ne_array[j][i])
        tmp_den_cl_array.append(SUM_ne_cl_array[j][i])
        tmp_sbp_array.append(SUM_sbp_fit[j][i])
        tmp_csbp_array.append(SUM_csbp_fit[j][i])
        tmp_kfit_array.append(SUM_kfit_array[j][i])
        tmp_m_array.append(SUM_mass_array[j][i])
        tmp_mfit_array.append(SUM_nfw_fit_array[j][i])
    tmp_t_array.sort()
    tmp_t2d_array.sort()
    tmp_den_array.sort()
    tmp_den_cl_array.sort()
    tmp_sbp_array.sort()
    tmp_csbp_array.sort()
    tmp_kfit_array.sort()
    tmp_m_array.sort()
    tmp_mfit_array.sort()
    M_10_ARRAY.append(tmp_m_array[IND_10])
    M_16_ARRAY.append(tmp_m_array[IND_16])
    M_50_ARRAY.append(tmp_m_array[IND_50])
    M_84_ARRAY.append(tmp_m_array[IND_84])
    M_90_ARRAY.append(tmp_m_array[IND_90])
    MFIT_10_ARRAY.append(tmp_mfit_array[IND_10])
    MFIT_16_ARRAY.append(tmp_mfit_array[IND_16])
    MFIT_50_ARRAY.append(tmp_mfit_array[IND_50])
    MFIT_84_ARRAY.append(tmp_mfit_array[IND_84])
    MFIT_90_ARRAY.append(tmp_mfit_array[IND_90])
    T_10_ARRAY.append(tmp_t_array[IND_10])
    T_16_ARRAY.append(tmp_t_array[IND_16])
    T_50_ARRAY.append(tmp_t_array[IND_50])
    T_84_ARRAY.append(tmp_t_array[IND_84])
    T_90_ARRAY.append(tmp_t_array[IND_90])
    TPROJ_10_ARRAY.append(tmp_t2d_array[IND_10])
    TPROJ_16_ARRAY.append(tmp_t2d_array[IND_16])
    TPROJ_50_ARRAY.append(tmp_t2d_array[IND_50])
    TPROJ_84_ARRAY.append(tmp_t2d_array[IND_84])
    TPROJ_90_ARRAY.append(tmp_t2d_array[IND_90])
    SBP_10_ARRAY.append(tmp_sbp_array[IND_10])
    SBP_16_ARRAY.append(tmp_sbp_array[IND_16])
    SBP_50_ARRAY.append(tmp_sbp_array[IND_50])
    SBP_84_ARRAY.append(tmp_sbp_array[IND_84])
    SBP_90_ARRAY.append(tmp_sbp_array[IND_90])
    CSBP_10_ARRAY.append(tmp_csbp_array[IND_10])
    CSBP_16_ARRAY.append(tmp_csbp_array[IND_16])
    CSBP_50_ARRAY.append(tmp_csbp_array[IND_50])
    CSBP_84_ARRAY.append(tmp_csbp_array[IND_84])
    CSBP_90_ARRAY.append(tmp_csbp_array[IND_90])
    KFIT_10_ARRAY.append(tmp_kfit_array[IND_10])
    KFIT_16_ARRAY.append(tmp_kfit_array[IND_16])
    KFIT_50_ARRAY.append(tmp_kfit_array[IND_50])
    KFIT_84_ARRAY.append(tmp_kfit_array[IND_84])
    KFIT_90_ARRAY.append(tmp_kfit_array[IND_90])
    DEN_10_ARRAY.append(tmp_den_array[IND_10])
    DEN_16_ARRAY.append(tmp_den_array[IND_16])
    DEN_50_ARRAY.append(tmp_den_array[IND_50])
    DEN_84_ARRAY.append(tmp_den_array[IND_84])
    DEN_90_ARRAY.append(tmp_den_array[IND_90])
    DEN_CL_10_ARRAY.append(tmp_den_cl_array[IND_10])
    DEN_CL_16_ARRAY.append(tmp_den_cl_array[IND_16])
    DEN_CL_50_ARRAY.append(tmp_den_cl_array[IND_50])
    DEN_CL_84_ARRAY.append(tmp_den_cl_array[IND_84])
    DEN_CL_90_ARRAY.append(tmp_den_cl_array[IND_90])
NAME=sys.argv[1][:-9]
tmp='array_plt.json'
fi=open(tmp,'w')
radius_model=[]
for i in rne_array:
    radius_model.append(str(i))
dat={
        "name": NAME,
        "radius": [r_array,re_array],
        "temperature": [t_array,te_array],
        "sbp": [sbp_array,sbpe_array],
        "radius_model": radius_model,
        "temperature_model": [T_50_ARRAY,T_16_ARRAY,T_84_ARRAY,T_10_ARRAY,T_90_ARRAY],
        "projected_temperature": [TPROJ_50_ARRAY,TPROJ_16_ARRAY,TPROJ_84_ARRAY,TPROJ_10_ARRAY,TPROJ_90_ARRAY],
        "sbp_model": [SBP_50_ARRAY,SBP_16_ARRAY,SBP_84_ARRAY,SBP_10_ARRAY,SBP_90_ARRAY],
        "csbp_model": [CSBP_50_ARRAY,CSBP_16_ARRAY,CSBP_84_ARRAY],
        "den_model": [DEN_50_ARRAY,DEN_16_ARRAY,DEN_84_ARRAY,DEN_10_ARRAY,DEN_90_ARRAY],
        "den_cl_model": [DEN_CL_50_ARRAY,DEN_CL_16_ARRAY,DEN_CL_84_ARRAY,DEN_CL_10_ARRAY,DEN_CL_90_ARRAY],
        "m_calc": [M_50_ARRAY,M_16_ARRAY,M_84_ARRAY,M_10_ARRAY,M_90_ARRAY],
        "m_fit": [MFIT_50_ARRAY,MFIT_16_ARRAY,MFIT_84_ARRAY,M_10_ARRAY,M_90_ARRAY],
        "k_fit": [KFIT_50_ARRAY,KFIT_16_ARRAY,KFIT_84_ARRAY,KFIT_10_ARRAY,KFIT_90_ARRAY],
        "r200": R200_0,
        }
json.dump(dat,fi,indent=2)
fi.close()
K_50_ARRAY=numpy.array(K_50_ARRAY)
sbpe_array=numpy.array(sbpe_array)
sbp_array=numpy.array(sbp_array)
csbp_array=np.array(csbp_array)
csbpe_array=np.array(csbpe_array)
t_array=numpy.array(t_array)
te_array=numpy.array(te_array)
sbpe_array=sbp_array*0.05+sbpe_array
te_array=t_array*0.12+te_array

plt.clf()
plt.loglog(rne_array,SBP_84_ARRAY,color='grey')
plt.loglog(rne_array,SBP_50_ARRAY,'b')
plt.loglog(rne_array,SBP_16_ARRAY,color='grey')
plt.fill_between(rne_array,SBP_16_ARRAY,SBP_84_ARRAY,color='grey')
plt.errorbar(rsbp_array,sbp_array,xerr=rsbpe_array,yerr=sbpe_array,color='k',linestyle='none')
plt.xlabel(name+'_Radius(kpc)')
plt.ylabel('SBP(cm^-2 pixel^-2 s^-1)')
plt.savefig('sbp.pdf',dpi=100)


plt.clf()
plt.loglog(rne_array,CSBP_84_ARRAY,color='grey')
plt.loglog(rne_array,CSBP_50_ARRAY,'b')
plt.loglog(rne_array,CSBP_16_ARRAY,color='grey')
plt.fill_between(rne_array,CSBP_16_ARRAY,CSBP_84_ARRAY,color='grey')
plt.errorbar(rcsbp_array,csbp_array,xerr=rcsbpe_array,yerr=csbpe_array,color='k',linestyle='none')
plt.xlabel(name+'_Radius(kpc)')
plt.ylabel('SBP(cm^-2 pixel^-2 s^-1)')
plt.savefig('csbp.pdf',dpi=100)
plt.clf()

plt.plot(rne_array,T_84_ARRAY,color='grey')
plt.plot(rne_array,T_50_ARRAY,'b')
plt.plot(rne_array,T_16_ARRAY,color='grey')
plt.plot(rne_array,TPROJ_50_ARRAY,color='yellow')
plt.fill_between(rne_array,T_16_ARRAY,T_84_ARRAY,color='grey')
plt.fill_between(rne_array,TPROJ_16_ARRAY,TPROJ_84_ARRAY,color='green',alpha=0.5)
r1_array=[]
r2_array=[]
t2_array=[]
t1_array=[]
re1_array=[]
te1_array=[]
re2_array=[]
te2_array=[]
for i in range(len(r_array)):
    if flag_tproj_array[i]=='1':
        r1_array.append(r_array[i])
        t1_array.append(t_array[i])
        re1_array.append(re_array[i])
        te1_array.append(te_array[i])
    else:
        r2_array.append(r_array[i])
        t2_array.append(t_array[i])
        re2_array.append(re_array[i])
        te2_array.append(te_array[i])

plt.errorbar(r1_array,t1_array,xerr=re1_array,yerr=te1_array,color='k',linestyle='none')
plt.errorbar(r2_array,t2_array,xerr=re2_array,yerr=te2_array,color='orange',linestyle='none')
plt.xlabel(name+'_Radius(kpc)')
plt.ylabel('Temperature(keV)')
plt.savefig('temperature.pdf',dpi=100)


plt.clf()
plt.loglog(rne_array,KFIT_50_ARRAY,'g')
plt.fill_between(rne_array,KFIT_16_ARRAY,KFIT_84_ARRAY,color='green',alpha=0.5)
plt.xlabel(name+'_Radius(kpc)')
plt.ylabel('Entropy(kev cm^2)')
plt.savefig('entropy.pdf',dpi=100)

plt.clf()
plt.loglog(rne_array,M_50_ARRAY,'b')
plt.fill_between(rne_array,M_16_ARRAY,M_84_ARRAY,color='grey')
plt.fill_between(rne_array,MFIT_16_ARRAY,MFIT_84_ARRAY,color='green',alpha=0.5)
plt.xlim(10,3000)
plt.ylim(1e11,1e16)
plt.xlabel(name+'_Radius(kpc)')
plt.ylabel('Mass(Msun))')
plt.savefig('mass.pdf',dpi=100)


###############
