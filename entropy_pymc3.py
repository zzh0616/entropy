#!/usr/bin/env python3
## $1 temperature data; radius in kpc
## $2 temperature parameter file
## $3 sbp data; 'radius err sbp err' radius in kpc
## $4 cooling function file , radius in kpc, radius should match density
## $5 cm_per_pixel for this source
## $6 fitted entropy profile for this source !!not used for this version
## $7 outfile name for fitted surface brightness
## $8 outfile name for fitted entropy
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
import matplotlib.pyplot as plt
import sys
import math
import numpy
import numpy as np
import scipy
import matplotlib
import re
import types
import pymc3 as pm
import json
import theano.tensor as tt
import theano
numpy.set_printoptions(linewidth=100000)
#theano.config.exception_verbosity='high'
pi=3.1415926
cm_per_pixel=np.float(sys.argv[5])
name=sys.argv[1][0:-4]
print(name)
print(cm_per_pixel)
G=6.67e-11 #m^3 kg^-1 s^-2
mp=1.67e-27  #kg
kev=1.6e-16 #J
kpc=3.086e19 #m
Msun=2e30 #kg
mu=0.61
##p=[rho,rs,delta,delta2]
def mod_nfw(r,p):
    rho=p[0]
    rs=p[1]
    delta=p[2]
    delta2=p[3]
    den=rho/(tt.power(r/rs,delta)*(1+tt.power(r/rs,delta2-delta)))*4*pi*r*r
    return den

def Modefied_Mnfw(r,p):
    x=tt.dscalar('x')
    x.get.test_value=1.0
    y=mod_nfw(x,p)
    f=theano.function([x],y)
    M=quad(f,0,r,p)
    return M[0]

def mod_nfw_divbyr(r,p):
    a=mod_nfw(r,p)/r
    return a
def nfw(r,p):
    new_p=[p[0],p[1],1.0,3.0]
    a=mod_nfw(r,new_p)
    return a
def nfw_divbyr(r,p):
    a=nfw(r,p)
    a=a/r
    return a


## NOT USED in this version
def dndr(ne,r0,p):
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
#    a1=p[12]
#    gamma1=p[13]
#    k1=p[14]
    nth_a=p[8]
    nth_b=p[9]
    nth_gamma=p[10]
#    delta=p[8]
    r200=R200_0
    r500=R500_0
#    a=p[9]
#    beta=BETA_0
#    nt=p[10]
#    nm=p[11]
    z=Z_0
    delta=p[11]
    delta2=p[12]
    T_array=p[13]
    i=numpy.int(numpy.round(r0/5))
    p_MMnfw=[rho,rs,delta,delta2]
    if i >=400:
#        print(i)
        i=399
    if i==0:
        i=1
    T=T_array[i]
    dTdr=(T_array[i]-T_array[i-1])/5
    Eg_ori=quad(nfw,0,r0,p_MMnfw)
    Eg_modified=quad(mod_nfw,0,r0,p_MMnfw)
    if Eg_ori[0]==0:
        m_factor=0
    else:
        m_factor=Eg_modified[0]/Eg_ori[0]
    z1=m_factor*-4*pi*rho*rs*rs*rs*(numpy.log((r0+rs)/rs)-r0/(rs+r0))*G*mu*mp/r0/r0\
            *Msun/kpc/kev   #unit in kev/kpc same as z2
#    m200=4*pi*rho*pow(rs,3)*(math.log((rs+r200)/rs)-r200/(rs+r200))
#    eta=1-a*numpy.power(1+z,beta)*numpy.power(r0/r500,nt)*\
#            numpy.power(m200/3/numpy.power(10,14),nm)
#    detadr=-a*numpy.power(1+z,beta)*numpy.power(m200/3/pow(10,14),nm)*nt*\
#            numpy.power(r0/r500,nt-1)/r500
#    eta=1
#    detadr=0
    eta=nth_a*(1+numpy.exp(-numpy.power(r0/R200M/nth_b,nth_gamma)))
    detadr=nth_a*numpy.exp(-numpy.power(r0/R200M/nth_b,nth_gamma))*\
            nth_gamma*-numpy.power(r0/R200M/nth_b,nth_gamma-1)/R200M/nth_b
    dnedr=(z1*eta+T/eta*detadr-dTdr)*ne/T
    return dnedr
### length of x(radius) and k_fit(entropy) must be the same
def calc_T(x,k_fit,m_factor,p):
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
#    a1=p[12]
#    gamma1=p[13]
#    k1=p[14]
#    delta=p[8]
    r200=R200_0
    r500=R500_0
#    a=p[9]
#    beta=BETA_0
#    nt=p[10]
#    nm=p[11]
    nth_a=p[8]
    nth_b=p[9]
    nth_gamma=p[10]
    delta=p[11]
    delta2=p[12]
    p_MMnfw=[rho,rs,delta,delta2]
    z=Z_0
#    m200=4*pi*rho*pow(rs,3)*(math.log((rs+r200)/rs)-r200/(rs+r200))
#    eta=1-a*numpy.power(1+z,beta)*numpy.power(r0/r500,nt)*\
#        numpy.power(m200/3/numpy.power(10,14),nm)
#    eta=1
#    m200=Modefied_Mnfw(r200,p_MMnfw)
#    Eg_ori=quad(nfw,0,r0,p_MMnfw)[0]/r+quad(nfw_divbyr,r0,numpy.inf,p_MMnfw)[0]
 #   Eg_modified=quad(mod_nfw,0,r0,p_MMnfw)[0]/r+quad(mod_nfw_divbyr,r0,numpy.inf,p_MMnfw)[0]
 #   if Eg_ori==0:
 #       m_factor=0
 #   else:
 #       m_factor=Eg_modified/Eg_ori
#    m_factor=1
#    print(m_factor)
    eta=nth_a*(1+numpy.exp(-numpy.power(r0/R200M/nth_b,nth_gamma)))
    Tx=(T0_0+(1-c3)*c1*rho*rs*rs*rs*numpy.log((rs+x)/rs)/x*m_factor)/(1/eta-c3-c2+c2*k/k_fit)
    return Tx

def entropy_model(x,a,b,c,k0):
    x=numpy.array(x)
    return a*tt.power(x,b)*numpy.exp(c*(-x/R200_0))+k0
def calc_mass(r_array,t_array,ne_array,p_eta):
    mass_array=[]
    a=p_eta[0]
    b=p_eta[1]
    c=p_eta[2]
    for i in range(len(r_array)):
        r0=r_array[i]
        eta=a*(1+numpy.exp(-numpy.power(r0/R200M/b,c)))
        detadr=a*numpy.exp(-numpy.power(r0/R200M/b,c))*\
                c*-numpy.power(r0/R200M/b,c-1)/R200M/b
#        eta=1
#        detaer=0
        if i==0:
            dTdr=(t_array[i+1]-t_array[i])/(r_array[i+1]-r_array[i])
            dnedr=(ne_array[i+1]-ne_array[i])/(r_array[i+1]-r_array[i])
            z1=(dnedr*t_array[i]/ne_array[i]+dTdr-t_array[i]/eta*detadr)/eta #kev/kpc
            mass=-z1/G/mu/mp*r0*r0/Msun*kpc*kev #Msun
            mass_array.append(mass)
        elif i==len(r_array)-1:
            dTdr=(t_array[i]-t_array[i-1])/(r_array[i]-r_array[i-1])
            dnedr=(ne_array[i]-ne_array[i-1])/(r_array[i]-r_array[i-1])
            z1=(dnedr*t_array[i]/ne_array[i]+dTdr-t_array[i]/eta*detadr)/eta #kev/kpc
            mass=-z1/G/mu/mp*r0*r0/Msun*kpc*kev #Msun
            mass_array.append(mass)
        else:
            dTdr_r=(t_array[i+1]-t_array[i])/(r_array[i+1]-r_array[i])
            dnedr_r=(ne_array[i+1]-ne_array[i])/(r_array[i+1]-r_array[i])
            dTdr_l=(t_array[i]-t_array[i-1])/(r_array[i]-r_array[i-1])
            dnedr_l=(ne_array[i]-ne_array[i-1])/(r_array[i]-r_array[i-1])
            dTdr=(dTdr_l+dTdr_r)/2
            dnedr=(dnedr_l+dnedr_r)/2
            z1=(dnedr*t_array[i]/ne_array[i]+dTdr-t_array[i]/eta*detadr)/eta #kev/kpc
            mass=-z1/G/mu/mp*r0*r0/Msun*kpc*kev #Msun
#        print(mass)
        if mass<=1e7:
            mass=1e7
        mass_array.append(mass)
    return mass_array

class FIT_MODEL(pm.distributions.distribution.Continuous):
    def __init__(self,n2,rs,a0,gamma0,delta,k0,n3,rho,T0,ne0,sbp_c,kmod_a,kmod_b,kmod_c,nth_a,nth_b,nth_gamma,kmod_k0,delta2,*args,**kwargs):
        super(FIT_MODEL,self).__init__(*args,**kwargs)
        self.n2=tt.as_tensor_variable(n2)
        self.rs=tt.as_tensor_variable(rs)
        self.a0=tt.as_tensor_variable(a0)
        self.gamma0=tt.as_tensor_variable(gamma0)
        self.delta=tt.as_tensor_variable(delta)
        self.k0=tt.as_tensor_variable(k0)
        self.n3=tt.as_tensor_variable(n3)
        self.rho=tt.as_tensor_variable(rho)
        self.T0=tt.as_tensor_variable(T0)
        self.ne0=tt.as_tensor_variable(ne0)
        self.sbp_c=tt.as_tensor_variable(sbp_c)
        self.kmod_a=tt.as_tensor_variable(kmod_a)
        self.kmod_b=tt.as_tensor_variable(kmod_b)
        self.kmod_c=tt.as_tensor_variable(kmod_c)
        self.nth_a=tt.as_tensor_variable(nth_a)
        self.nth_b=tt.as_tensor_variable(nth_b)
        self.nth_gamma=tt.as_tensor_variable(nth_gamma)
        self.kmod_k0=tt.as_tensor_variable(kmod_k0)
        self.delta2=tt.as_tensor_variable(delta2)
        self.testval=0.0

#    def logp (n2=self.n2,rs=self.rs,a0=self.a0,gamma0=self.gamma0,delta=self.delta,k0=self.k0,n3=self.n3,rho=self.rho,T0=self.T0,ne0=self.ne0,sbp_c=self.sbp_c,kmod_a=self.kmod_a,kmod_b=self.kmod_b,kmod_c=self.kmod_c,nth_a=self.nth_a,nth_b=self.nth_b,nth_gamma=self.nth_gamma,kmod_k0=self.kmod_k0,delta2=self.delta2):
    def logp(self,value):
        n2=self.n2;rs=self.rs;a0=self.a0;gamma0=self.gamma0;delta=self.delta;k0=self.k0;n3=self.n3;rho=self.rho;T0=self.T0;ne0=self.ne0;sbp_c=self.sbp_c;kmod_a=self.kmod_a;kmod_b=self.kmod_b;kmod_c=self.kmod_c;nth_a=self.nth_a;nth_b=self.nth_b;nth_gamma=self.nth_gamma;kmod_k0=self.kmod_k0;delta2=self.delta2
        lhood=0
        rne_array=numpy.array(range(1,2001,5))
        tmp_array=numpy.array(range(0,2005,5))
        y0=[ne0]
        p=[N1_0,n2,rs,a0,gamma0,k0,n3,rho,nth_a,nth_b,nth_gamma,delta,delta2]
#        print(p)
#        print([kmod_a,kmod_b,kmod_c,kmod_k0])
        T1_array=[]
        K_ARRAY_FIT=entropy_model(rne_array,kmod_a,kmod_b,kmod_c,kmod_k0)
        m_factor=[]
        p_MMnfw=[rho,rs,delta,delta2]
        m_factor_tmp=[]
        for j in range(1,2201,100):
            r=j
            v=tt.scalar('v')
            v.tag.test_value=1
            y1=nfw(v,p_MMnfw)
            y2=y1/v
            f1=theano.function([v],y1)
            f2=theano.function([v],y2)
            e_ori=quad(f1,0,r)[0]/r+quad(f2,r,numpy.inf)[0]
            y1=mod_nfw(v,p_MMnfw)
            y2=mod_nfw_divbyr(v,p_MMnfw)
            f1=theano.function([v],y1)
            f2=theano.function([v],y2)
            e_mod=quad(f1,0,r,p_MMnfw)[0]/r+quad(f2,r,numpy.inf,p_MMnfw)[0]
            m_tmp=e_mod/e_ori
            m_factor_tmp.append(m_tmp)
        for i in range(1,2001,5):
            if i==1:
                m_factor.append(m_factor_tmp[0])
                continue
            for j in range(1,2201,100):
                if j>=i:
                    j_this=numpy.int((j-1)/100)
                    diff=m_factor_tmp[j_this]-m_factor_tmp[j_this-1]
                    m_factor.append(m_factor_tmp[j_this-1]+diff*(i-j+100)/100)
                    break
        m_factor=numpy.array(m_factor)
        print(m_factor_tmp)
        print(m_factor)
        T1_array=calc_T(rne_array,K_ARRAY_FIT,m_factor,p)
        ne_array=numpy.power(K_ARRAY_FIT/T1_array,-1.5)
        T_array=T1_array
        p_nth=[nth_a,nth_b,nth_gamma]
        m_array=calc_mass(rne_array,T_array,ne_array,p_nth)
        ne_array=numpy.insert(ne_array,0,0.0)
        rtmp=numpy.array(list(range(41)))
        rmass_array=30*numpy.power(1.11,rtmp)
        for i in range(len(rmass_array)):
            rmass_array[i]=numpy.int(rmass_array[i])
        rmass_array=numpy.insert(rmass_array,0,10)
        rmass_array=numpy.insert(rmass_array,1,20)
        for i in range(len(rmass_array)):
            if rmass_array[i]>=10:
                r_this=numpy.int(numpy.round(rmass_array[i]/5))
                Mnfw_model=numpy.log(Modefied_Mnfw(rmass_array[i],[rho,rs,delta,delta2]))
                M_this=numpy.log(m_array[r_this])
                M_this_err=0.5
                print(rmass_array[i],Mnfw_model,M_this,M_this_err,Mnfw_model-M_this)
                if rmass_array[i]<=rsbp_array[-1]:
                    if rmass_array[i]>=50:
                        lhood=lhood+gpob(Mnfw_model,M_this,M_this_err)*1
                    elif rmass_array[i]>=30:
                        lhood=lhood+gpob(Mnfw_model,M_this,M_this_err)*0.5
                elif rmass_array[i]>R200_0:
                    lhood=lhood+gpob(Mnfw_model,M_this,M_this_err)*0.4
                elif rmass_array[i]>200:
                    lhood=lhood+gpob(Mnfw_model,M_this,M_this_err)*1.0
                else:
                    lhood=lhood+gpob(Mnfw_model,M_this,M_this_err)*1.0
        for i in range(len(r_array)):
            r_this=numpy.int(numpy.round(r_array[i]/5))
            t_this=T_array[r_this]
            te_tmp=te_array[i]+t_array[i]*0.0
            print(t_this,t_array[i])
            lhood=lhood+gpob(t_this,t_array[i],te_array[i])*1.2
        for i in range(len(rsbp_array)):
            sbp_this=deproject_model.calc_sb(rsbp_array[i],tmp_array,ne_array,cfunc_use_array,cm_per_pixel)
            sbp_this=abs(sbp_this)
            sbp_this=sbp_this+sbp_c
            sbp_this=abs(np.log(sbp_this))
            tmp_sbp=abs(np.log(sbp_array[i]))
            tmp_sbpe=abs(numpy.log(1+sbpe_array[i]/sbp_array[i]))+0.1
            print(rsbp_array[i],sbp_this,tmp_sbp,tmp_sbpe,sbp_this-tmp_sbp)
            if rsbp_array[i]>10:
                lhood=lhood+gpob(sbp_this,tmp_sbp,tmp_sbpe)*1
            else:
                lhood=lhood+gpob(sbp_this,tmp_sbp,tmp_sbpe)*0.5
        return lhood

r_array=[]
re_array=[]
t_array=[]
te_array=[]
rsbp_array=[]
rsbpe_array=[]
sbp_array=[]
sbpe_array=[]
for i in open(sys.argv[1]):
    r,rer,t,te=i.split()
    r=float(r)
    rer=float(rer)
    t=float(t)
    te=float(te)
    r_array.append(r)
    re_array.append(rer)
    t_array.append(t)
    te_array.append(te)


#r_array=numpy.array(r_array)
#re_array=numpy.array(re_array)
#t_array=numpy.array(t_array)
#te_array=numpy.array(te_array)

for i in open(sys.argv[3]):
    r,rer,sbp,sbpe=i.split()
    r=float(r)
    rer=float(rer)
    sbp=float(sbp)
    sbpe=float(sbpe)
    rsbp_array.append(r)
    rsbpe_array.append(rer)
    sbp_array.append(sbp)
    sbpe_array.append(sbpe)
#rsbp_array=numpy.array(rsbp_array)
#rsbpe_array=numpy.array(rsbpe_array)
#sbp_array=numpy.array(sbp_array)
#sbpe_array=numpy.array(sbpe_array)
FLAG_SBPC=0
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
            N3_MIN=float(i.split()[2])
            N3_MAX=float(i.split()[3])
            N3_ERR=N3_MAX-N3_MIN
        elif len(i.split())==5:
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
    elif re.match(r'^r200m\s',i):
        R200M=float(i.split()[1])
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
        FLAG_SBPC=1
if FLAG_SBPC==0:
    SBP_C0=1e-11
    SBP_CERR=1
    SBP_CMIN=0
    SBP_CMAX=sbp_array[-1]*0.9

cfunc_ori_array=[]
r_cfunc_array=[]
cfunc_use_array=[]
cfunc_ori_array.append(0)
r_cfunc_array.append(0)
cfunc_use_array.append(0)
for i in open(sys.argv[4]):
    r_cfunc_array.append(float(i.split()[0]))
    cfunc_ori_array.append(float(i.split()[1]))
for i in range(1,2001,5):
    if i==0:
        continue
    for j in range(len(r_cfunc_array)):
        if r_cfunc_array[j]>i:
            cfunc_use_array.append(cfunc_ori_array[j])
            break
cfunc_use_array=numpy.array(cfunc_use_array)


model=pm.Model()
with model:
    T0=pm.Uniform('T0',lower=T0_MIN,upper=T0_MIN,testval=T0_0)
    n2=pm.Bound(pm.Normal,lower=0.1,upper=10)('n2',mu=N2_0,sd=N2_ERR,testval=N2_0)
    rs=pm.Bound(pm.Normal,lower=1,upper=RS_0*3)('rs',mu=RS_0,sd=RS_ERR,testval=RS_0)
    a0=pm.Bound(pm.Normal,lower=0,upper=np.inf)('a0',mu=A0_0,sd=A0_ERR,testval=A0_0)
    gamma0=pm.Bound(pm.Normal,lower=0,upper=np.inf)('gamma0',mu=GAMMA0_0,sd=GAMMA0_ERR,testval=GAMMA0_0)
    k0=pm.Bound(pm.Normal,lower=K0_MIN,upper=K0_MAX)('k0',mu=K0_0,sd=K0_ERR,testval=K0_0)
    n3=pm.Bound(pm.Normal,lower=N3_MIN,upper=N3_MAX)('n3',mu=N3_0,sd=N3_ERR,testval=N3_0)
    rho=pm.Bound(pm.Normal,lower=0,upper=np.inf)('rho',mu=RHO_0,sd=RHO_ERR,testval=RHO_0)
    ne0=pm.Bound(pm.Normal,lower=1e-5,upper=0.5)('ne0',mu=NEC_0,sd=NEC_ERR,testval=NEC_0)
    kmod_a=pm.Bound(pm.Normal,lower=0,upper=np.inf)('kmod_a',mu=KMOD_A0,sd=KMOD_AERR,testval=KMOD_A0)
    kmod_b=pm.Bound(pm.Normal,lower=0.1,upper=10)('kmod_b',mu=KMOD_B0,sd=KMOD_BERR,testval=KMOD_B0)
    kmod_c=pm.Bound(pm.Normal,lower=-10,upper=np.inf)('kmod_c',mu=KMOD_C0,sd=KMOD_CERR,testval=KMOD_C0)
    kmod_k0=pm.Bound(pm.Normal,lower=0,upper=np.inf)('kmod_k0',mu=KMOD_K0,sd=KMOD_K0ERR,testval=KMOD_K0)
    sbp_c=pm.Bound(pm.Normal,lower=SBP_CMIN,upper=SBP_CMAX)('sbp_c',mu=SBP_C0,sd=SBP_CERR,testval=SBP_C0)
    delta=pm.Bound(pm.Normal,lower=0,upper=2.99)('delta',mu=1.0,sd=0.5,testval=1.0)
    delta2=pm.Bound(pm.Normal,lower=2,upper=10)('delta2',mu=3.0,sd=0.55,testval=3.0)
    nth_a=pm.Bound(pm.Normal,lower=0.3,upper=0.5)('nth_a',mu=NTH_A_0,sd=NTH_A_ERR,testval=NTH_A_0)
    nth_b=pm.Bound(pm.Normal,lower=0,upper=np.inf)('nth_b',mu=NTH_B_0,sd=NTH_B_ERR,testval=NTH_B_0)
    nth_gamma=pm.Bound(pm.Normal,lower=0,upper=np.inf)('nth_gamma',mu=NTH_GAMMA_0,sd=NTH_GAMMA_ERR,testval=NTH_GAMMA_0)
    lh=FIT_MODEL('lh',n2=n2,rs=rs,a0=a0,gamma0=gamma0,delta=delta,k0=k0,n3=n3,rho=rho,T0=T0,ne0=ne0,sbp_c=sbp_c,kmod_a=kmod_a,kmod_b=kmod_b,kmod_c=kmod_c,nth_a=nth_a,nth_b=nth_b,nth_gamma=nth_gamma,kmod_k0=kmod_k0,delta2=delta2)
map_estimate=pm.find_MAP(model=model,fmin=scipy.optimize.fmin_powell)
print(map_estimate)
'''
M2=pymc.MCMC(set([T0,n2,rs,a0,gamma0,k0,n3,rho,ne0,sbp_c,temp_ne2,nth_a,nth_b,nth_gamma,kmod_a,kmod_b,kmod_c,kmod_k0,delta,delta2]))
aa=500
bb=200
cc=3
M2.sample(iter=aa,burn=bb,thin=cc)
median=numpy.int((aa-bb)/cc/2)

M2.write_csv("result.csv")
#k1_f=M2.trace('k1')[:]
#a1_f=M2.trace('a1')[:]
#gamma1_f=M2.trace('gamma1')[:]
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
#a_f=M2.trace('ae')[:]
#nt_f=M2.trace('nt')[:]
#nm_f=M2.trace('nm')[:]
nth_a_f=M2.trace('nth_a')[:]
nth_b_f=M2.trace('nth_b')[:]
nth_gamma_f=M2.trace('nth_gamma')[:]
kmod_a_f=M2.trace('kmod_a')[:]
kmod_b_f=M2.trace('kmod_b')[:]
kmod_c_f=M2.trace('kmod_c')[:]
kmod_k0_f=M2.trace('kmod_k0')[:]
ne_array=[]
T_array=[]
##############################
#print result
rne_array=numpy.array(range(1,2001,5))
tmp_array=numpy.array(range(0,2005,5))

#tc0_f=(T0_f+(1-n3_f)*N1_0*rho_f*rs_f*rs_f*rs_f*numpy.log((rs_f+1)/rs_f)-\
#        (a0_f+k0_f)*n2_f*numpy.power(ne0_f,2/3))/(1-n3_f-n2_f)
SUM_T_array=[]
SUM_ne_array=[]
SUM_ne_odearray=[]
SUM_sbp_odearray=[]
SUM_sbp_fit=[]
SUM_k_array=[]
SUM_kfit_array=[]
SUM_nfw_fit_array=[]
SUM_mass_array=[]
for i in range(len(ne0_f)):
#    y0=ne0_f[i]
    print(i)
    p=[n1_f,n2_f[i],rs_f[i],a0_f[i],gamma0_f[i],k0_f[i],n3_f[i],\
            rho_f[i],\
            nth_a_f[i],nth_b_f[i],nth_gamma_f[i],delta_f[i],delta2_f[i]]
#    print(p)
#    print(y0)
    K_ARRAY_FIT=entropy_model(rne_array,kmod_a_f[i],kmod_b_f[i],kmod_c_f[i],kmod_k0_f[i])
    p_nth=[nth_a_f[i],nth_b_f[i],nth_gamma_f[i]]
#    m_factor=1
    m_factor=[]
    p_MMnfw=[rho_f[i],rs_f[i],delta_f[i],delta2_f[i]]
    for j in range(len(rne_array)):
        r=rne_array[j]
        e_ori=quad(nfw,0,r,p_MMnfw)[0]/r+quad(nfw_divbyr,r,numpy.inf,p_MMnfw)[0]
        e_mod=quad(mod_nfw,0,r,p_MMnfw)[0]/r+quad(mod_nfw_divbyr,r,numpy.inf,p_MMnfw)[0]
        m_tmp=e_mod/e_ori
        m_factor.append(m_tmp)
    m_factor=numpy.array(m_factor)
    T_array=calc_T(rne_array,K_ARRAY_FIT,m_factor,p)
    SUM_T_array.append(T_array)
    p1=[*p,T_array]
    y0=np.power(K_ARRAY_FIT[0]/T_array[0],-1.5)
    y1=odeint(dndr,y0,rne_array,args=(p1,))
    ne_array=[]
    ne_odearray=[]
    for ii in range(0,len(y1)):
        ne_odearray.append(y1[ii][0])
    ne_odearray=numpy.array(ne_odearray)
    SUM_ne_odearray.append(ne_odearray)
    ne_array=numpy.power(K_ARRAY_FIT/T_array,-1.5)
    SUM_ne_array.append(ne_array)
    k_array=T_array*numpy.power(ne_odearray,-2/3)
    SUM_k_array.append(k_array)
    mass_array=calc_mass(rne_array,T_array,ne_array,p_nth)
    SUM_mass_array.append(mass_array)
    nfw_fitted_array=[]
    for j in range(len(rne_array)):
        nfw_fitted=Modefied_Mnfw(rne_array[j],[rho_f[i],rs_f[i],delta_f[i],delta2_f[i]])
        nfw_fitted_array.append(nfw_fitted)
    SUM_nfw_fit_array.append(nfw_fitted_array)
#NE_ARRAY_FIT=ne_array
    ne_array=numpy.insert(ne_array,0,0.0)
    ne_odearray=numpy.insert(ne_odearray,0,0.0)
    sbp_fit=[]
    sbp_odefit=[]
    tmp_r_use=range(1,2001,5)
    for j in range(1,2001,5):
        a=deproject_model.calc_sb(j,tmp_array,ne_array,cfunc_use_array,cm_per_pixel)
        a=a+sbp_c_f[i]
        b=deproject_model.calc_sb(j,tmp_array,ne_odearray,cfunc_use_array,cm_per_pixel)
        b=b+sbp_c_f[i]
        sbp_fit.append(a)
        sbp_odefit.append(b)
    SUM_sbp_fit.append(sbp_fit)
    SUM_sbp_odearray.append(sbp_odefit)
    k_array_fitted=entropy_model(rne_array,kmod_a_f[i],kmod_b_f[i],kmod_c_f[i],kmod_k0_f[i])
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
DEN_10_ARRAY=[]
DEN_16_ARRAY=[]
DEN_50_ARRAY=[]
DEN_84_ARRAY=[]
DEN_90_ARRAY=[]
SBPODE_10_ARRAY=[]
SBPODE_16_ARRAY=[]
SBPODE_50_ARRAY=[]
SBPODE_84_ARRAY=[]
SBPODE_90_ARRAY=[]
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
for i in range(400):
    tmp_t_array=[]
    tmp_k_array=[]
    tmp_den_array=[]
    tmp_sbpode_array=[]
    tmp_sbp_array=[]
    tmp_kfit_array=[]
    tmp_m_array=[]
    tmp_mfit_array=[]
    for j in range(len(ne0_f)):
        tmp_t_array.append(SUM_T_array[j][i])
        tmp_k_array.append(SUM_k_array[j][i])
        tmp_den_array.append(SUM_ne_array[j][i])
        tmp_sbpode_array.append(SUM_sbp_odearray[j][i])
        tmp_sbp_array.append(SUM_sbp_fit[j][i])
        tmp_kfit_array.append(SUM_kfit_array[j][i])
        tmp_m_array.append(SUM_mass_array[j][i])
        tmp_mfit_array.append(SUM_nfw_fit_array[j][i])
#    print(tmp_t_array)
#    print(tmp_k_array)
#    print(tmp_den_array)
#    print(tmp_sbp_array)
    tmp_t_array.sort()
    tmp_k_array.sort()
    tmp_den_array.sort()
    tmp_sbpode_array.sort()
    tmp_sbp_array.sort()
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
    SBP_10_ARRAY.append(tmp_sbp_array[IND_10])
    SBP_16_ARRAY.append(tmp_sbp_array[IND_16])
    SBP_50_ARRAY.append(tmp_sbp_array[IND_50])
    SBP_84_ARRAY.append(tmp_sbp_array[IND_84])
    SBP_90_ARRAY.append(tmp_sbp_array[IND_90])
    K_10_ARRAY.append(tmp_k_array[IND_10])
    K_16_ARRAY.append(tmp_k_array[IND_16])
    K_50_ARRAY.append(tmp_k_array[IND_50])
    K_84_ARRAY.append(tmp_k_array[IND_84])
    K_90_ARRAY.append(tmp_k_array[IND_90])
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
    SBPODE_10_ARRAY.append(tmp_sbpode_array[IND_10])
    SBPODE_16_ARRAY.append(tmp_sbpode_array[IND_16])
    SBPODE_50_ARRAY.append(tmp_sbpode_array[IND_50])
    SBPODE_84_ARRAY.append(tmp_sbpode_array[IND_84])
    SBPODE_90_ARRAY.append(tmp_sbpode_array[IND_90])
NAME=sys.argv[1][:-4]
tmp='array_plt.json'
fi=open(tmp,'w')
dat={
        "name": NAME,
        "radius": [r_array,re_array],
        "temperature": [t_array,te_array],
        "sbp": [sbp_array,sbpe_array],
        "radius_model": list(range(1,2001,5)),
        "temperature_model": [T_50_ARRAY,T_16_ARRAY,T_84_ARRAY,T_10_ARRAY,T_90_ARRAY],
        "sbp_model": [SBP_50_ARRAY,SBP_16_ARRAY,SBP_84_ARRAY,SBP_10_ARRAY,SBP_90_ARRAY],
        "sbp_ode": [SBPODE_50_ARRAY,SBPODE_16_ARRAY,SBPODE_84_ARRAY],
        "k_model": [K_50_ARRAY,K_16_ARRAY,K_84_ARRAY,K_10_ARRAY,K_90_ARRAY],
        "den_model": [DEN_50_ARRAY,DEN_16_ARRAY,DEN_84_ARRAY,DEN_10_ARRAY,DEN_90_ARRAY],
        "m_calc": [M_50_ARRAY,M_16_ARRAY,M_84_ARRAY,M_10_ARRAY,M_90_ARRAY],
        "m_fit": [MFIT_50_ARRAY,MFIT_16_ARRAY,MFIT_84_ARRAY,M_10_ARRAY,M_90_ARRAY],
        "k_fit": [KFIT_50_ARRAY,KFIT_16_ARRAY,KFIT_84_ARRAY,KFIT_10_ARRAY,KFIT_90_ARRAY],
        "r200": R200_0,
        }
json.dump(dat,fi,indent=2)
fi.close()
K_50_ARRAY=numpy.array(K_50_ARRAY)
#fi=open("entropy_fitted.txt",'w')
#print(K_50_ARRAY,file=fi)
#fi.close()
sbpe_array=numpy.array(sbpe_array)
sbp_array=numpy.array(sbp_array)
t_array=numpy.array(t_array)
te_array=numpy.array(te_array)
sbpe_array=sbp_array*0.001+sbpe_array
te_array=t_array*0.05+te_array
#print(sbp_fit)

plt.hold(False)
plt.loglog(tmp_r_use,SBP_84_ARRAY,color='grey')
plt.hold(True)
plt.loglog(tmp_r_use,SBP_50_ARRAY,'b')
plt.loglog(tmp_r_use,SBP_16_ARRAY,color='grey')
plt.fill_between(range(1,2001,5),SBP_16_ARRAY,SBP_84_ARRAY,color='grey')
plt.errorbar(rsbp_array,sbp_array,xerr=rsbpe_array,yerr=sbpe_array,color='k',linestyle='none')
plt.xlabel(name+'_Radius(kpc)')
plt.ylabel('SBP(cm^-2 pixel^-2 s^-1)')
plt.savefig('sbp.pdf',dpi=100)

plt.hold(False)
plt.loglog(tmp_r_use,SBPODE_84_ARRAY,color='grey')
plt.hold(True)
plt.loglog(tmp_r_use,SBPODE_50_ARRAY,'b')
plt.loglog(tmp_r_use,SBPODE_16_ARRAY,color='grey')
plt.fill_between(range(1,2001,5),SBPODE_16_ARRAY,SBPODE_84_ARRAY,color='grey')
plt.errorbar(rsbp_array,sbp_array,xerr=rsbpe_array,yerr=sbpe_array,color='k',linestyle='none')
plt.xlabel(name+'_Radius(kpc)')
plt.ylabel('SBP(cm^-2 pixel^-2 s^-1)')
plt.savefig('sbpode.pdf',dpi=100)

plt.hold(False)
plt.plot(rne_array,T_84_ARRAY,color='grey')
plt.hold(True)
plt.plot(rne_array,T_50_ARRAY,'b')
plt.plot(rne_array,T_16_ARRAY,color='grey')
plt.fill_between(range(1,2001,5),T_16_ARRAY,T_84_ARRAY,color='grey')
plt.errorbar(r_array,t_array,xerr=re_array,yerr=te_array,color='k',linestyle='none')
plt.xlabel(name+'_Radius(kpc)')
plt.ylabel('Temperature(keV)')
plt.savefig('temperature.pdf',dpi=100)


plt.hold(False)
plt.loglog(rne_array,K_84_ARRAY,color='grey')
plt.hold(True)
plt.loglog(rne_array,K_50_ARRAY,'b')
plt.loglog(rne_array,K_16_ARRAY,color='grey')
plt.fill_between(range(1,2001,5),K_16_ARRAY,K_84_ARRAY,color='grey')
plt.loglog(rne_array,KFIT_50_ARRAY,'g')
plt.fill_between(range(1,2001,5),KFIT_16_ARRAY,KFIT_84_ARRAY,color='green')
plt.xlabel(name+'_Radius(kpc)')
plt.ylabel('Entropy(kev cm^2)')
plt.savefig('entropy.pdf',dpi=100)

plt.hold(False)
plt.loglog(rne_array,M_50_ARRAY,'b')
plt.hold(True)
plt.fill_between(range(1,2001,5),M_16_ARRAY,M_84_ARRAY,color='grey')
plt.fill_between(range(1,2001,5),MFIT_16_ARRAY,MFIT_84_ARRAY,color='green')
plt.xlim(10,3000)
plt.ylim(1e11,1e16)
plt.xlabel(name+'_Radius(kpc)')
plt.ylabel('Mass(Msun))')
plt.savefig('mass.pdf',dpi=100)
'''

###############
