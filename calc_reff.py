#!/usr/bin/env python3
#$1 name of temperature data file ($f.txt)
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
from astropy.cosmology import FlatLambdaCDM
cosmo=FlatLambdaCDM(H0=71,Om0=0.27,Tcmb0=2.725)

G=6.67e-11 #m^3 kg^-1 s^-2
mp=1.67e-27  #kg
kev=1.6e-16 #J
kpc=3.086e19 #m
Msun=2e30 #kg
mu=0.61
N1_0=1.40733e-10
pi=3.1415926
T0_0=0
script_dir=''
tmp1=list(range(1,11))
tmp2=list(range(11,41,3))
tmp3=list(range(41,3001,5))
tmp1.extend(tmp2)
tmp1.extend(tmp3)
rne_array=numpy.array(tmp1)
tmp=sys.argv[0].split('/')
for i in range(len(tmp)-1):
    script_dir=script_dir+tmp[i]
    script_dir=script_dir+'/'
ta1=numpy.load(script_dir+'lrs_ori.npy')
ta2=numpy.load(script_dir+'lrs_dvr.npy')
ta3=numpy.load(script_dir+'hrs_ori.npy')
ta4=numpy.load(script_dir+'hrs_dvr.npy')
t_total=[ta1,ta2,ta3,ta4]

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
##p=[rho,rs,delta,delta2]
def mod_nfw(r,p):
    rho=p[0]
    rs=p[1]
    delta=p[2]
    delta2=p[3]
    den=rho/(np.power(r/rs,delta)*(numpy.power(1+r/rs,delta2-delta)))*4*pi*r*r
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
def calc_T(x,k_fit,m_factor,p):
    T0_0=0
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
    nth_a=p[8]
    nth_b=p[9]
    nth_gamma=p[10]
    delta=p[11]
    delta2=p[12]
    c4=p[13]
    p_MMnfw=[rho,rs,delta,delta2]
    eta=nth_a*(1+numpy.exp(-numpy.power(r0/R200M/nth_b,nth_gamma)))
    Tx=(T0_0+(1-c3*c4)*c1*rho*rs*rs*rs*numpy.log((rs+x)/rs)/x*m_factor)/(1/eta-c3-c2+c2*k/k_fit)
    return Tx

def entropy_model(x,a,b,c,k0):
    x=numpy.array(x)
    return a*numpy.power(x,b)*numpy.exp(c*(-x/R200_0))+k0
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
        elif i==len(r_array)-1:
            dTdr=(t_array[i]-t_array[i-1])/(r_array[i]-r_array[i-1])
            dnedr=(ne_array[i]-ne_array[i-1])/(r_array[i]-r_array[i-1])
            z1=(dnedr*t_array[i]/ne_array[i]+dTdr-t_array[i]/eta*detadr)/eta #kev/kpc
            mass=-z1/G/mu/mp*r0*r0/Msun*kpc*kev #Msun
        else:
            dTdr_r=(t_array[i+1]-t_array[i])/(r_array[i+1]-r_array[i])
            dnedr_r=(ne_array[i+1]-ne_array[i])/(r_array[i+1]-r_array[i])
            dTdr_l=(t_array[i]-t_array[i-1])/(r_array[i]-r_array[i-1])
            dnedr_l=(ne_array[i]-ne_array[i-1])/(r_array[i]-r_array[i-1])
            dTdr=(dTdr_l+dTdr_r)/2
            dnedr=(dnedr_l+dnedr_r)/2
            z1=(dnedr*t_array[i]/ne_array[i]+dTdr-t_array[i]/eta*detadr)/eta #kev/kpc
            mass=-z1/G/mu/mp*r0*r0/Msun*kpc*kev #Msun
        if mass<=1e7:
            mass=1e7
        mass_array.append(mass)
    return mass_array
def clumping_model(x,n1,n2,n3,n4,n5):
    return numpy.power(1+x,n1)*numpy.exp(x*n2)+n3*numpy.exp(-(x-n4)*(x-n4)/n5)

for i in open('param_zzh_for_py.txt'):
    if re.match(r'^R200\s',i):
        R200_0=float(i.split()[1])
        R200M=R200_0
        r200=R200_0
    if re.match(r'z\s',i):
        z=float(i.split()[1])
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

name=sys.argv[1]
M2=pymc.database.pickle.load('sampled.pickle')
aa=50000
bb=30000
cc=50
for i in open('global.cfg'):
    if re.match(r'^sbp_data_file',i):
        sbp_data_file=i.split()[1]
    if re.match(r'^sbp_eninfo',i):
        sbp_type=i.split()[3]
    if re.match(r'^temp_data_file',i):
        temp_data_file=i.split()[1]
if sbp_type=='CNT':
    cfunc_file='cfunc_for_density_fit_cnt.txt'
elif sbp_type=='ERG':
    cfunc_file='cfunc_for_density_fit_erg.txt'
cm_per_pixel=cosmo.kpc_proper_per_arcmin(z).value/60*0.492*kpc*100
for i in open(temp_data_file):
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
for i in open(sbp_data_file):
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
for i in open(cfunc_file):
    r_cfunc_array.append(float(i.split()[0]))
    cfunc_ori_array.append(float(i.split()[1]))
for i in open('cfunc_for_chandra_density_fit_cnt.txt'):
    r_cfunc_c_array.append(float(i.split()[0]))
    cfunc_cori_array.append(float(i.split()[1]))
for i in rne_array:
    if i==0:
        continue
    for j in range(len(r_cfunc_array)):
        if r_cfunc_array[j]>i:
            cfunc_use_array.append(cfunc_ori_array[j])
            cfunc_cuse_array.append(cfunc_cori_array[j])
            break
cfunc_use_array=numpy.array(cfunc_use_array)
cfunc_cuse_array=np.array(cfunc_cuse_array)
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
ne_array=[]
T_array=[]
p_all_f=numpy.array([n2_f,rs_f,a0_f,gamma0_f,k0_f,n3_f,rho_f,ne0_f,T0_f,sbp_c_f,delta_f,delta2_f,nth_a_f,nth_b_f,nth_gamma_f,kmod_a_f,kmod_b_f,kmod_c_f,kmod_k0_f,cp_p_f,cp_e_f,cp_g0_f,cp_x0_f,cp_sigma_f,c4_f])
numpy.save('p_all',p_all_f)
##############################

SUM_T_array=[]
SUM_ne_array=[]
SUM_ne_odearray=[]
SUM_sbp_odearray=[]
SUM_sbp_fit=[]
SUM_k_array=[]
SUM_kfit_array=[]
SUM_nfw_fit_array=[]
SUM_mass_array=[]
SUM_ne_cl_array=[]
p=[n1_f,n2_f.mean(),rs_f.mean(),a0_f.mean(),gamma0_f.mean(),k0_f.mean(),n3_f.mean(),rho_f.mean(),nth_a_f.mean(),nth_b_f.mean(),nth_gamma_f.mean(),delta_f.mean(),delta2_f.mean(),c4_f.mean()]
K_ARRAY_FIT=entropy_model(rne_array,kmod_a_f.mean(),kmod_b_f.mean(),kmod_c_f.mean(),kmod_k0_f.mean())
p_nth=[nth_a_f.mean(),nth_b_f.mean(),nth_gamma_f.mean()]
m_factor=[]
p_MMnfw=[rho_f.mean(),rs_f.mean(),delta_f.mean(),delta2_f.mean()]
for j in range(len(rne_array)):
    r=rne_array[j]
    e_ori=4*pi*rho_f.mean()*np.power(rs_f.mean(),3)/r*np.log((rs_f.mean()+r)/rs_f.mean())
    e_mod=modnfw_readarray.calc('n',r,p_MMnfw,t_total)/r+modnfw_readarray.calc('r',r,p_MMnfw,t_total)
    m_tmp=e_mod/e_ori
    m_factor.append(m_tmp)
m_factor=numpy.array(m_factor)
T_array=calc_T(rne_array,K_ARRAY_FIT,m_factor,p)
#    SUM_T_array.append(T_array)
p1=[*p,T_array]
ne_array=[]
ne_array=numpy.power(K_ARRAY_FIT/T_array,-1.5)
#    SUM_ne_array.append(ne_array)
mass_array=calc_mass(rne_array,T_array,ne_array,p_nth)
#    SUM_mass_array.append(mass_array)
nfw_fitted_array=[]
for j in range(len(rne_array)):
    nfw_fitted=Modefied_Mnfw(rne_array[j],[rho_f.mean(),rs_f.mean(),delta_f.mean(),delta2_f.mean()])
    nfw_fitted_array.append(nfw_fitted)
#    SUM_nfw_fit_array.append(nfw_fitted_array)
ne_cl_array=clumping_model(rne_array/R200_0,cp_p_f.mean(),cp_e_f.mean(),cp_g0_f.mean(),cp_x0_f.mean(),cp_sigma_f.mean())*ne_array
for j in range(len(ne_cl_array)):
    if ne_cl_array[j]<ne_array[j]:
        ne_cl_array[j]=ne_array[j]
ne_cl_array=numpy.insert(ne_cl_array,0,0.0)
#    SUM_ne_cl_array.append(ne_cl_array)
ne_array=numpy.insert(ne_array,0,0.0)
sbp_fit=[]
tmp_r_use=rne_array
sbp_model_array=[]
for j in range(len(rsbp_array)):
    a=deproject_model.calc_sb(rsbp_array[j],tmp_array,ne_cl_array,cfunc_use_array,cm_per_pixel)
    a=a+sbp_c_f.mean()
    sbp_model_array.append(a)
sbp_model_array=np.array(sbp_model_array)
t_model_array=[]
for i in range(len(r_array)):
    for j in range(len(rne_array)):
        if rne_array[j]>r_array[i]:
            r_this=j
            break
    t_this=T_array[r_this]
    t_model_array.append(t_this)
t_model=np.array(t_model_array)
rtmp=numpy.array(list(range(45)))
rmass_array=30*numpy.power(1.11,rtmp)
for j in range(len(rmass_array)):
    rmass_array[j]=numpy.int(rmass_array[j])
rmass_array=numpy.insert(rmass_array,0,10)
rmass_array=numpy.insert(rmass_array,1,20)
mass_model=[]
for j in range(len(rmass_array)):
    Mnfw_model=numpy.log(Modefied_Mnfw(rmass_array[j],[rho_f.mean(),rs_f.mean(),delta_f.mean(),delta2_f.mean()]))
    mass_model.append(Mnfw_model)
mass_model=np.array(mass_model)
nbin=0
reffi_1=np.zeros(len(rsbp_array)+len(r_array)+len(rmass_array))
reffi_2=np.zeros(len(rsbp_array)+len(r_array)+len(rmass_array))
for i in range(len(ne0_f)):
#    print(i)
    p=[n1_f,n2_f[i],rs_f[i],a0_f[i],gamma0_f[i],k0_f[i],n3_f[i],\
            rho_f[i],\
            nth_a_f[i],nth_b_f[i],nth_gamma_f[i],delta_f[i],delta2_f[i],c4_f[i]]
    K_ARRAY_FIT=entropy_model(rne_array,kmod_a_f[i],kmod_b_f[i],kmod_c_f[i],kmod_k0_f[i])
    p_nth=[nth_a_f[i],nth_b_f[i],nth_gamma_f[i]]
    m_factor=[]
    p_MMnfw=[rho_f[i],rs_f[i],delta_f[i],delta2_f[i]]
    for j in range(len(rne_array)):
        r=rne_array[j]
        e_ori=4*pi*rho_f[i]*np.power(rs_f[i],3)/r*np.log((rs_f[i]+r)/rs_f[i])
        e_mod=modnfw_readarray.calc('n',r,p_MMnfw,t_total)/r+modnfw_readarray.calc('r',r,p_MMnfw,t_total)
        m_tmp=e_mod/e_ori
        m_factor.append(m_tmp)
    m_factor=numpy.array(m_factor)
    T_array=calc_T(rne_array,K_ARRAY_FIT,m_factor,p)
#    SUM_T_array.append(T_array)
    p1=[*p,T_array]
    ne_array=[]
    ne_array=numpy.power(K_ARRAY_FIT/T_array,-1.5)
#    SUM_ne_array.append(ne_array)
    mass_array=calc_mass(rne_array,T_array,ne_array,p_nth)
#    SUM_mass_array.append(mass_array)
    nfw_fitted_array=[]
    for j in range(len(rne_array)):
        nfw_fitted=Modefied_Mnfw(rne_array[j],[rho_f[i],rs_f[i],delta_f[i],delta2_f[i]])
        nfw_fitted_array.append(nfw_fitted)
#    SUM_nfw_fit_array.append(nfw_fitted_array)
    ne_cl_array=clumping_model(rne_array/R200_0,cp_p_f[i],cp_e_f[i],cp_g0_f[i],cp_x0_f[i],cp_sigma_f[i])*ne_array
    for j in range(len(ne_cl_array)):
        if ne_cl_array[j]<ne_array[j]:
            ne_cl_array[j]=ne_array[j]
    ne_cl_array=numpy.insert(ne_cl_array,0,0.0)
#    SUM_ne_cl_array.append(ne_cl_array)
    ne_array=numpy.insert(ne_array,0,0.0)
    sbp_fit=[]
    tmp_r_use=rne_array
    for j in range(len(rsbp_array)):
        a=deproject_model.calc_sb(rsbp_array[j],tmp_array,ne_cl_array,cfunc_use_array,cm_per_pixel)
        a=a+sbp_c_f[i]
        p_this=j
        reffi_1[p_this]=reffi_1[p_this]+np.power(a-sbp_model_array[j],2)
        reffi_2[p_this]=reffi_2[p_this]+np.power(a-sbp_array[j],2)
    for j in range(len(r_array)):
        for k in range(len(rne_array)):
            if rne_array[k]>r_array[j]:
                r_this=k
                break
        t_this=T_array[r_this]
        p_this=len(rsbp_array)+j
        reffi_1[p_this]=reffi_1[p_this]+np.power(t_this-t_model[j],2)
        reffi_2[p_this]=reffi_2[p_this]+np.power(t_this-t_array[j],2)
    rtmp=numpy.array(list(range(45)))
    rmass_array=30*numpy.power(1.11,rtmp)
    for j in range(len(rmass_array)):
        rmass_array[j]=numpy.int(rmass_array[j])
    rmass_array=numpy.insert(rmass_array,0,10)
    rmass_array=numpy.insert(rmass_array,1,20)
    for j in range(len(rmass_array)):
        if rmass_array[j]>=10:
            for tmp in range(len(rne_array)):
                if rne_array[tmp]>=rmass_array[j]:
                    r_this=tmp
                    break
            Mnfw_model=numpy.log(Modefied_Mnfw(rmass_array[j],[rho_f[i],rs_f[i],delta_f[i],delta2_f[i]]))
            M_this=numpy.log(mass_array[r_this])
            p_this=len(rsbp_array)+len(r_array)+j
            reffi_1[p_this]=reffi_1[p_this]+np.power(Mnfw_model-mass_model[j],2)
            reffi_2[p_this]=reffi_2[p_this]+np.power(Mnfw_model-M_this,2)

reffi_1=np.append(np.append(np.array(sbpe_array)+np.array(sbp_array)*0.03,np.array(te_array)+np.array(t_array)*0.05),np.array(mass_model)*0.08)

reffi_1=reffi_1*reffi_1*400
#print(reffi_1,reffi_2)
reffi=1-reffi_2/reffi_1
#print(reffi)
reff=reffi.sum()/(len(t_array)+len(rsbp_array)+len(rmass_array))
print(reff)
#NAME=sys.argv[1][:-4]
'''

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
DEN_CL_10_ARRAY=[]
DEN_CL_16_ARRAY=[]
DEN_CL_50_ARRAY=[]
DEN_CL_84_ARRAY=[]
DEN_CL_90_ARRAY=[]
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
for i in range(len(rne_array)):
    tmp_t_array=[]
    tmp_k_array=[]
    tmp_den_array=[]
    tmp_den_cl_array=[]
    tmp_sbpode_array=[]
    tmp_sbp_array=[]
    tmp_kfit_array=[]
    tmp_m_array=[]
    tmp_mfit_array=[]
    for j in range(len(ne0_f)):
        tmp_t_array.append(SUM_T_array[j][i])
#        tmp_k_array.append(SUM_k_array[j][i])
        tmp_den_array.append(SUM_ne_array[j][i])
        tmp_den_cl_array.append(SUM_ne_cl_array[j][i])
#        tmp_sbpode_array.append(SUM_sbp_odearray[j][i])
        tmp_sbp_array.append(SUM_sbp_fit[j][i])
        tmp_kfit_array.append(SUM_kfit_array[j][i])
        tmp_m_array.append(SUM_mass_array[j][i])
        tmp_mfit_array.append(SUM_nfw_fit_array[j][i])
    tmp_t_array.sort()
#    tmp_k_array.sort()
    tmp_den_array.sort()
    tmp_den_cl_array.sort()
#    tmp_sbpode_array.sort()
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
#    K_10_ARRAY.append(tmp_k_array[IND_10])
#    K_16_ARRAY.append(tmp_k_array[IND_16])
#    K_50_ARRAY.append(tmp_k_array[IND_50])
#    K_84_ARRAY.append(tmp_k_array[IND_84])
#    K_90_ARRAY.append(tmp_k_array[IND_90])
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
#    SBPODE_10_ARRAY.append(tmp_sbpode_array[IND_10])
#    SBPODE_16_ARRAY.append(tmp_sbpode_array[IND_16])
#    SBPODE_50_ARRAY.append(tmp_sbpode_array[IND_50])
#    SBPODE_84_ARRAY.append(tmp_sbpode_array[IND_84])
#    SBPODE_90_ARRAY.append(tmp_sbpode_array[IND_90])
NAME=sys.argv[1][:-4]
'''
