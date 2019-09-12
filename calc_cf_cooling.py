#!/usr/bin/env python3
##############
#USEAGE:
#$1: parameter result file(csv) from pymc
#$2: saved array_data_file from entopry_mcmc.py
#$3: cooling function file used to calculate lx
#$4: luminosity distance in unit cm
#
# PURPOSE:
# 1 calculate mass deposition rate (within cooling radius)
# 2 calcuate "extra energy loss" vs radius
# 3 calcualte lx vs radius distribution
# 4 plot necessary fit results !!!to be added
###################
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
kpc_cm=3.086e21   #cm
#dr=5*kpc
pi=3.1415926
kev_to_erg=1.6e-9
Gyr=1e9*365.24*24*3600 #s
M2=pymc.database.pickle.load('sampled.pickle')

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
dl=np.float(sys.argv[4]) #cm calcuted by Gu script
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
#    dkdr=a0*gamma0*numpy.power(r0,gamma0-1)
    r200=R200_0
    nth_a=p[8]
    nth_b=p[9]
    nth_gamma=p[10]
    delta=p[11]
    delta2=p[12]
    p_MMnfw=[rho,rs,delta,delta2]
    eta=nth_a*(1+numpy.exp(-numpy.power(r0/R200M/nth_b,nth_gamma)))
    Tx=(T0_0+(1-c3)*c1*rho*rs*rs*rs*numpy.log((rs+x)/rs)/x*m_factor)/(1/eta-c3-c2+c2*k/k_fit)
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
        eta=0.6+0.4*a*(1+numpy.exp(-numpy.power(r0/R200M/b,c)))
        detadr=a*numpy.exp(-numpy.power(r0/R200M/b,c))*\
                c*-numpy.power(r0/R200M/b,c-1)/R200M/b
        detadr=0.4*detadr
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

# read cooling function file
cfunc_file=sys.argv[3]
r_cfunc_array=[]
cfunc_ori_array=[]
cfunc_lx_use_array=[]
for i in open(cfunc_file,'r'):
    r_cfunc_array.append(float(i.split()[0]))
    cfunc_ori_array.append(float(i.split()[1]))

for i in rne_array:
    i=float(i)
    if i==0:
        continue
    for j in range(len(r_cfunc_array)):
        if r_cfunc_array[j]>i:
            cfunc_lx_use_array.append(cfunc_ori_array[j])
            break
cfunc_lx_use_array=np.array(cfunc_lx_use_array)
for i in open('param_zzh_for_py.txt'):
    if re.match(r'^R200\s',i):
        R200_0=float(i.split()[1])
        R200M=R200_0
        r200=R200_0
r_array=[]
re_array=[]
t_array=[]
te_array=[]
rsbp_array=[]
rsbpe_array=[]
sbp_array=[]
sbpe_array=[]
M2=pymc.database.pickle.load('sampled.pickle')
aa=30000
bb=10000
cc=50
for i in open('global.cfg'):
    if re.match(r'^sbp_cfg',i):
        sbp_cfg=i.split()[1]
    if re.match(r'^radius_sbp_file',i):
        sbp_data_file=i.split()[1]
for i in open(sbp_cfg):
    if re.match(r'^cm_per_pixel',i):
        cm_per_pixel=float(i.split()[1])
for i in open(sbp_data_file):
    r,rer,sbp,sbpe=i.split()
    r=float(r)
    rer=float(rer)
    sbp=float(sbp)
    sbpe=float(sbpe)
    rsbp_array.append(r)
    rsbpe_array.append(rer)
    sbp_array.append(sbp)
    sbpe_array.append(sbpe)
#M2.write_csv(name+"result.csv")
cfunc_ori_array=[]
r_cfunc_array=[]
cfunc_use_array=[]
cfunc_ori_array.append(0)
r_cfunc_array.append(0)
cfunc_use_array.append(0)
for i in open('cfunc_for_density_fit.txt'):
    r_cfunc_array.append(float(i.split()[0]))
    cfunc_ori_array.append(float(i.split()[1]))
for i in rne_array:
    if i==0:
        continue
    for j in range(len(r_cfunc_array)):
        if r_cfunc_array[j]>i:
            cfunc_use_array.append(cfunc_ori_array[j])
            break
cfunc_use_array=numpy.array(cfunc_use_array)
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
ne_array=[]
T_array=[]
p_all_f=numpy.array([n2_f,rs_f,a0_f,gamma0_f,k0_f,n3_f,rho_f,ne0_f,T0_f,sbp_c_f,delta_f,delta2_f,nth_a_f,nth_b_f,nth_gamma_f,kmod_a_f,kmod_b_f,kmod_c_f,kmod_k0_f,cp_p_f,cp_e_f,cp_g0_f,cp_x0_f,cp_sigma_f])
numpy.save('p_all',p_all_f)
##############################
#print result

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
SUM_ne_cl_array=[]
SUM_lx_array=[]
SUM_E_loss_this_all=[]
SUM_lx_sum=[]
for i in range(len(ne0_f)):
#    y0=ne0_f[i]
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
#        e_ori=quad(nfw,0,r,p_MMnfw)[0]/r+quad(nfw_divbyr,r,numpy.inf,p_MMnfw)[0]
        e_ori=4*pi*rho_f[i]*np.power(rs_f[i],3)/r*np.log((rs_f[i]+r)/rs_f[i])
        e_mod=modnfw_readarray.calc('n',r,p_MMnfw,t_total)/r+modnfw_readarray.calc('r',r,p_MMnfw,t_total)
#        e_mod=quad(mod_nfw,0,r,p_MMnfw)[0]/r+quad(mod_nfw_divbyr,r,numpy.inf,p_MMnfw)[0]
        m_tmp=e_mod/e_ori
        m_factor.append(m_tmp)
    m_factor=numpy.array(m_factor)
    T_array=calc_T(rne_array,K_ARRAY_FIT,m_factor,p)
    SUM_T_array.append(T_array)
    p1=[*p,T_array]
    y0=np.power(K_ARRAY_FIT[0]/T_array[0],-1.5)
#    y1=odeint(dndr,y0,rne_array,args=(p1,))
    ne_array=[]
    ne_odearray=[]
#    for ii in range(0,len(y1)):
#        ne_odearray.append(y1[ii][0])
#    ne_odearray=numpy.array(ne_odearray)
#    SUM_ne_odearray.append(ne_odearray)
    ne_array=numpy.power(K_ARRAY_FIT/T_array,-1.5)
    SUM_ne_array.append(ne_array)
#    k_array=T_array*numpy.power(ne_odearray,-2/3)
#    SUM_k_array.append(k_array)
    mass_array=calc_mass(rne_array,T_array,ne_array,p_nth)
    SUM_mass_array.append(mass_array)
    nfw_fitted_array=[]
    for j in range(len(rne_array)):
        nfw_fitted=Modefied_Mnfw(rne_array[j],[rho_f[i],rs_f[i],delta_f[i],delta2_f[i]])
        nfw_fitted_array.append(nfw_fitted)
    SUM_nfw_fit_array.append(nfw_fitted_array)
#NE_ARRAY_FIT=ne_array
    ne_cl_array=clumping_model(rne_array/R200_0,cp_p_f[i],cp_e_f[i],cp_g0_f[i],cp_x0_f[i],cp_sigma_f[i])*ne_array
    for j in range(len(ne_cl_array)):
        if ne_cl_array[j]<ne_array[j]:
            ne_cl_array[j]=ne_array[j]
    ne_cl_array=numpy.insert(ne_cl_array,0,0.0)
    SUM_ne_cl_array.append(ne_cl_array)
    ne_array=numpy.insert(ne_array,0,0.0)
#    ne_odearray=numpy.insert(ne_odearray,0,0.0)
    sbp_fit=[]
    sbp_odefit=[]
    tmp_r_use=rne_array
#    for j in rne_array:
#        a=deproject_model.calc_sb(j,tmp_array,ne_cl_array,cfunc_use_array,cm_per_pixel)
#        a=a+sbp_c_f[i]
#        b=deproject_model.calc_sb(j,tmp_array,ne_odearray,cfunc_use_array,cm_per_pixel)
#        b=b+sbp_c_f[i]
#        sbp_fit.append(a)
#        sbp_odefit.append(b)
#    SUM_sbp_fit.append(sbp_fit)
#    SUM_sbp_odearray.append(sbp_odefit)
    k_array_fitted=entropy_model(rne_array,kmod_a_f[i],kmod_b_f[i],kmod_c_f[i],kmod_k0_f[i])
    SUM_kfit_array.append(k_array_fitted)
    lx_array=[]
    E_loss_this_all=[]
    lx_1=0
    lx_2=0
    lx_3=0
    lx_4=0
    lx_5=0
    eloss_1=0
    eloss_2=0
    eloss_3=0
    eloss_4=0
    eloss_5=0
    for ii in range(len(rne_array)):
        r=rne_array[ii]
        ne_this=ne_array[ii] #cm^-3
        ne_cl_this=ne_cl_array[ii] #cm^-3
        t_this=T_array[ii]   #kev
        if ii==0:
            dr=rne_array[1]-rne_array[0]
        elif ii==len(rne_array)-1:
            dr=rne_array[ii]-rne_array[ii-1]
        else:
            dr=(rne_array[ii+1]-rne_array[ii-1])/2
        dr=dr*kpc_cm #cm
        lx_this=ne_cl_this*ne_cl_this/1.2*4*pi*r*r*dr*cfunc_lx_use_array[ii]*kpc_cm*kpc_cm*4*pi*dl*dl
        k_this=k_array_fitted[ii]
        k_ori_this=k0_f[i]+a0_f[i]*np.power(r,gamma0_f[i])
        E_loss_factor=1.5*n2_f[i]*t_this*\
            4*pi*r*r*kpc_cm*kpc_cm*dr*ne_this*(2.2-1.2)*kev_to_erg/10/Gyr #erg/s
        E_loss_sub=E_loss_factor*k_ori_this/k_this
#        E_loss_this=E_loss_factor*(1-E_loss_sub)
        E_loss_this=np.array([E_loss_factor-E_loss_sub,0])
        if r<r200:
            if r<0.1*r200:
                lx_1=lx_1+lx_this
                eloss_1=eloss_1+E_loss_this
            elif r<0.3*r200:
                lx_2=lx_2+lx_this
                eloss_2=eloss_2+E_loss_this
            elif r<0.5*r200:
                lx_3=lx_3+lx_this
                eloss_3=eloss_3+E_loss_this
            elif r<0.8*r200:
                lx_4=lx_4+lx_this
                eloss_4=eloss_4+E_loss_this
            else:
                lx_5=lx_5+lx_this
                eloss_5=eloss_5+E_loss_this
    E_loss_this_all=[eloss_1,eloss_2,eloss_3,eloss_4,eloss_5]
    lx_array=[lx_1,lx_2,lx_3,lx_4,lx_5]
    lx_sum=sum(lx_array)
    SUM_E_loss_this_all.append(E_loss_this_all)
    SUM_lx_array.append(lx_array)
    SUM_lx_sum.append(lx_sum)

lx_array_down=[]
lx_array_up=[]
E_loss_array_1_up=[]
E_loss_array_1_down=[]
E_loss_array_2_up=[]
E_loss_array_2_down=[]
IND_10=numpy.int(numpy.round((aa-bb)/cc*0.1))
IND_16=numpy.int(numpy.round((aa-bb)/cc*0.16))
IND_50=numpy.int(numpy.round((aa-bb)/cc*0.5))
IND_84=numpy.int(numpy.round((aa-bb)/cc*0.84))
IND_90=numpy.int(numpy.round((aa-bb)/cc*0.9))
for i in range(len(lx_array)):
    tmp_lx_array=[]
    tmp_eloss_array_1=[]
    tmp_eloss_array_2=[]
    tmp_other_array=[]
    for j in range(len(ne0_f)):
        tmp_lx_array.append(SUM_lx_array[j][i])
        tmp_eloss_array_1.append(SUM_E_loss_this_all[j][i][0])
        tmp_eloss_array_2.append(SUM_E_loss_this_all[j][i][1])
#        tmp_other_array.append(SUM_lx_array[j][i]+SUM_E_loss_this_all[j][i])
    tmp_lx_array.sort()
    tmp_eloss_array_1.sort()
    tmp_eloss_array_2.sort()
    lx_array_down.append(tmp_lx_array[IND_10])
    lx_array_up.append(tmp_lx_array[IND_90])
    E_loss_array_1_down.append(tmp_eloss_array_1[IND_10])
    E_loss_array_1_up.append(tmp_eloss_array_1[IND_90])
    E_loss_array_2_down.append(tmp_eloss_array_2[IND_10])
    E_loss_array_2_up.append(tmp_eloss_array_2[IND_90])

E_loss_array_down=np.array(E_loss_array_1_down)-np.array(E_loss_array_2_up)
E_loss_array_up=np.array(E_loss_array_1_up)-np.array(E_loss_array_2_down)
SUM_lx_sum.sort()
l200_up=SUM_lx_sum[IND_90]
l200=SUM_lx_sum[IND_50]
l200_down=SUM_lx_sum[IND_10]
print('l200:',l200,l200_down,l200_up)
print('lx_array:',lx_array_down[:],lx_array_up[:])
print('eloss_array:',E_loss_array_down,E_loss_array_up)
#print('other_energy:',other_array_down[:],other_array_up[:])


#read fitted paramter value
#param_file=sys.argv[1]
#result_file=sys.argv[2]



# !!! HOW TO CONSIDER THE EFFECT OF E(Z) HERE ?

# !!! EXACT PRINCIPLE TO CALCULATE THE MASS DEPOSITION RATE ?
'''
for i in range(len(rne_array)):
    r=rne_array[i]
    if r>R200:
        break
    ne_this=ne_array[i]   #cm^-3
    t_this=t_array[i]   #kev
    if i==0:
        dr=rne_array[1]-rne_array[0]
    elif i==len(rne_array)-1:
        dr=rne_array[i]-rne_array[i-1]
    else:
        dr=(rne_array[i+1]-rne_array[i-1])/2
    dr=dr*kpc_cm #cm
    lx_this=ne_this*ne_this*1.2*4*pi*r*r*dr*cfunc_use_array[i]*kpc_cm*kpc_cm*4*pi*dl*dl
    E_loss_this_all=[]
    for j in range(len(a0_f)):
        lx_this=ne_this*ne_this*1.2*4*pi*r*r*dr*cfunc_use_array[i]*kpc_cm*kpc_cm*4*pi*dl*dl
        k_this=kmod_a_f[j]*np.power(r,kmod_b_f[j])*np.exp(kmod_c_f[j]*(-r/R200)*1)+kmod_k0_f[j]
        k_ori_this=k0_f[j]+a0_f[j]*np.power(r,gamma0_f[j])
        E_loss_this=n2_f[j]*t_this*(k_this-k_ori_this)/k_this*\
            4*pi*r*r*kpc_cm*kpc_cm*dr*ne_this*2.2*kev_to_erg/10/Gyr #erg/s
        E_loss_this_all.append(E_loss_this)
    E_loss_this_all.sort()
    E_loss_up_this=E_loss_this_all[0]
    E_loss_down_this=E_loss_this_all[-1]
    lx_array.append(lx_this)
    E_loss_array_up.append(E_loss_up_this)
    E_loss_array_down.append(E_loss_down_this)
'''
#r500=   # !!!
#r200=   # !!!


