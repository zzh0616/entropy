#!/usr/bin/env python3
from zzh_model import gpob
from numpy.random import random
import zzh_model
import deproject_model
import modnfw_readarray
from modnfw_readarray import lintp
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
tmp1=list(range(1,11))
tmp2=list(range(11,41,3))
tmp3=list(range(41,3001,5))
tmp1.extend(tmp2)
tmp1.extend(tmp3)
rne_array=numpy.array(tmp1)

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
def main(name):
    for i in open('param_zzh_for_py.txt'):
	if re.match(r'^R200\s',i):
	    R200_0=float(i.split()[1])
	if re.match(r'z\s',i):
        z=float(i.split()[1])
	if re.match(r'aa\s',i):
		aa=float(i.split()[1])
		bb=float(i.split()[2])
		cc=float(i.split()[3])
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
	M2=pymc.database.pickle.load('sampled.pickle')
	for i in open('global.cfg'):
    	if re.match(r'^sbp_data_file',i):
        	sbp_data_file=i.split()[1]
    	if re.match(r'^temp_data_file',i):
        	temp_data_file=i.split()[1]
#cm_per_pixel=cosmo.kpc_proper_per_arcmin(z).value/60*0.492*kpc*100
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
	t_array=list(t_array)
	te_array=list(te_array)
	sum_array=np.load('sum_array.npy')

sbp_c_f=M2.trace('sbp_c')[:]
ne_array=[]
T_array=[]
#SUM_T_array=[]
#SUM_ne_array=[]
#SUM_ne_odearray=[]
#SUM_sbp_odearray=[]
#SUM_sbp_fit=[]
#SUM_k_array=[]
#SUM_kfit_array=[]
#SUM_nfw_fit_array=[]
#SUM_mass_array=[]
#SUM_ne_cl_array=[]
sbp_fit=[]
tmp_r_use=rne_array
sbp_model_array=[]
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
