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
def reassign(rout_array,rin_array,yin_array):
    if len(rin_array) != len(yin_array):
        print('warning: please check input of reassign')
        return 0
    yout_array=[]
    for i in range(len(rout_array)):
        if rout_array[i]<=rin_array[0]:
            yout=yin_array[0]
        elif rout_array[i]>=rin_array[-1]:
            yout=yin_array[-1]
        else:
            for j in range(len(rin_array)):
                if j==len(rin_array)-1:
                    print('warning: please check the result of the function reassign')
                    return -1
                if rout_array[i]<rin_array[j]:
                    yout=lintp(rout_array[i],rin_array[j-1],rin_array[j],yin_array[j-1],yin_array[j])
                    break
        yout_array.append(yout)
    yout_array=np.array(yout_array)
    return yout_array

def main(name,sum_array=[]):
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
    for i in open('param_zzh_for_py.txt'):
        if re.match(r'^R200\s',i):
            R200_0=float(i.split()[1])
        if re.match(r'z\s',i):
            z=float(i.split()[1])
        if re.match(r'aa\s',i):
            aa=float(i.split()[1])
            bb=float(i.split()[2])
            cc=float(i.split()[3])
    ind_50=int((aa-bb)/cc*0.50)
    ind_84=int((aa-bb)/cc*0.84)
    ind_16=int((aa-bb)/cc*0.16)
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
    t_array=np.array(t_array)
    te_array=np.array(te_array)
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
    if len(rcsbp_array)>0:
        flag_csbp=True
    else:
        flag_csbp=False
    if sum_array==[]:
        sum_array=np.load('sum_array.npy',allow_pickle=True)
    sum_t_array=np.array(list(sum_array[1]),dtype=float)
    sum_sbp_array=np.array(list(sum_array[6]),dtype=float)
    sum_csbp_array=np.array(list(sum_array[7]),dtype=float)
    sum_tproj_array=np.array(list(sum_array[8]),dtype=float)
    sum_mcalc_array=np.array(list(sum_array[4]),dtype=float)
    sum_mfit_array=np.array(list(sum_array[5]),dtype=float)
    t_array_center=np.sort(sum_t_array,0)[ind_50]
    sbp_array_center=np.sort(sum_sbp_array,0)[ind_50]
    csbp_array_center=np.sort(sum_csbp_array,0)[ind_50]
    tproj_array_center=np.sort(sum_tproj_array,0)[ind_50]
    mcalc_array_center=np.sort(sum_mcalc_array,0)[ind_50]
    mfit_array_center=np.sort(sum_mfit_array,0)[ind_50]
    t_model=reassign(r_array,rne_array,t_array_center)
    tproj_model_array=reassign(r_array,rne_array,tproj_array_center)
    sbp_model_array=reassign(rsbp_array,rne_array,sbp_array_center)
    if flag_csbp:
        csbp_model_array=reassign(rcsbp_array,rne_array,csbp_array_center)
    rtmp=numpy.array(list(range(45)))
    rmass_array=30*numpy.power(1.11,rtmp)
    for j in range(len(rmass_array)):
        rmass_array[j]=numpy.int(rmass_array[j])
    rmass_array=numpy.insert(rmass_array,0,10)
    rmass_array=numpy.insert(rmass_array,1,20)
    mass_model=reassign(rmass_array,rne_array,mfit_array_center)
    nbin=0
    reffi_1=np.zeros(len(rsbp_array)+len(r_array)+len(rmass_array)+len(rcsbp_array))
    reffi_2=np.zeros(len(rsbp_array)+len(r_array)+len(rmass_array)+len(rcsbp_array))
    for i in range(len(sum_t_array)):
#    print(i)
        t_fit_array=reassign(r_array,rne_array, sum_t_array[i])
        tproj_fit_array=reassign(r_array,rne_array,sum_tproj_array[i])
        mass_array=reassign(rmass_array,rne_array,sum_mcalc_array[i])
        nfw_fitted_array=reassign(rmass_array,rne_array,sum_mfit_array[i])
        tmp_r_use=rne_array
        sbp_fit_array=reassign(rsbp_array,rne_array,sum_sbp_array[i])
        if flag_csbp:
            csbp_fit_array=reassign(rcsbp_array,rne_array,sum_csbp_array[i])
        for j in range(len(rsbp_array)):
            a=sbp_fit_array[j]
            p_this=j
            reffi_1[p_this]=reffi_1[p_this]+np.power(a-sbp_model_array[j],2)
            reffi_2[p_this]=reffi_2[p_this]+np.power(a-sbp_array[j],2)
        for j in range(len(r_array)):
            if flag_tproj_array[j] == "1":
                t_this=t_fit_array[j]
                t_model_this=t_model[j]
            else:
                t_this=tproj_fit_array[j]
                t_model_this=tproj_model_array[j]
            p_this=len(rsbp_array)+j
            reffi_1[p_this]=reffi_1[p_this]+np.power(t_this-t_model_this,2)
            reffi_2[p_this]=reffi_2[p_this]+np.power(t_this-t_array[j],2)
        for j in range(len(rmass_array)):
            if rmass_array[j]>=30:
                Mnfw_model=nfw_fitted_array[j]
                M_this=mass_array[j]
                p_this=len(rsbp_array)+len(r_array)+j
                reffi_1[p_this]=reffi_1[p_this]+np.power(Mnfw_model-mass_model[j],2)
                reffi_2[p_this]=reffi_2[p_this]+np.power(Mnfw_model-M_this,2)
        if flag_csbp:
            for j in range(len(rcsbp_array)):
                csbp_fit=csbp_fit_array[j]
                p_this=len(rsbp_array)+len(r_array)+len(rmass_array)+j
                reffi_1[p_this]=reffi_1[p_this]+np.power(csbp_fit-csbp_model_array[j],2)
                reffi_2[p_this]=reffi_2[p_this]+np.power(csbp_fit-csbp_array[j],2)
    merr=0.15+te_array.sum()/t_array.sum()
    reffi_1=np.append(np.append(np.append(np.array(sbpe_array)+np.array(sbp_array)*0.10,np.array(te_array)+np.array(t_array)*0.08),np.array(mass_model)*0.20),np.array(csbpe_array)+np.array(csbp_array)*0.10)

    reffi_1=reffi_1*reffi_1*400
#print(reffi_1,reffi_2)
    reffi=1-reffi_2/reffi_1
#    print(reffi)
    reff=reffi.sum()/(len(t_array)+len(rsbp_array)+len(rmass_array)+len(rcsbp_array))
#    print(reff)
    return reff
if __name__=='__main__':
    name=sys.argv[1]
    main(name)
