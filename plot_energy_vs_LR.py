#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import re
import json
import analyze_db
from astropy.cosmology import FlatLambdaCDM
from calc_reff import reassign
cosmo=FlatLambdaCDM(H0=71,Om0=0.27,Tcmb0=2.725)
kpc=3.086e19 #m
pi=3.1415926
LR_array=[]
energy_array=[]
Tave_array=[]
numtot_array=[]
csb_array=[]
csbe_array=[]
tcool_array=[]
m200_array=[]
for f in open(sys.argv[1],'r'):
    f=f[0:-1]
    filename=f+'/'+f+'_result.csv'
#    for i in open(filename,'r'):
#        if re.match(r'^k0,',i):
#            k0=np.float(i.split(',')[1])
    filename=f+'/'+f+'_suminfo.txt'
    for i in open(filename,'r'):
        if re.match(r'^Efeed',i):
            energy=np.float(i.split()[1])
        if re.match(r'^Tave',i):
            T_ave=np.float(i.split()[1])
        if re.match(r'^r200',i):
            r200=np.float(i.split()[1])
        if re.match(r'^csb',i):
            csb=np.float(i.split()[1])
            csbe=(np.float(i.split()[3])-np.float(i.split()[2]))/2
        if re.match(r'^tcool',i):
            tcool=(np.float(i.split()[1]))
        if re.match(r'^m200',i):
            m200=np.float(i.split()[1])
    tcool_array.append(tcool)
    csb_array.append(csb)
    csbe_array.append(csbe)
    m200_array.append(m200)
    filename=f+'/'+f+'_L1.4G.txt'
    for i in open(filename,'r'):
        LR=np.float(i)
    filename=f+'/param_zzh_for_py.txt'
    for i in open(filename,'r'):
        if re.match(r'^z',i):
            z=np.float(i.split()[1])
    dl=cosmo.luminosity_distance(z).value*kpc*1000 #m
    LR=LR*4*pi*dl*dl*1e-29 #W Hz^-1
    filename=f+'/gasnumber.npy'
    numtot=np.load(filename)
    numtot_array.append(numtot)
    LR_array.append(LR)
    energy_array.append(energy/np.power(T_ave,0.0))
    Tave_array.append(T_ave)
Tave_array=np.array(Tave_array)
numtot_array=np.array(numtot_array)
energy_array=np.array(energy_array)
avefeed_array=energy_array/numtot_array
plt.clf()
plt.loglog(energy_array,LR_array,'+')
plt.xlabel('feedback energy (kev)')
plt.ylabel('L(1.4GHz) (W Hz^-1)')
plt.savefig('energy_vs_LR.pdf')
plt.clf()
plt.loglog(Tave_array,energy_array,'+')
plt.xlabel('average temperature(kev)')
plt.ylabel('total feedback energy (kev)')
plt.savefig('T_vs_feedback.pdf')
plt.clf()
plt.loglog(m200_array,energy_array,'+')
plt.xlabel('m200 (Msun)')
plt.ylabel('feedback energy (keV)')
plt.savefig('m200_vs_feedback.pdf')
plt.clf()
SUM_Efeed_array=[]
SUM_Efeed_array_up=[]
SUM_Efeed_array_down=[]
SUM_CC_Efeed_array=[]
SUM_CC_Efeed_array_up=[]
SUM_CC_Efeed_array_down=[]
SUM_NCC_Efeed_array=[]
SUM_NCC_Efeed_array_up=[]
SUM_NCC_Efeed_array_down=[]
SUM_EICM_array=[]
SUM_CC_EICM_array=[]
SUM_NCC_EICM_array=[]
i=0
kc_cc_array=[]
kc_ncc_array=[]
LR_cc_array=[]
LR_ncc_array=[]
kc_array=[]
r_scaled_array=np.arange(0.001,1.4,0.001)
r_new_array=r_scaled_array*r200
sum_sbpe_array=[]
sum_kce_array=[]
for f in open(sys.argv[1],'r'):
    f=f[0:-1]
    print(f)
    filename=f+'/'+f+'_plt.json'
    fi=open(filename)
    dat=json.load(fi)
    r_array=np.array(dat['radius_model'],dtype=float)
    R200_0=np.float(dat['r200'])
    kc=np.float(dat['k_fit'][0][1])
    kc_array.append(kc)
    kce=(-np.float(dat['k_fit'][1][1])+np.float(dat['k_fit'][2][1]))/2
    sum_kce_array.append(kce)
    sbpe=(-np.float(dat['sbp_model'][1][1])+np.float(dat['sbp_model'][2][1]))/2/np.float(dat['sbp_model'][0][1])
    sum_sbpe_array.append(sbpe)
    filename=f+'/'+'sum_array.npy'
    sum_array=np.load(filename,allow_pickle=True)
    kfit_array=np.array(list(sum_array[0]),dtype=float)
    T_array=np.array(list(sum_array[1]),dtype=float)
    filename=f+'/p_all.npy'
    p=np.load(filename)
    for j in range(len(kfit_array)):
        j=int(len(kfit_array)/2)
        a0=p[2][j]
        k0=p[4][j]
        gamma0=p[3][j]
        n2=p[0][j]
        kori_array=a0*np.power(r_array,gamma0)+k0
        kfit=np.delete(kfit_array[j],0)
        EICM_array=n2*T_array[j]*(kfit-kori_array)/kfit
#        print(Efeed_array)
        EICM_new_array=reassign(r_new_array,r_array,EICM_array)
        SUM_EICM_array.append(EICM_new_array)
        filename=f+'/'+'Efeed_scaled_array.npy'
        [Efeed_array,Efeed_array_down,Efeed_array_up]=np.load(filename)
        SUM_Efeed_array.append(Efeed_array)
        SUM_Efeed_array_up.append(Efeed_array_up)
        SUM_Efeed_array_down.append(Efeed_array_down)
        plt.semilogx(LR_array[i],Efeed_array[0],'+')
        if tcool_array[i]>7.7:
#        plt.semilogx(r_scaled_array,Efeed_array/np.power(Tave_array[i],1),'r')
            SUM_NCC_Efeed_array.append(Efeed_array)
            SUM_NCC_Efeed_array_up.append(Efeed_array_up)
            SUM_NCC_Efeed_array_down.append(Efeed_array_down)
        else:
#        plt.semilogx(r_scaled_array,Efeed_array/np.power(Tave_array[i],1),'b')
            SUM_CC_Efeed_array.append(Efeed_array)
            SUM_CC_Efeed_array_up.append(Efeed_array_up)
            SUM_CC_Efeed_array_down.append(Efeed_array_down)
        break
    if tcool_array[i]<=7.7:
        LR_cc_array.append(LR_array[i])
        kc_cc_array.append(kc)
    else:
        LR_ncc_array.append(LR_array[i])
        kc_ncc_array.append(kc)
    i=i+1
#plt.savefig('LR_vs_Efeed0.pdf')
plt.clf()
plt.loglog(kc_cc_array,LR_cc_array,'b+')
plt.loglog(kc_ncc_array,LR_ncc_array,'r+')
plt.xlabel(r'central entropy (kev $\rm cm^2$)')
plt.ylabel(r'L(1.4 GHz) (W $\rm Hz^{-1}$)')
plt.savefig('LR_vs_kc.pdf')
plt.clf()
csb_array=np.array(csb_array)
csbe_array=np.array(csbe_array)
SUM_Efeed_array=np.array(SUM_Efeed_array)
#SUM_Efeed_array.sort(0)
#Efeed_center=SUM_Efeed_array[int(len(SUM_Efeed_array)/2)]
#Efeed_up=SUM_Efeed_array[int(len(SUM_Efeed_array)*0.84)]
#Efeed_down=SUM_Efeed_array[int(len(SUM_Efeed_array)*0.16)]
[Efeed_center,Efeed_down,Efeed_up]=analyze_db.sum_error_calc(SUM_Efeed_array,SUM_Efeed_array_down,SUM_Efeed_array_up)
print(Efeed_center)
#plt.semilogx(r_scaled_array,Efeed_center)
SUM_CC_Efeed_array=np.array(SUM_CC_Efeed_array)
#SUM_CC_Efeed_array.sort(0)
#Efeed_CC_center=SUM_CC_Efeed_array[int(len(SUM_CC_Efeed_array)/2)]
#Efeed_CC_up=SUM_CC_Efeed_array[int(len(SUM_CC_Efeed_array)*0.84)]
#Efeed_CC_down=SUM_CC_Efeed_array[int(len(SUM_CC_Efeed_array)*0.16)]
[Efeed_CC_center,Efeed_CC_down,Efeed_CC_up]=analyze_db.sum_error_calc(SUM_CC_Efeed_array,SUM_CC_Efeed_array_down,SUM_CC_Efeed_array_up)
plt.semilogx(r_scaled_array,Efeed_CC_center,'b')
SUM_NCC_Efeed_array=np.array(SUM_NCC_Efeed_array)
#SUM_NCC_Efeed_array.sort(0)
#Efeed_NCC_center=SUM_NCC_Efeed_array[int(len(SUM_NCC_Efeed_array)/2)]
#Efeed_NCC_up=SUM_NCC_Efeed_array[int(len(SUM_NCC_Efeed_array)*0.84)]
#Efeed_NCC_down=SUM_NCC_Efeed_array[int(len(SUM_NCC_Efeed_array)*0.16)]
[Efeed_NCC_center,Efeed_NCC_down,Efeed_NCC_up]=analyze_db.sum_error_calc(SUM_NCC_Efeed_array,SUM_NCC_Efeed_array_down,SUM_NCC_Efeed_array_up)
plt.loglog(r_scaled_array,Efeed_NCC_center,'r')
plt.xlim(0.03,2)
plt.ylim(0.05,50)
plt.fill_between(r_scaled_array,Efeed_CC_down,Efeed_CC_up+0.1,color='blue',alpha=0.7)
plt.fill_between(r_scaled_array,Efeed_NCC_down,Efeed_NCC_up+0.1,color='red',alpha=0.3)
plt.xlabel('Radius r/r200')
plt.ylabel('feedback energy per particle (kev)')
plt.savefig('Efeedback_vs_Radius.pdf')

plt.clf()
for i in range(len(Tave_array)):
    if tcool_array[i]>7.7:
        plt.loglog(Tave_array[i],energy_array[i],'r+')
    else:
        plt.loglog(Tave_array[i],energy_array[i],'b+')
plt.xlabel('Tave (kev)')
plt.ylabel('total feedback energy (kev)')
plt.savefig('csb_vs_feedback.pdf')
plt.clf()
f_array=[0,0,0,0]
for i in range(len(Tave_array)):
    if tcool_array[i]>7.7:
        if LR_array[i]>1e24:
            if f_array[0]==0:
                plt.loglog(kc_array[i],csb_array[i],'r*',label='NCC, AGN',markersize=8)
                f_array[0]=1
            else:
                plt.loglog(kc_array[i],csb_array[i],'r*',markersize=8)
            plt.errorbar(kc_array[i],csb_array[i],xerr=sum_kce_array[i],yerr=csbe_array[i],color='r')
        else:
            if f_array[1]==0:
                plt.loglog(kc_array[i],csb_array[i],'r+',label='NCC, normal',markersize=8)
                f_array[1]=1
            else:
                plt.loglog(kc_array[i],csb_array[i],'r+',markersize=8)
            plt.errorbar(kc_array[i],csb_array[i],xerr=sum_kce_array[i],yerr=csbe_array[i],color='r')
    else:
        if LR_array[i]>1e24:
            if f_array[2]==0:
                plt.loglog(kc_array[i],csb_array[i],'b*',label='CC, AGN',markersize=8)
                f_array[2]=1
            else:
                plt.loglog(kc_array[i],csb_array[i],'b*',markersize=8)
            plt.errorbar(kc_array[i],csb_array[i],xerr=sum_kce_array[i],yerr=csbe_array[i],color='b')
        else:
            if f_array[3]==0:
                plt.loglog(kc_array[i],csb_array[i],'b+',label='CC, normal',markersize=8)
                f_array[3]=1
            else:
                plt.loglog(kc_array[i],csb_array[i],'b+',markersize=8)
            plt.errorbar(kc_array[i],csb_array[i],xerr=sum_kce_array[i],yerr=csbe_array[i],color='b')
plt.xlabel('central entropy (kev cm^2)')
plt.ylabel('S ( < 40 kpc)/S ( < 400 kpc)')
plt.legend()
plt.savefig('kc_vs_csb.pdf')
plt.clf()
for f in open(sys.argv[1],'r'):
    f=f[0:-1]
    filename=f+'/'+f+'_result.csv'
    for i in open(filename):
        if re.match(r'^nth_a',i):
            nth_a=np.float(i.split(',')[1])
        if re.match(r'^nth_b',i):
            nth_b=np.float(i.split(',')[1])
        if re.match(r'nth_gamma',i):
            nth_gamma=np.float(i.split(',')[1])
    eta=1-nth_a*(1+np.exp(-np.power(r_array/R200_0/nth_b,nth_gamma)))
    plt.semilogx(r_array/r200,eta)
plt.xlabel('Radius (r/r200)')
plt.ylabel('none-thermal pressure fraction')
plt.savefig('nth_profile.pdf')
