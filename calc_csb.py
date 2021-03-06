#!/usr/bin/env python3

import sys
import numpy as np
import re
import json
from astropy.cosmology import FlatLambdaCDM
cosmo=FlatLambdaCDM(H0=71,Om0=0.27,Tcmb0=2.725)
pi=3.1415926
mpc=3.086e24 #cm
def main(name='',iput='no',sbp_array=[],t_array=[],ne_array=[],sbp_c=0,r500=-1,calcerr=False):
    json_file=name+'_plt.json'
    fi=open(json_file,'r')
    dat=json.load(fi)
    fi.close()
    r_array=np.array(dat['radius_model'],dtype=float)
    if iput=='no':
        sbp_array=np.array(dat['sbp_model'][0],dtype=float)
        t_array=np.array(dat['temperature_model'][0],dtype=float)
        ne_array=np.array(dat['den_cl_model'][0],dtype=float)
        result_file=name+'_result.csv'
        for i in open(result_file):
            if re.match(r'^sbp_c,',i):
                sbp_c=float(i.split(',')[1])
        summary_file=name+'_suminfo.txt'
        for i in open(summary_file):
            if re.match(r'^r500',i):
                r500=float(i.split()[1])
            if re.match(r'^r200',i):
                r200=float(i.split()[1])
            if re.match(r'^csb',i):
                csberr=(float(i.split()[3])-float(i.split()[2]))/2/float(i.split()[1])
            if re.match(r'^lx200',i):
                fxerr=(float(i.split()[3])-float(i.split()[2]))/2/float(i.split()[1])
    if r500 == -1:
        print('input error: please check the input of r500')
        return -1
    if iput=='yes':
        csberr=-1
        fxerr=-1
    hlim=1.0*r500
    llim=0.7*r500
    ex_reg=0*r500
    if iput=='yes':
        llim=0.2*r500
    sbp_array=sbp_array-sbp_c*0
    param_file='param_zzh_for_py.txt'
    for i in open(param_file):
        if re.match(r'^z',i):
            z=float(i.split()[1])
    dl=cosmo.luminosity_distance(z).value #Mpc
    dm=cosmo.comoving_distance(z).value  #Mpc
    cfunc_file='cfunc_for_lxcalc.txt'
    rcfunc_array=[]
    cfunc_array=[]
    for i in open(cfunc_file):
        rcfunc_array.append(float(i.split()[0]))
        cfunc_array.append(float(i.split()[1]))
    summary_file=name+'_suminfo.txt'
    c1=0
    c2=0
    if len(rcfunc_array)-len(r_array) != 10 or r_array[10] != rcfunc_array[10]:
        print('warning: please check the cooling function file')
        print(len(rcfunc_array),len(r_array),r_array[10],rcfunc_array[10])
    for i in range(len(r_array)):
        if r_array[i] <ex_reg:
            continue
        if i==0:
            c1=c1+pi*sbp_array[0]
            c2=c1
        elif i==len(r_array)-2:
            print("something wrong")
            return -1
        else:
            rlow=r_array[i]-(r_array[i]-r_array[i-1])/2
            rhigh=r_array[i]+(r_array[i+1]-r_array[i])/2
            if rhigh > hlim and rlow <= hlim:
                c2=c2+sbp_array[i]*pi*(hlim*hlim-rlow*rlow)
                break
            if rhigh <= hlim:
                c2=c2+sbp_array[i]*pi*(rhigh*rhigh-rlow*rlow)
                if rhigh > llim and rlow <= llim:
                    c1=c1+sbp_array[i]*pi*(llim*llim-rlow*rlow)
                if rhigh <= llim:
                    c1=c1+sbp_array[i]*pi*(rhigh*rhigh-rlow*rlow)
    csb=1-c1/c2
    if iput == 'no':
        if ex_reg != -1:
            outfile=name+'_csb.txt'
            fi=open(outfile,'w')
            print('csb=',csb,'csberr=',csberr,'ex=',ex_reg,'in=',llim,'out=',hlim,file=fi)
    ex_reg=0.0*r500
    hlim=1.0*r500
    lx=0
    for i in range(len(r_array)):
        if r_array[i]<ex_reg:
            continue
        rhigh=r_array[i]+(r_array[i+1]-r_array[i])/2
        rlow=r_array[i]-(r_array[i]-r_array[i-1])/2
        if rhigh > hlim :
            break
        lx=lx+cfunc_array[i]*ne_array[i]*ne_array[i]/1.2*4/3*3.1416*(rhigh**3-rlow**3)*(mpc/1000)**3*4*pi*dm*dm
    fx=lx/4/pi/dl/dl  #erg/cm^2/s
    etot=0
    erad=0
    for i in range(len(r_array)):
        if i==0:
            continue
        if r_array[i]<0.048*r500:
            T_this=t_array[i]
            ne_this=ne_array[i]
            cfunc_this=cfunc_array[i]
            etot=etot+ne_this*1.93*1.5*T_this*1.602e-9 #erg
            erad=erad+ne_array[i]*ne_array[i]/1.18*cfunc_this*4*pi*dl*dl*mpc*mpc #erg/s
    tcool=etot/erad/3600/24/365.2/1e9
    for i in range(len(r_array)):
        if r_array[i]>0.04*r500:
            cusp=-(np.log(ne_array[i])-np.log(ne_array[i-1]))/(np.log(r_array[i]-np.log(r_array[i-1])))
            break
    if iput == 'no':
        print(csb,csberr,tcool,fx,fxerr)
        if ex_reg != -1:
            print('tcool=',tcool,file=fi)
            print('fx=',fx,'err=',fxerr,'ex=',ex_reg,'hlim=',hlim,file=fi)
    return [csb,tcool]
if __name__=='__main__':
    name=sys.argv[1]
    if len(sys.argv)==3:
        calcerr=sys.argv[3]
    else:
        calcerr=False
    main(name=name,sbp_array=[],t_array=[],ne_array=[],sbp_c=0,iput='no',calcerr=calcerr)





