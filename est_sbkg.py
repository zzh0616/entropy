#!/usr/bin/env python3

import numpy as np
import json
import sys
import re
from astropy.cosmology import FlatLambdaCDM
import warnings
cosmo=FlatLambdaCDM(H0=71,Om0=0.27)
def main(name):
    result_file=name+'_result.csv'
    for i in open(result_file,'r'):
        if re.match('sbp_c',i):
            sbp_c=np.float(i.split(',')[1])
    for i in open('param_zzh_for_py.txt'):
        if re.match(r'^z',i):
            z=np.float(i.split()[1])
    for i in open('aic_flag.txt'):
        if re.match(r'^aic-flag',i):
            aic_flag=i.split()[1]
        if re.match(r'^aic-diff',i):
            aic_diff=i.split()[1]
    for i in open('global.cfg'):
        if re.match(r'^sbp_data_file',i):
            sbp_data_file=i.split()[1]
    exp_file='../exposure.txt'
    expo=0
    if '+' in name:
        name_use=name.replace('+','\+')
    else:
        name_use=name
    for i in open(exp_file):
        if re.match(name_use+'\s',i):
            expo=np.float(i.split()[1])
    if expo==0:
        warnings.warn('check the exposure time for '+name )
        expo=30
    rsbp_array=[]
    rsbpe_array=[]
    sbp_array=[]
    sbpe_array=[]
    rcsbp_array=[]
    rcsbpe_array=[]
    csbp_array=[]
    csbpe_array=[]
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
    json_file=name+'_plt.json'
    fi=open(json_file,'r')
    dat=json.load(fi)
    rne_array=np.array(dat['radius_model'],dtype=float)
    sbp_model=np.array(dat['sbp_model'][0],dtype=float)
    sum_file=name+'_suminfo.txt'
    for i in open(sum_file,'r'):
        if re.match('r500:',i):
            r500=np.float(i.split()[1])
        if re.match('r200:',i):
            r200=np.float(i.split()[1])
        if re.match('Tave',i):
            T=np.float(i.split()[1])
        if re.match('csb',i):
            csb=np.float(i.split()[1])
        if re.match('tcool',i):
            tcool=np.round(np.float(i.split()[1]),3)
        if re.match('lx200:',i):
            lx200=np.float(i.split()[1])
        if re.match('m500:',i):
            m500=np.float(i.split()[1])
    dl=cosmo.luminosity_distance(z).value #Mpc
    kpc_per_arcmin=cosmo.kpc_proper_per_arcmin(z).value #kpc/arcmin
    angle_size=r200/kpc_per_arcmin
    fx=np.round(lx200/dl/dl,3)
    for i in open('pcount.txt'):
        pcount=i[0:-1]
    for i in range(len(rne_array)):
        if rne_array[i]>=r200:
            ind=i
            break
#    print(sbp_c,rne_array[400],sbp_model[400])
    snr=(sbp_model[ind]-sbp_c)/sbp_c
    snr=np.round(snr,3)
    snr_m=np.round((sbp_array[-1]-sbp_c)/sbp_c,3)
    for i in open(name+'_csb.txt'):
        if re.match('csb',i):
            ocsb=np.round(float(i.split()[1]),3)
            ocsberr=np.round(float(i.split()[3]),3)
        if re.match('fx',i):
            ofx=float(i.split()[1])
            ofxerr=float(i.split()[3])
    print(name,rne_array[ind],snr,z,snr_m,T,csb,tcool,fx,pcount,aic_flag,aic_diff,m500,expo,angle_size,lx200,ocsb,ofx,ocsberr,ofxerr)
#    print(name,ofx,ocsb,aic_flag)

    return snr

if __name__=='__main__':
    name=sys.argv[1]
    main(name)
