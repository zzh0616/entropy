#!/usr/bin/env python3
import numpy as np
import json
import sys
import re
from calc_reff import reassign
def main(k,sum_array=[],mte=0.05,mse=0.05,mcse=0.05,ex_reg=[]):
    pi=3.1415926
    tmp1=list(range(1,11))
    tmp2=list(range(11,41,3))
    tmp3=list(range(41,3001,5))
    tmp1.extend(tmp2)
    tmp1.extend(tmp3)
    rne_array=np.array(tmp1)
    for i in open('param_zzh_for_py.txt'):
        if re.match(r'aa\s',i):
            aa=float(i.split()[1])
            bb=float(i.split()[2])
            cc=float(i.split()[3])
        if re.match(r'^R500',i):
            r500=float(i.split()[1])
        if re.match(r'^R200',i):
            r200=float(i.split()[1])
    ind_50=int((aa-bb)/cc*0.50)
    ind_84=int((aa-bb)/cc*0.84)
    ind_16=int((aa-bb)/cc*0.16)
    if ex_reg == []:
        ex_reg=0.1*r200
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
    sbp_array=np.array(sbp_array)
    sbpe_array=np.array(sbpe_array)
    csbp_array=np.array(csbp_array)
    csbpe_array=np.array(csbpe_array)
    te_array=te_array+t_array*mte
    csbpe_array=csbpe_array+csbp_array*mcse
    sbpe_array=sbpe_array+sbp_array*mse
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
    t_model_array=reassign(r_array,rne_array,t_array_center)
    tproj_model_array=reassign(r_array,rne_array,tproj_array_center)
    sbp_model_array=reassign(rsbp_array,rne_array,sbp_array_center)
    if flag_csbp:
        csbp_model_array=reassign(rcsbp_array,rne_array,csbp_array_center)
    rtmp=np.array(list(range(45)))
    rmass_array=30*np.power(1.11,rtmp)
    for j in range(len(rmass_array)):
        rmass_array[j]=np.int(rmass_array[j])
    rmass_array=np.insert(rmass_array,0,10)
    rmass_array=np.insert(rmass_array,1,20)
    mass_model=reassign(rmass_array,rne_array,mfit_array_center)
    lhood=0
    n=0
    print(ex_reg,mte,mse,mcse)
    for i in range(len(r_array)):
        if flag_tproj_array[i]=='1':
            t_model_this=tproj_model_array[i]
        else:
            t_model_this=t_model_array[i]
        lthis=1/np.sqrt(2*pi)/te_array[i]*np.exp(-(t_array[i]-t_model_this)**2/2/te_array[i]**2)
        if r_array[i]>=ex_reg:
            lhood=lhood+np.log(lthis)
#        print(t_array[i],te_array[i],t_model_this)
            n=n+1
    for i in range(len(sbp_array)):
        lthis=1/np.sqrt(2*pi)/sbpe_array[i]*np.exp(-(sbp_array[i]-sbp_model_array[i])**2/2/sbpe_array[i]**2)
        if rsbp_array[i]>ex_reg:
            lhood=lhood+np.log(lthis)
#        print(sbp_array[i],sbpe_array[i],sbp_model_array[i])
            n=n+1
    if flag_csbp:
        for i in range(len(csbp_array)):
            lthis=1/np.sqrt(2*pi)/csbpe_array[i]*np.exp(-(csbp_array[i]-csbp_model_array[i])**2/2/csbpe_array[i]**2)
            if rcsbp_array[i]>ex_reg:
                lhood=lhood+np.log(lthis)
#            print(csbp_array[i],csbpe_array[i],csbp_model_array[i])
                n=n+1
    aic=2*k-2*lhood
#    aicc=aic+(2*k*k+2*k)/(n-k-1)
    filename='aic_info.txt'
    fi=open(filename,'w')
    print('aic=',aic,file=fi)
 #   print('aicc=',aicc,file=fi)
    print(aic)
    return aic

if __name__ == '__main__':
    model_type=sys.argv[1]
    if model_type=='clump':
        num_mp=16
    elif model_type=='ori':
        num_mp=21
    elif model_type=='plain':
        num_mp=13
    elif model_type=='hydro':
        num_mp=18
    else:
        print('model type error')
        exit
    if len(sys.argv)==5:
        a=float(sys.argv[2])
        b=float(sys.argv[3])
        c=float(sys.argv[4])
        ex_reg=[]
    elif len(sys.argv)==6:
        a=float(sys.argv[2])
        b=float(sys.argv[3])
        c=float(sys.argv[4])
        ex_reg=float(sys.argv[5])
    else:
        a=0.03
        b=0.01
        c=0.01
        ex_reg=[]
    main(num_mp,[],a,b,c,ex_reg)



