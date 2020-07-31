#!/usr/bin/env python3

import numpy as np
import sys
import matplotlib.pyplot as plt
import deproject_model
import json
import re
from calc_reff import reassign
pi=3.1416
def c_mod(r_array,ne_array,rind1,rind2):
    ne_array=r_array
    r1=r_array[rind1]
    r2=r_array[rind2]
    print(r1,r2)
    print(ne_array[rind1],ne_array[rind2])
#    vtot=np.abs(4/3*pi*(np.power(r2,3)-np.power(r1,3)))
    gnsum=0
    emsum=0
    vtot=0
    for i in range(rind1,rind2+1):
        tmp=ne_array[i]*(np.power(r_array[i+1],3)-np.power(r_array[i-1],3))
        gnsum=gnsum+tmp
        emsum=emsum+tmp*ne_array[i]
        vtot=vtot+tmp/ne_array[i]
    gndave=gnsum/vtot
    gndcl=np.sqrt(emsum/vtot)
    return gndcl/gndave
def main(listfile,flag_ori_only):
    rwalker=np.array([0.202,0.266,0.385,0.474,0.559,0.774,0.984,1.094])
    swalker=np.array([0.703,0.934,1.188,1.397,1.642,1.954,1.693,1.872])
    swalker_low=np.array([0.873,1.103,1.388,1.794,2.085,3.194,2.236,2.671])
    swalker_up=np.array([0.456,0.677,0.876,1.100,1.197,1.026,1.364,1.269])
#    rwalker=np.linspace(0.1,2,1000)
#    swalker=4.4*np.power(rwalker,1.1)*np.exp(-rwalker**2)
    first=0
    rscale_array=np.linspace(0.1,2,1000)
    sum_kori_array=[]
    sum_kbeta_array=[]
    sum_kcl_array=[]
    for f in open(listfile):
        f=f.split()[0]
        if flag_ori_only=='0':
            filename=f+'/beta_fit_result_all.npy'
        else:
            filename=f+'/beta_fit_result_all_orionly.npy'
        p_all=np.load(filename)
        p_center=p_all[0]
        p_use=p_center[0:6]
        json_file=f+'/'+f+'_plt.json'
        fi=open(json_file)
        dat=json.load(fi)
        rne_array=np.array(dat['radius_model'],dtype=float)
        t_array=np.array(dat['temperature_model'][0],dtype=float)
        ne_array=np.array(dat['den_model'][0],dtype=float)
        necl_array=np.array(dat['den_cl_model'][0],dtype=float)
        cl_add_fac=np.power(necl_array/ne_array,0.2)
#        test=c_mod(rne_array,ne_array,400,600)
#        print(test)
        kori_array=np.array(dat['k_fit'][0])
        kori_array=t_array*np.power(ne_array,-2/3)
        kcl_array=t_array*np.power(necl_array,-2/3)/cl_add_fac
        den=deproject_model.dbeta(rne_array,*p_use)
        kbeta_array=t_array*np.power(den,-2/3)
        info_file=f+'/'+f+'_suminfo.txt'
        for i in open(info_file):
            if re.match(r'^r200',i):
                r200=float(i.split()[1])
        rnorm=0.3*r200
        for i in range(len(rne_array)):
            if rne_array[i]>=rnorm:
                kori_norm=kori_array[i]
                kbeta_norm=kbeta_array[i]
                break
        rnew_array=rscale_array*r200
        kori_new_array=reassign(rnew_array,rne_array,kori_array/kori_norm)
        kbeta_new_array=reassign(rnew_array,rne_array,kbeta_array/kbeta_norm)
        kcl_new_array=reassign(rnew_array,rne_array,kcl_array/kori_norm)
        sum_kori_array.append(kori_new_array)
        sum_kbeta_array.append(kbeta_new_array)
        sum_kcl_array.append(kcl_new_array)
#        if first==0:
#            plt.loglog(rne_array/r200,kori_array/kori_norm,'b',label='new-model')
#            plt.loglog(rne_array/r200,kbeta_array/kbeta_norm,'g',label='beta')
#            first=1
#        else:
#            plt.loglog(rne_array/r200,kori_array/kori_norm,'b')
#            plt.loglog(rne_array/r200,kbeta_array/kbeta_norm,'g')
    plt.loglog(rwalker,swalker,'k',label='Walker et al. 2012a')
    plt.fill_between(rwalker,swalker_low,swalker_up,color='k',alpha=0.3)
    sum_kori_array=np.array(sum_kori_array)
    sum_kbeta_array=np.array(sum_kbeta_array)
    sum_kcl_array=np.array(sum_kcl_array)
    print(sum_kori_array.shape)
    sum_kori_array.sort(0)
    sum_kbeta_array.sort(0)
    sum_kcl_array.sort(0)
    ind50=int(len(sum_kbeta_array)*0.5)
    ind16=int(len(sum_kbeta_array)*0.16)
    ind84=int(len(sum_kbeta_array)*0.84)
#    plt.loglog(rscale_array,sum_kori_array[ind50],'b',label='TSM model')
#    plt.loglog(rscale_array,sum_kbeta_array[ind50],'g',label=r'Double-$\beta$ model')
#    plt.fill_between(rscale_array,sum_kori_array[ind16],sum_kori_array[ind84],color='b',alpha=0.3)
#    plt.fill_between(rscale_array,sum_kbeta_array[ind16],sum_kbeta_array[ind84],color='g',alpha=0.3)
    plt.loglog(rscale_array,sum_kcl_array[ind50]*1.1,'g',label='RTI model prediction (gas clumping switched off)')
    plt.fill_between(rscale_array,sum_kcl_array[ind16]*1.1,sum_kori_array[ind84]*1.1,color='g',alpha=0.3)
    plt.legend(loc='upper left')
    plt.xlim(0.1,1.09)
    plt.ylim(0.4,5)
    plt.xlabel(r'Radius ($r/r_{200}$)')
    plt.ylabel(r'Entropy $K/K(0.3r_{200})$')
    if flag_ori_only=='0':
        plt.savefig('compare_with_beta.pdf')
    else:
        plt.savefig('compare_with_beta_orionly.pdf')
    return 0

if __name__=='__main__':
    listfile=sys.argv[1]
    flag_ori_only=sys.argv[2]
    main(listfile,flag_ori_only)
