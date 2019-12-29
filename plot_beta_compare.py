#!/usr/bin/env python3

import numpy as np
import sys
import matplotlib.pyplot as plt
import deproject_model
import json
import re
from calc_reff import reassign
def main(listfile):
    rwalker=np.array([0.204,0.367,0.520,0.649,0.776,0.979,1.086])
    swalker=np.array([0.68,1.142,1.556,1.758,1.974,1.683,1.872])
#    rwalker=np.linspace(0.1,2,1000)
#    swalker=4.4*np.power(rwalker,1.1)*np.exp(-rwalker**2)
    first=0
    rscale_array=np.linspace(0.1,2,1000)
    sum_kori_array=[]
    sum_kbeta_array=[]
    for f in open(listfile):
        f=f.split()[0]
        filename=f+'/beta_fit_result_all.npy'
        p_all=np.load(filename)
        p_center=p_all[0]
        p_use=p_center[0:6]
        json_file=f+'/'+f+'_plt.json'
        fi=open(json_file)
        dat=json.load(fi)
        rne_array=np.array(dat['radius_model'],dtype=float)
        t_array=np.array(dat['temperature_model'][0],dtype=float)
        kori_array=np.array(dat['k_fit'][0])
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
        sum_kori_array.append(kori_new_array)
        sum_kbeta_array.append(kbeta_new_array)
#        if first==0:
#            plt.loglog(rne_array/r200,kori_array/kori_norm,'b',label='new-model')
#            plt.loglog(rne_array/r200,kbeta_array/kbeta_norm,'g',label='beta')
#            first=1
#        else:
#            plt.loglog(rne_array/r200,kori_array/kori_norm,'b')
#            plt.loglog(rne_array/r200,kbeta_array/kbeta_norm,'g')
    plt.loglog(rwalker,swalker,'r',label='average profile from Walker')
    plt.fill_between(rwalker,swalker*0.7,swalker*1.4,color='k',alpha=0.3)
    sum_kori_array=np.array(sum_kori_array)
    sum_kbeta_array=np.array(sum_kbeta_array)
    print(sum_kori_array.shape)
    sum_kori_array.sort(0)
    sum_kbeta_array.sort(0)
    ind50=int(len(sum_kbeta_array)*0.5)
    ind16=int(len(sum_kbeta_array)*0.16)
    ind84=int(len(sum_kbeta_array)*0.84)
    plt.loglog(rscale_array,sum_kori_array[ind50],'b',label='ori model')
    plt.loglog(rscale_array,sum_kbeta_array[ind50],'g',label='beta model')
    plt.fill_between(rscale_array,sum_kori_array[ind16],sum_kori_array[ind84],color='b',alpha=0.3)
    plt.fill_between(rscale_array,sum_kbeta_array[ind16],sum_kbeta_array[ind84],color='g',alpha=0.3)
    plt.legend()
    plt.xlim(0.1,1.2)
    plt.ylim(0.4,5)
    plt.xlabel('radius (r/r200)')
    plt.ylabel('scaled entropy K/K(0.3r200)')
    plt.savefig('compare_with_beta.pdf')
    return 0

if __name__=='__main__':
    listfile=sys.argv[1]
    main(listfile)
