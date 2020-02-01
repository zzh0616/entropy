#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
import sys
import json
import re
from calc_reff import reassign
def main(dir_data="../data",dir_plain="../plain",dir_clump="../clump",listfile="list.txt"):
    r_scaled_array=np.arange(0.001,1.4,0.001)
    sum_kclump=[]
    sum_kdata=[]
    sum_kplain=[]
    for f in open(listfile):
        name=f.split(',')[0]
        info_file=dir_clump+'/'+name+'/'+name+"_suminfo.txt"
        for i in open(info_file):
            if re.match(r'^r200',i):
                r200=float(i.split()[1])
            if re.match(r'^k200',i):
                k200=float(i.split()[1])
        knorm=k200
        rnorm=r200
        for d in [dir_clump,dir_data,dir_plain]:
            info_file=d+'/'+name+'/'+name+"_plt.json"
            fi=open(info_file,'r')
            dat=json.load(fi)
            kfit=np.array(dat["k_fit"][0],dtype=float)
            kfit_down=np.array(dat["k_fit"][1],dtype=float)
            kfit_up=np.array(dat["k_fit"][2],dtype=float)
            rmodel=np.array(dat["radius_model"],dtype=float)
            k_scaled=kfit/knorm
            rnew=r_scaled_array*r200
            knew=reassign(rnew,rmodel,k_scaled)
            if d==dir_clump:
                sum_kclump.append(knew)
            elif d==dir_data:
                sum_kdata.append(knew)
            elif d==dir_plain:
                sum_kplain.append(knew)
    sum_kdata=np.array(sum_kdata)
    sum_kclump=np.array(sum_kclump)
    sum_kplain=np.array(sum_kplain)
    sum_kclump.sort(0)
    sum_kdata.sort(0)
    sum_kplain.sort(0)
    ind50=int(len(sum_kclump)*0.5)
    ind16=int(len(sum_kclump)*0.16)
    ind84=int(len(sum_kclump)*0.84)
    plt.loglog(r_scaled_array,sum_kplain[ind50],'g',label='Case 1')
    plt.fill_between(r_scaled_array,sum_kplain[ind16],sum_kplain[ind84],color='lightgreen',alpha=1.0)
    plt.loglog(r_scaled_array,sum_kclump[ind50],'k',label='Case 2')
    plt.fill_between(r_scaled_array,sum_kclump[ind16],sum_kclump[ind84],color='lightgrey',alpha=1.0)
    plt.loglog(r_scaled_array,sum_kdata[ind50],'b',label='Case 3')
    plt.fill_between(r_scaled_array,sum_kdata[ind16],sum_kdata[ind84],color='lightblue',alpha=1.0)
    plt.xlabel(r'Radius ($r/r_{200}$)')
    plt.ylabel(r'Entropy ($K/K_{200}$)')
    plt.legend()
    plt.savefig('entropy_model_compare.pdf')
    return 0

if __name__=="__main__":
    main()
