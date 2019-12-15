#!/usr/bin/env python3
import numpy as np
import sys
import json
import re

def main():
    name=sys.argv[1]
    json_file=name+'_plt.json'
    fi=open(json_file,'r')
    dat=json.load(fi)
    T_array=np.array(dat['temperature_model'][0],dtype=float)
    r_array=np.array(dat['radius_model'],dtype=float)
    fi.close()
    fi=open('T_dump_tmp.dat','w')
    for i in range(len(r_array)):
        if T_array[i]<=0.08:
            T_array[i]=0.08
        if T_array[i]>20:
            T_array[i]=20
        print(r_array[i],T_array[i],file=fi)
    for i in range(10):
        print(r_array[-1]+5+5*i,T_array[-1],file=fi)
    fi.close()
    return 0

def calc_a0():
    name=sys.argv[1]
    info_file=name+'_suminfo.txt'
    for f in open(info_file,'r'):
        if re.match(r'^r200',i):
            r200=float(i.split()[1])
        if re.match(r'^k200',i):
            k200=float(i.split()[1])
    a0=k200*1.43/(r200**1.1)
    print(a0)

if __name__=='__main__':
    main()
