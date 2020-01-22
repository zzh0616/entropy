#!/usr/bin/env python3
import numpy as np
import sys
import json
import re

def main(name):
    json_file=name+'_plt.json'
    fi=open(json_file,'r')
    dat=json.load(fi)
    T_array=np.array(dat['temperature_model'][0],dtype=float)
    r_array=np.array(dat['radius_model'],dtype=float)
    fi.close()
    fi=open('T_dump_tmp.dat','w')
    ra_array=[]
    va_array=[]
    rnh_array=[]
    vnh_array=[]
    fia=open('A_dump_tmp.dat','w')
    finh=open('NH_dump_tmp.dat','w')
    abund_file=name+'_abund.txt'
    nh_file=name+'_nh.txt'
    for i in open(abund_file):
        r_a,v_a=np.array(i.split(),dtype=float)
        ra_array.append(r_a)
        va_array.append(v_a)
    for i in open(nh_file):
        r_nh,v_nh=np.array(i.split(),dtype=float)
        rnh_array.append(r_nh)
        vnh_array.append(v_nh)
# start to output
    for i in range(len(r_array)):
        if T_array[i]<=0.08:
            T_array[i]=0.08
        if T_array[i]>20:
            T_array[i]=20
        if r_array[i]<=ra_array[0]:
            abund=va_array[0]
        elif r_array[i]>=ra_array[-1]:
            abund=va_array[-1]
        else:
            for j in range(len(ra_array)):
                if r_array[i]>ra_array[j]:
                    abund=va_array[j]+(va_array[j+1]-va_array[j])/(ra_array[j+1]-ra_array[j])*(r_array[i]-ra_array[j])
        print('%.4e %.4e'%(r_array[i],abund),file=fia)
        if r_array[i]<=rnh_array[0]:
            nh=vnh_array[0]
        elif r_array[i]>=rnh_array[-1]:
            nh=vnh_array[-1]
        else:
            for j in range(len(rnh_array)):
                if r_array[i]>rnh_array[j]:
                    nh=vnh_array[j]+(vnh_array[j+1]-vnh_array[j])/(rnh_array[j+1]-rnh_array[j])*(r_array[i]-rnh_array[j])
        print('%.4e %.4e'%(r_array[i],nh),file=finh)
        print(r_array[i],T_array[i],file=fi)
    for i in range(10):
        print(r_array[-1]+5+5*i,T_array[-1],file=fi)
        print(r_array[-1]+5+5*i,va_array[-1],file=fia)
        print(r_array[-1]+5+5*i,vnh_array[-1],file=finh)
    fia.close()
    fi.close()
    finh.close()

    return 0

def calc_a0(name):
    info_file=name+'_suminfo.txt'
    for f in open(info_file,'r'):
        if re.match(r'^r200',i):
            r200=float(i.split()[1])
        if re.match(r'^k200',i):
            k200=float(i.split()[1])
    a0=k200*1.43/(r200**1.1)
    print(a0)

if __name__=='__main__':
    name=sys.argv[1]
    mode=sys.argv[2]
    if mode == "calc_a0":
        calc_a0(name)
    else:
        main(name)
