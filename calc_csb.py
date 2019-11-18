#!/usr/bin/env python3

import sys
import numpy as np
import re
import json

pi=3.1415926
def main(name):
#    for i in open('global.cfg'):
#        if re.match(r'^sbp_data_file',i):
#            sbp_data_file=i.split()[1]
#    r_array=[]
#    sbp_array=[]
#    for i in open(sbp_data_file):
#        if len(i)!=5:
#            print('format error for the sbp_data_file')
#            return -1
#        elif i.split()[4]='o':
    json_file=name+'_plt.json'
    fi=open(json_file,'r')
    dat=json.load(fi)
    r_array=np.array(dat['radius_model'],dtype=float)
    sbp_array=np.array(dat['sbp_model'][0],dtype=float)
    result_file=name+'_result.csv'
    for i in open(result_file):
        if re.match(r'^sbp_c,',i):
            sbp_c=float(i.split(',')[1])
    sbp_array=sbp_array-sbp_c
    c1=0
    c2=0
    for i in range(len(r_array)):
        if i==0:
            c1=c1+pi*sbp_array[0]
            c2=c1
        elif i==len(r_array)-2:
            print("something wrong")
            return -1
        else:
            rlow=r_array[i]-(r_array[i]-r_array[i-1])/2
            rhigh=r_array[i]+(r_array[i+1]-r_array[i])/2
            if rhigh > 400 and rlow <= 400:
                c2=c2+sbp_array[i]*pi*(400*400-rlow*rlow)
                break
            if rhigh <= 400:
                c2=c2+sbp_array[i]*pi*(rhigh*rhigh-rlow*rlow)
                if rhigh > 40 and rlow <= 40:
                    c1=c1+sbp_array[i]*pi*(40*40-rlow*rlow)
                if rhigh <= 40:
                    c1=c1+sbp_array[i]*pi*(rhigh*rhigh-rlow*rlow)
    csb=c1/c2
    print(csb)

if __name__=='__main__':
    name=sys.argv[1]
    main(name)





