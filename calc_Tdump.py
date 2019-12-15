#!/usr/bin/env python3
import numpy as np
import sys
import json

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
        print(r_array[i],T_array[i],file=fi)
    for i in range(10):
        print(r_array[-1]+5+5*i,T_array[-1],file=fi)
    fi.close()
    return 0
if __name__=='__main__':
    main()
