#!/usr/bin/env python3

import numpy as np
from scipy.integrate import quad
from time import time
#from multiprocessing import Process
#import multiprocessing as mp

def f(x,p):
    return np.power(x,1)/np.power(x,p[0])/np.power(1+x,p[1]-p[0])

def main():
    tot=np.zeros([120,320,1600])
    a1=time()
    for i in range(tot.shape[0]):
        i_use=i*0.025+0.025
        a2=time()
        print(a2-a1)
        for j in range(tot.shape[1]):
            j_use=j*0.025+2.025
            print(j/320,i/120)
            for k in range(tot.shape[2]):
                if k<99:
                    k_use=k*0.0002+0.0002
                elif k<1349:
                    k_use=(k-99)*0.004+0.0200    # for rs>200
                else:
                    k_use=5.0200+(k-1349)*0.05
                tot[i,j,k]=quad(f,k_use,np.inf,[i_use,j_use])[0]
#                print(i_use,j_use,k_use)

    fi=np.save('hrs_dvr',tot)


if __name__=='__main__':
    main()
