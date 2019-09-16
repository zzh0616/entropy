#!/usr/bin/env python3

import numpy as np
import numpy
from scipy.integrate import quad
from time import time
import sys
#linear interpolation
pi=3.1415926
def lintp(x,x1,x2,y1,y2):
    if x>=x1-1e-10 and x<=x2+1e-10:
        return y1+(y2-y1)/(x2-x1)*(x-x1)
    else :
        return np.nan

def calc(a,r,p,ttal):
#    t1=time()
    rs=p[1]
    rho=p[0]
    delta1=p[2]
    delta2=p[3]
# verify if the parameters are in the reasonable range
    if delta1<=0 or delta1>3 :
        return np.nan
    if delta2<=2 or delta2>10 :
        return np.nan
    r_use=r/rs
#    if r/r_use<=0 or r/r_use>= :
#        return np.nan
    flag=0
    if a=='r':
        name2='_dvr.npy'
        flag=flag+1
    elif a=='n':
        name2='_ori.npy'
    if rs<20:
        return np.nan
    elif rs<200:
        name1='lrs'
        if r_use<0.0005 or r_use>161:
            return np.nan
        elif r_use<0.05:
            r_left=np.int(r_use/0.0005-1)
            r_right=r_left+1
            r_use_left=r_left*0.0005+0.0005
            r_use_right=r_use_left+0.0005
        elif r_use<10.05:
            r_left=np.int((r_use-0.05)/0.01+99)
            r_right=r_left+1
            r_use_left=(r_left-99)*0.01+0.05
            r_use_right=r_use_left+0.01
        else:
            r_left=np.int((r_use-10.05)/0.3+1099)
            r_right=r_left+1
            r_use_left=10.05+(r_left-1099)*0.3
            r_use_right=r_use_left+0.3
    else:
        name1='hrs'
        flag=flag+2
        if r_use<0.0002 or r_use>17.5:
            return np.nan
        elif r_use<0.02:
            r_left=np.int(r_use/0.0002-1)
            r_right=r_left+1
            r_use_left=0.0002*r_left+0.0002
            r_use_right=r_use_left+0.0002
        elif r_use<5.02:
            r_left=np.int((r_use-0.02)/0.004+99)
            r_right=r_left+1
            r_use_left=0.004*(r_left-99)+0.02
            r_use_right=r_use_left+0.004
        else:
            r_left=np.int((r_use-5.02)/0.05+1349)
            r_right=r_left+1
            r_use_left=5.02+0.05*(r_left-1349)
            r_use_right=r_use_left+0.05
    name=name1+name2
    t=ttal[flag]
    delta1_left=np.int((delta1-0.025)/0.025)
    delta1_right=delta1_left+1
    delta1_use_left=delta1_left*0.025+0.025
    delta1_use_right=delta1_use_left+0.025
    delta2_left=np.int((delta2-2.025)/0.025)
    delta2_right=delta2_left+1
    delta2_use_left=delta2_left*0.025+2.025
    delta2_use_right=delta2_use_left+0.025
    c_d1l_d2l_rl=t[delta1_left,delta2_left,r_left]
    c_d1r_d2l_rl=t[delta1_right,delta2_left,r_left]
    c_d1l_d2r_rl=t[delta1_left,delta2_right,r_left]
    c_d1r_d2r_rl=t[delta1_right,delta2_right,r_left]
    c_d1l_d2l_rr=t[delta1_left,delta2_left,r_right]
    c_d1r_d2l_rr=t[delta1_right,delta2_left,r_right]
    c_d1l_d2r_rr=t[delta1_left,delta2_right,r_right]
    c_d1r_d2r_rr=t[delta1_right,delta2_right,r_right]
    c_d2l_rl=lintp(delta1,delta1_use_left,delta1_use_right,c_d1l_d2l_rl,c_d1r_d2l_rl)
    c_d2r_rl=lintp(delta1,delta1_use_left,delta1_use_right,c_d1l_d2r_rl,c_d1r_d2r_rl)
    c_d2l_rr=lintp(delta1,delta1_use_left,delta1_use_right,c_d1l_d2l_rr,c_d1r_d2l_rr)
    c_d2r_rr=lintp(delta1,delta1_use_left,delta1_use_right,c_d1l_d2r_rr,c_d1r_d2r_rr)
    c_rl=lintp(delta2,delta2_use_left,delta2_use_right,c_d2l_rl,c_d2r_rl)
    c_rr=lintp(delta2,delta2_use_left,delta2_use_right,c_d2l_rr,c_d2r_rr)
    c=lintp(r_use,r_use_left,r_use_right,c_rl,c_rr)
#    print(r_use_left,r_use_right)
#    print(c_rl,c_rr,c)
    val=4*pi*rs*rs*rs*rho*c
    if a=='r':
        val=val/rs
#    t2=time()
#    print(t2-t1)
    return val

def mod_nfw(r,p):
    rho=p[0]
    rs=p[1]
    delta=p[2]
    delta2=p[3]
    den=rho/(np.power(r/rs,delta)*(np.power(1+r/rs,delta2-delta)))*4*pi*r*r
    return den

def mod_nfw_divbyr(r,p):
    a=mod_nfw(r,p)/r
    return a

def ori(a,r,p):
    if a=='r':
        M=quad(mod_nfw_divbyr,r,np.inf,p)[0]
    elif a=='n':
        M=quad(mod_nfw,0,r,p)[0]
    else:
        M=np.nan
    return M

def main():
    a=sys.argv[1]
    r=np.float(sys.argv[2])
    p=np.array(sys.argv[3].split(','),dtype=float)
    ta1=numpy.load('lrs_ori.npy')
    ta2=numpy.load('lrs_dvr.npy')
    ta3=numpy.load('hrs_ori.npy')
    ta4=numpy.load('hrs_dvr.npy')
    t_total=[ta1,ta2,ta3,ta4]
    print(a,r,p)
    t1=time()
    M_true=ori(a,r,p)
    t2=time()
    M_use=calc(a,r,p,t_total)
    t3=time()
    print('M_true:',M_true,' M_use:',M_use)
    print('T_int:',t2-t1,' T_appr:',t3-t2)

if __name__=='__main__':
    main()


