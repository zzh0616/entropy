#!/usr/bin/env python3

import sys
import numpy
import scipy
import re
import numpy as np
import math
pi=3.1416
cm_per_kpc=3.0857e21
## calc_v_ring: calc the volume:
## a sphere(rout) sbutracted by a cycline(rin)
## input in unit cm output in unit cm^3
def calc_v_ring(rout,rin):
    if rin<rout:
        tmp=rout*rout-rin*rin
        return 4.18879*math.sqrt(tmp*tmp*tmp) #4.18879=4/3*pi
    return 0

## calculate the surface brightness on the given radius,
## with given density profile and cooling function
## r is the radius in unit kpc to calculate
## rlist gives the outer raidus of each annuli
## the density and cooling function is the average of each annuli
## rlist[0]=den[0]=cool_func[0]=0
def calc_sb(r,r_list,den,cool_func,cm_per_pixel):
    ne_np_ratio=1.2
    tot=len(r_list)
#    kpc_per_pixel=cm_per_pixel/cm_per_kpc
    r_use=-1
    r_list=numpy.array(r_list,dtype=float)
    den=numpy.array(den,dtype=float)
    cool_func=numpy.array(cool_func,dtype=float)
    r_list_ori=r_list
    r_list=r_list*cm_per_kpc
    tmp_array=den*den*cool_func/ne_np_ratio
    for i in range(tot):
        if r<r_list_ori[i]:
            r_use=i
            break
    if r_use==-1:
        return numpy.NaN
    if tot>=r_use:
#        vol_total=0
        sb_tmp=0
        for i in range(tot)[r_use:tot]:
            vol=calc_v_ring(r_list[i],r_list[r_use-1])-calc_v_ring(\
                    r_list[i-1],r_list[r_use-1])-calc_v_ring(r_list[i],\
                    r_list[r_use])+calc_v_ring(r_list[i-1],r_list[r_use])
#            vol_total=vol_total+vol
            sb_tmp=sb_tmp+tmp_array[i]*vol
        pixel_out=r_list[r_use]/cm_per_pixel
        pixel_in=r_list[r_use-1]/cm_per_pixel
        area=pi*(pixel_out*pixel_out-pixel_in*pixel_in)
        sb_tmp=sb_tmp/area
#        print(r_use,vol_total)
        return sb_tmp

#calculate 2d temperature porfile
def calc_projT(r,r_list,T_array,den_array):
    tot=len(r_list)
    i_use=-1
    T_array=np.array(T_array,dtype=float)
    den_array=np.array(den_array,dtype=float)
    tmp=den_array*den_array
    tmp1=tmp*np.power(T_array,0.25)
    tmp2=tmp*np.power(T_array,-0.75)
    for i in range(len(tmp1)):
        if np.isnan(tmp1[i]) or np.isinf(tmp1[i]):
            tmp1[i]=0
        if np.isnan(tmp2[i]) or np.isinf(tmp2[i]):
            tmp2[i]=0
    for i in range(tot):
        if r<r_list[i]:
            i_use=i
            break
    if i_use==-1:
        return np.NaN
    if tot>i_use:
        w1=0
        w2=0
        for i in range(tot)[i_use:tot]:
            vol=calc_v_ring(r_list[i],r_list[i_use-1])-calc_v_ring(\
                    r_list[i-1],r_list[i_use-1])-calc_v_ring(r_list[i],\
                    r_list[i_use])+calc_v_ring(r_list[i-1],r_list[i_use])
#            w1=w1+den_array[i]*den_array[i]*math.pow(T_array,0.25)
#            w2=w2+den_array[i]*den_array[i]*math.pow(T_array,-0.75)
            w1=w1+tmp1[i]*vol
            w2=w2+tmp2[i]*vol
        T2d=w1/w2
        return T2d
## beta model
def beta(r,a,beta,rc):
    return a*numpy.power(1+r*r/rc/rc,-3.0/2.0*beta)
## dbeta model
def dbeta(r,a0,beta0,rc0,a1,beta1,rc1):
    return a0*numpy.power(1+r*r/rc0/rc0,-3.0/2.0*beta0)+\
            a1*numpy.power(1+r*r/rc1/rc1,-3.0/2.0*beta1)
## nfw model output:msun
def Mnfw(r,rho0,rs):
    return 4*pi*rho0*rs*rs*rs*(numpy.log((rs+r)/rs)-r/(rs+r))
