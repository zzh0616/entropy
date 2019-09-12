#!/usr/bin/env python3

import sys
import numpy
import scipy
import re

pi=3.1416
cm_per_kpc=3.0857e21
## calc_v_ring: calc the volume:
## a sphere(rout) sbutracted by a cycline(rin)
## input in unit kpc output in unit cm^3
def calc_v_ring(rout,rin):
    rin=rin*cm_per_kpc
    rout=rout*cm_per_kpc
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
    kpc_per_pixel=cm_per_pixel/cm_per_kpc
    r_use=-1
    for i in range(tot):
        if r<r_list[i]:
            r_use=i
            break
    if r_use==-1:
        return numpy.NaN
#   len(r_list)=r:
#        vol=calc_v_ring(r_list[r-2],r_list[r-1])
#        return vol*cool_func[r-1]*den[r-1]*den[r-1]/ne_np_ratio
    if len(r_list)>=r_use:
        vol_total=0
        sb_tmp=0
        for i in range(tot)[r_use:tot]:
            vol=calc_v_ring(r_list[i],r_list[r_use-1])-calc_v_ring(\
                    r_list[i-1],r_list[r_use-1])-calc_v_ring(r_list[i],\
                    r_list[r_use])+calc_v_ring(r_list[i-1],r_list[r_use])
            vol_total=vol_total+vol
            sb_tmp=sb_tmp+den[i]*den[i]*cool_func[i]/ne_np_ratio*vol
        pixel_out=r_list[r_use]/kpc_per_pixel
        pixel_in=r_list[r_use-1]/kpc_per_pixel
        area=pi*(pixel_out*pixel_out-pixel_in*pixel_in)
        sb_tmp=sb_tmp/area
#        print(r_use,vol_total)
        return sb_tmp

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
