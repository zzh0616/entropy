#!/usr/bin/env python3

import numpy as np
import numpy
from modnfw_readarray import mod_nfw,calc
import sys
import re
import json
import matplotlib.pyplot as plt
pi=3.1415926
mp=1.67e-27  #kg
kev=1.6e-16 #J
kpc=3.086e21 #cm
Msun=2e30 #kg

def critical_radius(p_od,p_Mnfw,z,t_total):
    od=p_od[0]
    od_type=p_od[1]
    Ez=(0.27*np.power(1+z,3)+0.73)**0.5
    H0=2.3e-18
    H=H0*Ez
    G=6.673e-8 #cm^3 g^-1 s^-2
    cri_density=3*H*H/8/pi/G
    kpc=3.086e21 #cm
    Msun=1.99e33 #g
    radius=np.nan
    for r in range(1,8000):
        if od_type=='m':
            mass=calc('n',r,p_Mnfw,t_total)
            area=4/3*pi*r*r*r
            den_this=mass/area
        elif od_type=='c':
            print('WARNING: **NOT** A CORRECT WAY TO CALCULATE ')
            den_this=mod_nfw(r,p_Mnfw)/4/pi/r/r
        den_this=den_this*Msun/kpc/kpc/kpc
        od_this=den_this/cri_density
        if od_this<=od:
            radius=r
            break
    return radius

def cal_mass(r,p_mnfw,t_total):
    mass=calc('n',r,p_mnfw,t_total)
    return mass

def main():
    if len(sys.argv)==1:
        print('useage: calc_mass.py <"cfg/input"> <does not matter for cfg/parameters for input> <param_zzh_for_py/z> <json_plt_file>')
        return 1
    if sys.argv[1]=='input':
        p_mnfw=np.array(sys.argv[2].split(','),dtype=float)
        z=float(sys.argv[3])
    elif sys.argv[1]=='cfg':
        p_all=np.load('p_all.npy')
        rho_all=p_all[6]
        rs_all=p_all[1]
        delta_all=p_all[10]
        delta2_all=p_all[11]
        try :
            print(sys.argv[3])
        except NameError:
            param_file='param_zzh_for_py.txt'
        except IndexError:
            param_file='param_zzh_for_py.txt'
        else:
            param_file=sys.argv[3]
        for i in open(param_file,'r'):
            if re.match(r'^z\s',i):
                z=float(i.split()[1])
    script_dir=''
    tmp=sys.argv[0].split('/')
    for i in range(len(tmp)-1):
        script_dir=script_dir+tmp[i]
        script_dir=script_dir+'/'
#    print(script_dir)
    ta1=numpy.load(script_dir+'lrs_ori.npy')
    ta2=numpy.load(script_dir+'lrs_dvr.npy')
    ta3=numpy.load(script_dir+'hrs_ori.npy')
    ta4=numpy.load(script_dir+'hrs_dvr.npy')
    t_total=[ta1,ta2,ta3,ta4]
    fi=open(sys.argv[4],'r')
    dat=json.load(fi)
    den_array=numpy.array(dat['den_model'][0])
    den_up=numpy.array(dat['den_model'][2])
    den_down=numpy.array(dat['den_model'][1])
    rne_array=numpy.array(dat['radius_model'],dtype=float)
    r200_array=[]
    r500_array=[]
    m500_array=[]
    m200_array=[]
    m3000k_array=[]
    SUM_den_tot=[]
    for i in range(len(rho_all)):
        p_mnfw=[rho_all[i],rs_all[i],delta_all[i],delta2_all[i]]
        den_tot=mod_nfw(rne_array,p_mnfw)/4/pi/rne_array/rne_array
        SUM_den_tot.append(den_tot)
        r200m=critical_radius([200,'m'],p_mnfw,z,t_total)
        r500m=critical_radius([500,'m'],p_mnfw,z,t_total)
        m200m=cal_mass(r200m,p_mnfw,t_total)
        m500m=cal_mass(r500m,p_mnfw,t_total)
        m3000k=cal_mass(3000,p_mnfw,t_total)-cal_mass(2995,p_mnfw,t_total)
        r200_array.append(r200m)
        r500_array.append(r500m)
        m200_array.append(m200m)
        m500_array.append(m500m)
        m3000k_array.append(m3000k)
    id_up=int(len(rho_all)*0.83)
    id_down=int(len(rho_all)*0.17)
    id_mid=int(len(rho_all)*0.5)
    r200_array=numpy.array(r200_array)
    r500_array=numpy.array(r500_array)
    m200_array=numpy.array(m200_array)
    m500_array=numpy.array(m500_array)
    m3000k_array=numpy.array(m3000k_array)
    r200_array.sort()
    r200_c=r200_array[id_mid]
    r200_u=r200_array[id_up]
    r200_d=r200_array[id_down]
    print(r200_c,r200_d,r200_u)
    m200_array.sort()
    m200_c=m200_array[id_mid]
    m200_u=m200_array[id_up]
    m200_d=m200_array[id_down]
    m3000k_array.sort()
    print(m200_c,m200_d,m200_u)
    print(m3000k_array[id_mid])
    for i in range(len(rne_array)):
        r=rne_array[i]
        if r>=r200_c:
            rne_id=i
            break
    gmas=0
    gmas_up=0
    gmas_low=0
    gmas_array=[]
    gmas_uparray=[]
    gmas_lowarray=[]
#    den_array=den_up
#    print(den_array)
    for i in range(len(rne_array)):
        r=rne_array[i]
        if i==0:
            gmas=gmas+4*pi*r*r*rne_array[i]*den_array[i]*1.93*0.61*mp/Msun*kpc*kpc*kpc
            gmas_up=gmas_up+4*pi*r*r*rne_array[i]*den_up[i]*1.93*0.61*mp/Msun*kpc*kpc*kpc
            gmas_low=gmas_low+4*pi*r*r*rne_array[i]*den_down[i]*1.93*0.61*mp/Msun*kpc*kpc*kpc
        else:
            gmas=gmas+4*pi*r*r*(rne_array[i]-rne_array[i-1])*den_array[i]*1.93*0.61*mp/Msun*kpc*kpc*kpc
            gmas_up=gmas_up+4*pi*r*r*(rne_array[i]-rne_array[i-1])*den_up[i]*1.93*0.61*mp/Msun*kpc*kpc*kpc
            gmas_low=gmas_low+4*pi*r*r*(rne_array[i]-rne_array[i-1])*den_down[i]*1.93*0.61*mp/Msun*kpc*kpc*kpc
#        print(r,gmas)
        if i==rne_id+1:
            gmas_200=gmas
            gmas200_low=gmas_low
            gmas200_up=gmas_up
        gmas_array.append(gmas)
        gmas_uparray.append(gmas_up)
        gmas_lowarray.append(gmas_low)
    print(gmas_200,gmas200_low,gmas200_up)
    print(gmas_array[-1]-gmas_array[-2],rne_array[-1]-rne_array[-2])
    SUM_den_tot=np.array(SUM_den_tot)
    SUM_den_tot.sort(0)
    den_tot=SUM_den_tot[id_mid]
    den_totlow=SUM_den_tot[id_down]
    den_totup=SUM_den_tot[id_up]
    plt.loglog(rne_array,gmas_array,color='b')
    plt.fill_between(rne_array,gmas_lowarray,gmas_uparray,color='grey')
    plt.xlabel('Radius (kpc)')
    plt.xlim(30,2000)
    plt.ylim(1e10,2e14)
    plt.ylabel('Gas Mass ($M_{\odot}$)')
    plt.savefig('gmas.pdf')
    plt.clf()
    plt.loglog(rne_array,den_tot,color='b')
    plt.fill_between(rne_array,den_totlow,den_totup,color='grey')
    plt.xlim(30,2000)
    plt.xlabel('Radius (kpc)')
    plt.ylabel('Total Density ($M_{\odot}/kpc^{3}$)')
    plt.savefig('tden.pdf')
    return 0

if __name__=='__main__':
    main()
