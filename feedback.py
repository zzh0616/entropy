#!/usr/bin/env python3

import numpy as np
from scipy.optimize import minimize
import json
import sys
import re
from matplotlib import pyplot as plt
##calcuate the radius before non-gravitational heating
def radius_calc(r_array,m_array,kmodel_array,fb,ne_array,nth_array,plot_flag=0):
    mp=1.627e-27 #kg
    mu=0.6
    mue=1.16
    pi=3.1415926
    msun=2e30 #kg
    kpc_per_cm=3.086e21
    G=6.67e-11 #SI unit
    kev_to_j=1.6e-16
    pgold_array=np.copy(r_array)   #initialize, to be calculated
    #calculate the pressure profile if only gravity is considered using kmodel_array
    for i in range(len(r_array)):
        if i==0:
            dr=r_array[-1]-r_array[-2]
            ngout=(m_array[-1]-m_array[-2])*fb/(4/3*pi*(r_array[-1]**3-r_array[-2]**3))/mu/mp*msun/kpc_per_cm**3  #cm ^-3
#            print(ngout,ne_array[-1]*1.93)
            Pgout=ngout**(5/3)*kmodel_array[-1]/(1.93)**(2/3)/nth_array[-1] #kev*cm^-3
            pgold_array[-1]=Pgout
        elif i==len(r_array)-1:
            rindex=len(r_array)-i-1
            p_last=pgold_array[rindex+1]
            dr=(r_array[rindex+1]-r_array[rindex])/1 #kpc
            dptot=-dr*kpc_per_cm/100*(mp*mue**0.4*mu**0.6*G*m_array[rindex]*msun/(r_array[rindex]*(kpc_per_cm/100))**2)/kev_to_j*(pgold_array[rindex+1]*nth_array[rindex+1]**1/kmodel_array[rindex+1])**0.6 # kev*cm^-3
#            dp=dptot*nth_array[rindex+1]+pgold_array[rindex+1]*(nth_array[rindex+1]-nth_array[rindex])/nth_array[rindex+1]
            pthis=p_last-dptot
            pgold_array[rindex]=pthis
        else:
            rindex=len(r_array)-i-1
            p_last=pgold_array[rindex+1]
            dr=(r_array[rindex+1]-r_array[rindex])/1 #kpc
            dptot=-dr*kpc_per_cm/100*(mp*mue**0.4*mu**0.6*G*m_array[rindex]*msun/(r_array[rindex]*(kpc_per_cm/100))**2)/kev_to_j*(pgold_array[rindex+1]*nth_array[rindex+1]**1/kmodel_array[rindex+1])**0.6 # kev*cm^-3
#            dp=dptot*nth_array[rindex+1]+pgold_array[rindex+1]*(nth_array[rindex+1]-nth_array[rindex])/nth_array[rindex+1]
#            rindex=rindex-1
#            p_last=p_last-dptot
#            dr=(r_array[rindex+1]-r_array[rindex])/1 #kpc
#            dptot=0.5*(dptot-dr*kpc_per_cm/100*(mp*mue**0.4*mu**0.6*G*m_array[rindex]*msun/(r_array[rindex]*(kpc_per_cm/100))**2)/kev_to_j*(p_last*nth_array[rindex+1]**1/kmodel_array[rindex+1])**0.6) # kev*cm^-3
#            dp=dptot*nth_array[rindex+1]+pgold_array[rindex+1]*(nth_array[rindex+1]-nth_array[rindex])/nth_array[rindex+1]
#            rindex=rindex+1
            pthis=pgold_array[rindex+1]-dptot*1.0
            pgold_array[rindex]=pthis
    pgold_array=pgold_array*nth_array
    ne_old_array=(pgold_array/kmodel_array)**0.6*(mue/mu)**-0.6
    if plot_flag ==1 :
        plt.loglog(r_array,pgold_array,'k')
        plt.loglog(r_array,p_ref,'b')
        plt.show()
#    ne_old_array=np.copy(ne_array)
        plt.plot(r_array,ne_old_array/ne_array,'k')
        plt.loglog(r_array,ne_array,'b')
        plt.show()
    rold_array=[]
    rold_c=0
    mg_old_array=[]
    mg_new_array=[]
    count=0
    count_old=1
    mg_old=0
    iold=0
    flag=0
    mgnew_this=0
    rold_this=0
    #calculate the radius correspond to r_array
    for i in r_array:
        mgnew_this=mgnew_this+4*pi/3*(i**3-iold**3)*ne_array[count]*mue*mp*kpc_per_cm**3/msun
        mg_new_array.append(mgnew_this)
        count=count+1
        iold=i
        if flag==1:
            rold_array.append(np.NAN)
            mg_old_array.append(np.NAN)
            continue
#        print(mg_new_array)
        while mg_old<=mgnew_this:
            rold_this=rold_this+1
            if rold_this>=r_array[-1]:
                flag=1
                break
            if rold_this>=r_array[count_old]:
                count_old=count_old+1
            ne_old=ne_old_array[count_old]+(rold_this-r_array[count_old])*(ne_old_array[count_old]-ne_old_array[count_old-1])/(r_array[count_old]-r_array[count_old-1])
            mg_old=mg_old+4/3*pi*(rold_this**3-(rold_this-1)**3)*ne_old*mue*mp*kpc_per_cm**3/msun
        mg_old_array.append(mg_old)
        rold_array.append(rold_this)
    rold_array=np.array(rold_array)
    if plot_flag == 1:
        plt.plot(rold_array,mg_old_array,'k')
        plt.plot(r_array,mg_new_array,'b')
        plt.show()
    return rold_array

def kmodel(r,a,b,c):
    k_this=a*np.power(r,b)+c
    return k_this

def nth_calc(r,a,b,gamma,r200):
    return a*(1+np.exp(-np.power(r/r200/b,gamma)))

def kref(r,a,b,c,k0,r200):
    return a*np.power(r,b)*np.exp(-r/c/r200)+k0


def test():
    fi=open(sys.argv[1])
    dat=json.load(fi)
    for i in open(sys.argv[2]):
        if re.match(r'^n2,',i):
            n2=float(i.split(',')[1])
        if re.match(r'^nth_a,',i):
            nth_a=float(i.split(',')[1])
        if re.match(r'^nth_b,',i):
            nth_b=float(i.split(',')[1])
        if re.match(r'^nth_gamma,',i):
            nth_gamma=float(i.split(',')[1])
        if re.match(r'^a0,',i):
            a0=float(i.split(',')[1])
        if re.match(r'^gamma0,',i):
            gamma0=float(i.split(',')[1])
        if re.match(r'^k0,',i):
            k0=float(i.split(',')[1])
    global r_array,m_array,fb,ne_array,nth_array,kth_array
    r_array=np.array(dat['radius_model'],dtype=float)
    m_array=np.array(dat['m_fit'][0],dtype=float)
    t_array=np.array(dat['temperature_model'][0],dtype=float)
    ne_array=np.array(dat['den_model'][0],dtype=float)
    global p_ref,k_obs
    p_ref=t_array*ne_array*1.93
    k_obs=np.array(dat['k_fit'][0],dtype=float)
    k_obs=t_array*ne_array**(-2/3)
    kth_array=kmodel(r_array,a0,gamma0,k0)
    fb=0.115
    r200=float(sys.argv[3])
    nth_array=nth_calc(r_array,nth_a,nth_b,nth_gamma,r200)
    radius_calc(r_array,m_array,k_obs,fb,ne_array,nth_array)

def main():
    fi=open(sys.argv[1])
    dat=json.load(fi)
    for i in open(sys.argv[2]):
        if re.match(r'^n2,',i):
            n2=float(i.split(',')[1])
        if re.match(r'^nth_a,',i):
            nth_a=float(i.split(',')[1])
        if re.match(r'^nth_b,',i):
            nth_b=float(i.split(',')[1])
        if re.match(r'^nth_gamma,',i):
            nth_gamma=float(i.split(',')[1])
        if re.match(r'^a0,',i):
            a0=float(i.split(',')[1])
        if re.match(r'^gamma0,',i):
            gamma0=float(i.split(',')[1])
        if re.match(r'^k0,',i):
            k0=float(i.split(',')[1])
    global r_array,m_array,fb,nth_array,kth_array,ne_array
    r_array=np.array(dat['radius_model'],dtype=float)
    m_array=np.array(dat['m_fit'][0],dtype=float)
    r200=np.float(dat['r200'])
#    a0=1.6
#    gamma0=0.97
#    k0=100
    kth_array=kmodel(r_array,a0,gamma0,k0)
    fb=0.16
    ne_array=np.array(dat['den_model'][0],dtype=float)
    t_array=np.array(dat['temperature_model'][0],dtype=float)
#    r200=float(sys.argv[3])
    nth_array=nth_calc(r_array,nth_a,nth_b,nth_gamma,r200)
    global p_ref,k_obs
    p_ref=t_array*ne_array*1.93
    p=[a0,gamma0,k0]
    f = lambda p1: poss(*p1)
    result=minimize(f,p,method='Powell')
#    print(result.x)
#    print(poss(*p),poss(*result.x),len(r_array))
    a0f,gamma0f,k0f=np.abs(result.x)
  #  a1,gamma1,k1=[1.3,1.025,220]
  #  p=[a1,gamma1,k1]
    kmodel_array=kmodel(r_array,a0f,gamma0f,k0f)
    r_old_array=radius_calc(r_array,m_array,kmodel_array,fb,ne_array,nth_array,0)
    kold_array=kmodel(r_old_array,a0f,gamma0f,k0f)
    kcompare_array=kmodel(r_array,a0,gamma0,k0)
#    plt.plot(kold_array/kcompare_array)
#    plt.show()
    kobs_array=np.array(dat['k_fit'][0],dtype=float)
    count=0
    flag=0
    dq_array=[]
    Eg_array=[]
    pi=3.1415926
    Msun=2e30
    kpc=3.086e19
    kpc_per_cm=3.086e21
    kev=1.6e-16
    G=6.67e-11
    mu=0.6
    mp=1.027e-27
    for i in range(len(r_array)):
        kold=kold_array[i]
        kobs=kobs_array[i]
        if flag==1:
            dq_array.append(np.NAN)
            Eg_array.append(np.NAN)
            continue
        while r_old_array[count]<=r_array[i] or r_old_array[count]==np.NAN:
            if count>=len(r_old_array)-1:
                flag=1
                break
            count=count+1
        Mold=m_array[count]
        Mnow=m_array[i]
        dq=n2*1.5*(kobs-kold)/kobs
        Eg=G*mu*mp*(-Mnow/r_array[i]+Mold/r_old_array[i])*Msun/kpc/kev
        dq_array.append(dq)
        Eg_array.append(Eg)
    dq_array=np.array(dq_array)
    Eg_array=np.array(Eg_array)
#    plt.plot(r_array,r_old_array)
#    plt.show()
    plt.plot(r_array,dq_array,'r',label='extra heating(dQ)')
    plt.plot(r_array,Eg_array,'g',label='gravity')
#    plt.savefig('E_ICM.pdf')
    #calculate Lx
    cfunc_file=sys.argv[3]
    dl=np.float(sys.argv[4])
    r_cfunc_array=[]
    cfunc_ori_array=[]
    cfunc_lx_use_array=[]
    for i in open(cfunc_file,'r'):
        r_cfunc_array.append(float(i.split()[0]))
        cfunc_ori_array.append(float(i.split()[1]))

    for i in r_array:
        i=float(i)
        if i==0:
            continue
        for j in range(len(r_cfunc_array)):
            if r_cfunc_array[j]>i:
                cfunc_lx_use_array.append(cfunc_ori_array[j])
                break
    cfunc_lx_use_array=np.array(cfunc_lx_use_array)
    lx_array=[]
    EL_array=[]
    Gyr=365.24*3600*24*1e9 #s
    erg_to_kev=6.24e8
    for i in range(len(r_array)):
        if i == 0:
            V_this=4/3*pi*r_array[0]**3*kpc_per_cm**3
        else:
            V_this=4/3*pi*(r_array[i]**3-r_array[i-1]**3)*kpc_per_cm**3
        lx_this=ne_array[i]**2/1.2*V_this*cfunc_lx_use_array[i]*4*pi*dl*dl
        E_loss=lx_this*5*Gyr*erg_to_kev/(ne_array[i]*V_this*1.93)
        lx_array.append(lx_this)
        EL_array.append(E_loss)
    lx_array=np.array(lx_array)
#    print(lx_array.sum())
    EL_array=np.array(EL_array)
    Efeed_array=EL_array+Eg_array+dq_array
    plt.plot(r_array,EL_array,'b',label='radiation')
    plt.plot(r_array,Efeed_array,'k',label='total feedback')
    plt.xlabel('radius (kpc)')
    plt.ylabel('energy per particle (kev)')
    name=sys.argv[1][0:-9]+'_feedback.pdf'
    plt.xlim(1,2000)
    plt.legend()
    plt.savefig(name)
    Efeed_tot=0
    for i in range(len(r_array)):
        if r_array[i]<= 50:
            if i==0:
                V_this=4/3*pi*r_array[0]**3*kpc_per_cm**3
            else:
                V_this=4/3*pi*(r_array[i]**3-r_array[i-1]**3)*kpc_per_cm**3
            Efeed_tot=Efeed_tot+Efeed_array[i]*ne_array[i]*V_this*1.93
    print(Efeed_tot)
    return 0

def poss(a0,gamma0,k0):
#   global r_array
    k0=np.abs(k0)
    a0=np.abs(a0)
    gamma0=np.abs(gamma0)
    kmodel_array=kmodel(r_array,a0,gamma0,k0)
    r_old_array=radius_calc(r_array,m_array,kmodel_array,fb,ne_array,nth_array)
    count=0
    rnow=0
    like=0
    flag=0
    kmodel_array=kmodel(r_old_array,a0,gamma0,k0)
    for i in range(600):
        if kmodel_array[i] == kmodel_array[i]:
            like=like+(kmodel_array[i]-kth_array[i])**2/kth_array[i]/kth_array[i]
#    like=np.sum((kth_array-kmodel_array)**2)
#    for i in range(len(r_old_array)):
#        if flag==1:
#            break
#        while r_array[count]<=r_old_array[i]:
#            count=count+1
#            if count>=len(r_array)-1 or r_old_array[i]==np.NAN:
#                flag=1
#                break
#        kold=kmodel_array[i]
#        kth=kth_array[i]
#        like=like+(kth-kold)*(kth-kold)/kth/kth*10*10
#    like=like+(gamma0-1.1)*(gamma0-1.1)*2000
#    print(like)
    return like
if __name__ == "__main__":
    main()

