#!/usr/bin/env python3

import numpy as np
from scipy.optimize import minimize
import json
import sys
import re
from matplotlib import pyplot as plt
import analyze_db
import calc_csb
import calc_reff
from modnfw_readarray import mod_nfw,calc,lintp
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import z_at_value
import astropy.units as u
cosmo=FlatLambdaCDM(H0=71,Om0=0.27,Tcmb0=2.725)

def mz(t,a0,a1,a2):
    ttot=cosmo.age(0).value
    m=a0*np.power(t,a1-a2*np.log(t/ttot)) #10^14 Msun
    return m
def lm(x,t,a0,a2):
    z=z_at_value(cosmo.age,t*u.Gyr)
    ez=cosmo.efunc(z)
    return(a0+7/3*np.log10(ez)+a2*np.log10(x/1e14))
def calc_teff(m500,z):
    tnow=cosmo.age(z).value
    tfin=cosmo.age(3).value
    pm1=[0.06265,1.94,0.55]
    pm2=[0.05633,1.10,0.88]
    pm3=[0.01876,0.64,0.96]
    mref1=mz(tnow,*pm1)
    mref2=mz(tnow,*pm2)
    mref3=mz(tnow,*pm3)
    if m500>np.sqrt(mref1*mref2):
        pm_use=pm1
    elif m500>np.sqrt(mref2*mref3):
        pm_use=pm2
    else:
        pm_use=pm3
    mnorm=m500/mz(tnow,*pm_use)
    t_cp=[cosmo.age(0.0).value,cosmo.age(0.25).value,cosmo.age(0.5).value,cosmo.age(0.6).value,cosmo.age(0.8).value,cosmo.age(1.0).value,cosmo.age(1.5).value,cosmo.age(2.0).value]
    pl_tot=[[0.294,1.345],[0.332,1.320],[0.359,1.312],[0.348,1.380],[0.400,1.375],[0.388,1.526],[0.370,1.360],[0.434,1.509]]
    time_array=np.linspace(tnow,tfin,200)
    dt=time_array[0]-time_array[1]
    teff=0
    pl_tot=np.array(pl_tot)
    for i in range(len(time_array)):
        if time_array[i]>t_cp[1]:
            pl_use=lintp(time_array[i],t_cp[1],t_cp[0],pl_tot[1],pl_tot[0])
        elif time_array[i]>t_cp[2]:
            pl_use=lintp(time_array[i],t_cp[2],t_cp[1],pl_tot[2],pl_tot[1])
        elif time_array[i]>t_cp[3]:
            pl_use=lintp(time_array[i],t_cp[3],t_cp[2],pl_tot[3],pl_tot[2])
        elif time_array[i]>t_cp[4]:
            pl_use=lintp(time_array[i],t_cp[4],t_cp[3],pl_tot[4],pl_tot[3])
        elif time_array[i]>t_cp[5]:
            pl_use=lintp(time_array[i],t_cp[5],t_cp[4],pl_tot[5],pl_tot[4])
        elif time_array[i]>t_cp[6]:
            pl_use=lintp(time_array[i],t_cp[6],t_cp[5],pl_tot[6],pl_tot[5])
        elif time_array[i]>t_cp[7]:
            pl_use=lintp(time_array[i],t_cp[7],t_cp[6],pl_tot[7],pl_tot[6])
        else:
            pl_use=pl_tot[7]
        m_this=mnorm*mz(time_array[i],*pm_use)
        l_this=lm(m_this,time_array[i],*pl_use)
        l_this=np.power(10,l_this)
        if i==0:
            l_norm=l_this
        teff=teff+l_this/l_norm*dt
    return teff
# calculate the value (in kpc) of the r_delta(delta=200,etc.)
def critical_radius(od,p_Mnfw,z,t_total):
    Ez=(0.27*np.power(1+z,3)+0.73)**0.5
    pi=3.1415926
    H0=2.3e-18
    H=H0*Ez
    G=6.673e-8 #cm^3 g^-1 s^-2
    cri_density=3*H*H/8/pi/G
    kpc=3.086e21 #cm
    Msun=1.99e33 #g
    radius=np.nan
    for r in range(1,8000):
        mass=calc('n',r,p_Mnfw,t_total)
        area=4/3*pi*r*r*r
        den_this=mass/area
        den_this=den_this*Msun/kpc/kpc/kpc
        od_this=den_this/cri_density
        if od_this<=od:
            radius=r
            break
    return radius

def kmodel(r,a,b,c):
    k_this=a*np.power(r,b)+c
    return k_this

def nth_calc(r,a,b,gamma,r200):
    return a*(1+np.exp(-np.power(r/r200/b,gamma)))

def main():
    pi=3.1416
    kpc=3.086e19
    kpc_per_cm=3.086e21
    kev=1.6e-16
    G=6.67e-11
    mu=0.61
    mp=1.672e-27
    Gyr=365.24*3600*24*1e9 #s
    erg_to_kev=6.24e8
    Msun=1.99e30
    path=sys.argv[0].split('/')
    name=sys.argv[1]
    cfunc_file=sys.argv[2]
    dm=np.float(sys.argv[3])
    mode_sum=sys.argv[4]
    script_dir=''
    for i in range(len(path)-1):
        script_dir=script_dir+path[i]
        script_dir=script_dir+'/'
    ta1=np.load(script_dir+'lrs_ori.npy')
    ta2=np.load(script_dir+'lrs_dvr.npy')
    ta3=np.load(script_dir+'hrs_ori.npy')
    ta4=np.load(script_dir+'hrs_dvr.npy')
    for i in open('param_zzh_for_py.txt','r'):
        if re.match(r'z\s',i):
            z=float(i.split()[1])
    json_file=name+'_plt.json'
    t_total=[ta1,ta2,ta3,ta4]
    fi=open(json_file)
    dat=json.load(fi)
    r_array=np.array(dat['radius_model'],dtype=float)
#    r200=np.float(dat['r200'])
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
    param_file="p_all.npy"
    p=np.load(param_file)
    if mode_sum == "calc":
        SUM_array=analyze_db.main(t_total,name,True,'ktdcmnr',False)
        np.save('sum_array',SUM_array)
        sum_k_array=np.array(SUM_array[0])
        sum_T_array=np.array(SUM_array[1])
        sum_ne_array=np.array(SUM_array[2])
        sum_ne_cl_array=np.array(SUM_array[3])
        sum_mass_array=np.array(SUM_array[4])
        sum_mnfw_array=np.array(SUM_array[5])
        sum_sbp_array=np.array(SUM_array[6])
        sum_csb_array=np.array(SUM_array[7])
        sum_tproj_array=np.array(SUM_array[8])
    elif mode_sum == "read":
        sumarray_file='sum_array.npy'
        SUM_array=np.load(sumarray_file,allow_pickle=True)
        sum_k_array=np.array(list(SUM_array[0]),dtype=float)
        sum_T_array=np.array(list(SUM_array[1]),dtype=float)
        sum_ne_array=np.array(list(SUM_array[2]),dtype=float)
        sum_ne_cl_array=np.array(list(SUM_array[3]),dtype=float)
        sum_mass_array=np.array(list(SUM_array[4]),dtype=float)
        sum_mnfw_array=np.array(list(SUM_array[5]),dtype=float)
        sum_sbp_array=np.array(list(SUM_array[6]),dtype=float)
        sum_csb_array=np.array(list(SUM_array[7]),dtype=float)
        sum_tproj_array=np.array(list(SUM_array[8]),dtype=float)
    else:
        print("please check the mode_sum parameter")
        return -1
    sum_feedback_array=[]
    sum_dq_array=[]
    sum_EL_array=[]
    ind_50=int(len(p[0])*0.5)
    ind_84=int(len(p[0])*0.84)
    ind_16=int(len(p[0])*0.16)
    sum_Efeed=[]
    sum_dq=[]
    sum_r200_array=[]
    sum_r500_array=[]
    sum_m200_array=[]
    sum_m500_array=[]
    sum_lx200_array=[]
    sum_lx500_array=[]
    sum_gm200_array=[]
    sum_gm500_array=[]
    sum_fg500_array=[]
    sum_fg200_array=[]
    sum_k200_array=[]
    sum_Tave_array=[]
    sum_num_tot=[]
    sum_m3000k_array=[]
    sum_gm3000k_array=[]
    sum_fg3000k_array=[]
    sum_csb_array=[]
    sum_tcool_array=[]
    for i in range(len(p[0])):
        print(name,i)
        lx_array=[]
        EL_array=[]
        a0=p[2][i]
        k0=p[4][i]
        gamma0=p[3][i]
        n2=p[0][i]
        sbp_c=p[9][i]
        kobs=sum_k_array[i]
        sbp_array=sum_sbp_array[i]
        kth=a0*np.power(r_array,gamma0)+k0
        T_array=sum_T_array[i]
        dq_array=n2*T_array*(kobs-kth)/kobs
        ne_array=sum_ne_array[i]
        ne_cl_array=sum_ne_cl_array[i]
        p_mnfw=[p[6][i],p[1][i],p[10][i],p[11][i]]
        r500=critical_radius(500,p_mnfw,z,t_total)
        r200=critical_radius(200,p_mnfw,z,t_total)
        sum_r200_array.append(r200)
        sum_r500_array.append(r500)
        m200=calc('n',r200,p_mnfw,t_total)
        m500=calc('n',r500,p_mnfw,t_total)
        m3000k=calc('n',3000,p_mnfw,t_total)
        sum_m200_array.append(m200)
        sum_m500_array.append(m500)
        sum_m3000k_array.append(m3000k)
        flag_r200=0
        flag_r500=0
        lx200=-1
        lx500=-1
        teff=calc_teff(m500,z)
        for j in range(len(r_array)):
            if j == 0:
                V_this=4/3*pi*r_array[0]**3*kpc_per_cm**3
            else:
                V_this=4/3*pi*(r_array[j]**3-r_array[j-1]**3)*kpc_per_cm**3
            lx_this=ne_cl_array[j]**2/1.2*V_this*cfunc_lx_use_array[j]*4*pi*dm*dm
            E_loss=lx_this*teff*Gyr*erg_to_kev/(ne_array[j]*V_this*1.93)
            lx_array.append(lx_this)
            EL_array.append(E_loss)
            if flag_r500==0:
                if r_array[j]>r500:
                    lx500=np.array(lx_array).sum()
                    flag_r500=1
            if flag_r200==0:
                if r_array[j]>r200:
                    lx200=np.array(lx_array).sum()
                    flag_r200=1
        Efeedback_array=EL_array+dq_array
        sum_feedback_array.append(Efeedback_array)
        sum_dq_array.append(dq_array)
        sum_EL_array.append(EL_array)
        sum_lx200_array.append(lx200)
        sum_lx500_array.append(lx500)
        Efeed_tot=0
        dq_tot=0
        num_tot=0
        for j in range(len(r_array)):
            if r_array[j]<= r500 and r_array[j]>=0.00*r200:
                if j==0:
                    V_this=4/3*pi*r_array[0]**3*kpc_per_cm**3
                else:
                    V_this=4/3*pi*(r_array[j]**3-r_array[j-1]**3)*kpc_per_cm**3
                Efeed_tot=Efeed_tot+Efeedback_array[j]*ne_array[j]*V_this*1.93
                dq_tot=dq_tot+dq_array[j]*ne_array[j]*V_this*1.93
                num_tot=num_tot+ne_array[j]*V_this*1.93
        sum_Efeed.append(Efeed_tot)
        sum_dq.append(dq_tot)
        sum_num_tot.append(num_tot)
        flag_r200=0
        flag_r500=0
        gmas=0
        for j in range(len(r_array)):
            if j==0:
                V_this=4/3*pi*r_array[0]**3*kpc_per_cm**3
            else:
                V_this=4/3*pi*(r_array[j]**3-r_array[j-1]**3)*kpc_per_cm**3
            gmas=gmas+ne_array[j]*1.93*0.61*V_this*mp/Msun
            if flag_r500==0:
                if r_array[j]>r500:
                    gmas500=gmas
                    flag_r500=1
            if flag_r200==0:
                if r_array[j]>r200:
                    gmas200=gmas
                    flag_r200=1
        sum_gm3000k_array.append(gmas)
        fg200=gmas200/m200
        fg500=gmas500/m500
        fg3000k=gmas/m3000k
        sum_fg3000k_array.append(fg3000k)
        sum_gm200_array.append(gmas200)
        sum_gm500_array.append(gmas500)
        sum_fg200_array.append(fg200)
        sum_fg500_array.append(fg500)
        T200=G*m200*Msun*mu*mp/2/r200/kpc/kev
        Ez=(0.27*(1+z)**3+0.73)**0.5
        k200=362*T200*Ez**(-4/3)*(0.27/0.3)**(-4/3)
        sum_k200_array.append(k200)
        flag1=0
        flag2=0
        count=0
        t_tmp=0
        T_array=sum_T_array[i]
        for j in range(len(r_array)):
            if flag1==0:
                if r_array[j]>=0.2*r500:
                    flag1=1
            elif flag1==1:
                if flag2==0:
                    if r_array[j]<=0.5*r500:
                        t_tmp=t_tmp+T_array[j]
                        count=count+1
                    else:
                        flag2=1
        Tave=t_tmp/count
        sum_Tave_array.append(Tave)
        csb,tcool=calc_csb.main(name=name,iput='yes',t_array=T_array,sbp_array=sbp_array,ne_array=ne_array,sbp_c=sbp_c,r500=r500)
        sum_csb_array.append(csb)
        sum_tcool_array.append(tcool)
    out_file=name+'_suminfo.txt'
    fi=open(out_file,'w')
    sum_Efeed=np.sort(sum_Efeed)
    sum_dq=np.sort(sum_dq)
    sum_num_tot=np.sort(sum_num_tot)
    sum_r200_array=np.sort(sum_r200_array)
    sum_r500_array=np.sort(sum_r500_array)
    sum_m200_array=np.sort(sum_m200_array)
    sum_m500_array=np.sort(sum_m500_array)
    sum_lx200_array=np.sort(sum_lx200_array)
    sum_lx500_array=np.sort(sum_lx500_array)
    sum_gm200_array=np.sort(sum_gm200_array)
    sum_gm500_array=np.sort(sum_gm500_array)
    sum_fg200_array=np.sort(sum_fg200_array)
    sum_fg500_array=np.sort(sum_fg500_array)
    sum_k200_array=np.sort(sum_k200_array)
    sum_Tave_array=np.sort(sum_Tave_array)
    sum_m3000k_array=np.sort(sum_m3000k_array)
    sum_gm3000k_array=np.sort(sum_gm3000k_array)
    sum_fg3000k_array=np.sort(sum_fg3000k_array)
    sum_csb_array=np.sort(sum_csb_array)
    sum_tcool_array=np.sort(sum_tcool_array)
    reff=calc_reff.main(name,SUM_array)
    print('r200:',sum_r200_array[ind_50],sum_r200_array[ind_16],sum_r200_array[ind_84],file=fi)
    print('r500:',sum_r500_array[ind_50],sum_r500_array[ind_16],sum_r500_array[ind_84],file=fi)
    print('m200:',sum_m200_array[ind_50],sum_m200_array[ind_16],sum_m200_array[ind_84],file=fi)
    print('m500:',sum_m500_array[ind_50],sum_m500_array[ind_16],sum_m500_array[ind_84],file=fi)
    print('lx200:',sum_lx200_array[ind_50],sum_lx200_array[ind_16],sum_lx200_array[ind_84],file=fi)
    print('lx500:',sum_lx500_array[ind_50],sum_lx500_array[ind_16],sum_lx500_array[ind_84],file=fi)
    print('gm200:',sum_gm200_array[ind_50],sum_gm200_array[ind_16],sum_gm200_array[ind_84],file=fi)
    print('gm500:',sum_gm500_array[ind_50],sum_gm500_array[ind_16],sum_gm500_array[ind_84],file=fi)
    print('fg200:',sum_fg200_array[ind_50],sum_fg200_array[ind_16],sum_fg200_array[ind_84],file=fi)
    print('fg500:',sum_fg500_array[ind_50],sum_fg500_array[ind_16],sum_fg500_array[ind_84],file=fi)
    print('k200:',sum_k200_array[ind_50],sum_k200_array[ind_16],sum_k200_array[ind_84],file=fi)
    print('Tave(0.2-0.5r500):',sum_Tave_array[ind_50],sum_Tave_array[ind_16],sum_Tave_array[ind_84],file=fi)
    print('Efeed(0-1r500):',sum_Efeed[ind_50],sum_Efeed[ind_16],sum_Efeed[ind_84],file=fi)
    print('gnum(0-1r500):',sum_num_tot[ind_50],sum_num_tot[ind_16],sum_num_tot[ind_84],file=fi)
    print('dq(0-1r500):',sum_dq[ind_50],sum_dq[ind_16],sum_dq[ind_84],file=fi)
    print('m3000k:',sum_m3000k_array[ind_50],sum_m3000k_array[ind_16],sum_m3000k_array[ind_84],file=fi)
    print('gm3000k:',sum_gm3000k_array[ind_50],sum_gm3000k_array[ind_16],sum_gm3000k_array[ind_84],file=fi)
    print('fg3000k:',sum_fg3000k_array[ind_50],sum_fg3000k_array[ind_16],sum_fg3000k_array[ind_84],file=fi)
    print('tcool:', sum_tcool_array[ind_50],sum_tcool_array[ind_16],sum_tcool_array[ind_84],file=fi)
    print('csb:',sum_csb_array[ind_50],sum_csb_array[ind_16],sum_csb_array[ind_84],file=fi)
    print('reff:',reff,file=fi)
    fi.close()
    sum_feedback_array=np.array(sum_feedback_array)
    sum_feedback_array.sort(0)
    sum_EL_array=np.array(sum_EL_array)
    sum_EL_array.sort(0)
    sum_dq_array=np.array(sum_dq_array)
    sum_dq_array.sort(0)
    dq_array_center=sum_dq_array[ind_50]
    EL_array_center=sum_EL_array[ind_50]
    feedback_array_center=sum_feedback_array[ind_50]
    feedback_array_up=sum_feedback_array[ind_84]
    feedback_array_down=sum_feedback_array[ind_16]
    np.save('Efeedback',feedback_array_center)
    plt.clf()
    plt.plot(r_array,dq_array_center,label='extra energy change')
    plt.fill_between(r_array,sum_dq_array[ind_16],sum_dq_array[ind_84],alpha=0.3)
    plt.plot(r_array,EL_array_center,'k',label='radiation')
    plt.fill_between(r_array,sum_EL_array[ind_16],sum_EL_array[ind_84],color='k',alpha=0.3)
    plt.plot(r_array,feedback_array_center,label='total feedback')
    plt.fill_between(r_array,feedback_array_down,feedback_array_up,alpha=0.3)
    plt.xlabel('radius (kpc)')
    plt.ylabel('energy per particle (kev)')
    feedback_file=name+'_feedback.pdf'
    plt.xlim(1,2000)
    plt.ylim(-5,5)
    plt.legend()
    plt.savefig(feedback_file)
    Efeed_scaled_array=[]
    Efeed_scaled_array_up=[]
    Efeed_scaled_array_down=[]
    scaled_array=np.arange(0.001,1.4,0.001)
    for i in range(len(scaled_array)):
        for j in range(len(r_array)):
            if scaled_array[i]*r200<r_array[j]:
                Efeed=(feedback_array_center[j])
                Efeedup=feedback_array_up[j]
                Efeeddown=feedback_array_down[j]
                Efeed_scaled_array.append(Efeed)
                Efeed_scaled_array_up.append(Efeedup)
                Efeed_scaled_array_down.append(Efeeddown)
                break
            if j==len(r_array)-1:
                Efeed_scaled_array.append(feedback_array_center[j])
                Efeed_scaled_array_up.append(feedback_array_up[j])
                Efeed_scaled_array_down.append(feedback_array_down[j])
    Efeed_scaled_array_tot=np.array([Efeed_scaled_array,Efeed_scaled_array_down,Efeed_scaled_array_up])
    rnew_array=scaled_array*r200
    sum_dq_scaled_array=[]
    for i in range(len(sum_dq_array)):
        dq_scaled_array_this=calc_reff.reassign(rnew_array,r_array,sum_dq_array[i])
        sum_dq_scaled_array.append(dq_scaled_array_this)
    sum_dq_scaled_array=np.sort(sum_dq_scaled_array,0)
    dq_scaled_array_tot=np.array([sum_dq_scaled_array[ind_50],sum_dq_scaled_array[ind_16],sum_dq_scaled_array[ind_84]])
    np.save('Efeed_scaled_array',Efeed_scaled_array_tot)
    np.save('dq_scaled_array',dq_scaled_array_tot)


if __name__ == "__main__":
    main()

