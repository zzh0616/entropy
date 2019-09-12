#!/usr/bin/env python3
## $1 temperature data; radius in kpc
## $2 temperature parameter file
## $3 sbp data; 'radius err sbp err' radius in kpc
## $4 cooling function file , radius in kpc, radius should match density
## $5 outfile name for fitted parameter, if exist
## $6 outfile name for fitted temperature
## $7 outfile name for fitted surface brightness
## $8 outfile name for fitted entropy
####
#tempory var define area
cm_per_pixel=1.80e21
####
from scipy.optimize import minimize
from scipy.optimize import fsolve
from scipy.integrate import odeint
from scipy.integrate import ode
from zzh_model import gpob
import zzh_model
import entropy_model
import deproject_model
import matplotlib.pyplot as plt
import sys
import math
import numpy
import scipy
import matplotlib
import re
import numpy as np
import types
global T0_0,T0_MIN,T0_MAX,N1_0,N1_MIN,N1_MAX,N2_0,N2_ERR
global RS_0,RS_ERR,A0_0,A0_ERR,GAMMA0_0,GAMMA0_ERR,K0_0,K0_MIN,K0_MAX
global A1_0,A1_ERR,GAMMA1_0,GAMMA1_ERR,K1_0,K1_ERR
global N3_0,N3_MIN,N3_MAX,RHO_0,RHO_ERR,Z_0,Z_MIN,Z_MAX
global A_0,A_ERR,BETA_0,BETA_MIN,BETA_MAX,NT_0,NT_ERR
global NM_0,NM_MIN,NM_MAX,X_0,X_MIN,X_MAX
global r_array
global re_array
global t_array
global te_array
global ORI_POST
global P0
global counter_tm
#global ne_global_array
numpy.set_printoptions(linewidth=1000)
pi=3.1415926
## y=[T,ne]
def dydr(y,r0,p):
    G=6.67e-11 #m^3 kg^-1 s^-2
    mp=1.67e-27  #kg
    kev=1.6e-16 #J
    kpc=3.086e19 #m
    Msun=2e30 #kg
    mu=0.61
    rs=p[2]
    rho=p[7]
    c2=p[1]
    c3=p[6]
    c4=1-c3-c2
    a0=p[3]
    gamma0=p[4]
    k0=p[5]
    c1=p[0]
    k=a0*numpy.power(r0,gamma0)+k0
    dkdr=a0*gamma0*numpy.power(r0,gamma0-1)
    z1=-4*pi*rho*rs*rs*rs*(numpy.log((r0+rs)/rs)-r0/(rs+r0))*G*mu*mp/r0/r0\
            *Msun/kpc/kev   #unit in kev/kpc same as z2
    z2=(1-c3)*c1*rho*rs*rs*rs*(-numpy.log((rs+r0)/rs)/r0/r0+1/(r0*(rs+r0)))
    ne=y[1]
    T=y[0]
#    dTdr=z1-(z2-c4*z1-c2*dkdr*numpy.power(ne,2/3))/(2/3*k*c2/T*numpy.power(\
#            ne,2/3)-c4)
#    dnedr=(z2-c4*z1-c2*dkdr*numpy.power(ne,2/3))/(2*k*c2/3/numpy.power(\
#            ne,1/3)-c4*T/ne)
    tmp_n0=0.0371
    tmp_beta=0.5
    tmp_rc=32.7
    dnedr=10*(deproject_model.beta(r0+0.05,tmp_n0,tmp_beta,tmp_rc)-\
            deproject_model.beta(r0-0.05,tmp_n0,tmp_beta,tmp_rc))
    dTdr=z1-T/ne*dnedr
#    dTdr=1/c4*(z2-c2*k*2/3*numpy.power(ne,-1/3)*dnedr-c2*dkdr*numpy.power(ne,2/3))
    return [dTdr,dnedr]

#    T=entropy_model.main(r0,R500_0,R200_0,*p)
#    if T<=0:
#        return [0
#    left=4*pi*rho*math.pow(rs,3)*(math.log((rs+r0)/rs)-r0/(rs+r0))
#    dTdr=(zzh_model.main(r0+10,R500_0,R200_0,*p)-T)/10
#    return (-left*msun*G*mu*mp/(r0*kpc)/(r0*kpc)/(T*kev)*kpc-dTdr/T)*ne
counter_tm=0
def postpob(n1,n2,rs,a0,gamma0,k0,n3,rho,tc0,ne0):
    global counter_tm
#    global ne_global_array
    counter_tm=counter_tm+1
    print(counter_tm)
    z=Z_0
    a=0
    beta=BETA_0
    n1=N1_0
    if n3 < N3_MIN:
        n3=N3_MIN
    if n3>N3_MAX:
        n3=N3_MAX
#    if nm<NM_MIN:
#        nm=NM_MIN
#    if nm>NM_MAX:
#        nm=NM_MAX
    if k0>K0_MAX:
        k0=K0_MAX
    if k0<K0_MIN:
        k0=K0_MIN
    n2=abs(n2)
#    nt=abs(nt)
#    a=abs(a)
    rs=abs(rs)
    rho=abs(rho)
    ne0=abs(ne0)
#    t0=abs(t0)
    k0_err=K0_0*1/2
    prior=gpob(n2,N2_0,N2_ERR)+gpob(rs,RS_0,RS_ERR)+gpob(a0,A0_0,A0_ERR)+\
            gpob(gamma0,GAMMA0_0,GAMMA0_ERR)+gpob(rho,RHO_0,RHO_ERR)+\
            gpob(k0,K0_0,k0_err)+gpob(ne0,NE_0,NE_0_ERR)+gpob(tc0,TC,TC_ERR)
#gpob(a,A_0,A_ERR)+gpob(nt,NT_0,NT_ERR)+gpob(nm,0.2,0.3)+\
    likehood=0
    ne_array=[]
    T_array=[]
    rne_array=range(5,3005,10)
    tmp_array=range(0,3010,10)
    y0=[tc0,ne0]
    p=[n1,n2,rs,a0,gamma0,k0,n3,rho]
    y=odeint(dydr,y0,rne_array,args=(p,),Dfun=None,h0=1e-5,hmax=1,hmin=1e-6,\
            rtol=[0.01,0.01],atol=[1e-3,1e-6],mxstep=1000000)
#    ne_global_array=ne_array
    ne_array.append(0)
    for i in range(0,len(y)):
        ne_array.append(y[i][1])
        T_array.append(y[i][0])
#    print(ne_array)
#    T_this=
#    kobs=pow(ne_this,-2/3)*T_this
#    p=[t0,n1,n2,rs,a0,gamma0,k0,kobs,n3,z,rho,a,beta,nt,nm,x]
    for i in range(len(r_array)):
        r_this=numpy.int(r_array[i]/10-1)
        t_this=T_array[r_this]
        likehood=likehood+gpob(t_this,t_array[i],te_array[i])
    for i in range(len(rsbp_array)):
        sbp_this=deproject_model.calc_sb(rsbp_array[i],tmp_array,ne_array,cfunc_use_array,cm_per_pixel)
        print(sbp_this)
        sbp_this=abs(math.log(sbp_this))
        tmp_sbp=abs(math.log(sbp_array[i]))
        print(sbp_this,tmp_sbp)
        tmp_sbpe=abs(math.log(sbpe_array[i]))
#        sbp_this=abs(sbp_this)
#        tmp_sbp=abs(sbp_array[i])
#        tmp_sbpe=abs(sbpe_array[i])
        likehood=likehood+0.3*gpob(sbp_this,tmp_sbp,tmp_sbpe)
    return -(prior+likehood)

def adjust(p):
    if p[0]<T0_MIN:
        p[0]=T0_MIN
    if p[0]>T0_MAX:
        p[0]=T0_MAX
    p[1]=N1_0
    if p[6]<K0_MIN:
        p[6]=K0_MIN
    if p[6]>K0_MAX:
        p[6]=K0_MAX
    if p[10]<N3_MIN:
        p[10]=N3_MIN
    if p[10]>N3_MAX:
        p[10]=N3_MAX
    p[11]=Z_0
    p[14]=BETA_0
    if p[16]<NM_MIN:
        p[16]=NM_MIN
    if p[16]>NM_MAX:
        p[16]=NM_MAX
    if p[17]<X_MIN:
        p[17]=X_MIN
    if p[17]>X_MAX:
        p[17]=X_MAX
    p[2]=abs(p[2])
    p[9]=abs(p[9])
    p[7]=abs(p[7])
    p[15]=abs(p[15])
    p[13]=abs(p[13])
    p[18]=abs(p[18])
    p[3]=abs(p[3])
    p[12]=abs(p[12])
    p[13]=0
    return 0
def print_result(p):
    print('t0\t%.4e'%(p[0]))
    print('n1\t%.4e'%(p[1]))
    print('n2\t%.4e'%(p[2]))
    print('rs\t%.4e'%(p[3]))
    print('a0\t%.4e'%(p[4]))
    print('gamma0\t%.4e'%(p[5]))
    print('k0\t%.4e'%(p[6]))
    print('a1\t%.4e'%(p[7]))
    print('gamma1\t%.4e'%(p[8]))
    print('k1\t%.4e'%(p[9]))
    print('n3\t%.4e'%(p[10]))
    print('z\t%.4e'%(p[11]))
    print('rho\t%.4e'%(p[12]))
    print('a\t%.4e'%(p[13]))
    print('beta\t%.4e'%(p[14]))
    print('nt\t%.4e'%(p[15]))
    print('nm\t%.4e'%(p[16]))
    print('x\t%.4e'%(p[17]))
    print('ne0\t%.4e'%(p[18]))
    return 0

def plot_result(ori_param,fin_param):
    t_fit=[]
    t_ori=[]
    ori_param[18:100]=[]
    for i in range(2000):
        if i == 0:
            continue
        tmp=zzh_model.main(i,R500_0,R200_0,*fin_param)
        t_fit.append(tmp)
        tmp=zzh_model.main(i,R500_0,R200_0,*ori_param)
        t_ori.append(tmp)
    plt.plot(t_fit,'b')
#    plt.plot(t_ori,'r')
    plt.plot(r_array,t_array,'g*')
    plt.show()
    return 0
ORI_POST=1.0
r_array=[]
re_array=[]
t_array=[]
te_array=[]
rsbp_array=[]
rsbpe_array=[]
sbp_array=[]
sbpe_array=[]
for i in open(sys.argv[1]):
    r,rer,t,te=i.split()
    r=float(r)
    rer=float(rer)
    t=float(t)
    te=float(te)
    r_array.append(r)
    re_array.append(rer)
    t_array.append(t)
    te_array.append(te)


r_array=numpy.array(r_array)
re_array=numpy.array(re_array)
t_array=numpy.array(t_array)
te_array=numpy.array(te_array)

for i in open(sys.argv[3]):
    r,rer,sbp,sbpe=i.split()
    r=float(r)
    rer=float(rer)
    sbp=float(sbp)
    sbpe=float(sbpe)
    rsbp_array.append(r)
    rsbpe_array.append(rer)
    sbp_array.append(sbp)
    sbpe_array.append(sbpe)
rsbp_array=numpy.array(rsbp_array)
rsbpe_array=numpy.array(rsbpe_array)
sbp_array=numpy.array(sbp_array)
sbpe_array=numpy.array(sbpe_array)

for i in open(sys.argv[2],'r'):
    if re.match(r'^T0\s',i):
        T0_0=float(i.split()[1])
        T0_MIN=float(i.split()[2])
        T0_MAX=float(i.split()[3])
    if re.match(r'^N1\s',i):
        N1_0=float(i.split()[1])
        N1_MIN=float(i.split()[2])
        N1_MAX=float(i.split()[3])
    elif re.match(r'^N2\s',i):
        N2_0=float(i.split()[1])
        N2_ERR=float(i.split()[2])
    elif re.match(r'^rs\s',i):
        RS_0=float(i.split()[1])
        RS_ERR=float(i.split()[2])
    elif re.match(r'^a0\s',i):
        A0_0=float(i.split()[1])
        A0_ERR=float(i.split()[2])
    elif re.match(r'^gamma0\s',i):
        GAMMA0_0=float(i.split()[1])
        GAMMA0_ERR=float(i.split()[2])
    elif re.match(r'k0\s',i):
        K0_0=float(i.split()[1])
        K0_MIN=float(i.split()[2])
        K0_MAX=float(i.split()[3])
    elif re.match(r'^a1\s',i):
        A1_0=float(i.split()[1])
        A1_ERR=float(i.split()[2])
    elif re.match(r'^gamma1\s',i):
        GAMMA1_0=float(i.split()[1])
        GAMMA1_ERR=float(i.split()[2])
    elif re.match(r'^k1\s',i):
        K1_0=float(i.split()[1])
        K1_ERR=float(i.split()[2])
    elif re.match(r'^N3\s',i):
        N3_0=float(i.split()[1])
        N3_MIN=float(i.split()[2])
        N3_MAX=float(i.split()[3])
    elif re.match(r'^z\s',i):
        Z_0=float(i.split()[1])
        Z_MIN=float(i.split()[2])
        Z_MAX=float(i.split()[3])
    elif re.match(r'^rho\s',i):
        RHO_0=float(i.split()[1])
        RHO_ERR=float(i.split()[2])
    elif re.match(r'^a\s',i):
        A_0=float(i.split()[1])
        A_ERR=float(i.split()[2])
    elif re.match(r'^beta\s',i):
        BETA_0=float(i.split()[1])
        BETA_MIN=float(i.split()[2])
        BETA_MAX=float(i.split()[3])
    elif re.match(r'^nt\s',i):
        NT_0=float(i.split()[1])
        NT_ERR=float(i.split()[2])
    elif re.match(r'^nm\s',i):
        NM_0=float(i.split()[1])
        NM_MIN=float(i.split()[2])
        NM_MAX=float(i.split()[3])
    elif re.match(r'^R200\s',i):
        R200_0=float(i.split()[1])
        R200_MIN=float(i.split()[2])
        R200_MAX=float(i.split()[3])
    elif re.match(r'^R500\s',i):
        R500_0=float(i.split()[1])
        R500_MIN=float(i.split()[2])
        R500_MAX=float(i.split()[3])
    elif re.match(r'^x\s',i):
        X_0=float(i.split()[1])
        X_MIN=float(i.split()[2])
        X_MAX=float(i.split()[3])
    elif re.match(r'^alpha\s',i):
        alpha_0=float(i.split()[1])
        alpha_err=float(i.split()[2])
    elif re.match(r'^nec\s',i):
        NE_0=float(i.split()[1])
        NE_0_ERR=float(i.split()[2])
    elif re.match(r'^TC\s',i):
        TC=float(i.split()[1])
        TC_ERR=float(i.split()[2])
cfunc_ori_array=[]
r_cfunc_array=[]
cfunc_use_array=[]
cfunc_ori_array.append(0)
r_cfunc_array.append(0)
cfunc_use_array.append(0)
for i in open(sys.argv[4]):
    r_cfunc_array.append(float(i.split()[0]))
    cfunc_ori_array.append(float(i.split()[1]))

for i in range(5,3005,10):
    if i==0:
        continue
    for j in range(len(r_cfunc_array)):
        if r_cfunc_array[j]>i:
            cfunc_use_array.append(cfunc_ori_array[j])
            break

#P0=[T0_0, N1_0, N2_0, RS_0, A0_0, GAMMA0_0, K0_0, A1_0, GAMMA1_0, \
#        K1_0, N3_0, Z_0, RHO_0, A_0, BETA_0, NT_0, NM_0, X_0, NE_0,TC_0]
P0=[N1_0,N2_0,RS_0,A0_0,GAMMA0_0,K0_0,N3_0,RHO_0,TC,NE_0]
ORI_POST=postpob(*P0)
p_tmp=P0
#bnds=((T0_MIN,T0_MAX),\
#        (N1_MIN,N1_MAX),(None,None),(None,None),(None,None),(None,None),\
#        (K0_MIN,K0_MAX),(None,None),(None,None),(None,None),(N3_MIN,N3_MAX),\
#        (Z_MIN,Z_MAX),(None,None),(None,None),(BETA_MIN,BETA_MAX),(None,None),\
#        (NM_MIN,NM_MAX))
f = lambda p1: postpob(*p1)
result=scipy.optimize.minimize(f,p_tmp,method='Powell', options={'maxiter':300,'disp':True})
print(result.x)
n1_f=N1_0
n2_f=result.x[1]
rs_f=result.x[2]
a0_f=result.x[3]
gamma0_f=result.x[4]
k0_f=result.x[5]
n3_f=result.x[6]
if n3_f > N3_MAX:
    n3_f=N3_MAX
if n3_f <N3_MIN:
    n3_f=N3_MIN
rho_f=result.x[7]
ne0_f=result.x[9]
tc0_f=result.x[8]
ne_array=[]
T_array=[]
rne_array=numpy.array(range(5,3005,10))
tmp_array=numpy.array(range(0,3010,10))
y0=[tc0_f,ne0_f]
p=[n1_f,n2_f,rs_f,a0_f,gamma0_f,k0_f,n3_f,rho_f,R500_0]
#y=odeint(dydr,y0,rne_array,args=(p,),Dfun=dfdy)
y=odeint(dydr,y0,rne_array,args=(p,),Dfun=None,h0=1e-5,hmax=1,hmin=1e-6,\
        rtol=[0.01,0.01],atol=[1e-3,1e-6])
ne_array.append(0)
for i in range(len(y)):
    ne_array.append(y[i][1])
    T_array.append(y[i][0])
ne_array=numpy.array(ne_array)
T_array=numpy.array(T_array)
ne_array=numpy.insert(ne_array,0,0.0)
sbp_fit=[]
tmp_r_use=range(5,3000,10)
for i in range(5,3000,10):
    a=deproject_model.calc_sb(i,tmp_array,ne_array,cfunc_use_array,cm_per_pixel)
    sbp_fit.append(a)
print(sbp_fit)
#plt.loglog(tmp_r_use,sbp_fit,'b-')
plt.hold
#plt.loglog(rsbp_array,sbp_array,'g*')
plt.show()
plt.plot(rne_array,T_array,'b-')
plt.hold
plt.plot(r_array,t_array,'g*')
plt.show()




#print(result.message)
#print(result.x)
#adjust(result.x)
#print(result.x)
#final_pos=postpob(*result.x)
#print(ORI_POST,final_pos)
#print_result(result.x)
#result_temp=[result.x[0],result.x[1],result.x[2],result.x[3],result.x[4]\
#        ,result.x[5],result.x[6],result.x[7],result.x[8],result.x[9],
#        result.x[10],result.x[11],result.x[12],result.x[13],result.x[14]
#        ,result.x[15],result.x[16],result.x[17]]
#plot_result(P0,result_temp)
#ne_array=[]
#rne_array=range(5,3005,10)
#ne_array=scipy.integrate.odeint(dnedr,result.x[18],rne_array,args=(result_temp,))
#print(ne_array)
#sbp_fit=[]
#tmp_array=range(0,3010,10)
#ne_array=numpy.insert(ne_array,0,0.0)
#tmp_r_use=range(5,1000,10)
#for i in range(5,1000,10):
#    a=deproject_model.calc_sb(i,tmp_array,ne_array,cfunc_use_array,cm_per_pixel)
#    sbp_fit.append(a)
#print(sbp_fit)
#plt.loglog(tmp_r_use,sbp_fit,'b-')
#plt.hold
#plt.loglog(rsbp_array,sbp_array,'g*')
#plt.show()



if len(sys.argv) >= 10:
    if final_pos>-1e-100:
        fi=open(sys.argv[4],'a')
        print('%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e'\
            %(R500_0,R200_0,result.x[0],result.x[1],result.x[2],result.x[3],\
            result.x[4],result.x[5],result.x[6],result.x[7],result.x[8],result.x[9],\
            result.x[10],result.x[11],result.x[12],result.x[13],result.x[14],\
            result.x[15],result.x[16],result.x[17],result.x[18]),file=fi)
        fi.close()
if len(sys.argv)>=10:
    if final_pos>-1e-100:
        fi=open(sys.argv[5],'a')
        print(t_array,file=fi)
        fi.close()

if len(sys.argv)>=10:
    rkout_array=[]
    rkouterr_array=[]
    kout_array=[]
    if final_pos>-1e-100:
        for i in open(sys.argv[6]):
            tmp=float(i.split()[0])
            tmperr=float(0)
            rkout_array.append(tmp)
            rkouterr_array.append(tmperr)
            kout=entropy_model.entropy(tmp,R500_0,result.x[7],result.x[8],result.x[9],result.x[18])
            kout_array.append(kout)
        fi=open(sys.argv[7],'a')
        for i in range(len(rkout_array)):
            print('%.4e %.4e %.4e'%(rkout_array[i],rkouterr_array[i],kout_array[i]),file=fi)
        fi.close()

if len(sys.argv)>=9:
    fi=open(sys.argv[8],'a')
    for i in range(3000):
        rtmp=float(i)
        kout=entropy_model.entropy(rtmp,R500_0,result.x[7],result.x[8],result.x[9],result.x[18])
        print('%.4e %.4e'%(rtmp,kout),file=fi)
    fi.close()


