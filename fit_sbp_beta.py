#!/usr/bin/env python3
##
# $1 sbp_data_file: radius,radius_err,sb,sb_err, radius in kpc
# $2 sbp_cfg_file: as GU
# $3 cooling function data file: radius in kpc
# $4 name of output parameter file
# $5 name of output sb_data file
import sys
import numpy
import numpy as np
import scipy
import math
import re
from scipy.optimize import fsolve
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from zzh_model import gpob
import deproject_model
import json
import random
numpy.set_printoptions(linewidth=4900)
def swap(a,b):
    return [b,a]

def shuffle(a_array,err_array):
    if len(a_array)!=len(err_array):
        print('error from shuffle in fit_sbp_beta.py: please check the length of input arrays')
        return -1
    a_array=np.array(a_array)
    err_array=np.array(err_array)
    anew=random.gauss(a_array,err_array)
    return anew

def adjust(par,p_min_array,p_max_array):
    for i in range(len(par)):
        if par[i]<p_min_array[i]:
            par[i]=p_min_array[i]
        if par[i]>p_max_array[i]:
            par[i]=p_max_array[i]
#    if model=='dbeta':
#        if par[1]<par[4]:
#            par[1]=par[4]
#        if par[2]>par[5]:
#            [par[2],par[5]]=swap(par[2],par[5])
    return par

def likehood(par,p_min_array,p_max_array,rmodel_array,r_array,sb_array,sbe_array,cfunc_use_array,cfunc_cuse_array,flag_array,model):
    pro=0
    par=adjust(par,p_min_array,p_max_array)
    den_array=[]
    den_array.append(0)
    if model=='beta':
        bkg=par[3]
        par_use=par[0:3]
        den_array=deproject_model.beta(rmodel_array,*par_use)
        pro=pro+gpob(par[1],0.66,0.3)
    elif model=='dbeta':
        bkg=par[6]
        par_use=par[0:6]
#        pro=pro+gpob(par_use[4],0.66,0.05)
        den_array=deproject_model.dbeta(rmodel_array,*par_use)
    else:
        print('error:model does not exist yet')
        return 'NaN'
    for i in range(len(r_array)):
        if r_array[i]==0:
            continue
        if flag_array[i]=='o':
            sb_m=deproject_model.calc_sb(r_array[i],rmodel_array,den_array,cfunc_use_array)
            sb_m=sb_m+bkg
            pro=pro+gpob(sb_m,sb_array[i],sbe_array[i])
        if flag_array[i]=='c':
            sb_m=deproject_model.calc_sb(r_array[i],rmodel_array,den_array,cfunc_cuse_array)
            pro=pro+gpob(sb_m,sb_array[i],sbe_array[i])
    return -pro

def plot_result(par,p_min_array,p_max_array,rne_array,r_array,sb_array,sbe_array,cfunc_use_array,cfunc_cuse_array,flag_array,model,re_array):
    den_sum_array=[]
    sbp_sum_array=[]
    sbpc_sum_array=[]
    par=np.array(par)
    if len(par.shape) == 1:
        print('waring: make sure it is a center only plot')
        ind50=0
        ind16=0
        ind84=0
    elif len(par.shape) == 2:
        if len(par)< 6:
            print('error from plot_result in fit_sbp_beta: need more monte-carlo times')
            return -1
        else:
            ind50=int(len(par)*0.5)
            ind16=int(len(par)*0.16)
            ind84=int(len(par)*0.84)
    for i in range(len(par)):
        plot_y=[]
        plot_yc=[]
        den=[]
        den.append(0)
        plot_y.append(0)
        plot_yc.append(0)
        if model=='beta':
            bkg=par[i][3]
            par_use=par[i][0:3]
            den=deproject_model.beta(rne_array,*par_use)
        if model=='dbeta':
            bkg=par[i][6]
            par_use=par[i][0:6]
            den=deproject_model.dbeta(rne_array,*par_use)
        den_sum_array.append(den)
        for i in rne_array:
            if i==0:
                continue
            plot_y.append(deproject_model.calc_sb(i,rne_array,den,cfunc_use_array)+bkg)
            plot_yc.append(deproject_model.calc_sb(i,rne_array,den,cfunc_cuse_array))
        sbp_sum_array.append(plot_y)
        sbpc_sum_array.append(plot_yc)
    sbp_sum_array=np.sort(sbp_sum_array,0)
    sbpc_sum_array=np.sort(sbpc_sum_array,0)
    den_sum_array=np.sort(den_sum_array,0)
    plt.loglog(rne_array,sbp_sum_array[ind50],'r-',lw=0.5,label='data from others')
    plt.loglog(rne_array,sbpc_sum_array[ind50],'b-',lw=0.5,label='our chandra data')
    plt.xlim(20,2000)
    plt.errorbar(r_array,sb_array,xerr=re_array,yerr=sbe_array,color='k',linestyle='none')
#    print(den)
    plt.fill_between(rne_array,sbp_sum_array[ind16],sbp_sum_array[ind84],color='r',alpha=0.3)
    plt.fill_between(rne_array,sbpc_sum_array[ind16],sbpc_sum_array[ind84],color='b',alpha=0.3)
    plt.legend()
    plt.savefig('sbp_beta_fit.pdf',dpi=100)
    return [den_sum_array[ind50],den_sum_array[ind16],den_sum_array[ind84]]

def main(sbp_data_file,param_file,cfunc_file,cfunc_cfile,json_file):
    tmp1=list(range(0,11))
    tmp2=list(range(11,41,3))
    tmp3=list(range(41,3001,5))
    tmp1.extend(tmp2)
    tmp1.extend(tmp3)
    rne_array=numpy.array(tmp1)
    r_array=[]
    re_array=[]
    sb_array=[]
    sbe_array=[]
    flag_array=[]
    cfunc_ori_array=[]
    r_cfunc_array=[]
    cfunc_use_array=[]
    cfunc_cori_array=[]
    rc_cfunc_array=[]
    cfunc_cuse_array=[]
    p0=[]
    p_min_array=[]
    p_max_array=[]
    flag_array.append('0')
    r_array.append(0)
    re_array.append(0)
    sb_array.append(0)
    sbe_array.append(0)
    cfunc_ori_array.append(0)
    cfunc_cuse_array.append(0)
    r_cfunc_array.append(0)
    cfunc_use_array.append(0)
    for i in open(sbp_data_file):
        tmp=i.split()
        r_array.append(float(tmp[0]))
        re_array.append(float(tmp[1]))
        sb_array.append(float(tmp[2]))
        sbe_array.append(float(tmp[3]))
        flag_array.append(tmp[4])
    r_array=numpy.array(r_array)
    re_array=numpy.array(re_array)
    sb_array=numpy.array(sb_array)
    sbe_array=numpy.array(sbe_array)
    sbe_array=sbe_array+sb_array*0.1

    for i in open(param_file,'r'):
        if re.match(r'^n01\s',i):
            model='dbeta'
            n01=float(i.split()[1])
            if len(i.split())==4:
                n01_min=float(i.split()[2])
                n01_max=float(i.split()[3])
        if re.match(r'^n02\s',i):
            n02=float(i.split()[1])
            if len(i.split())==4:
                n02_min=float(i.split()[2])
                n02_max=float(i.split()[3])
        if re.match(r'^rc1\s',i):
            rc1=float(i.split()[1])
            if len(i.split())==4:
                rc1_min=float(i.split()[2])
                rc1_max=float(i.split()[3])
        if re.match(r'^rc2\s',i):
            rc2=float(i.split()[1])
            if len(i.split())==4:
                rc2_min=float(i.split()[2])
                rc2_max=float(i.split()[3])
        if re.match(r'^beta1\s',i):
            beta1=float(i.split()[1])
            if len(i.split())==4:
                beta1_min=float(i.split()[2])
                beta1_max=float(i.split()[3])
        if re.match(r'^beta2\s',i):
            beta2=float(i.split()[1])
            if len(i.split())==4:
                beta2_min=float(i.split()[2])
                beta2_max=float(i.split()[3])
        if re.match(r'^n0\s',i):
            model='beta'
            n0=float(i.split()[1])
            if len(i.split())==4:
                n0_min=float(i.split()[2])
                n0_max=float(i.split()[3])
        if re.match(r'^rc\s',i):
            rc=float(i.split()[1])
            if len(i.split())==4:
                rc_min=float(i.split()[2])
                rc_max=float(i.split()[3])
        if re.match(r'^beta\s',i):
            beta=float(i.split()[1])
            if len(i.split())==4:
                beta_min=float(i.split()[2])
                beta_max=float(i.split()[3])
        if re.match(r'^bkg\s',i):
            bkg=float(i.split()[1])
            if bkg==0:
                bkg=1e-10
            if len(i.split())==4:
                bkg_min=float(i.split()[2])
                bkg_max=float(i.split()[3])

    for i in open(cfunc_file):
        r_cfunc_array.append(float(i.split()[0]))
        cfunc_ori_array.append(float(i.split()[1]))
    for i in open(cfunc_cfile):
        rc_cfunc_array.append(float(i.split()[0]))
        cfunc_cori_array.append(float(i.split()[1]))

#    r_max=int(r_array[len(r_array)-1]+re_array[len(re_array)-1])

    for i in rne_array:
        if i==0:
            continue
        for j in range(len(r_cfunc_array)):
            if r_cfunc_array[j]>i:
                cfunc_use_array.append(cfunc_ori_array[j])
                break
        for j in range(len(rc_cfunc_array)):
            if rc_cfunc_array[j]>i:
                cfunc_cuse_array.append(cfunc_cori_array[j])
                break
    if model=='beta':
        p0=[n0,beta,rc,bkg]
        p_min_array=[0,0,0,0]
        p_max_array=[9999,9999,9999,9999]
        if 'n0_min' in locals().keys():
            p_min_array[0]=n0_min
            p_max_array[0]=n0_max
        if 'beta_min' in locals().keys():
            p_min_array[1]=beta_min
            p_max_array[1]=beta_max
        if 'rc_min' in locals().keys():
            p_min_array[2]=rc_min
            p_max_array[2]=rc_max
        if 'bkg_min' in locals().keys():
            p_max_array[3]=bkg_min
            p_min_array[3]=bkg_max
    if model=='dbeta':
        p0=[n01,beta1,rc1,n02,beta2,rc2,bkg]
        p_min_array=[0,0,0,0,0,0,0]
        p_max_array=[9999,9999,9999,9999,9999,9999,9999]
        if 'n01_min' in locals().keys():
            p_min_array[0]=n01_min
            p_max_array[0]=n01_max
        if 'beta1_min' in locals().keys():
            p_min_array[1]=beta1_min
            p_max_array[1]=beta1_max
        if 'rc1_min' in locals().keys():
            p_min_array[2]=rc1_min
            p_max_array[2]=rc1_max
        if 'n02_min' in locals().keys():
            p_min_array[3]=n02_min
            p_max_array[3]=n02_max
        if 'beta2_min' in locals().keys():
            p_min_array[4]=beta2_min
            p_max_array[4]=beta2_max
        if 'rc2_min' in locals().keys():
            p_min_array[5]=rc2_min
            p_max_array[5]=rc2_max
        if 'bkg_min' in locals().keys():
            p_min_array[6]=bkg_min
            p_max_array[6]=bkg_max
    p_all=[]
    for i in range(100):
        if i==0:
            sb_array_this=sb_array
        else:
            sb_array_this=shuffle(sb_array,sbe_array)
        result=minimize(likehood,p0,args=(p_min_array,p_max_array,rne_array,r_array,sb_array_this,sbe_array,cfunc_use_array,cfunc_cuse_array,flag_array,model),method='Powell')
        adjust(result.x,p_min_array,p_max_array)
        print(result.x)
        init_pro=likehood(p0,p_min_array,p_max_array,rne_array,r_array,sb_array,sbe_array,cfunc_use_array,cfunc_cuse_array,flag_array,model)
        final_pro=likehood(result.x,p_min_array,p_max_array,rne_array,r_array,sb_array,sbe_array,cfunc_use_array,cfunc_cuse_array,flag_array,model)
        print(init_pro,final_pro)
        p_all.append(result.x)
    sum_ne_array=plot_result(p_all,p_min_array,p_max_array,rne_array,r_array,sb_array,sbe_array,cfunc_use_array,cfunc_cuse_array,flag_array,model,re_array)
#    if len(sys.argv)>4:
#        fi=open(sys.argv[4],'a')
#        print(result.x,file=fi)
#        fi.close()
#    if len(sys.argv)>5:
#        fi=open(sys.argv[5],'a')
#        print(sb_array,file=fi)
#        fi.close()
#    return 0
    fi=open(json_file)
    dat=json.load(fi)
    t_array=np.array(dat['temperature_model'][0],dtype=float)
    t_array=np.insert(t_array,0,0)
    t_array_up=np.array(dat['temperature_model'][2],dtype=np.float)
    t_array_down=np.array(dat['temperature_model'][1],dtype=np.float)
    t_array_up=np.insert(t_array_up,0,0)
    t_array_down=np.insert(t_array_down,0,0)
    ne_array=sum_ne_array[0]
    ne_array_down=sum_ne_array[1]
    ne_array_up=sum_ne_array[2]
    te_array=(t_array_up-t_array_down)/2/t_array
    nee_array=(ne_array_up-ne_array_down)/2/ne_array
    k_array=t_array*np.power(ne_array,-2/3)
    ke_array=np.sqrt(np.square(te_array)+np.square(nee_array))*k_array
    kori_array=np.array(dat['k_fit'][0],dtype=np.float)
    kori_array=np.insert(kori_array,0,0)
    kori_array_up=np.array(dat['k_fit'][2],dtype=np.float)
    kori_array_down=np.array(dat['k_fit'][1],dtype=np.float)
    kori_array_down=np.insert(kori_array_down,0,0)
    kori_array_up=np.insert(kori_array_up,0,0)
    plt.clf()
    plt.loglog(rne_array,k_array,label='beta model fit')
    plt.fill_between(rne_array,k_array-ke_array,k_array+ke_array,alpha=0.3)
    plt.loglog(rne_array,kori_array,label='ori model')
    plt.fill_between(rne_array,kori_array_down,kori_array_up,alpha=0.3)
    plt.legend()
    plt.xlim(10,2000)
    plt.savefig('entropy_beta_compare.pdf')
    np.save('beta_fit_result_all',p_all)

if __name__=='__main__':
    param_file=sys.argv[2]
    sbp_data_file=sys.argv[1]
    cfunc_file=sys.argv[3]
    cfunc_cfile=sys.argv[4]
    json_file=sys.argv[5]
    main(sbp_data_file,param_file,cfunc_file,cfunc_cfile,json_file)

