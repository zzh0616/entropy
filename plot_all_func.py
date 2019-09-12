#!/usr/bin/env python3

# plot pdf from json file
# plot merged entropy from all json file

# $1  mode  1:plot temp,sbp,mass fit and entropy profile in local dir.
#           2:plot entropy profiles for all sources in the list file
#           3:compare entropy profiles with/without extra heating
# $2 in  list of the name of the source
#
# in mode 2 this program will try to read extra data to compare the entropy in outer region, which is named as a1795_compare.txt
#the compare file should have 4 rows, first for radius,2-4 for entropy,lowerlimit,upperlimit,both should be scaled.
import matplotlib.pyplot as plt
import sys
import json
import numpy as np
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import pylab
import re as the_re
import os
import pymc
#mode = sys.argv[1]
def main(mode,name):
    if mode[0] == "1" or mode == "2" or mode =="3" :
#    if mode == "2":
#        pylab.figure('k_global')
#    if mode == "1":
#        pp=PdfPages('all.pdf')
#    for name in open(sys.argv[2]):
        name=name[0:-1]
        print(name)
        #name=sys.argv[2]
        json_file = name+"/"+name+"_plt.json"
        global_file= name + '/global.cfg'
        for i in open(global_file,'r'):
            if the_re.match(r'^radius_sbp_file',i):
                sbp_file_tmp = i.split()[1]
        sbp_file= name +"/" + sbp_file_tmp
        fi=open(json_file,'r')
        dat = json.load(fi)
        sbp_array = []
        sbpe_array= []
        rsbp_array = []
        rsbpe_array = []
        for i in open(sbp_file):
            [r, re, sbp, sbpe] = i.split()
            sbp_array.append(np.float(sbp))
            sbpe_array.append(np.float(sbpe))
            rsbp_array.append(np.float(r))
            rsbpe_array.append(np.float(re))
        r200 = np.float(dat["r200"])
        sbp_array = np.array(sbp_array)
        sbpe_array = np.array(sbpe_array)
        rsbp_array = np.array(rsbp_array)
        rsbpe_array = np.array(rsbpe_array)
        if rsbp_array[0] == rsbpe_array[0]:
            rsbpe_array[0] = rsbpe_array[0]*0.99
        sbpe_array=sbpe_array+sbp_array*0.03
        t_array = []
        te_array = []
        rt_array = []
        rte_array = []
        temp_file=name+"/"+name+".txt"
        for i in open(temp_file):
            [rt,rte,t,te] = i.split()
            t_array.append(np.float(t))
            te_array.append(np.float(te))
            rt_array.append(np.float(rt))
            rte_array.append(np.float(rte))
        t_array = np.array(t_array)
        te_array = np.array(te_array)
        rt_array = np.array(rt_array)
        rte_array = np.array(rte_array)
        te_array=te_array+t_array*0.03
        if rt_array[0] == rte_array[0]:
            rte_array[0] = rte_array[0]*0.99
        r_model = dat["radius_model"]
        r_model = np.array(r_model,dtype=float)
        sbp_fit = dat["sbp_model"][0]
        sbp_fit = np.array(sbp_fit)
        sbp_fit_down = dat["sbp_model"][1]
        sbp_fit_down = np.array(sbp_fit_down)
        sbp_fit_up = dat["sbp_model"][2]
        sbp_fit_up = np.array(sbp_fit_up)
        k_fit = dat["k_fit"][0]
        k_fit = np.array(k_fit)
        k_fit_down = dat["k_fit"][1]
        k_fit_down = np.array(k_fit_down)
        k_fit_up = dat["k_fit"][2]
        k_fit_up = np.array(k_fit_up)
        ne_array = np.array(dat["den_model"][0],dtype=float)
        ne_cl_array = np.array(dat["den_cl_model"][0],dtype=float)
        m_calc = dat["m_calc"][0]
        m_calc = np.array(m_calc)
        m_calc_down = dat["m_calc"][1]
        m_calc_down = np.array(m_calc_down)
        m_calc_up = dat["m_calc"][2]
        m_calc_up = np.array(m_calc_up)
        m_fit = dat["m_fit"][0]
        m_fit=np.array(m_fit)
        m_fit_down = dat["m_fit"][1]
        m_fit_down = np.array(m_fit_down)
        m_fit_up = dat["m_fit"][2]
        m_fit_up = np.array(m_fit_up)
        t_fit = dat["temperature_model"][0]
        t_fit = np.array(t_fit)
        t_fit_down = dat["temperature_model"][1]
        t_fit_down = np.array(t_fit_down)
        t_fit_up = dat["temperature_model"][2]
        t_fit_up = np.array(t_fit_up)
        r200= np.float(dat["r200"])
        fitting_param_file=name+"/"+"param_zzh_for_py.txt"
        for i in open(fitting_param_file):
            if the_re.match(r'^R500',i):
                r500=np.float(i.split()[1])
        param_file=name+"/"+name+"_result.csv"
        compare_file=name+"/"+name+"_compare_mod.txt"
        for i in open(param_file):
            if the_re.match(r'^kmod_a',i):
                kmod_a=np.float(i.split(",")[1])
            if the_re.match(r'^kmod_b',i):
                kmod_b=np.float(i.split(',')[1])
            if the_re.match(r'^kmod_c',i):
                kmod_c=np.float(i.split(',')[1])
            if the_re.match(r'kmod_k0',i):
                kmod_k0=np.float(i.split(',')[1])
            if the_re.match(r'^cp_e',i):
                cp_e=np.float(i.split(',')[1])
            if the_re.match(r'^cp_p',i):
                cp_p=np.float(i.split(',')[1])
            if the_re.match(r'^cp_g0',i):
                cp_g0=np.float(i.split(',')[1])
            if the_re.match(r'^cp_x0',i):
                cp_x0=np.float(i.split(',')[1])
            if the_re.match(r'^cp_sigma',i):
                cp_sigma=np.float(i.split(',')[1])
            if the_re.match(r'^a0',i):
                a0=np.float(i.split(',')[1])
            if the_re.match(r'^gamma0',i):
                gamma0=np.float(i.split(',')[1])
            if the_re.match(r'^k0',i):
                k0=np.float(i.split(',')[1])
            if the_re.match(r'^sbp_c',i):
                sbp_c=np.float(i.split(',')[1])
                sbp_c_down=sbp_c-np.float(i.split(',')[2])
            if the_re.match(r'^ref_sbp_c',i):
                ref_sbp_c=np.float(i.split(',')[1])
                ref_sbp_c_down=ref_sbp_c-np.float(i.split(',')[2])
        k_norm=kmod_a*np.power(0.3*r200,kmod_b)*np.exp(kmod_c*-0.3)+kmod_k0
        c_factor=np.power(1+0.3,cp_p)*np.exp(0.3*(cp_e))+0+cp_g0*np.exp(-(0.3-cp_x0)*(0.3-cp_x0)/cp_sigma)
#            k_norm=k_norm*np.power(c_factor,-2/3)
        r_norm=r200
        name_tmp=list(name)
        if name[0:5]=='abell':
            name_tmp[1:5]=''
        for i in range(len(name)):
            if name[i]=='_':
                name_tmp[i]=" "
            if 'a'<=name[i]<='z':
                name_tmp[i]=name_tmp[i].upper()
            if name[0:4]=='400d':
                name_tmp[3]='d'
            if name[0:5]=='hydra':
                name_tmp[0:5]='Hydra'
        name_use=''.join(name_tmp)
        if mode == '1t':
            matplotlib.rcParams['xtick.direction'] = 'in'
#            fig,axarr=plt.subplots(2,2,sharex='col')
#            fig.subplots_adjust(hspace=0,wspace=0,left=0.2,right=0.8)
#            ax=axarr[0][1]
#            ax.loglog(r_model,m_fit,'b')
#            ax.fill_between(r_model,m_calc_down,m_calc_up,color='green')
#            ax.fill_between(r_model,m_fit_down,m_fit_up,color='grey')
#            ax.yaxis.set_ticks_position('right')
#            ax.set_xlim(10,2000)
#            ax.text(1.25,0.5,'Total Mass ($M_{\odot}$)',rotation=90,transform=ax.transAxes,verticalalignment="center",horizontalalignment="left")
#            ymin,ymax=plt.ylim()
#            print(ymin,ymax)
#            ax.loglog([r200-1,r200],[1e12,1e15])
#            ax=axarr[1][1]
#            ax.loglog(r_model,k_fit,'b')
#            ax.fill_between(r_model,k_fit_down,k_fit_up,color='grey')
#            ax.yaxis.set_ticks_position('right')
#            ax.text(1.25,0.5,'Entropy(keV cm^2)',rotation=90,transform=ax.transAxes,verticalalignment="center",horizontalalignment="left")
#            ax.set_xlabel('radius (kpc)')
#            fig.suptitle(name+' (r200= ' + np.str(r200)+' kpc)')
#            ymin,ymax=plt.ylim()
#            ax.vlines(r200,ymin,ymax)
#        ax.set_xlim(10,2000)
#            ax=axarr[0][0]
            rt_array_new=[]
            t_array_new=[]
            te_array_new=[]
            rte_array_new=[]
            i_cuse=-1
            for i in range(len(rt_array)):
                if rt_array[i]+rte_array[i]>=0.3*r500 or i==len(rt_array)-1:
                    rt_array_new.append(rt_array[i])
                    t_array_new.append(t_array[i])
                    te_array_new.append(te_array[i])
                    rte_array_new.append(rte_array[i])
                    if i_cuse == -1:
                        i_cuse=i
            plt.plot(r_model,t_fit)
            plt.errorbar(rt_array[0:i_cuse],t_array[0:i_cuse],xerr=rte_array[0:i_cuse],yerr=te_array[0:i_cuse],color='k',linestyle='none')
            plt.errorbar(rt_array_new,t_array_new,xerr=rte_array_new,yerr=te_array_new,color='c',linestyle='none')
#            plt.plot(r_model,t_fit_down)
#            plt.plot(r_model,t_fit_up)
            plt.fill_between(r_model,t_fit_down,t_fit_up,color='grey')
#            ax.set_ylabel('Temperature (keV)')
            plt.xlim(10,2000)
            plt.ylim(2,21)
            plt.text(20,3,name_use,fontsize=5)
#            ymin,ymax=plt.ylim()
#            ax.loglog([r200,r200],[1,20])i
        if mode =='1s':
 #           ax=axarr[1][0]
            rsbp_array_new=[]
            sbp_array_new=[]
            sbpe_array_new=[]
            rsbpe_array_new=[]
            i_cuse=-1
            for i in range(len(rsbp_array)):
                if rsbp_array[i]>=0.3*r500:
                    rsbp_array_new.append(rsbp_array[i])
                    sbp_array_new.append(sbp_array[i])
                    sbpe_array_new.append(sbpe_array[i])
                    rsbpe_array_new.append(rsbpe_array[i])
                    if i_cuse == -1:
                        i_cuse=i
            plt.errorbar(rsbp_array[0:i_cuse],sbp_array[0:i_cuse],xerr=rsbpe_array[0:i_cuse],yerr=sbpe_array[0:i_cuse],color='k',linestyle='none')
            plt.errorbar(rsbp_array_new,sbp_array_new,xerr=rsbpe_array_new,yerr=sbpe_array_new,color='c',linestyle='none')
#            plt.errorbar(rsbp_array,sbp_array,xerr=rsbpe_array,yerr=sbpe_array,color='k',linestyle='none')
            plt.loglog(r_model,sbp_fit+ref_sbp_c)
            plt.vlines(0.3*r500,1e-11,1e-6,linestyles='dashed')
#            plt.loglog(r_model,sbp_fit_up)
#            plt.loglog(r_model,sbp_fit_down)
            sbp_fit_down=sbp_fit_down+ref_sbp_c_down
#            sbp_fit_down=sbp_fit_down-sbp_fit_down[0:600].min()+ref_sbp_c_down
            sbp_fit_up=sbp_fit_up+ref_sbp_c*2-ref_sbp_c_down
            plt.fill_between(r_model[0:500],sbp_fit_down[0:500],sbp_fit_up[0:500],color='grey')
            plt.xlim(10,2000)
            plt.ylim(3e-10,5e-7)
            plt.text(15,1e-9,name_use,fontsize=5)
#            ax.set_ylabel('SBP (cm^-2 pixel^-2 s^-1)')
#            ax.set_xlabel('radius (kpc)')
#            ymin,ymax=plt.ylim()
 #           ax.vlines(r200,ymin,ymax)
#        ax.set_xlim(10,2000)
#            plt.savefig('pdf/'+name+'.pdf',dpi=10)
#            plt.savefig(name+"/"+name+'.pdf')
#            pp.savefig(fig)
        if mode == "2" :
            r_entropy=r_model/r_norm
#            matplotlib.rcParams['xtick.direction'] = 'in'
            k_entropy=k_fit/k_norm
            k_entropy_up=k_fit_up/k_norm
            k_entropy_down=k_fit_down/k_norm
 #           pylab.figure('k_local')
#            plt.clf()
            plt.loglog(r_entropy,k_entropy)
#            plt.text(0.1,3.5,'with clumping correction',color='lightskyblue')
            plt.fill_between(r_entropy,k_entropy_down,k_entropy_up,color='grey')
#            plt.xlabel('Radius ($r/r_{200}$)')
#            plt.ylabel('Entropy ($k/k(0.3r_{200})$)')
            plt.xlim(0.05,2)
            plt.ylim(0.1,4)
            plt.plot([0.01,1],[0.0237,4])
#            plt.tick_params(axis='x',direction='in')
            k_entropy_uncorr=k_entropy*np.power(ne_cl_array/ne_array,-2/3)
            k_entropy_uncorr_up=k_entropy_up*np.power(ne_cl_array/ne_array,-2/3)
            k_entropy_uncorr_down=k_entropy_down*np.power(ne_cl_array/ne_array,-2/3)
            plt.loglog(r_entropy,k_entropy_uncorr,'g')
#            plt.text(0.1,3,'without clumping correction',color='g')
            plt.fill_between(r_entropy,k_entropy_uncorr_down,k_entropy_uncorr_up,color='palegreen',alpha=1)
            plt.text(0.07,2.5,name_use,fontsize=6)
            if os.path.isfile(compare_file):
                print("compare file found")
                r_compare=[]
                k_compare=[]
                k_compare_low=[]
                k_compare_up=[]
                for i in open(compare_file):
                    r_compare.append(float(i.split()[0]))
                    k_compare.append(float(i.split()[1]))
                    k_compare_low.append(-float(i.split()[2])+float(i.split()[1]))
                    k_compare_up.append(float(i.split()[3])-float(i.split()[1]))
                r_compare=np.array(r_compare)/r_norm
                k_compare=np.array(k_compare)/k_norm
                k_compare_low=np.array(k_compare_low)/k_norm
                k_compare_up=np.array(k_compare_up)/k_norm
                plt.errorbar(r_compare,k_compare,yerr=[k_compare_low,k_compare_up],color='k',linestyle='none',marker="+")
#                plt.text(0.1,2.6,'other\'s result',color='k')
#            plt.savefig(name+'/'+name+'_k_scaled.pdf')
#            pylab.figure('k_global')
#            plt.loglog(r_entropy, k_entropy)
        if mode == "3":
            pylab.figure('k_compare_all')
            k_ori=a0*np.power(r_model,gamma0)+k0
            k_mod=kmod_a*np.power(r_model,kmod_b)*np.exp(-kmod_c*r_model/r200)+kmod_k0
            plt.loglog(r_model/r200,k_ori/k_norm,'b',alpha=0.1)
            plt.loglog(r_model/r200,k_fit/k_norm,'r',alpha=0.3)
'''    if mode == "2":
        pylab.figure('k_global')
        plt.xlabel('Radius ($r/r_{200}$)')
        plt.ylabel('Entropy ($k/k(0.3r_{200})$)')
#        plt.ylabel('Entropy')
        std_x=[0.123,1]
        std_y=[0.4,4]
        plt.xlim(0.1,2.0)
        plt.ylim(0.4,4)
        plt.loglog(std_x,std_y,'k')
        plt.savefig('entropy.pdf')
    if mode == '1':
        pp.close()
    if mode == '3':
        plt.xlabel('Radius ($r/r_{200}$)')
        plt.ylabel('Entropy ($k/k(0.3r_{200})$)')
        plt.savefig('entropy_totalvsgravity.pdf')
if mode == 'c':
    for name in sys.argv[2]:
        json_file = name+"/"+name+"_plt.json"
        db_file= name+"/"+'sample.pickle'
        r_model=np.array(json_file["radius_model"],dtype=float)
        M2=pymc.database.pickle.load(db_file)
        a0_f=M2.trace('a0')[:]
        r200=np.float(json_file['r200'])
        gamma0_f=M2.trace('gamma0')[:]
        k0_f=M2.trace('k0')[:]
        kmod_a_f=M2.trace('kmod_a')[:]
        kmod_b_f=M2.trace('kmod_b')[:]
        kmod_c_f=M2.trace('kmod_c')[:]
        kmod_k0_f=M2.trace('kmod_k0')[:]
        k_ori_total=[]
        k_mod_total=[]
        for i in range(len(a0_f)):
            k_ori=a0_f[i]*np.power(r_model,gamma0_f[i])+k0_f[i]
            k_mod=kmod_a_f[i]*np.power(r_model,kmod_b_f[i])*np.exp(kmod_c_f[i]*-r_model/R200)+kmod_k0_f[i]
            k_ori_total.append(k_ori)
            k_mod_total.append(k_mod)
'''

