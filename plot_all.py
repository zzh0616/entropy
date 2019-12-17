#!/usr/bin/env python3

# plot pdf from json file
# plot merged entropy from all json file

# $1  mode  1:plot temp,sbp,mass fit and entropy profile in local dir.
#           2:plot entropy profiles for all sources in the list file
#           3:compare entropy profiles with/without extra heating, now abandoned
#           4:plot the clumping factor porfiles for all clusters
#           5:plot 3 sample entropy profiles in one figure with errorbar for all sample cluster
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
mode = sys.argv[1]
if mode == "1" or mode == "2" or mode =="3" or mode == "4" or mode == "5":
    if mode == "2":
        pylab.figure('k_global')
    if mode == "1":
        pp=PdfPages('all.pdf')
    if mode == "5":
        count5=0
        tmp_name=1
    for name in open(sys.argv[2]):
        name=name[0:-1]
#        print(name)
        #name=sys.argv[2]
        json_file = name+"/"+name+"_plt.json"
        global_file= name + '/global.cfg'
        for i in open(global_file,'r'):
            if the_re.match(r'^sbp_data_file',i):
                sbp_file_tmp = i.split()[1]
        sbp_file= name +"/" + sbp_file_tmp
        fi=open(json_file,'r')
        dat = json.load(fi)
        sbp_array = []
        sbpe_array= []
        rsbp_array = []
        rsbpe_array = []
        for i in open(sbp_file):
            [r, re, sbp, sbpe,flag_proj] = i.split()
            sbp_array.append(np.float(sbp))
            sbpe_array.append(np.float(sbpe))
            rsbp_array.append(np.float(r))
            rsbpe_array.append(np.float(re))
        info_file=name+"/"+name+"_suminfo.txt"
        for i in open(info_file):
            if the_re.match(r'^k200',i):
                k200=float(i.split()[1])
        r200 = np.float(dat["r200"])
        sbp_array = np.array(sbp_array)
        sbpe_array = np.array(sbpe_array)
        rsbp_array = np.array(rsbp_array)
        rsbpe_array = np.array(rsbpe_array)
        if rsbp_array[0] == rsbpe_array[0]:
            rsbpe_array[0] = rsbpe_array[0]*0.99
        t_array = dat["temperature"][0]
        t_array = np.array(t_array)
        t_ave=t_array.mean()
        te_array = dat["temperature"][1]
        te_array = np.array(te_array)
        te_array=te_array+t_array*0.12
        rt_array = dat["radius"][0]
        rt_array = np.array(rt_array)
        rte_array = dat["radius"][1]
        rte_array = np.array(rte_array)
        t_ave=t_array.mean()
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
        csbp_fit = dat["csbp_model"][0]
        csbp_fit = np.array(csbp_fit)
        csbp_fit_down = dat["csbp_model"][1]
        csbp_fit_down = np.array(csbp_fit_down)
        csbp_fit_up = dat["csbp_model"][2]
        csbp_fit_up = np.array(csbp_fit_up)
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
        param_file=name+"/"+name+"_result.csv"
        fitting_param_file=name+"/"+"param_zzh_for_py.txt"
        for i in open(fitting_param_file):
            if the_re.match(r'^R500',i):
                r500=np.float(i.split()[1])
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
        knorm_index=0
        for i in range(len(r_model)):
            if 1.0*r500<r_model[i]:
                knorm_index=i
                break
            if i==len(r_model)-1:
                print('something wrong with 0.3r500')
                knorm_index=0
#        k_norm=k_fit[knorm_index]
        k_norm=k200
        c_factor=np.power(1+0.3,cp_p)*np.exp(0.3*(cp_e))+0+cp_g0*np.exp(-(0.3-cp_x0)*(0.3-cp_x0)/cp_sigma)
#            k_norm=k_norm*np.power(c_factor,-2/3)
        r_norm=r200
        if mode == "1":
            matplotlib.rcParams['xtick.direction'] = 'in'
            fig,axarr=plt.subplots(2,2,sharex='col')
            fig.subplots_adjust(hspace=0,wspace=0,left=0.2,right=0.8)
            ax=axarr[0][1]
            ax.loglog(r_model,m_fit,'b')
            ax.fill_between(r_model,m_calc_down,m_calc_up,color='green')
            ax.fill_between(r_model,m_fit_down,m_fit_up,color='grey')
            ax.yaxis.set_ticks_position('right')
            ax.set_xlim(10,2000)
            ax.text(1.25,0.5,'Total Mass ($M_{\odot}$)',rotation=90,transform=ax.transAxes,verticalalignment="center",horizontalalignment="left")
#            ymin,ymax=plt.ylim()
#            print(ymin,ymax)
#            ax.loglog([r200-1,r200],[1e12,1e15])
            ax=axarr[1][1]
            ax.loglog(r_model,k_fit/k_norm,'b')
            ax.fill_between(r_model,k_fit_down/k_norm,k_fit_up/k_norm,color='grey')
            ax.yaxis.set_ticks_position('right')
            ax.text(1.25,0.5,'Entropy(keV cm^2)',rotation=90,transform=ax.transAxes,verticalalignment="center",horizontalalignment="left")
            ax.set_xlabel('radius (kpc)')
            fig.suptitle(name+' (r200= ' + np.str(r200)+' kpc)')
#            ymin,ymax=plt.ylim()
#            ax.vlines(r200,ymin,ymax)
#        ax.set_xlim(10,2000)
            ax=axarr[0][0]
            ax.loglog(r_model,t_fit)
            ax.errorbar(rt_array,t_array,xerr=rte_array,yerr=te_array,color='k',linestyle='none')
            ax.loglog(r_model,t_fit_down)
            ax.loglog(r_model,t_fit_up)
            ax.fill_between(r_model,t_fit_down,t_fit_up,color='grey')
            ax.set_ylabel('Temperature (keV)')
            ax.set_xlim(10,2000)
            ax.set_ylim(1,20)
#            ymin,ymax=plt.ylim()
#            ax.loglog([r200,r200],[1,20])
            ax=axarr[1][0]
            ax.errorbar(rsbp_array,sbp_array,xerr=rsbpe_array,yerr=sbpe_array,color='k',linestyle='none')
            ax.loglog(r_model,sbp_fit)
            ax.loglog(r_model,sbp_fit_up)
            ax.loglog(r_model,sbp_fit_down)
            ax.fill_between(r_model,sbp_fit_down,sbp_fit_up,color='grey')
            ax.loglog(r_model,csbp_fit)
            ax.loglog(r_model,csbp_fit_up)
            ax.loglog(r_model,csbp_fit_down)
            ax.fill_between(r_model,csbp_fit_down,csbp_fit_up,color='grey')
            ax.set_ylabel('SBP') # (cm^-2 pixel^-2 s^-1)')
            ax.set_xlabel('radius (kpc)')
#            ymin,ymax=plt.ylim()
 #           ax.vlines(r200,ymin,ymax)
#        ax.set_xlim(10,2000)

            plt.savefig('pdf/'+name+'.pdf',dpi=10)
            plt.savefig(name+"/"+name+'.pdf')
            pp.savefig(fig)
        if mode == "2" or mode == "5":
            r_entropy=r_model/r_norm
            k_entropy=k_fit/k_norm
            k_entropy_up=k_fit_up/k_norm
            k_entropy_down=k_fit_down/k_norm
        if mode == "2":
            pylab.figure('k_local')
            plt.clf()
            plt.loglog(r_entropy,k_entropy)
            if k_entropy[-1] < 2:
                print(1)
                continue
            plt.text(0.1,3.5,'with clumping correction',color='lightskyblue')
            plt.fill_between(r_entropy,k_entropy_down,k_entropy_up,color='grey')
            plt.xlabel('Radius ($r/r_{200}$)')
            plt.ylabel('Entropy ($k/k(0.3r_{500})$)')
            plt.xlim(0.01,1.5)
            plt.ylim(0.04,4)
            plt.plot([0.01,0.3,1],[0.0237/3,1/3,4/3])
            k_entropy_uncorr=k_entropy*np.power(ne_cl_array/ne_array,-2/3)
            k_entropy_uncorr_up=k_entropy_up*np.power(ne_cl_array/ne_array,-2/3)
            k_entropy_uncorr_down=k_entropy_down*np.power(ne_cl_array/ne_array,-2/3)
            plt.loglog(r_entropy,k_entropy_uncorr,'g')
            plt.text(0.1,3,'without clumping correction',color='g')
            plt.fill_between(r_entropy,k_entropy_uncorr_down,k_entropy_uncorr_up,color='g')
            plt.title(name)
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
                plt.text(0.1,2.6,'other\'s result',color='k')
            plt.savefig(name+'/'+name+'_k_scaled.pdf')
            pylab.figure('k_global')
            info_file=name+'/'+name+'_suminfo.txt'
            for i in open(info_file,'r'):
                if the_re.match(r'^tcool',i):
                    tcool=float(i.split()[1])
            if tcool <= 7.7:
                color='b'
            elif tcool > 7.7  :
                color='r'
            plt.loglog(r_entropy, k_entropy,color,linewidth=0.5)
            plt.fill_between(r_entropy,k_entropy_down,k_entropy_up,alpha=0.3,color='grey')
        if mode == "5":
            if count5==0:
                plt.fill_between(r_entropy,k_entropy_down,k_entropy_up,alpha=0.5,color='grey')
                plt.loglog(r_entropy, k_entropy,label=name)
            if count5==1:
                plt.fill_between(r_entropy,k_entropy_down,k_entropy_up,alpha=0.5,color='green')
                plt.loglog(r_entropy, k_entropy,label=name)
            if count5==2:
                plt.fill_between(r_entropy,k_entropy_down,k_entropy_up,alpha=0.5,color='grey')
                plt.loglog(r_entropy, k_entropy,label=name,color='yellow')
            if count5==3:
                plt.fill_between(r_entropy,k_entropy_down,k_entropy_up,alpha=0.5,color='green')
                plt.loglog(r_entropy, k_entropy,label=name,color='red')
            plt.legend()
            count5=count5+1
            if count5==4:
                plt.xlabel(r'Radius (r/${\rm r_{200}}$)')
                plt.ylabel(r'Entropy (k/${\rm k(0.3r_{500})}$)')
                std_x=[0.1,1.5]
                std_y=[0.48*0.95,9.36*0.95]
                plt.xlim(0.07,1.5)
                plt.ylim(0.5,9)
                plt.loglog(std_x,std_y,'k',linewidth=1.2)
                plt.savefig('pdf/'+str(tmp_name)+'.pdf')
                plt.clf()
                tmp_name=tmp_name+1
                count5=0
        if mode == "3":
            pylab.figure('k_compare_all')
            k_ori=a0*np.power(r_model,gamma0)+k0
            plt.loglog(r_model/r200,k_ori/k_norm,'b',alpha=0.1)
            plt.loglog(r_model/r200,k_fit/k_norm,'r',alpha=0.3)
        if mode == "4":
            c_factor_array=np.power(1+r_model/r200,cp_p)*np.exp(r_model/r200*(cp_e))+0+cp_g0*np.exp(-(r_model/r200-cp_x0)*(r_model/r200-cp_x0)/cp_sigma)
            if c_factor_array.min()<0.95:
                print(name)
                for ii in range(len(c_factor_array)):
                    if c_factor_array[ii]<=1:
                        c_factor_array[ii]=1
#                continue
            c_factor_array=np.sqrt(c_factor_array)
            if t_ave <= 4:
                color='b'
            elif t_ave <=7 :
                color='g'
            else:
                color='r'
            plt.plot(r_model/r200,c_factor_array,color='grey',linewidth=0.4)
    if mode == "2":
        pylab.figure('k_global')
        plt.xlabel(r'Radius (r/${\rm r_{200}}$)')
        plt.ylabel(r'Entropy (k/${\rm k_{200}}$)')
#        plt.ylabel('Entropy')
        std_x=[0.1,1.5]
        std_y=[0.112,2.202]
        plt.xlim(0.07,1.5)
        plt.ylim(0.1,9)
        plt.loglog(std_x,std_y,'k',linewidth=1.2)
        plt.plot([0.08,0.11],[7,7],'k')
        plt.plot([0.08,0.11],[5,5],'b')
        plt.plot([0.08,0.11],[3.6,3.6],'r')
        plt.text(0.12,6.85,'Baseline entropy profile from Voit+05',fontsize=8)
        plt.text(0.12,4.9,'Scaled entropy porfiles of cool core clusters',fontsize=8,color='k')
        plt.text(0.12,3.5,'Scaled entropy profiles of none cool core clusters',fontsize=8,color='k')
        plt.savefig('entropy.pdf')
    if mode == "1":
        pp.close()
    if mode == "3":
        plt.xlabel('Radius ($r/r_{200}$)')
        plt.ylabel('Entropy ($k/k(0.3r_{200})$)')
        plt.savefig('entropy_totalvsgravity.pdf')
    if mode == "4":
        plt.xlim(0,1)
#        plt.ylim(0.5,3)
        plt.xlabel(r'Radius (r/${\rm r_{200}}$)')
        plt.ylabel(r'$(<\rho^2>/<\rho>^2)^{0.5}$')
        reference_factor=np.sqrt(np.power(1+r_model/r200,-3.7)*np.exp(r_model/r200*(3.7))+2.0*np.exp(-(r_model/r200-0.018)*(r_model/r200-0.018)/0.0002))
        plt.plot(r_model/r200,reference_factor,'b')
        ax=plt.gca()
#        ax.set_yticks([1,2,3,4,5,6,7,8,9,10])
#        ax.set_yticklabels(('1','2','3','4','','6','','','','','10'))
        ax.set_ylim(0.5,3)
        plt.plot([0.03,0.13],[2.75,2.75],'b')
        plt.plot([0.03,0.13],[2.50,2.50],'grey')
        plt.text(0.15,2.73,'Reference clumping factor profile from Vazza+13',color='b')
        plt.text(0.15,2.48,'Clumping factor profiles of sample clusters')
        plt.savefig('clumping_factor.pdf')
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
            k_mod=k_fit #temp
            k_ori_total.append(k_ori)
            k_mod_total.append(k_mod)


