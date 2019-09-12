#!/usr/bin/env python3
## $1 fitted parameter file 
## $2 ori temperature data file
## $3 sbprofile file to provide radiu for output entropy file
## $4 output entropy data file
## $5 output temperature fitting plot
import entropy_model


import matplotlib.pyplot as plt
import sys
import math
import numpy
import scipy
import matplotlib
import re
import types
import json
t_tot=[]
k_tot=[]
for f in open(sys.argv[1]):
    p=f.split()
    for i in range(len(p)):
        p[i]=float(p[i])
    t_array=[]
    k_array=[]
    for r in range(5000):
        if r == 0:
            continue
        t_array.append(entropy_model.main(r,*p))
        k_array.append(entropy_model.entropy(r,p[0],p[6],p[7],p[8],p[17]))
    k_tot.append(k_array)
    t_tot.append(t_array)
r500=p[0]
t_up=[]
t_middle=[]
t_down=[]
k_up=[]
k_middle=[]
k_down=[]
for r in range(5000):
    if r == 0:
        continue
    tmp=[]
    num=len(t_tot)
    for i in range(num):
        tmp.append(t_tot[i][r-1])
    num_down=int(num*0.16)
    num_middle=int(num*0.5)
    num_up=int(num*0.84)
    tmp=sorted(tmp)
    t_up.append(tmp[num_up])
    t_middle.append(tmp[num_middle])
    t_down.append(tmp[num_down])
    tmp=[]
    num=len(k_tot)
    for i in range(num):
        tmp.append(k_tot[i][r-1])
    num_down=int(num*0.16)
    num_middle=int(num*0.5)
    num_up=int(num*0.84)
    tmp=sorted(tmp)
    k_up.append(tmp[num_up])
    k_middle.append(tmp[num_middle])
    k_down.append(tmp[num_down])

#print(t_up)
#print(t_middle)
#print(t_down)
r_array=[]
re_array=[]
t_array=[]
te_array=[]
for i in open(sys.argv[2]):
    r,rer,t,te=i.split()
    r=float(r)
    rer=float(rer)
    t=float(t)
    te=float(te)
    r_array.append(r)
    re_array.append(rer)
    t_array.append(t)
    te_array.append(te)
r500_x=[]
r500_y=[]
tmp=1.5*r500
_1_5r500 = 1.5*r500
for i in range(2000):
    r500_x.append(tmp)
    r500_y.append(i/100)
fi=open(sys.argv[4],'w')
for i in open(sys.argv[3]):
    r=float(i.split()[0])
    ri=int(r)
    tmp_down=k_down[ri]-k_middle[ri]
    tmp_up=k_up[ri]-k_middle[ri]
    print('%.4e %.4e %.4e %4e'%(r,k_middle[ri],tmp_down,tmp_up),file=fi)
fi.close()
model_radius=list(range(4999))
re_array[0]=re_array[0]-1
plt.semilogx(t_up,color='grey')
plt.semilogx(t_middle,'b')
plt.semilogx(t_down,color='grey')
plt.semilogx(r500_x,r500_y,'m--')
plt.errorbar(r_array,t_array,xerr=re_array,yerr=te_array,color='k',linestyle='none')
dat ={
        "name": sys.argv[3],
        "radius": [r_array,re_array],
        "temperature": [t_array,te_array],
        "radius_model": model_radius,
        "temperature_model": [t_middle,t_down,t_up],
        "1_5r500": _1_5r500,
        }


tmp=sys.argv[5]+'_plt.json'
fi=open(tmp,'w')
json.dump(dat,fi,indent=2)
fi.close()
plt.xlim(1,5000)
plt.ylim(0.5,8.6)
plt.xlabel('Radius (kpc)')
plt.ylabel('Temperature (keV)')
#plt.text(0.8*r500,13,r'1.5r_{500}')
plt.fill_between(range(4999),t_down,t_up,color='grey')
#plt.title(sys.argv[3])
plt.text(0.03*r500,12,sys.argv[3] ,color="black",
        fontsize=15)
tmp=sys.argv[5]+'.pdf'
plt.savefig(tmp,dpi=100)
#plt.show()
