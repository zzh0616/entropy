#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
import sys
p=np.load('p_all.npy')
array_total=[]
name=sys.argv[1]

fi=open(name+'_result_from_p_all.csv','w')

print('n2,',str(p[0].mean())+',',p[0].std(),file=fi)
print('rs,',str(p[1].mean())+',',p[1].std(),file=fi)
print('a0,',str(p[2].mean())+',',p[2].std(),file=fi)
print('gamma0,',str(p[3].mean())+',',p[3].std(),file=fi)
print('k0,',str(p[4].mean())+',',p[4].std(),file=fi)
print('n3,',str(p[5].mean())+',',p[5].std(),file=fi)
print('rho,',str(p[6].mean())+',',p[6].std(),file=fi)
print('ne0,',str(p[7].mean())+',',p[7].std(),file=fi)
print('T0,','0,','0',file=fi)
print('sbp_c,',str(p[9].mean())+',',p[9].std(),file=fi)
print('delta,',str(p[10].mean())+',',p[10].std(),file=fi)
print('delta2,',str(p[11].mean())+',',p[11].std(),file=fi)
print('nth_a,',str(p[12].mean())+',',p[12].std(),file=fi)
print('nth_b,',str(p[13].mean())+',',p[13].std(),file=fi)
print('nth_gamma,',str(p[14].mean())+',',p[14].std(),file=fi)
print('kmod_a,',str(p[15].mean())+',',p[15].std(),file=fi)
print('kmod_b,',str(p[16].mean())+',',p[16].std(),file=fi)
print('kmod_c,',str(p[17].mean())+',',p[17].std(),file=fi)
print('kmod_k0,',str(p[18].mean())+',',p[18].std(),file=fi)
print('cp_p,',str(p[19].mean())+',',p[19].std(),file=fi)
print('cp_e,',str(p[20].mean())+',',p[20].std(),file=fi)
print('cp_g0,',str(p[21].mean())+',',p[21].std(),file=fi)
print('cp_x0,',str(p[22].mean())+',',p[22].std(),file=fi)
print('cp_sigma,',str(p[23].mean())+',',p[23].std(),file=fi)

fi.close()
