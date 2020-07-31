#!/usr/bin/env python3

import numpy as np
import sys
import re
from astropy.cosmology import FlatLambdaCDM
cosmo=FlatLambdaCDM(H0=71,Om0=0.27)

def calc(ra1,dec1,ra2,dec2,z,h0=71,Om0=0.27):
    cosmo=FlatLambdaCDM(H0=h0,Om0=Om0)
    da=cosmo.angular_diameter_distance(z).value #Mpc<F8>
    dist=((ra1-ra2)**2+(dec1-dec2)**2)**0.5*da*1000/180*3.1416 #kpc, approximation for small angle
    return dist

def main(name,xcent_file,bcg_file_flag=True):
    param_file='param_zzh_for_py.txt'
    for i in open(param_file):
        if re.match(r'z\s',i):
            z=float(i.split()[2])
    ra_bcg=input('input the BCG RA for '+name+' in format ##:##:## or in degree\n')
    dec_bcg=input('input the BCG DEC for '+name+' in format ##:##:## or in degree\n')
    if len(ra_bcg.split(':'))==3:
        ra_tmp=np.array(ra_bcg.split(':'),dtype=float)
        ra_bcg=ra_tmp[0]*15+ra_tmp[1]*15/60+ra_tmp[2]*15/3600
    else:
        ra_bcg=float(ra_bcg)
    if len(dec_bcg.split(':'))==3:
        dec_tmp=np.array(dec_bcg.split(':'),dtype=float)
        if dec_tmp[0]<0:
            dec_tmp[1]=-dec_tmp[1]
            dec_tmp[2]=-dec_tmp[2]
        dec_bcg=dec_tmp[0]+dec_tmp[1]/60+dec_tmp[2]/3600
    else:
        dec_bcg=float(dec_bcg)
    for i in open(xcent_file):
            xc=i
    ra_tmp=np.array(xc.split(',')[0].split('(')[1].split(':'),dtype=float)
    dec_tmp=np.array(xc.split(',')[1].split(')')[0].split(':'),dtype=float)
    ra_xc=ra_tmp[0]*15+ra_tmp[1]/4+ra_tmp[2]*15/3600
    if dec_tmp[0]<0:
        dec_tmp[1]=-dec_tmp[1]
        dec_tmp[2]=-dec_tmp[2]
    dec_xc=dec_tmp[0]+dec_tmp[1]/60+dec_tmp[2]/3600
    print(ra_xc,dec_xc,ra_bcg,dec_bcg)
    dist=calc(ra_bcg,dec_bcg,ra_xc,dec_xc,z)
    print(dist)
    if bcg_file_flag:
        fi=open(name+'_bcg_offset.txt','w')
        print('BCG_RA:',ra_bcg,file=fi)
        print('BCG_DEC:',dec_bcg,file=fi)
        print('XRAY_RA:',ra_xc,file=fi)
        print('XRAY_DEC:',dec_xc,file=fi)
        print('OFFSET(kpc):',dist,file=fi)

if __name__=='__main__':
    name=sys.argv[1]
    xcent_path=sys.argv[2]
    xcent_file=xcent_path+'/centroid_wcs.reg'
    main(name,xcent_file)

