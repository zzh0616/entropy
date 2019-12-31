#!/usr/bin/env python3

import sys
import numpy as np
import json
import os
from astropy.cosmology import FlatLambdaCDM
cosmo=FlatLambdaCDM(H0=71,Om0=0.27,Tcmb0=2.73)
pi=3.1416
def main():
    listfile=sys.argv[1]
    csvfile=sys.argv[2]
    filename_array=[]
    std_name_array=[]
    dist=-1
    for i in open(listfile):
        filename_array.append(i.split(',')[0])
        std_name_array.append(i.split(',')[1][0:-1])
    for i in open(csvfile):
        stdname=i.split(',')[0]
        if ';' in stdname:
            stdname=stdname.split(';')[0]
#        print(stdname)
        if stdname not in std_name_array:
            print('not found',stdname)
        elif stdname == "Abell 1246" or stdname == "Abell 2163" or stdname =="ZWCL 3146":
            print(stdname)
            j=std_name_array.index(stdname)
            nameuse=filename_array[j]+'/bcg_dist.txt'
            fi=open(nameuse,'w')
            dist=-1
            print(dist,file=fi)
            fi.close()
        else:
            bcg_ra=float(i.split(',')[6])
            bcg_dec=float(i.split(',')[7])
            cl_ra=float(i.split(',')[3])
            cl_dec=float(i.split(',')[4])
            z=float(i.split(',')[5])
            da=cosmo.angular_diameter_distance(z).value
            arc=np.sqrt((bcg_ra-cl_ra)**2+(bcg_dec-cl_dec)**2)/180*pi
            dist=arc*da
            j=std_name_array.index(stdname)
            nameuse=filename_array[j]+'/'+'bcg_dist.txt'
            print(stdname,nameuse)
            fi=open(nameuse,'w')
            print(dist,file=fi)
            fi.close()
if __name__=='__main__':
    main()


'''
fi=open(sys.argv[1],'r')
dat=json.load(fi)
name=dat['Source Name'].upper()
obsid=dat['Obs. ID']
z=dat['redshift']
if 'lyt' in sys.argv[1]:
    ra=dat['R. A.']
    dec=dat['Dec.']
    tave=dat['T(0.2-0.5 R500)']
    terr=dat['T_err(0.2-0.5 R500)']
#    std_name=dat['Unified Name']
    ra=ra.replace('h',':')
    ra=ra.replace('m',':')
    ra=ra.replace('s','')
    dec=dec.replace('d',':')
    dec=dec.replace('m',':')
    dec=dec.replace('s','')
else:
    ra=dat['RA']
    dec=dat['DEC']
    if len(sys.argv) >= 3:
        a=open(sys.argv[2],'r')
        for i in a:
            tave=round(np.float(i.split()[0]),3)
            terr=round(np.float(i.split()[2])-np.float(i.split()[1]),3)
    else:
        tave=''
        terr=''
    std_name=''


print(name,',',obsid,',',ra,',',dec,',',z,',',tave,',',terr)
'''
