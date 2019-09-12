#!/usr/bin/env python3

import sys
import numpy as np
import json
import os

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

