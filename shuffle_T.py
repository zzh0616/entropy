#!/usr/bin/python

import sys
import scipy

output_file=open(sys.argv[2],'w')
for i in open(sys.argv[1]):
    r,re,c,s=i.strip().split()
    c=float(c)
    s=float(s)
    
    c1=-1

    while c1<=0:
        c1=scipy.random.normal(0,1)*s+c
    
    output_file.write("%s\t%s\t%s\t%s\n"%(r,re,c1,s))

