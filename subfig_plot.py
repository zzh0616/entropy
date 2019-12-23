#!/usr/bin/env python3

#$1 : mode 1t for temperature compare vs chadnra
#          1s for surface brightness compare vs chandra
#          2 for entropy compare vs others
#$2 : list for path of the sources to be drawn
#$3,$4, row and column for the subfigures

import matplotlib.pyplot as plt
import sys
import plot_all
import matplotlib
import numpy as np
mode=sys.argv[1]
x=int(sys.argv[3])
y=int(sys.argv[4])
if x*y != len(open(sys.argv[2]).readlines()):
    print('lines not match, please check')
    exit
matplotlib.rcParams['xtick.direction']='in'
fig, axarr=plt.subplots(x,y,sharex='all',sharey='all')

fig.subplots_adjust(hspace=0,wspace=0,bottom=0.15,top=None)
index=0
for f in open(sys.argv[2]):
    index=index+1
#    if mode=='2':
    ax_this=plt.subplot(x,y,index)
    if np.mod(index-1,y) != 0:
        plt.setp([ax_this.get_yticklabels()],visible=False)
    if np.floor((index-1)/y) != x-1:
        plt.setp([ax_this.get_xticklabels()],visible=False)
    fi=open('_tmp.tmp','w')
    print(f,end='',file=fi)
    fi.close()
    plot_all.main('3','_tmp.tmp')
    if mode=='1s':
        ax_this.set_xticks([10,100,1000])
        ax_this.set_xticklabels(['10','$10^2$','$10^3$'],fontsize=6)
    if mode=='1t':
        ax_this.set_xscale("log")
        ax_this.set_yscale("log")
        ax_this.set_xticks([10,100,1000])
        ax_this.set_xticklabels(['10','$10^2$','$10^3$'],fontsize=6)

#        plt.xticklabels=None
#plt.setp([ax.get_xticklabels() for ax in axarr[:-1,:].reshape(-1)],visiable=False)
#plt.setp([ax.get_yticklabels() for ax in axarr[:, 1:].reshape(-1)],visiable=False)
if mode == '2':
    fig.text(0.45,0.08,r'Radius (${\rm r/r_{200}}$)')
    fig.text(0.02,0.4,r'Entropy (${\rm K/K_{200}}$)',rotation=90)
    blue_line=matplotlib.lines.Line2D([],[])
    orange_line=matplotlib.lines.Line2D([],[],color='orange')
    black_line=matplotlib.lines.Line2D([],[],color='k')
    green_line=matplotlib.lines.Line2D([],[],color='green',linestyle='--',alpha=0.5)
    green_patch=matplotlib.patches.Patch(color='palegreen',alpha=1,lw=2)
#    test=matplotlib.patches.Rectangle([0,0],1,1,color='green',lw=2,alpha=0.7)
    grey_patch=matplotlib.patches.Patch(color='grey',lw=2)
    z=np.random.randn(10)
    black_cross=matplotlib.lines.Line2D([],[],color='k',marker='+',markersize=5,markerfacecolor='k',linestyle='none')
    plt.legend(handles=[(grey_patch,blue_line),black_line,green_line],labels=['fitted entropy profile','baseline entropy','outermost radius of FoV'],fontsize=6,loc='upper center',bbox_to_anchor=(-1.5,-0.6),ncol=6)
    fig.savefig('entropy_subfig.pdf')
if mode == '1t':
    fig.text(0.45,0.03,'Radius (kpc)')
    fig.text(0.02,0.62,'Temperature (keV)',rotation=90)
    fig.savefig('temperature_compare_subfig.pdf')
if mode == '1s':
    fig.text(0.45,0.04,'Radius (kpc)')
    fig.text(0.02,0.73,r'Surface Brightness (${\rm cm^{-2}\ pixel^{-2}\ s^{-1}}$)',rotation=90)
    fig.savefig('sbp_compare_subfig.pdf')

