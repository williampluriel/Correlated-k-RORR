#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 11:40:14 2018

@author: pluriel
"""

import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy.optimize import curve_fit

nb = 15
nT = 1
nP = 0

k_rebin = np.load('/home/pluriel/data_dir/make_corrk_save/xsec_combined/k_corr_conv_rebin.npy')
#k_rebin = k_rebin[:,:,:,::-1]

B = np.load('/home/pluriel/data_dir/make_corrk_save/xsec_combined/k_corr_conv.npy')
B_arr = np.load("/home/pluriel/data_dir/make_corrk_save/xsec_combined/bande_sample_HR.npy")
G_arr2 = np.load('/home/pluriel/data_dir/make_corrk_save/xsec_combined/weights_conv.npy')
G_arr = np.array([0.252362754191094,0.255197877582394,0.191542751879209,     
                         0.126425173849644,7.735535321334741E-002,4.489907382534229E-002,
                         2.501283456640186E-002,1.345856933556737E-002,7.016455573275246E-003,
                         3.547640305300775E-003,1.737978015237359E-003,8.223551811163857E-004,
                         3.734677496729497E-004,1.608929328184595E-004,6.425542517147573E-005,
                         2.256637440673820E-005,0.000000 ])

sum_garr = 1-np.cumsum(G_arr)
sum_garr2 = 1-np.cumsum(G_arr2[nT,nP,nb,:])

k_dist = np.load('/home/pluriel/data_dir/make_corrk_save/xsec_combined/output_CO2/k_corr_CO2_200_1.0000e-05.npy')
k_dist2 = np.load('/home/pluriel/data_dir/make_corrk_save/xsec_combined/output_H2O/k_corr_H2O_200_1.0000e-05.npy')

print B[0,0,7,273:289]
print np.shape(k_rebin)
print k_rebin[0,0,7,-1]

print sum_garr
print sum_garr2
print np.shape(B)
print "Convolution :"
print B[0,0,nb,:]
print "k_corr H2O :"
print k_dist2[nb,:]
print "k_corr CO2 :"
print k_dist[nb,:]

#plt.plot(sum_garr2[:],B[0,0,nb,:],'m--')
#plt.plot(sum_garr[:-1],k_dist[nb,1:],'r')
#plt.plot(sum_garr[:-1],k_dist2[nb,1:],'b')

plt.hlines(k_dist[nb,1:], xmin=sum_garr[:-1],xmax=sum_garr[1:],color='g',linestyles='--',label='k-corr CO$_2$',linewidth=2)
plt.vlines(sum_garr[:-1], ymin=k_dist[nb,1:],ymax=k_dist[nb,:-1],color='g',linestyles='--',linewidth=2)

plt.hlines(k_dist2[nb,1:], xmin=sum_garr[:-1],xmax=sum_garr[1:],color='b',label='k-corr H$_2$O',linewidth=2)
plt.vlines(sum_garr[:], ymin=k_dist2[nb,1:],ymax=k_dist2[nb,:-1],color='b',linewidth=2)

plt.hlines(B[nT,nP,nb,1:], xmin=sum_garr2[:-1],xmax=sum_garr2[1:],color='k', linestyles='-.',label='k-corr convolution',linewidth=2)
plt.vlines(sum_garr2[:], ymin=B[nT,nP,nb,1:],ymax=B[nT,nP,nb,:-1],color='k',linestyles="-.",linewidth=2)

plt.hlines(k_rebin[nT,nP,nb,1:], xmin=sum_garr[:-1],xmax=sum_garr[1:],color='m', linestyles='-.',label='k-corr convolution rebin',linewidth=2)
plt.vlines(sum_garr[:], ymin=k_rebin[nT,nP,nb,1:],ymax=k_rebin[nT,nP,nb,:-1],color='m',linestyles="-.",linewidth=2)

plt.xscale('log')
plt.yscale('log')
plt.xlim([1e-16,1])
plt.ylim([1e-40,1e-15])
plt.xlabel('g',fontsize=24)
plt.ylabel('Absorption coeff',fontsize=24)
plt.title("Atmospheric inventory: H$_2$O(20%) CO$_2$(80%); T = 200 K; P = 1 Pa",fontsize=24)
plt.gca().tick_params(axis='both',which='major',labelsize=24)
plt.legend(loc="lower left",fontsize=24)
plt.show()

#plt.hlines(A[0,0,:,7], xmin=B_arr[:-1],xmax=B_arr[1:])
#plt.vlines(B_arr[1:], ymin=A[0,0,1:,7],ymax=A[0,0,:-1,7])
#
#plt.hlines(B[:,7], xmin=B_arr[:-1],xmax=B_arr[1:])
#plt.vlines(B_arr[1:], ymin=B[1:,7],ymax=B[:-1,7])
#
#plt.xscale('log')
#plt.yscale('log')
#plt.legend(loc='best')
#plt.show()