#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
import pickle
from Kcorr_convert_library import *

adress = "/data1/pluriel/make_corrk/xsec_combined"
directory = "/data1/pluriel/make_corrk/xsec_combined/kcorr"
data = "EXOMOL"
bande = "bands"
typeband = "IR" # several choices TBD
typegauss = "gj1214b" # several choices TBD

T_sample_HR = np.load("%s/%s/T_sample.npy"%(adress,data))
print 'High resolution k temperature sample :'
print T_sample_HR
gauss_sample = np.array("%s/%s/gauss_sample_HR.npy"%(adress,data))
print 'Gauss coefficients :'
print gauss_sample
bande_sample = np.load("%s/%s/bande_sample_HR.npy"%(adress,bande))

T_pre = 200.
P_pre = 1e-5
n_species_pre = np.array(["H2O","CO2"])

data1 = pickle.load(open("%s/%s/%s/%s_T%i_P%.4e.xsec.pickle"%(adress,data,n_species_pre[0],n_species_pre[0],T_pre,P_pre)))
bande_sample_HR1 = data1[:,0] # [cm-1]
crossection1 = data1[:,1] # [cm^2/molec]
data2 = pickle.load(open("%s/%s/%s/%s_T%i_P%.4e.xsec.pickle"%(adress,data,n_species_pre[1],n_species_pre[1],T_pre,P_pre)))
crossection2 = data2[:,1] # [cm^2/molec]
bande_sample_HR2 = data2[:,0] # [cm-1]


print np.shape(bande_sample_HR1),bande_sample_HR1[0:11466971]
print np.shape(bande_sample_HR2),bande_sample_HR2


#crossection = crossection1[0:11466971] + crossection2

#np.save("%s/%s/%s%s_T%i_P%.4e.npy"%(adress,data,n_species_pre[0],n_species_pre[1],T_pre,P_pre),crossection)
