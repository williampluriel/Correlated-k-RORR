#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
import pickle
from Kcorr_convert_library import *

adress = "/data1/pluriel/make_corrk/xsec_combined"
directory = "/data1/pluriel/make_corrk/xsec_combined/kcorr"
data = "EXOMOL"
bande = "bands"
typeband = "VI" # several choices TBD
typegauss = "gj1214b" # several choices TBD

#toto = pickle.load(open('%s/data_caldas/NA.db'%(adress)))
#T_sample_HR = toto['t']
T_sample_HR = np.load("%s/%s/T_sample.npy"%(adress,data))
print 'High resolution k temperature sample :'
print T_sample_HR
gauss_sample = np.load("%s/%s/gauss_sample_HR.npy"%(adress,bande))
print 'Gauss coefficients :'
print gauss_sample
#bande_sample = np.load("%s/%s/%s/bande_sample.npy"%(adress,typegauss,typeband))
bande_sample = np.load("%s/%s/gj1214b/bande_sample_%s.npy"%(adress,bande,typeband)) #[cm-1]
print 'Bande sample for GCM simulation :'
for i in range(bande_sample.size - 1) :
    print '   %.8E cm-1   ----------------->  %.8E cm-1 '%(bande_sample[i],bande_sample[i+1])

Nb = bande_sample.size-1
Nq = gauss_sample.size-1
Ngres = 200
NT = T_sample_HR.size
n_species_abs = np.array(["H2O","CO2"])
#n_species_abs = np.array(["K","NA"])
Nspe = n_species_abs.size

# Loop on the species
for ispe in range(0,1):
    n_species_pre = n_species_abs[ispe]
    print n_species_pre
    #P_sample_HR_bar = toto['p']
    P_sample_HR_bar = np.load("%s/%s/P_sample.npy" % (adress,data))  ## [bar]
    P_sample_HR = P_sample_HR_bar*1.0e5  ## [Pa]
    NP = P_sample_HR.size
    k_corr = np.zeros((NT, NP, Nb, Nq+1), dtype=np.float64)

    # Loop on the temperatures
    for iT in range(0,1):
        T_pre = T_sample_HR[iT] # [K]
	T = iT

        # Loop on the pressures
        for iP in range(0,1):
            P_pre = P_sample_HR[iP] # [Pa]
	    P = iP
            k_corr[iT,iP,:,:] = KCross_convert(n_species_pre,T_pre,P_pre,bande_sample,Nb,Nq,Ngres,adress,directory,data,bande,T,P) # [cm^2/molecule]
            #print "successfully saved"

    np.save("%s/k_corr_H2O_CO2_test.npy"%(adress),k_corr)
    #np.save("%s/%s/k_corr_%s_%s.npy"%(adress,typegauss,n_species_pre,typeband),k_corr)
    print n_species_abs[ispe], " successfully saved"
