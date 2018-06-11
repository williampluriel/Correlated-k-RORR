#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import pickle

adress = "/data1/pluriel/make_corrk/xsec_combined"
directory = "/data1/pluriel/make_corrk/xsec_combined/kcorr"
data = "EXOMOL"
bande = "bands"
typeband = "IR" # several choices TBD
typegauss = "gauss_compo0" # several choices TBD
typeplanet = "earth" # planetary name

k_corr = np.load('%s/k_corr_RORR_gcm_%s_%s_%s.npy'%(directory,typeband,typegauss,typeplanet))
k_corr = k_corr.flat[:]

k_rebin = np.zeros((1,k_corr.size))

k_rebin[0,:] = k_corr[:]

print np.shape(k_corr)

np.savetxt('%s/%s/corrk_gcm_%s.dat'%(adress,typeplanet,typeband),k_rebin)

print "data saved"