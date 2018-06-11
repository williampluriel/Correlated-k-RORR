#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#"""
#Created on Wed Dec 20 14:30:17 2017
#
#@author: pluriel
#"""

import numpy as np
import pickle

def KCross_convert(n_species_pre,T_pre,P_pre,bande_sample,Nb,Nq,Ngres,adress,directory,data,bande,T,P):

    #Nb = bande_sample.size-1

    print 'Number of k bands : ', Nb
    print 'Number of Gauss points : ', Nq

    #toto = pickle.load(open('%s/data_caldas/K.db'%(adress)))
    #data = toto['xsecarr']
    #crossection = np.load('%s/data_caldas/%s_extrapoled.npy'%(adress,n_species_pre))
    #crossection = crossection[T,P,:]
    #bande_sample_HR = toto['wno']
    
    #data = pickle.load(open("%s/%s/%s/%s_T%i_P%.4e.xsec.pickle"%(adress,data,n_species_pre,n_species_pre,T_pre,P_pre/(1.0e5))))
    #bande_sample_HR = data[:,0] # [cm-1]
    #crossection = data[:,1] # [cm^2/molec]

    data = np.load("%s/cross_mixt_H2O_CO2_100_1e-5.npy"%(adress))
    bande_sample_HR = data[:,0] # [cm-1]
    crossection = data[:,1] # [cm^2/molec]

    # Here we eliminate bands where the wavelenght is >= or <= than bands in low resolution
    #zone, = np.where((bande_sample_HR >= bande_sample[0])*(bande_sample_HR <= bande_sample[bande_sample.size - 1]))

    #ind_min = np.amin(zone)
    #ind_max = np.amax(zone)

    #bande_sample_HR = bande_sample_HR[ind_min:ind_max+1]
    
    NbHR = bande_sample_HR.size
    print 'Size of the high resolution source : ', NbHR
    
    zone, = np.where(crossection[:]!=0)
    zone2, = np.where(crossection[:]==0)
    val_min = np.amin(crossection[zone])
    crossection[zone2] = np.ones(np.size(zone2))*val_min*1.0e-10

    # Creating Gauss coefficient 
  
    y = np.zeros(Nq+1)
    w = np.zeros(Nq+1)

#    w = np.array([ 0.008807 ,   0.02030071 , 0.03133602 , 0.04163837 , 0.05096506 , 0.05909727,
#      0.06584432 , 0.07104805 , 0.07458649 , 0.07637669 , 0.07637669 , 0.07458649,
#      0.07104805,  0.06584432,  0.05909727,  0.05096506 , 0.04163837 , 0.03133602,
#      0.02030071 , 0.008807])

    for i in range(Nq+1):

        x = i/float(Nq+2)

        y[i] = np.sin(np.pi*x)*np.exp(-12*x)

    for i in range(Nq):

        w[i] = y[i+1]

    w_tot = np.nansum(w)
    w = w/w_tot

    gw = np.zeros(Nq+1)

    for i in range(Nq):

        if i == 0 :

            gw[0] = w[0]

        else :

            gw[i] = gw[i-1] + w[i]

    print 'Gauss coefficients : '
    print w

    k_dist = np.zeros((Nb,Nq+1))
            
    k_spec = np.zeros(NbHR)

    print "Considered temperature : ", T_pre
    print "Considered pressure : ", P_pre
               
    Cn = P_pre/(8.314*T_pre)*6.022e+23 ## [molecule/m^3] absorbant density
    k_spec = 0.0001*Cn*crossection ## [m-1] opacity

    overflow = 0

    ind = 0

    # Main loop on the bands
    for i_bande in range(Nb):
        #print "bande : ",i_bande+1
        
        # if : permet de completer les bandes ou il n y a pas de donnees d absorption par 0
        if  ((bande_sample[i_bande] < bande_sample_HR[-1])*
	    (bande_sample[i_bande] > bande_sample_HR[0])) :

            wh, = np.where((bande_sample_HR[ind:] >= bande_sample[i_bande])*(bande_sample_HR[ind:] < bande_sample[i_bande+1]))
            nu_bande = bande_sample_HR[wh+ind] # decoupage en bande
            nu_size = nu_bande.size
            k_bande = k_spec[wh+ind] # opacity

            ind += wh.size

            nuw_bande = np.zeros(nu_size)

            for nu in range(1,nu_size-1) :

                nuw_bande[nu] = (nu_bande[nu+1] - nu_bande[nu-1])/2.

            nuw_bande[0] = nuw_bande[1]
            nuw_bande[nu_size-1] = nuw_bande[nu_size - 2]
            nuw_tot = np.nansum(nuw_bande)

            k_min = np.amin(k_bande)

#            if k_min < 1.0e-40 :
#
#                k_min = 1.0e-40

            k_max = np.amax(k_bande)

            lkmax = np.log10(k_max)
            lkmin = np.log10(k_min)

            dlogK = (lkmax - lkmin)/float(Ngres - 1)

            g_fine = np.zeros(Ngres)
            k_levs = np.zeros(Ngres)

            # Loop on Ngres
            for N in range(Ngres):

                k_lev = 10**(lkmin + (N+1)*dlogK)
                k_levs[N] = k_lev

                #print "k_lev = ", k_lev

                addline, = np.where(k_bande < k_lev)
                g_fine[N] = np.nansum(nuw_bande[addline])/float(nuw_tot)

            k_dist_temp = np.zeros(Nq+1)

            ik = 0

            # Loop on the Gauss points
            for i_q in range(Nq):

                if i_q == 0:
                    wh, = np.where((g_fine < gw[i_q])*(g_fine != 1.)) # wh == weight

                    if wh.size == Ngres:
                        overflow += 1
                        print('g-inversion overflow')
                    else :
                        if wh.size != 0 : k_dist_temp[i_q] = np.nansum(k_levs[wh])/(float(wh.size))
                        else : k_dist_temp[i_q] = k_min
                else :
                    ik += wh.size
                    wh, = np.where((g_fine[ik:] < gw[i_q])*(g_fine[ik:] != 1.))

                    if wh.size == Ngres:
                        overflow += 1
                        print('g-inversion overflow')
                    else :
                        if wh.size != 0 : k_dist_temp[i_q] = np.nansum(k_levs[wh+ik])/(float(wh.size))
                        else : k_dist_temp[i_q] = k_dist_temp[i_q - 1] ## [m-1] homogenous to an opacity
                
        else :
            k_dist_temp = np.zeros(Nq+1)
            for i_q in range(Nq):
                k_dist_temp[i_q] = 1.0e-55 # carefull, random value

        ## k_dist [cm^2/molecule] homogenous to a cross section :
        #k_dist[i_bande,:] = k_dist_temp # [m-1] absorption coeff
        k_dist[i_bande,:] = 1.3807e-19*k_dist_temp*T_pre/P_pre # [cm^2/molecule]
        #print "calcul ok"
    print(overflow)

    #np.save("%s/k_corr_%s_%i_%.4e.npy"%(directory,n_species_pre,T_pre,P_pre/(1.0e5)),k_dist)

    return k_dist #[cm^2/molecule]

## dec_first = 3
## dec_final = 4


#--------------------------------------------------------------------------------------------------------#
        
def k_corr_data_write(k_corr_data,k_correlated_file,dec_first,dec_final,P_dim,T_dim,dim_bande,dim_gauss):

    data = open(k_correlated_file,'w')
    first = ''
    for df in range(dec_first):
        first += ' '
    final = ''
    for dfl in range(dec_final):
        final += ' '

    #bar = ProgressBar(dim_gauss*dim_bande,'Kcorr record')

    for m in range(dim_gauss):
        for l in range(dim_bande):
                for j in range(P_dim):
                    for i in range(T_dim):

                        if i == 0 and j == 0 and l == 0 and m == 0 :
                            data.write(first)

                        data.write('%.16E'%(k_corr_data[i,j,l,m]))
                        data.write(final)

#---------------------------------------------------------------------------------------------------------#