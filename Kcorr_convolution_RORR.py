#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# author : William PLURIEL
# date : 22/01/2018

# --------------------------------------------------------------------------
# Treatment of overlapping gaseous absorption with the correlated-k method
# Random overlap method with resorting and rebinning (RORR) ; Admundsen 2017
# --------------------------------------------------------------------------

import numpy as np
import pickle

adress = "/data1/pluriel/make_corrk/xsec_combined"
directory = "/data1/pluriel/make_corrk/xsec_combined/kcorr"
data = "EXOMOL"
bande = "bands"
typeband = "VI" # several choices TBD
typegauss = "gauss_compo0" # several choices TBD
typeplanet = "gj1214b" # planetary name

Composition = np.array(1) # 0 : read a file with the compositions (T,P,nSpecies)
                          # 1 : imposes a fixed composition whatever is T and P
		          # 2 : one tracer's specie and fixed composition for the other species

T_sample_HR = np.load("%s/%s/T_sample.npy"%(adress,data)) # [K]
print "k-correlated temperature sample: "
print T_sample_HR
P_sample_HR = np.load("%s/%s/P_sample.npy"%(adress,data)) # [bar]
print "k-correlated pressure sample: "
print P_sample_HR
bande_sample = np.load("%s/%s/%s/bande_sample_%s.npy"%(adress,bande,typeplanet,typeband)) #[cm-1]
print 'Bande sample for GCM simulation :'
print bande_sample
gauss_sample = np.load("%s/%s/gauss_sample_HR.npy"%(adress,bande))
print "Gauss coefficients: "
print gauss_sample
Q_ratio = np.loadtxt('/home/pluriel/datagcm/corrk_data/gj1214b/Q_var.dat') # mixing ratio
T_compo = np.loadtxt("%s/earth/T.dat"%(adress))
print "Composition temperature sample: "
print T_compo # [K]
P_compo = np.loadtxt("%s/earth/p.dat"%(adress)) # [mbar]
P_compo = (10**(P_compo[:]))/1.e3 # [bar]
print "Composition pressure sample: "
print P_compo


#n_species = np.array(["H2","He","H2O","CH4","N2","NH3","CO","CO2"])
#n_species_abs = np.array(["H2O","CH4","NH3","CO","CO2"]) #,"K","NA"])
n_species = np.array(["N2","H2O","CO2"])
n_species_abs = np.array(["H2O","CO2"])

NT = T_sample_HR.size
print "NT =", NT
NP = P_sample_HR.size
print "NP =", NP
Nb = bande_sample.size
print "Nb =", Nb-1
Nq = gauss_sample.size
print "Nq =", Nq
Nspe = n_species_abs.size
print "Nspe =", Nspe
Ngres = 200
print "Ngres =", Ngres
NQ = Q_ratio.size
print "NQ =", NQ
NT_comp = T_compo.size
print "NT_comp =", NT_comp
NP_comp = P_compo.size
print "NP_comp =", NP_comp

# Choice of the composition
if Composition == 0 :
    compo = np.zeros((Nspe, NP_comp, NT_comp, NQ), dtype=np.float64)
    for iQ in range(NQ):
        compo[:,:,:,iQ] = np.load('%s/compo/compo_earth.npy'%(adress))
    
elif Composition == 1 :
    compo = np.zeros((Nspe,NP,NT,NQ), dtype=np.float64)
    compo_temp = np.array([0.5,0.5]) # manually enter the desired values [spe1,spe2,...]
    for ispe in range(Nspe):
        for iQ in range(NQ):
	    compo[ispe,:,:,iQ] = compo_temp[ispe]
	
elif Composition == 2 :
    # Put tracer's specie in last position
    compo = np.zeros((Nspe,NP,NT,NQ), dtype=np.float64)
    compo_spe_var = Q_ratio # Carreful, you must create a file Q_var.dat : 
                            # 1 column with mixing ratio of tracer's specie
    compo_temp = np.array([0.0035,0.965]) # manually enter the desired values [spe1,spe2,...]
    for ispe in range(Nspe-1):
         for iQ in range(NQ):
	     compo[ispe,:,:,iQ] = compo_temp[ispe]
	     compo[-1,:,:,iQ] = compo_spe_var[iQ]
    
else :
    print "Error during composition data load, check composition choice."

# Inialisation of all table
k_corr = np.zeros((Nspe, NT, NP, NQ, Nb-1, Nq), dtype=np.float64)
k_corr_temp = np.zeros((NT,NP,NQ,Nb-1, Nq, Nq), dtype=np.float64)
k_corr_tot = np.zeros((NT,NP,NQ,Nb-1, Nq*Nq), dtype=np.float64)
k_corr_rebin_gcm = np.zeros((Nq,Nb-1,NQ,NP,NT), dtype=np.float64)
k_corr_rebin = np.zeros((NT,NP,NQ,Nb-1,Nq), dtype=np.float64)
w = np.zeros((Nq,Nq),dtype=np.float64)
w_temp = np.zeros((Nq*Nq),dtype=np.float64)
w_tot = np.zeros((NT,NP,NQ,Nb-1, Nq*Nq), dtype=np.float64)
w_sum = np.zeros((NT,NP,NQ,Nb-1, Nq*Nq), dtype=np.float64)
w_rebin = np.zeros(Nq,dtype=np.float64)
composition = np.zeros((Nspe,NP,NT,NQ), dtype=np.float64)

# Read correlated-k data
for ispe in range(Nspe):
    n_spe = n_species_abs[ispe]
    print n_spe
    for iQ in range(NQ):
        #k_corr[ispe,:,:,:,:] = np.load("%s/gj1214b/output_%s/k_corr_%s.npy"%(adress,n_spe,n_spe))
        k_corr[ispe,:,:,iQ,:,:] = np.load("%s/%s/k_corr_%s_%s.npy"%(adress,typeplanet,n_spe,typeband))

k_corr_rebin[:,:,:,:,:] = k_corr[0,:,:,:,:,:]

# Beginning of the RORR convolution method
for ispe in range(1,Nspe): # begin at 1 because the first step must mixt 2 gases
    
    print "----------------------" 
    print "----------------------" 
    print "previous gas mixture",n_species_abs[0:ispe]
    print "added gas now", n_species_abs[ispe]
   
    for ib in range(Nb-1):
        print "---------"
        print "bande ",ib+1
        print "---------"
        for iT in range(0,1):
            
            T_ref = T_sample_HR[iT] # [K]
	
            if T_ref > T_compo[0] :

                wh, = np.where(T_compo > T_ref)

                if wh.size != 0 :

                    i_Tub = wh[0]
                    i_Tdb = i_Tub - 1
                    T_ub = T_compo[i_Tub]
                    T_db = T_compo[i_Tdb]
                    coeff_1 = (T_ref - T_db)/(T_ub - T_db)
                    coeff_1b = 1 - coeff_1

                else :

                    i_Tub, i_Tdb = T_compo.size - 1, T_compo.size - 1
                    coeff_1 = 0.
                    coeff_1b = 1.

            else :

                i_Tub, i_Tdb = 0, 0
                coeff_1 = 0.
                coeff_1b = 1.

            for iP in range(0,1) :

                P_ref = P_sample_HR[iP] # [bar]

                if P_ref > P_compo[0] :

                    wh, = np.where(P_compo > P_ref)

                    if wh.size != 0 :

                        i_Pub = wh[0]
                        i_Pdb = i_Pub - 1
                        P_ub = P_compo[i_Pub]
                        P_db = P_compo[i_Pdb]
                        coeff_2 = (P_ref - P_db)/(P_ub - P_db)
                        coeff_2b = 1 - coeff_2

                    else :

                        i_Pub, i_Pdb = P_compo.size - 1, P_compo.size - 1
                        coeff_2 = 0.
                        coeff_2b = 1.

                else :

                    i_Pub, i_Pdb = 0, 0
                    coeff_2 = 0.
                    coeff_2b = 1.
                
		for iQ in range(NQ):
		
                    compo_u_u = compo[ispe,i_Pub,i_Tub,iQ]
                    compo_u_d = compo[ispe,i_Pub,i_Tdb,iQ]
                    compo_d_u = compo[ispe,i_Pdb,i_Tub,iQ]
                    compo_d_d = compo[ispe,i_Pdb,i_Tdb,iQ]

                    c12, c1b2, c12b, c1b2b = coeff_1*coeff_2, coeff_1b*coeff_2, coeff_1*coeff_2b, coeff_1b*coeff_2b
		
		    compo1 = compo_u_u*c12 + compo_u_d*c1b2
                    compo2 = compo_d_u*c12b + compo_d_d*c1b2b

                    composition[ispe,iP,iT,iQ] = compo1 + compo2
                        
                    if Composition == 1 : composition[0,:,:,:] = compo[0]
		    
		    P_pre = P_sample_HR[iP]*1.0e5  # [Pa]
    
                    for iq1 in range(Nq):
                        for iq2 in range(Nq):
                            w[iq1,iq2] = gauss_sample[iq1]*gauss_sample[iq2] # weights
                            w_temp = w.flat[:] # reshape on 1 ligne
			
                            # Je me base sur l'équation (11) de Amundsen 2017, qui ne prend en compte que les espèces absorbantes.
		    	    # Or ,ici, je veux que kcorr comptabilise aussi les espèces non-absorpbantes, d'où une multiplication
			    # de k_corr_temp par la somme des mixing ratio (cf démonstration dans mon cahier de thèse).
                            k_corr_temp[iT,iP,iQ,ib,iq1,iq2] = (k_corr_rebin[iT,iP,iQ,ib,iq1]*np.sum(composition[0:ispe,iP,iT,iQ]) + 
                                                              k_corr[ispe,iT,iP,iQ,ib,iq2]*composition[ispe,iP,iT,iQ])
                            k_corr_tot = k_corr_temp.reshape(NT,NP,NQ,Nb-1,Nq*Nq)
                
                    # Resorting
                    nn = np.argsort(k_corr_tot[iT,iP,iQ,ib,:])
                    k_corr_tot[iT,iP,iQ,ib,:] = k_corr_tot[iT,iP,iQ,ib,nn]
                    w_tot[iT,iP,iQ,ib,:] = w_temp[nn]
                
                    # Rebinning
                    w_rebin = np.cumsum(gauss_sample)
                    w_sum[iT,iP,iQ,ib,:] = np.cumsum(w_tot[iT,iP,iQ,ib,:])

                    for iq1 in range(1,Nq):
                        wh, = np.where((w_sum[iT,iP,iQ,ib,:] < w_rebin[iq1])*(w_sum[iT,iP,iQ,ib,:] > w_rebin[iq1-1]))  
						 
		        k_corr_rebin_gcm[iq1,ib,iQ,iP,iT] = (np.nansum(k_corr_tot[iT,iP,iQ,ib,wh]*(w_tot[iT,iP,iQ,ib,wh]))
		                                         /(np.nansum(w_tot[iT,iP,iQ,ib,wh])))
                        k_corr_rebin[iT,iP,iQ,ib,iq1] = k_corr_rebin_gcm[iq1,ib,iQ,iP,iT]    
			
                    
                            
    # Saving data
    
    #np.save("%s/k_corr_RORR_%s_%s_%s.npy"%(directory,typeband,typegauss,typeplanet),k_corr_tot)
    #np.save("%s/weights_RORR_%s_%s_%s.npy"%(directory,typeband,typegauss,typeplanet),w_tot)
    #np.save("%s/k_corr_RORR_%s_%s_%s.npy"%(directory,typeband,typegauss,typeplanet),k_corr_rebin)
    
    np.save("%s/k_corr_RORR_H2O_CO2.npy"%(adress),k_corr_rebin)
    
    #np.save("%s/k_corr_RORR_gcm_%s_%s_%s.npy"%(directory,typeband,typegauss,typeplanet),k_corr_rebin_gcm)
    #np.save("%s/weights_RORR_%s_%s_%s.npy"%(directory,typeband,typegauss,typeplanet),w_rebin)
    
    print "Data saved"

########################################################################################################    
    
#    genpos=np.maximum(gen,np.zeros(gen.size))
#    k_cloud_rmd = np.zeros((P.size,bande_sample.size))
#    rmin=r_cloud[0]
#    rmax=r_cloud[-1]
#    if r_eff>1 :
#        Ncloud=r_eff
#        fac=np.cbrt(3/(4*np.pi*rho_p*Ncloud))
#        r_effs=fac*np.cbrt(genpos[:])
#    else :
#        r_effs=r_eff*np.ones(P.size)        
#    r_effs=np.maximum(r_effs,rmin*np.ones(P.size))
#    r_effs=np.minimum(r_effs,rmax*np.ones(P.size))
#
#    if Script == True :
#        if rank == rank_ref :
#	    print('cloud_scatttering: min max particle radii= %f, %f micron'%(np.amin(r_effs)/1.e-6,np.amax(r_effs)/1.e-6))
#            bar = ProgressBar(bande_sample.size,'Clouds scattering computation progression')
#
#    local_nu=bande_sample
#    numin=bande_cloud[0]
#    numax=bande_cloud[-1]
#    local_nu=np.maximum(local_nu,numin*np.ones(bande_sample.size))
#    local_nu=np.minimum(local_nu,numax*np.ones(bande_sample.size))
#
#    i_r_up = np.searchsorted(r_cloud, r_effs)
#    coeff_r=(r_effs - r_cloud[i_r_up-1])/(r_cloud[i_r_up] - r_cloud[i_r_up-1])
#
#    i_b_up = np.searchsorted(bande_cloud, local_nu)
#    coeff_band=(local_nu - bande_cloud[i_b_up-1])/( bande_cloud[i_b_up] - bande_cloud[i_b_up-1])
#    for i_bande in range(bande_sample.size) :
#        Qint = (1.-coeff_band[i_bande])*Qext[:,i_b_up[i_bande]-1]+ coeff_band[i_bande]* Qext[:,i_b_up[i_bande]]
#        Q_fin = (1.-coeff_r) *Qint[i_r_up-1]+ coeff_r*Qint[i_r_up]
#        k_cloud_rmd[:,i_bande] = 3/4.*(Q_fin*genpos*P*M/(rho_p*r_effs*R_gp*T))
#
#        if Script == True :
#            if rank == rank_ref :
#                bar.animate(i_bande + 1)