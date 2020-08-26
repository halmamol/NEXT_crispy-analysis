"""
Set of functions for cut efficiency, background acceptand and f.o.m. evaluation
"""

import numpy as np
import matplotlib.pyplot as plt
import tables as tb
import pandas as pd
import math 


def blobthreshold(minE, maxE, stepsE):
    energy=[]
    nsteps = int(maxE/stepsE)
    
    for i in range(0,nsteps):
        energy.append(minE + i*stepsE)

    return energy

def nevents_afterthreshold(thresholdE, reco_events):
    ncount = 0
    nevents_afterthreshold_energy = []

    for cut_nevent in range(0,len(thresholdE)):

        for tot_nevent in range(0,len(reco_events[0])):
            #print(f'Cut with energy {blobcut_energy[cut_nevent]}, element {tot_nevent} from total events {len(reco_selectron_e)}')
            if (reco_events[2][tot_nevent] > thresholdE[cut_nevent]):
                ncount=ncount+1

        nevents_afterthreshold_energy.append(ncount)
        ncount=0
        
    return nevents_afterthreshold_energy

def sqrterror_array(values):
    values_error = []
    
    for i in range(0, len(values)):
        values_error.append(math.sqrt(values[i]))
    return values_error

def ratio_error(f, a, b, a_error, b_error):
    f_error = f*math.sqrt((a_error/a)**2
                                  +(b_error/b)**2)
    return f_error

def efficiencyterms(nevents_afterthreshold_signal, nevents_afterthreshold_bkg, 
                    nevents_afterthreshold_error_signal, nevents_afterthreshold_error_bkg,
                    reco_signal_e, reco_bkg_e):
    
    totalevents_signal = len(reco_signal_e[0])
    totalevents_bkg = len(reco_bkg_e[0])

    totalevents_signal_error=math.sqrt(totalevents_signal)
    totalevents_bkg_error=math.sqrt(totalevents_bkg)

    e = nevents_afterthreshold_signal/totalevents_signal
    b = nevents_afterthreshold_bkg/totalevents_bkg

    fom = []
    for i in range(0,len(e)):
        fom.append(e[i]/math.sqrt(b[i]))
        
    fom_error, e_error, b_error = [], [], []

    for i in range(0,len(e)):
        e_error.append(ratio_error(e[i],nevents_afterthreshold_signal[i], totalevents_signal, 
                                   nevents_afterthreshold_error_signal[i], totalevents_signal_error))
        b_error.append(ratio_error(b[i],nevents_afterthreshold_bkg[i], totalevents_bkg, 
                                   nevents_afterthreshold_error_bkg[i], totalevents_bkg_error))
    
    for j in range(0,len(e)):
        fom_error.append(math.sqrt((e_error[j]/b[j])**2
                                      +(b_error[j]*e[j]/(2*b[j]**(3/2)))**2))
        
    return e, b, fom, e_error, b_error, fom_error


def best_fom(e,b,fom,energy):
    best_fom = {max(fom)}
    best_e = 100*e[fom.index(max(fom))]
    best_b = 100*b[fom.index(max(fom))]
    best_E = energy[fom.index(max(fom))]
    print(f'Best fom is {best_fom}, that corresponds to:')
    print("- signal efficiency of {:.2f}%".format(best_e))
    print("- background acceptance of {:.2f}%".format(best_b))
    print(f'- energy threshold of {best_E} MeV')
