"""
Set of functions for reco vs true assignation
"""

import numpy as np
import matplotlib.pyplot as plt
import tables as tb
import pandas as pd
import math 

def blobassignation(true_info, reco_info, tag):
    diff_b1_x, diff_b1_y, diff_b1_z = [], [], []
    diff_b2_x, diff_b2_y, diff_b2_z = [], [], []

    diff_b1_x_test, diff_b1_y_test, diff_b1_z_test = [], [], []
    diff_b2_x_test, diff_b2_y_test, diff_b2_z_test = [], [], []

    #[0] e_reco, [1] e_blob1, [2] e_blob2
    reco_e, reco_e_track, reco_e_blob1, reco_e_blob2 = [], [], [], []

    totalevents, events12AB, events12BA, eventsNA = 0, 0, 0, 0

    direction12, direction21 = False, False

    for nevt in range(1,true_info.event_id.nunique()):
    
        if nevt not in reco_info.event.unique():
            continue
        if nevt not in true_info.event_id.unique():
            continue
        if (tag == 'signal'):
            true_x_pA = float(true_info[(true_info.event_id == nevt) &
                                               (true_info.primary == True) &
                                               (true_info.particle_id == 1)].final_x.values)
            true_y_pA = float(true_info[(true_info.event_id == nevt) &
                                               (true_info.primary == True) &
                                               (true_info.particle_id == 1)].final_y.values)
            true_z_pA = float(true_info[(true_info.event_id == nevt) &
                                               (true_info.primary == True) &
                                               (true_info.particle_id == 1)].final_z.values)

            true_x_pB = float(true_info[(true_info.event_id == nevt) &
                                               (true_info.primary == True) &
                                               (true_info.particle_id == 2)].final_x.values)
            true_y_pB = float(true_info[(true_info.event_id == nevt) &
                                               (true_info.primary == True) &
                                               (true_info.particle_id == 2)].final_y.values)
            true_z_pB = float(true_info[(true_info.event_id == nevt) &
                                               (true_info.primary == True) &
                                               (true_info.particle_id == 2)].final_z.values)
        elif (tag == 'bkg'):
            true_x_pA = float(true_info[(true_info.event_id == nevt) &
                                               (true_info.primary == True)].initial_x.values)
            true_y_pA = float(true_info[(true_info.event_id == nevt) &
                                               (true_info.primary == True)].initial_y.values)
            true_z_pA = float(true_info[(true_info.event_id == nevt) &
                                               (true_info.primary == True)].initial_z.values)

            true_x_pB = float(true_info[(true_info.event_id == nevt) &
                                               (true_info.primary == True)].final_x.values)
            true_y_pB = float(true_info[(true_info.event_id == nevt) &
                                               (true_info.primary == True)].final_y.values)
            true_z_pB = float(true_info[(true_info.event_id == nevt) &
                                               (true_info.primary == True)].final_z.values)
    
        reco_x_blobA = float(reco_info[reco_info.event == nevt].blob1_x.values)
        reco_y_blobA = float(reco_info[reco_info.event == nevt].blob1_y.values)
        reco_z_blobA = float(reco_info[reco_info.event == nevt].blob1_z.values)
    
        reco_x_blobB = float(reco_info[reco_info.event == nevt].blob2_x.values)
        reco_y_blobB = float(reco_info[reco_info.event == nevt].blob2_y.values)
        reco_z_blobB = float(reco_info[reco_info.event == nevt].blob2_z.values)
    
        reco_e_blobA = float(reco_info[reco_info.event == nevt].eblob1.values)
        reco_e_blobB = float(reco_info[reco_info.event == nevt].eblob2.values)
    
        reco_e_track.append(float(reco_info[reco_info.event == nevt].energy.values))
    
        R_12 = math.sqrt((true_x_pB-true_x_pA)**2+(true_y_pB-true_y_pA)**2+(true_z_pB-true_z_pA)**2)/2
        R_AB = math.sqrt((reco_x_blobA-reco_x_blobB)**2+(reco_y_blobA-reco_y_blobB)**2+(reco_z_blobA-reco_z_blobB)**2)/2
        d_A1 = math.sqrt((true_x_pA-reco_x_blobA)**2+(true_y_pA-reco_y_blobA)**2+(true_z_pA-reco_z_blobA)**2)
        d_B1 = math.sqrt((true_x_pA-reco_x_blobB)**2+(true_y_pA-reco_y_blobB)**2+(true_z_pA-reco_z_blobB)**2)
        d_A2 = math.sqrt((true_x_pB-reco_x_blobA)**2+(true_y_pB-reco_y_blobA)**2+(true_z_pB-reco_z_blobA)**2)
        d_B2 = math.sqrt((true_x_pB-reco_x_blobB)**2+(true_y_pB-reco_y_blobB)**2+(true_z_pB-reco_z_blobB)**2)
    
        diff_b1_x.append(reco_x_blobA - true_x_pA)
        diff_b1_y.append(reco_y_blobA - true_y_pA)
    
        diff_b2_x.append(reco_x_blobB - true_x_pB)
        diff_b2_y.append(reco_y_blobB - true_y_pB)
    

    
        totalevents = totalevents+1

    
        if d_A1<R_12 or d_B2<R_12 or (d_A1<d_B1 and  d_A1<d_A2 and d_A1<d_B2) or (d_B2<d_B1 and d_B2<d_A2 and d_B2<d_A1):
            events12AB = events12AB+1
            direction12 = True
            #print('direction A1 B2')
        elif d_A2<R_12 or d_B1<R_12 or (d_A2<d_B1 and d_A2<d_A1 and d_A2<d_B2) or (d_B1<d_B2 and d_B1<d_A1 and d_B1<d_A2):
            events12BA = events12BA+1
            #print('direction A2 B1')
            direction21 = True
        else:
            eventsNA = eventsNA+1
        
        if (direction12 == True) & (reco_e_blobA > reco_e_blobB):
            reco_x_blob1 = reco_x_blobA
            reco_y_blob1 = reco_y_blobA
            reco_z_blob1 = reco_z_blobA
        
            reco_x_blob2 = reco_x_blobB
            reco_y_blob2 = reco_y_blobB
            reco_z_blob2 = reco_z_blobB
    
            true_x_p1 = true_x_pA
            true_y_p1 = true_y_pA
            true_z_p1 = true_z_pA
        
            true_x_p2 = true_x_pB
            true_y_p2 = true_y_pB
            true_z_p2 = true_z_pB
        
            reco_e_blob1.append(reco_e_blobA)
            reco_e_blob2.append(reco_e_blobB)

        elif (direction12 == True) & (reco_e_blobA < reco_e_blobB):
            reco_x_blob2 = reco_x_blobA
            reco_y_blob2 = reco_y_blobA
            reco_z_blob2 = reco_z_blobA
        
            reco_x_blob1 = reco_x_blobB
            reco_y_blob1 = reco_y_blobB
            reco_z_blob1 = reco_z_blobB
    
            true_x_p2 = true_x_pA
            true_y_p2 = true_y_pA
            true_z_p2 = true_z_pA
        
            true_x_p1 = true_x_pB
            true_y_p1 = true_y_pB
            true_z_p1 = true_z_pB
        
            reco_e_blob1.append(reco_e_blobB)
            reco_e_blob2.append(reco_e_blobA)
        
        elif (direction21 == True) & (reco_e_blobA > reco_e_blobB):
            reco_x_blob2 = reco_x_blobB
            reco_y_blob2 = reco_y_blobB
            reco_z_blob2 = reco_z_blobB
        
            reco_x_blob1 = reco_x_blobA
            reco_y_blob1 = reco_y_blobA
            reco_z_blob1 = reco_z_blobA
        
            true_x_p2 = true_x_pA
            true_y_p2 = true_y_pA
            true_z_p2 = true_z_pA
        
            true_x_p1 = true_x_pB
            true_y_p1 = true_y_pB
            true_z_p1 = true_z_pB
        
            reco_e_blob1.append(reco_e_blobA)
            reco_e_blob2.append(reco_e_blobB)
        
        elif (direction21 == True) & (reco_e_blobA < reco_e_blobB):
            reco_x_blob1 = reco_x_blobB
            reco_y_blob1 = reco_y_blobB
            reco_z_blob1 = reco_z_blobB
        
            reco_x_blob2 = reco_x_blobA
            reco_y_blob2 = reco_y_blobA
            reco_z_blob2 = reco_z_blobA
        
            true_x_p1 = true_x_pA
            true_y_p1 = true_y_pA
            true_z_p1 = true_z_pA
        
            true_x_p2 = true_x_pB
            true_y_p2 = true_y_pB
            true_z_p2 = true_z_pB
        
            reco_e_blob1.append(reco_e_blobB)
            reco_e_blob2.append(reco_e_blobA)
        
        diff_b1_x.append(reco_x_blob1 - true_x_p1)
        diff_b1_y.append(reco_y_blob1 - true_y_p1)
    
        diff_b2_x.append(reco_x_blob2 - true_x_p2)
        diff_b2_y.append(reco_y_blob2 - true_y_p2)
        direction12, direction21 = False, False
    

        
    reco_e_blob1 = np.array(reco_e_blob1, dtype=np.float)
    reco_e_blob2 = np.array(reco_e_blob2, dtype=np.float)
    reco_e_track = np.array(reco_e_track, dtype=np.float)

    reco_e.append(reco_e_track)
    reco_e.append(reco_e_blob1)
    reco_e.append(reco_e_blob2)

    print(f'Total events = {totalevents}')
    print(f'Events 12 = {events12AB} ({100*events12AB/totalevents}%)')
    print(f'Events 21 = {events12BA} ({100*events12BA/totalevents}%)')
    print(f'Events NA = {eventsNA} ({100*eventsNA/totalevents}%)')
    
    return diff_b1_x, diff_b1_y, diff_b2_x, diff_b2_y, reco_e
