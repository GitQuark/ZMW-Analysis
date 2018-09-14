#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 21:31:22 2017

@author: Bobby
"""

import numpy as np
import pandas as pd
from peakdetect import peakdet

def findpeaks(data, thresh):
    baseline = np.median(data[data<np.median(data)])
#    print(baseline)
    noise = np.std(data[data<np.median(data)])
    limit = baseline+thresh*noise
    firing = np.where(data>limit)[0]
#    firing = np.where(data>1.4)[0]
    #### locate the points where the current crosses the threshold ####
    
    startandend = np.diff(firing)
    startpoints = np.insert(startandend, 0, 2)
    endpoints = np.insert(startandend, -1, 2)
    startpoints = np.where(startpoints>1)[0]
    endpoints = np.where(endpoints>1)[0]
    startpoints = firing[startpoints]
    endpoints = firing[endpoints]
    
    #### Eliminate events that start before file or end after file ####
    
    if startpoints[0] == 0:
        startpoints = np.delete(startpoints,0)
        endpoints = np.delete(endpoints,0)
    if endpoints [-1] == len(data):
        startpoints = np.delete(startpoints,-1)
        endpoints = np.delete(endpoints,-1)
    
    #### Track points back up to baseline to find true start and end ####
    
    numberofevents=len(startpoints)
    highthresh = baseline + 0*noise
#    highthresh = .8
    
    for j in range(numberofevents):
        stp = startpoints[j] #mark initial guess for starting point
        while data[stp] > highthresh and stp > 0:
            stp = stp-1 # track back until we return to baseline
        startpoints[j] = stp # mark true startpoint
    
        ep = endpoints[j] #repeat process for end point
        if ep == len(data) -1:  # sure that the current returns to baseline
            endpoints[j] = -1              # before file ends. If not, mark points for
            startpoints[j] = -1              # deletion and break from loop
            ep = 0
            break
        while data[ep] > highthresh:
            ep = ep+1
            if ep == len(data) -1:  # sure that the current returns to baseline
                endpoints[j] = -1              # before file ends. If not, mark points for
                startpoints[j] = -1              # deletion and break from loop
#                ep = 0
                break
#            else:
#                try:
#                    if ep < startpoints[j+1]: # if we hit the next startpoint before we
#                        startpoints[j+1] = 0    # return to baseline, mark for deletion
#                        endpoints[j] = 0                  # and break out of loop
#                        ep = 0
#                        break
#                except:
#                    IndexError
            endpoints[j] = ep
    
    startpoints = startpoints[startpoints!=-1] # delete those events marked for
    endpoints = endpoints[endpoints!=-1]       # deletion earlier
    numberofevents = len(startpoints)
    
    startpoints = np.unique(startpoints)
    endpoints = np.unique(endpoints)
        
    return startpoints, endpoints, limit

def cleanpeaks(red,blue,redstart,redend,bluestart,blueend,redlimit,bluelimit):
    i = -1
    j = 0
    while i < len(redstart)-1:
#        while i < 41:
        i += 1
        j = 0
        re = redend[i]
        rs = redstart[i]
        redset = list(range(rs,re)) 
        while j < len(bluestart):
            be = blueend[j]
            bs = bluestart[j]
            blueset = list(range(bs,be))
 ##############  Start by comparing overlaps  ################# 
############## If peaks start at same time, move one ##############              
            if rs == bs:
                if red[rs] == blue[bs]:
                    print('youre fucked')
                elif  red[rs] > blue[bs]: 
                    bluestart[j] += 1
                    bs +=1                                
                elif  red[rs] < blue[bs]: 
                    redstart[i] += 1
                    rs +=1   
                    redset = list(range(rs,re)) 
############## If peaks end at same time, move one ##############
            if re == be:
                if red[re] == blue[be]:
                    print('youre fucked')
                elif  red[re] > blue[be]: 
                    blueend[j] -= 1
                    be -= 1
                elif  red[re] < blue[be]: 
                    redend[i] -= 1 
                    re -= 1
                    redset = list(range(rs,re)) 
############## Now check if there is any overlap ############## 
            intersect = set(blueset).intersection(redset)
############## if no overlap is found move to next event ##############
            if intersect == set():
                j+=1
                continue
            elif len(intersect) == 1:
                j+=1
                continue
############## If an overlap is found check for crossing points ##############
            else:
                intersect = list(intersect)
                intersect.sort()
                diff = red[intersect] - blue[intersect]
                diff[diff > 0] = 1
                diff[diff < 0] = -1
                diff = np.diff(diff)
                crossing = [intersect[i] for i in np.where(diff != 0.0)[0]]
                diff = diff[diff != 0] 
#                print(i,rs,re,bs,be,crossing)                      
############## if no crossing points exist XXXXXXXXXXXXX  ##############            
                if crossing == []:
                    if rs > bs:
                        if re < be:
                            if (red[intersect] < blue[intersect]).all():
                                redstart = np.delete(redstart,i)
                                redend = np.delete(redend,i)
#                                print('red deleted1')
                            else:
                                bluestart[j] = redend[i]+1
                                j -= 1
                                continue
                        else:
#                                redend[i] = blueend[j]
                            if np.max(blue[rs:be]) < np.max(red[rs:be]):
                                bluestart = np.delete(bluestart,j)
                                blueend = np.delete(blueend,j)
#                                print('blue deleted1')
                            else:
                                redstart[i] = blueend[j]+1
                                rs = redstart[i]
                                redset = list(range(rs,re))
#                                print('red deleted2')
                    elif rs < bs: 
                        if re > be:
                            if (red[intersect] > blue[intersect]).all():
                                    bluestart = np.delete(bluestart,j)
                                    blueend = np.delete(blueend,j) 
#                                    print('blue deleted2')
                            else:
                                pass
#                                print('something else2')
                        else:
#                                blueend[j] = redend[i]
                            if np.max(red[bs:re]) < np.max(blue[bs:re]):
                                bluestart[j] = redend[i]+1
#                                print('blue deleted3')
                            else:
#                                    redend[i] = bluestart[j]-1
#                                    re = redend[i]
#                                    redset = list(range(rs,re))
                                bluestart[j] = redend[i]
                                i-=1
#                                print('red deleted3')
                                continue
                    else:
                        pass
#                        print('something else')
                    i -= 1
                    j = len(bluestart)
######## if crossing points are found, check which color comes first #########            
################################ Red Starts #########################  
                else:    
                    if rs < bs:
                        if red[bs] < blue[bs]:
                            bluestart[j] = redstart[i]-1
                            j -= 1                                                                                                  
                            continue
############# Red Starts and ends, we move blue start and end #################                            
                        if re > be:                                    
###### In some cases, red ends but does not cross, so we move blue end ##########
                            if len(diff) % 2 != 0:
                                if red[be] > blue[be]:
                                    redend = np.sort(np.insert(redend,i,blueend[j]-1))
                                    redstart = np.sort(np.insert(redstart,i,redend[i]))
                                    rs = redstart[i]
                                    re =  redend[i]
                                    redset = list(range(rs,re))
                                    j -= 1
#                                    print('check me1')
                                    continue
                                else:
                                    redend[i] = blueend[j]-1
                                    re =  redend[i]
                                    redset = list(range(rs,re))
                                    bluestart[j] = crossing[0]
                                    j -= 1
#                                    print('check me2')
                                    continue
#                            print('rr')
                            bluestart[j] = crossing[0]
                            blueend[j] = crossing[-1]-1
                            for n,c in enumerate(crossing):
                                if n % 2 == 0:
                                    redend = np.sort(np.insert(redend,i,c-1))
                                    if n == 0:
                                        pass
                                    else:
                                        bluestart = np.sort(np.insert(bluestart,j,c))
                                else:
                                    redstart = np.sort(np.insert(redstart,i,c))
                                    if n != len(crossing)-1:
                                        blueend = np.sort(np.insert(blueend,j,c-1))
        
############# Red Starts and blue ends, move red end, blue start ##############                                                       
                        else: 
                            bluestart[j] = crossing[0]
                            redend[i] = crossing[-1]-1
#                            print('rb')
                            for n,c in enumerate(crossing):
                                if n % 2 == 0:
                                    if n == 0:
                                        pass
                                    else:
                                        redend = np.sort(np.insert(redend,i,c-1))
                                        bluestart = np.sort(np.insert(bluestart,j,c))
                                else:
                                    if n != len(crossing)-1:
                                        blueend = np.sort(np.insert(blueend,j,c-1))
                                        redstart = np.sort(np.insert(redstart,i,c))                                
################################ Blue Starts #########################                                   
                    elif rs > bs:  
                        if red[rs] > blue[rs]:
                            redstart[i] = bluestart[j]-1
                            rs = redstart[i]
                            redset = list(range(rs,re))
                            j -= 1
                            continue
############# Blue Starts and ends, we move red start and end #################
                        if be > re:
###### In some cases, blue ends but does not cross, so we move red end ##########
                            if len(diff) % 2 != 0:
                                blueend = np.sort(np.insert(blueend,j,redend[i]))
                                bluestart = np.sort(np.insert(bluestart,j,blueend[j]))
                                j -= 1
#                                print('check me3')
                                continue
#                            print('bb')
                            redstart[i] = crossing[0]
                            redend[i] = crossing[-1]-1
                            for n,c in enumerate(crossing):
                                if n % 2 == 0:
                                    blueend = np.sort(np.insert(blueend,j,c-1))
                                    if n == 0:
                                        pass
                                    else:
                                        redstart = np.sort(np.insert(redstart,i,c))
                                else:
                                    bluestart = np.sort(np.insert(bluestart,j,c))
                                    if n != len(crossing)-1:
                                        redend = np.sort(np.insert(redend,i,c-1))
############# Blue Starts and red ends, move blue end, red start ##############                                                       
                        else: 
###### In some cases, red ends but does not cross, so we move blue end ##########
                            if len(diff) % 2 == 0:
                                blueend[j] = redend[i]+1
#                                print('check me4')
                                continue                                
#                            print('br')
                            redstart[i] = crossing[0]
                            blueend[j] = crossing[-1]-1 
                            for n,c in enumerate(crossing):
                                if n % 2 == 0:
                                    if n == 0:
                                        pass
                                    else:
                                        blueend = np.sort(np.insert(blueend,j,c-1))
                                        redstart = np.sort(np.insert(redstart,i,c))
                                else:
                                    if n != len(crossing)-1:
                                        redend = np.sort(np.insert(redend,i,c-1))
                                        bluestart = np.sort(np.insert(bluestart,j,c))  
                    i -= 1
                    j = len(bluestart)  

    keep = np.where(redstart < redend)
    redstart = redstart[keep]
    redend = redend[keep]

    keep = np.where(bluestart < blueend)
    bluestart = bluestart[keep]
    blueend = blueend[keep]
    
    for i,x in enumerate(redstart):
        if np.max(red[x:redend[i]]) < redlimit:
            redstart[i] = -1
            redend[i] = -1
                  
    redstart = redstart[redstart != -1]
    redend = redend[redend != -1]


    for i,x in enumerate(bluestart):
        if np.max(blue[x:blueend[i]]) < bluelimit:
            bluestart[i] = -1
            blueend[i] = -1
    bluestart = bluestart[bluestart != -1]
    blueend = blueend[blueend != -1]
    
    return redstart,redend,bluestart,blueend
    

def findSubPeaks(red, blue, redstart,redend,bluestart,blueend, redsubthresh, bluesubthresh):
    pulsecount = 0
    df = pd.DataFrame({"ident":[],'stimes':[],'etimes':[],
                       'peaks':[],'mins':[]})
    
    df['peaks'] = df['peaks'].astype(object)
    df['mins'] = df['mins'].astype(object)    
    
    for i,x in enumerate([[redstart,redend],[bluestart,blueend]]): 
        df = df.append(pd.DataFrame({"ident":[i]*len(x[0]),
        'stimes':x[0],'etimes':x[1]}),ignore_index=True)
    
    df = df.sort_values(by = 'stimes').reset_index()
    
    for i,x in enumerate(df.ident):
        sp = int(df.stimes[i])
        ep = int(df.etimes[i])
        if x == 0:
            intensity = red
            noise = redsubthresh
        else:
            intensity = blue
            noise = bluesubthresh
        peaks, mins = peakdet(v = intensity[sp:ep],delta = noise, x = np.arange(sp,ep))
        if len(peaks) == 0 or len(mins) == 0:
            pulsecount += 1
        else: 
            df.set_value(i,'peaks',peaks.astype(object))
            df.set_value(i,'mins',mins.astype(object))
            pulsecount += 1 + len(peaks)

    print(pulsecount)
        
    return df

#### Note: add deletion of low amplitude subpeaks #########
