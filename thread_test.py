#!/usr/bin/env python

#RMS 2018
#Test of processing mseed files in parallel.
#Use the bash call time ./thread_test.py to report how long this takes compared to the serial version 

from threading import Thread 
import obspy as op
import time 
import numpy as np
import glob


def main():

    Nthreads = 4 
    allfiles = glob.glob('waveforms2/*.mseed') #Get a list of all the available files

    #Run threaded processing
    threaded(Nthreads,allfiles)
    #Run the same processing in serial
    serial(allfiles)

def threaded(Nthreads,allfiles):

        chunksize = int(np.floor(len(allfiles)/(Nthreads)))
        indices = []

        #Generate a list that specifies which files each thread should work on 

        j = 0
        for i in range(0,len(allfiles),chunksize):
            if j >= Nthreads:
                break
            else:
                indices.append([i,i+chunksize])
            j += 1 
                       
        indices[-1][-1] = len(allfiles)

        print("Threads to ")
        print(indices)


        #Initiate threads

        threads = []    
        for i in range(Nthreads):
            t = Thread(target=SeismicProcess, args=(indices[i],allfiles,i,20,0.05,2))
            threads.append(t)
            t.start()

def serial(allfiles):

    fmin = 0.05
    fmax = 2.0
    samp = 20

    St = op.Stream()

    for fname in allfiles:

        tr=op.read(fname,format='mseed')
        tr.interpolate(samp)
        St += tr[0]

    St.detrend(type='demean')
    St.detrend(type='linear')
    St.filter(type='bandpass',freqmin=fmin,freqmax=fmax)

    k = 1
    for trace in St:
        fname = '%s_%i.mseed' %(baseout,k)
        trace.write(fname,format='mseed')
        k += 1
    

def SeismicProcess(indices,fileslist,jobid,samp=None,fmin=None,fmax=None):
    
    '''indices is a list of the start and stop index from the fileslist to process. fmin and fmax are the
    bandpass filter limits. This applies standard seismic processing in parallel. Will probably be useful in 
    the fiber project'''
    
    processfiles = fileslist[indices[0]:indices[1]]
    
    St = op.Stream()
    print("I am thread %i with %i files to process" %(jobid,len(processfiles)))
    
    baseout = 'test'
    
    if samp:
        
        for fname in processfiles:
            tr = op.read(fname,format='mseed') 
            tr.interpolate(samp)
            St += tr[0]
    
    else:
        
        for fname in processfiles:
            tr = op.read(fname,format='mseed')
            St += tr[0]
            
    
    St.detrend(type='demean')
    St.detrend(type='linear')
    
    if fmin:
        St.filter(type='bandpass',freqmin=fmin,freqmax=fmax)
        
    
    k = 1
    for trace in St:
        fname = '%s_%i.mseed' %(baseout,k)
        trace.write(fname,format='mseed')
        k += 1
    

if __name__ == '__main__':

    main()
        


