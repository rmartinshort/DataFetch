#!/usr/bin/env python
#RMS 2018

#Get a list of the events that might have been recorded on the Fairbanks DAS array

from General_Fetch import Fetch
from obspy import UTCDateTime

def main():

    starttime = '2016-07-29'
    endtime = '2016-10-04'

    params = Fetch(starttime=UTCDateTime(starttime),endtime=UTCDateTime(endtime),\
    		minlatitude=55,maxlatitude=70,minlongitude=-160,maxlongitude=-140)

    #Get all local events
    centercoords = [-147.6732,64.8752]
    #params.fetchEvents(minmag=0.0)
    #params.writeEvents(centercoords=centercoords)


    #Teleseismic/large events
    centercoords = [64.8752,-147.6732]
    minradius = 0
    maxradius = 180
    minmag = 6
    params.fetchEvents(centercoords=centercoords,minradius=minradius,maxradius=maxradius,minmag=minmag)
    params.writeEvents(centercoords=centercoords)


if __name__ == "__main__":

    main()
