#!/usr/bin/env python
#RMS 2018

#Example showing how to use General_Fetch to download continuous data

#Download some continuous data and event database associated with that timeframe. Ready for attempt to
#split the dataset into events vs noise

from General_Fetch import Fetch
from obspy import UTCDateTime


def main():

	network = 'TA'
	station = 'TOLK'
	starttime = "2017-08-01"
	endtime = "2017-08-03"
	centercoords = [58, -145] #Center of region for teleseismic request
	minradius = 0
	maxradius = 180

	params = Fetch(network=network,station=station,starttime=UTCDateTime(starttime),endtime=UTCDateTime(endtime),\
		minlatitude=55,maxlatitude=70,minlongitude=-160,maxlongitude=-140)

	#Get all global teleseisms 
	params.fetchEvents(centercoords=centercoords,minradius=minradius,maxradius=maxradius,minmag=5)

	params.writeEvents(centercoords=[-149.57,68.64]) #includes the coordinates of TOLK station 

	#Get all local events
	params.fetchEvents(minmag=0.0)
	params.writeEvents(centercoords=[-149.57,68.64])

	#Download the data 
	params.GetData(datadirpath='Continuous_req_example')

if __name__ == '__main__':

	main()
