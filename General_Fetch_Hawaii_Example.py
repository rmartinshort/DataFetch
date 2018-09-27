#!/usr/bin/env python 

from General_Fetch import Fetch
from obspy import UTCDateTime


def main():

	network='HV'
	station=None
	channel='*'
	starttime = "2018-06-04"
	endtime = "2018-06-06"

	minlongitude=-156
	maxlongitude=-154.8
	minlatitude=19
	maxlatitude=20
	minmag = 3.0

	#Set up test instance. Ensure that the stations you request (can leave this blank to request all stations)
	#are within the boundary box given

	test = Fetch(network=network,channel=channel,station=station,starttime=UTCDateTime(starttime),endtime=UTCDateTime(endtime),\
		minlatitude=minlatitude,maxlatitude=maxlatitude,minlongitude=minlongitude,maxlongitude=maxlongitude)

	#Fetch the details of the events that are within the request region 
	test.fetchEvents(minmag=minmag)

	#Write the event details to a file
	test.writeEvents()

	#Get all the station information
	test.fetchInventory()

	#Write the stations to a file
	test.writeStations()

	#Get the data and store in a file called waveforms_dir. You can give it any path
	#Note that the data are downloaded as mseed files, one per component
	print("Getting data")
	test.GetData(req_type='event',datadirpath='hawaii_test')

	#Remove instrument response. Default is to displacement
	#test.CorrectResponse()

if __name__ == '__main__': 

	main()