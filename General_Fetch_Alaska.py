#!/usr/bin/env python 

from General_Fetch import Fetch
from obspy import UTCDateTime

def main():

	network='TA'
	station=None
	starttime = "2014-01-01"
	endtime = "2014-02-01"
	#centercoords = [64, -149]
	minmag = 3.0
	maxmag = 9.0
	mindepth = 10 #km
	maxdepth = None #km
	minlat=51
	maxlat=71
	minlon=-166
	maxlon=-132
	station_autoselect_flag = True

	#Set up test instance. Ensure that the stations you request (can leave this blank to request all stations)
	#are within the boundary box given

	test = Fetch(network=network,station=station,starttime=UTCDateTime(starttime),endtime=UTCDateTime(endtime),
		minlatitude=minlat,maxlatitude=maxlat,minlongitude=minlon,maxlongitude=maxlon,maxdepth=maxdepth,\
		mindepth=mindepth, vmodel='ak135',station_autoselect=station_autoselect_flag)

	#Fetch the details of the events that are within the request region 
	test.fetchEvents(minmag=minmag, mindepth=mindepth,maxdepth=maxdepth,maxmag=maxmag)

	#Get all the station information
	test.fetchInventory()

	#Write the event details to a file
	test.writeEvents()

	#Write the stations to a file
	test.writeStations()
	
	#Get the data and store in a file called waveforms_dir. You can give it any path
	#Note that the data are downloaded as mseed files, one per component
	print("Getting data")
	test.GetData(req_type='event',datadirpath='alaska_data_'+str(starttime)+'_'+str(endtime)+'_'+\
		str(minlat)+'_'+str(maxlat)+'_'+str(minlon)+'_'+str(maxlon)+'_'+str(minmag)+'_'+str(mindepth)+\
		'_km', vmodel='ak135')

	#Remove instrument response. Default is to displacement
	#test.CorrectResponse()

if __name__ == '__main__': 

	main()
