#!/usr/bin/env python

#RMS 2018
#Easy way to fetch data and metadata for events and/or stations. This was writen for the fiber seismology project
#but should be useful for all applications. Also contains functionality to write information to pandas dataframe
#Uses the obspy mass downloader to obtain event information

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth, locations2degrees
import obspy as op
import datetime
import glob
import sys


class Fetch:
	def __init__(self,network=None,station=None,level='channel',channel='BH*',starttime=None,endtime=None,\
		minlongitude=None,maxlongitude=None,minlatitude=None,maxlatitude=None,clientname="IRIS"):

		self.client = Client(clientname) 
		self.clientname = clientname

		self.network = network
		self.station = station
		self.level = level
		self.channel = channel

		if endtime == 'today':
			endtime = str(datetime.datetime.today())

		self.starttime = UTCDateTime(starttime)
		self.endtime = UTCDateTime(endtime)

		self.minlatitude = minlatitude
		self.minlongitude = minlongitude
		self.maxlatitude = maxlatitude
		self.maxlongitude = maxlongitude

	    #Quake catalog object 
		self.quake_cat = None

	def fetchInventory(self):

		'''Get an obspy inventory containing all the station information'''

		if self.station != 'None':
			self.inventory = self.client.get_stations(network=self.network,station=self.station,level=self.level,\
				channel=self.channel,starttime=self.starttime,endtime=self.endtime,minlongitude=self.minlongitude,\
				minlatitude=self.minlatitude,maxlongitude=self.maxlongitude,maxlatitude=self.maxlatitude)
		else:
			self.inventory = self.client.get_stations(network=self.network,station=None,level=self.level,\
				channel=self.channel,starttime=self.starttime,endtime=self.endtime,minlongitude=self.minlongitude,
				minlatitude=self.minlatitude,maxlongitude=self.maxlongitude,maxlatitude=self.maxlatitude)

	def fetchEvents(self,centercoords=None,minradius=None,maxradius=None,minmag=6,tofile=None,display=True):

		'''Get an obspy quake catalog containing the event information that was requested. If centercoords and min/max radius
		are set, then the program will use those to fetch. If not, it will use the user-supplied box coordinates. User supplied dates and
		times are also used'''

		self.minmag = minmag
		self.minradius = minradius
		self.maxradius = maxradius
		self.centercoords = centercoords 

		if centercoords:

			print "\nGathering earthquakes using center/radius info\n"

			self.quake_cat = self.client.get_events(starttime=self.starttime,endtime=self.endtime,latitude=self.centercoords[0],\
				longitude=self.centercoords[1],minradius=self.minradius,maxradius=self.maxradius,minmagnitude=minmag)

		else:

			print "\nGathering earthquakes within bounding box\n"

			self.quake_cat = self.client.get_events(starttime=self.starttime,endtime=self.endtime,minlatitude=self.minlatitude,\
				maxlatitude=self.maxlatitude,minlongitude=self.minlongitude,maxlongitude=self.maxlongitude,minmagnitude=self.minmag)

		if display == True:

			print "---------------------------------"
			print "Got the following events"
			print "---------------------------------"
			print self.quake_cat

	def writeEvents(self,catalog):

		'''Write event information to file, which can be loaded as a pandas dataframe'''

	def writeStations(self,catalog):

		'''Write station information to file, which can be loaded as a pandas dataframe'''

	def writeRays(self,catalog):

		'''Write station-event information to file, which can be loaded as a pandas dataframe'''

	def GetData(self,stationdirpath='stations',datadirpath='waveforms',req_type='continuous'):

		'''Call obspy mass downloader to get waveform data'''

		self.stationdirpath = stationdirpath
		self.datadirpath = datadirpath

		from obspy.clients.fdsn.mass_downloader import RectangularDomain, CircularDomain,\
		Restrictions, MassDownloader

		if self.quake_cat == None:

			print "Stop: Must call fetchEvents first to get event catalog to download from"

		if req_type == 'continuous':
			domain = RectangularDomain(minlatitude=self.minlatitude,maxlatitude=self.maxlatitude,\
				minlongitude=self.minlongitude,maxlongitude=self.maxlongitude)

			#Download data in daily segements - may want to change

			restrictions = Restrictions(\
                           starttime=starttime,endtime=endtime,\
                           chunklength_in_sec=86400,\
                           channel=self.channel,station=self.station,location="",\
                           reject_channels_with_gaps=False,\
                           minimum_length=0.0,minimum_interstation_distance_in_m=100.0)

			#Call mass downloader to get the waveform information

			mdl = MassDownloader(providers=[self.client])

			mdl.download(domain, restrictions, mseed_storage=datadirpath, stationxml_storage=stationdirpath)

		elif req_type == 'events':

			#Add option for non-continuous download - event/station pairing for example

			domain = RectangularDomain(minlatitude=self.minlatitude,maxlatitude=self.maxlatitude,\
				minlongitude=self.minlongitude,maxlongitude=self.maxlongitude)

			for event in self.quake_cat:

				print "Downloading data for event %s" %event

				#For each event, download the waveforms at all stations requested

				origin_time = event.origins[0].time

				if self.network:

					restrictions = Restrictions(starttime=origin_time - (5 * 60),endtime=origin_time + 3600,\
						reject_channels_with_gaps=False, minimum_length=0.95, minimum_interstation_distance_in_m=10E3,\
						channel=self.channel,location="",network=self.network,station=self.station)

				else:

					restrictions = Restrictions(starttime=origin_time - (5 * 60),endtime=origin_time + 3600,\
						reject_channels_with_gaps=False, minimum_length=0.95, minimum_interstation_distance_in_m=10E3,\
						channel=self.channel)

				mdl = MassDownloader(providers=[self.clientname])

				mdl.download(domain, restrictions, mseed_storage=datadirpath,\
					stationxml_storage=stationdirpath)

	def Set_datapaths(self,waveforms_path="waveforms",station_path="stations"):

		'''Set the directory names where downloaded data can be found'''

		self.stationdirpath = station_path 
		self.datadirpath = waveforms_path


	def CorrectResponse(self,resptype='displacement'):

		'''Correct downloaded data for insrument response'''

		if resptype == 'displacement':
			outtype = "DISP"
		elif resptype == 'velocity':
			outtype = "VEL"
		elif resptype == 'acceleration':
			outtype = "ACC"
		else:
			print "User input correction unit not valid: use displacement, velocity or acceleration"
			sys.exit(1)

		#get the station data
		station_path = '%s/*.xml' %self.stationdirpath
		stations = glob.glob(station_path)

		for station in stations:

			inv = op.read_inventory(station)

			stationname = station.split('/')[1][:-4]

			waveforms_path = '%s/%s*.mseed' %(self.datadirpath,stationname)
			waveforms = glob.glob(waveforms_path)
			
			stream = op.Stream()
			added_waveforms = []

			for waveform in waveforms:

				print "adding waveform %s to stream" %waveform

				try:
					st = op.read(waveform,format='mseed')
					stream += st[0]
					added_waveforms.append(waveform)
				except:
					print "Could not add %s to stream" %waveform
					continue

			print "\nCorrecting responses in stream\n"

			stream.remove_response(inventory=inv,output=outtype)

			#write corrected waveforms to mseed output

			i = 0
			for trace in stream:
				outname = '%s_%s.mseed' %(added_waveforms[i][:-6],outtype)
				trace.write(outname,format='mseed')
				i += 1









if __name__ == '__main__':

	network = "TA"
	station = "H22K"
	starttime = "2016-08-01"
	endtime = "2016-09-01"
	centercoords = [58, -145]
	minradius = 30
	maxradius = 120
	minmag = 6.0 

	test = Fetch(network=network,station=station,starttime=UTCDateTime(starttime),endtime=UTCDateTime(endtime),\
		minlatitude=55,maxlatitude=70,minlongitude=-160,maxlongitude=-140)
	#test.fetchEvents(centercoords=centercoords,minradius=minradius,maxradius=maxradius,minmag=minmag)
	#test.GetData(req_type='events')
	test.Set_datapaths()
	test.CorrectResponse()







