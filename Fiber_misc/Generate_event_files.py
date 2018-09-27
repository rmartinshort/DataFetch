#!/usr/bin/env python
#RMS 2018

import pandas as pd
import numpy as np
import obspy as op
import datetime

def main():

	DAS_df = pd.read_csv("Fnames_soilM_interp.dat")
	DAStimes = DAS_df['time'].apply(lambda x: op.UTCDateTime(x))

	##These files have come from Organize_events_for_fiber_compare.py
	LQ_df = pd.read_csv("Local_events_fb_fiber.dat")
	GB_df = pd.read_csv("Global_events_fb_fiber.dat")

	#determine the files that are closest in time to the start and end times of each event
	fnames = DAS_df['timestamp']
	isevent = np.zeros(len(fnames)) 
	eventid = np.empty(len(fnames))
	eventid[:] = np.nan

	#lengths of buffers on either side of the event arrival time, in seconds
	#Need this because the DAS timing is inaccurate
	bufferlen0 = 5*60
	bufferlen1 = 5*60

	events = LQ_df[:100]

	for index, row in events.iterrows():
	 	print("Working on local event %i" %index)
	 	stime = op.UTCDateTime(row['atime'])-bufferlen0
	 	etime = op.UTCDateTime(row['atime'])+bufferlen1

	 	#Get the start and end indices, corresponding to files that are cloest in time to the defined trim points 
	 	ind_start = np.argmin(np.array(DAS_df['time'].apply(lambda x: abs(stime - op.UTCDateTime(x)))))
	 	ind_end = np.argmin(np.array(DAS_df['time'].apply(lambda x: abs(etime - op.UTCDateTime(x)))))

	 	print(ind_start,ind_end)

	 	#Note that event is present
	 	isevent[ind_start:ind_end] = 1.0
	 	#Note that the event has this ID
	 	eventid[ind_start:ind_end] = index
		
        DAS_df['isevent'] = isevent
	DAS_df['event_id'] = eventid
	event_0 = DAS_df[DAS_df['event_id']==0]

	for eventID in events.index:
	 	ev = DAS_df[DAS_df['event_id']==eventID]
	 	ev_2EW = ev[ev['line']=='Line2EW']
	 	ev_CSN = ev[ev['line']=='LineCSN']
	 	name1 = "Local_event_%s_files_2EW.dat" %(eventID)
	 	name2 = "Local_event_%s_files_CSN.dat" %(eventID)
	 	ev_2EW['fname'].to_csv(name1,index=False)
	 	ev_CSN['fname'].to_csv(name2,index=False)

	#----------------------------------------------------------
	#Global events
	#----------------------------------------------------------

	#determine the files that are closest in time to the start and end times of each event
	#fnames = DAS_df['timestamp']
	#isevent = np.zeros(len(fnames))
	#eventid = np.empty(len(fnames))
	#eventid[:] = np.nan

	#for index, row in GB_df.iterrows():
	#	print("Working on global event %i" %index)
	#	stime = op.UTCDateTime(row['atime'])-bufferlen0
	#	etime = op.UTCDateTime(row['atime'])+bufferlen1

	#	print(op.UTCDateTime(row['atime']))

	#	#Get the start and end indices, corresponding to files 
	#	ind_start = np.argmin(np.array(DAS_df['time'].apply(lambda x: abs(stime - op.UTCDateTime(x)))))
	#	ind_end = np.argmin(np.array(DAS_df['time'].apply(lambda x: abs(etime - op.UTCDateTime(x)))))

	#	print(ind_start,ind_end)

	#	isevent[ind_start:ind_end] = 1.0
	#	eventid[ind_start:ind_end] = index

	#DAS_df['isevent'] = isevent
	#DAS_df['event_id'] = eventid

	#for eventID in GB_df.index:
	#	ev = DAS_df[DAS_df['event_id']==eventID]
	#	ev_2EW = ev[ev['line']=='Line2EW']
	#	ev_CSN = ev[ev['line']=='LineCSN']
	#	name1 = "Global_event_%s_files_2EW.dat" %(eventID)
	#	name2 = "Global_event_%s_files_CSN.dat" %(eventID)
	#	ev_2EW['fname'].to_csv(name1,index=False)
	#	ev_CSN['fname'].to_csv(name2,index=False)


if __name__ == '__main__':

	main()
