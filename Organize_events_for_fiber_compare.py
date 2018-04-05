#!/usr/bin/env python
#RMS 2018
#Open an event file made by FetchData and determine the fiber file name timestamps that correspond
#to the event arrivals. This will enable easy extraction of this data from the dataset

import pandas as pd
import obspy as op


def main():

    event_file_name = "Global_Events_2016-07-29T00:00:00.000000Z_2016-10-04T00:00:00.000000Z_55_-160_70_-140_6.dat"

    fb_df = pd.read_csv(event_file_name,
                   sep=' ',names=['lon','lat','dep','mag','otime','atime','aphase','dist'])

    fb_df.dropna(inplace=True)
    fb_df.sort_values(by='mag',ascending=False,inplace=True)
    fb_df.reset_index(drop=True,inplace=True)

    exp_start_time = op.UTCDateTime('20160813185229')

    stimes = []
    etimes = []

    #Determine the timestamps of the start and end files for each event
    for index, row in fb_df.iterrows():

       tsince_start = op.UTCDateTime(row['atime']) - exp_start_time
       (mins,seconds) = divmod(tsince_start,60)
       startftime = exp_start_time + (mins-10)*60
       endftime = startftime + 15*60
       sname = "%i%02d%02d%02d%02d%02d%02d" %(startftime.year,startftime.month,startftime.day,\
       startftime.hour,startftime.minute,startftime.second,startftime.microsecond)
       sname = sname[:-2]
       ename = "%i%02d%02d%02d%02d%02d%02d" %(endftime.year,endftime.month,endftime.day,\
       endftime.hour,endftime.minute,endftime.second,endftime.microsecond)
       ename = ename[:-2]
       stimes.append(sname)
       etimes.append(ename)

    fb_df['start_file'] = stimes
    fb_df['end_file'] = etimes

    fb_df.to_csv("Global_events_fb_fiber.dat",index=False)


if __name__ == "__main__":

    main()
