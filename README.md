# DataFetch

#### A tool to fetch continuous or event-based waveforms and their associated metadata. This should be generic enough to adapt to any seismology project. Makes used of obspy's mass downloader.

This can be used to make event-based or continuous requests for any combination of station, network and events. Optionally it will also download and remove instrument response. See the example scripts for more details of its use

It can also be used to download event-based data only from stations within the Shear Wave Window (SWW), which is necessary for local shear wave splitting studies. See the Alaska example for a use case her. 

