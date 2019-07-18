Author: Arshia Ruina
Date: 18 July 2019

These textfiles contain the names of the data files that can be used for running a script for analysis.
Creating these textfiles makes it easy to access all the files by just reading the text files line by line inside the script.
They were made in the following way:

In a python terminal,

	ruina@gridvm10:~/DAMPE/He_flux_analysis/datafiles$ python
	Python 2.7.14 |Anaconda, Inc.| (default, Nov 20 2017, 18:04:19) 
	[GCC 7.2.0] on linux2
	Type "help", "copyright", "credits" or "license" for more information.
	>>> outfile = open('20181019.txt','w')
	>>> from glob import glob
	>>> for f in glob("/beegfs/dampe/prod/FM/FlightData/2A/20181019/*/*.root"):
	...     outfile.write(f)
	... 
	>>> [Ctrl+d to exit]

After this, the filenames inside the text file created need to separated into new lines. To do this in vim,

1. Search for the pattern before which a new line has to be inserted.
For example, in command mode, type Â?/beegfsÂ» which will highlight these patterns.

2. Then type Â«:s//\r&/gÂ

Voila!
