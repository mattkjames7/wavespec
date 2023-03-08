import os

ModulePath = os.path.dirname(__file__) + "/__data/"

#this is where we will store the Lomb-Scargle stuff
libLSpath = os.path.dirname(__file__)+"/__data/liblombscargle/"
libLSfile = libLSpath + "liblombscargle.so"
libLS = None

#this is where we will store the irregular filter stuff
libFpath = os.path.dirname(__file__)+"/__data/libfilter/"
libFfile = libFpath + "libfilter.so"
libF = None
