import os

#this is where we will store the Lomb-Scargle stuff
libLSpath = os.path.dirname(__file__)+"/__data/liblombscargle/"
libLSfile = libLSpath + "liblombscargle.so"
libLS = None
