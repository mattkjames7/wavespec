import os
import numpy as np
from . import Globals

def _CheckFirstImport():
	#check if we need root or not!
	path = os.path.dirname(__file__)
	if '/usr/local/' in path:
		sudo = 'sudo '
	else:
		sudo = ''
	

	#first of all - check if the shared object exists
	if not os.path.isfile(Globals.libLSfile):
		print("liblombscargle.so not found - attempting compilation!")
		CWD = os.getcwd()
		os.chdir(Globals.libLSpath)
		os.system(sudo+'make')
		os.chdir(CWD)
