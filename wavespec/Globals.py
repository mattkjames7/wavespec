import os
import platform

ModulePath = os.path.dirname(__file__) + "/__data/"


def _lib_extension_preference():
	systype = platform.system()
	if systype == 'Darwin':
		return ['dylib', 'so', 'dll']
	if systype == 'Windows':
		return ['dll', 'pyd', 'so', 'dylib']
	return ['so', 'dylib', 'dll']


def _resolve_libfile(libdir, basename):
	# Prefer the platform-native extension, but fall back to others.
	for ext in _lib_extension_preference():
		candidate = os.path.join(libdir, f"{basename}.{ext}")
		if os.path.isfile(candidate):
			return candidate
	# Return the primary expected name even when missing, so warnings stay clear.
	return os.path.join(libdir, f"{basename}.{_lib_extension_preference()[0]}")

#this is where we will store the Lomb-Scargle stuff
libLSpath = os.path.dirname(__file__)+"/__data/liblombscargle/"
libLSfile = _resolve_libfile(libLSpath, "liblombscargle")
libLS = None

#this is where we will store the irregular filter stuff
libFpath = os.path.dirname(__file__)+"/__data/libfilter/"
libFfile = _resolve_libfile(libFpath, "libfilter")
libF = None
