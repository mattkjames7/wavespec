import os
from . import Globals


def _CheckFirstImport():
	"""Do not build on import. Only check presence of bundled shared libs.

	This function used to attempt to compile bundled C libraries on first
	import. Building is now performed during packaging via `setup.py`.
	If the shared libraries are missing at runtime, print a clear warning
	so users know the package was not installed/built correctly.
	"""

	missing = []
	if not os.path.isfile(Globals.libLSfile):
		missing.append(Globals.libLSfile)
	if not os.path.isfile(Globals.libFfile):
		missing.append(Globals.libFfile)

	if missing:
		print('WARNING: the following bundled shared libraries are missing:')
		for m in missing:
			print('  -', m)
		print('These libraries are built during packaging. Reinstall the package or run `python setup.py build` in the source tree to produce them.')
