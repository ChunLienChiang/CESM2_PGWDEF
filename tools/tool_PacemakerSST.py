"""
Create_PacemakerSST.py
==========================
Create the SST forcing (CTL clim. + 60-year warming footprint) files.
"""

import json

Config = json.load(open('../config.json'))


if (__name__ == '__main__'):

	"""
	When this file is called, create the SST forcing files.
	"""