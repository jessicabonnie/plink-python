#!/usr/bin/python2.7

import logging
import sys, os

debug = None
log_file = './log/test_logger.log'

log = logging.getLogger()
ch = logging.StreamHandler()

if os.path.exists(os.path.dirname(log_file)):
	fh = logging.FileHandler(log_file)
else:
	raise "log directory does not exist (" + os.path.dirname(log_file)+")"
	sys.exit(1)

log.addHandler(ch)
log.addHandler(fh)

ch_fmt = logging.Formatter("%(levelname)s\t: %(message)s")
fh_fmt = logging.Formatter("%(asctime)s %(process)d (%(levelname)s)\t: %(message)s")
ch.setFormatter(ch_fmt)
fh.setFormatter(fh_fmt)

if debug:
	log.setLevel(logging.DEBUG)
else:
	log.setLevel(logging.INFO)

log.debug('this is a debug message')
log.info('this is an informational message')
log.warning('this is a warning message')
log.error('this is an error message')
log.critical('this is a critical message')

