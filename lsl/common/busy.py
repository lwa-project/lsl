# -*- coding: utf-8 -*-

"""
Module to make a blinking ASCII busy indicator.
"""

import sys
import time
import threading

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['BusyIndicator', '__version__', '__revision__', '__all__']


class BusyIndicator(object):
	def __init__(self, message='Busy', interval=0.5, dots=3):
		self.message = message
		self.interval = interval
		self.dots = dots
		
		self.thread = None
		self.alive = threading.Event()
		
	def start(self):
		if self.thread is not None:
			self.stop()
			
		self.thread = threading.Thread(target=self.run, name='indicator')
		self.thread.setDaemon(1)
		self.alive.set()
		self.thread.start()
		
	def stop(self):
		if self.thread is not None:
			sys.stdout.write('%s\r' % (' '*(len(self.message)+self.dots)))
			sys.stdout.write("Done\n")
			
			self.alive.clear()
			self.thread.join()
			self.thread = None
			
	def run(self):
		i = 0
		while self.alive.isSet():
			sys.stdout.write('%s%s%s\r' % (self.message, '.'*i, ' '*(self.dots-i)))
			sys.stdout.flush()
			
			i += 1
			i %= (self.dots+1)
			time.sleep(self.interval)