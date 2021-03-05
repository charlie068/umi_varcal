try:
	import psutil
except:
	exit()
	
import time
import os
import sys


tmp_name = sys.argv[1]

def RAM(tmp_name):

	processID = os.getpid()

	start = round(psutil.virtual_memory()[3]/1000000000.0, 2)
	last = 0

	counter = 0
	rams = [0]

	done = False
	while not done:

		if os.path.isfile(tmp_name+'/.done'):
			done = True
			os.remove(tmp_name+'/.done')


		time.sleep(1)
		current = round(psutil.virtual_memory()[3]/1000000000.0, 2)
		current -= start

		if current >= last:
			
			try:
				os.remove(tmp_name+'/.'+str(last)+'.ram')
			except:
				pass

			try:
				f = open(tmp_name+'/.'+str(current)+'.ram', 'w')
				f.close()
			except:
				pass
			
			last = current




RAM(tmp_name)
