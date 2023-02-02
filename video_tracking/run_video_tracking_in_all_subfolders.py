
import sys # this is for system paths
import os
import time

# change -m 2, or -m 1
# 1 is for selecting and 2 for running


def scanfolder(startingDir):
	for path, dirs, files in os.walk(startingDir):
		for f in files:
			if f.endswith('.avi'):
				file_and_path = os.path.join(path,f)
				print("found video file with name: ", file_and_path)
				# for tracking
				os.system("python3 video_tracking_tetrahymena_with_shape.py" + " -m 2 -s 0 -i 0 -v " + '"' + file_and_path + '"')
				# for estimating densities
				# os.system("python3 video_tracking_tetrahymena.py" + " -s 20 -m 2 -i 1 -v " + '"' + file_and_path + '"')

		for d in dirs:
			print("now check dir: ", d)
			scanfolder(d)

def doingNothing():
	return True


# Here I specify the path to cv2.
# On my computer the actual name is cv2.cpython-35m-darwin.so
# and I create a link that I name cv2.so
sys.path.append('/usr/local/lib/python3.7/site-packages/')
sys.path.append('/home/aperna/.local/lib/python3.7/dist-packages/')


# file dialog 
import tkinter as tk
from tkinter import filedialog
	
root = tk.Tk()
root.withdraw()
root.update()

dir_path = filedialog.askdirectory(title='Please select a directory', initialdir=os.getcwd())
root.update()
root.destroy()
print(dir_path)

scanfolder(dir_path)


