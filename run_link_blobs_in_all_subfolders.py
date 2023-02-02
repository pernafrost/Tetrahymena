
import sys # this is for system paths
import os
import time



def scanfolder(startingDir):
	for path, dirs, files in os.walk(startingDir):
		for f in files:
			if f.endswith('blobs.csv'):
				file_and_path = os.path.join(path,f)
				print("found blob file with name: ", file_and_path)
				# for tracking
				os.system("python3 link_blobs_into_tracks.py" + " -f " + file_and_path)
				# os.system("python link_blobs_into_tracks.py" + " -f " + '"' + file_and_path + '"')

		for d in dirs:
			print("now check dir: ", d)
			scanfolder(d)

def doingNothing():
	return True

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


