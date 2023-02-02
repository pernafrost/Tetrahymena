

# import the necessary packages
import argparse
import sys # this is for system paths
import numpy as np
import os
import trackpy as tp
import pandas as pd

#import matplotlib as mpl
#import matplotlib.pyplot as plt


# file dialog 
import tkinter as tk
from tkinter import filedialog

# Here I specify the path to cv2.
# On my computer the actual name is cv2.cpython-35m-darwin.so
# and I create a link that I name cv2.so
sys.path.append('/usr/local/lib/python3.7/site-packages/')

# # Import cv2
# import cv2

ap = argparse.ArgumentParser()
# ap.add_argument("-m", "--mode", nargs='?',  default=0, type = int, help="running mode: 0 = do everything, 1 = only select ROI, 2 = only execute without selecting ROI")
ap.add_argument("-f", "--file", type = str, default = '', help="blob file name (including path)")
# ap.add_argument("-s", "--sample", default = 0, help="Flag indicating how many frames to sample use negative numbers for all frames")
# ap.add_argument("-i", "--imview", default = 1, help="Flag indicating whether to show tracked video frames during tracking")

args = vars(ap.parse_args())

# if args["mode"] != None:
	# runningMode = args["mode"]
# else:
	# runningMode = 0

# if args.get("sample", True):
	# maxNumFrames = int(args["sample"])

# if args.get("imview", True):
	# showImage = int(args["imview"])

# if runningMode == 1: # It must show the image in order to select the ROI
	# showImage = 1
	

# print("running mode is ", runningMode)

file_path = args["file"]

if os.path.isfile(file_path) == False:
	root = tk.Tk()
	root.withdraw()
	root.update()
	file_path = filedialog.askopenfilename()

pathName, fileName = os.path.split(file_path)

print(file_path)	
data = pd.read_csv(file_path, skipinitialspace = True)#, dtype = {"revenues" : "float64"})
print(data.columns)

print(data.head())

data.rename(columns={'xCentroid':'x', 'yCentroid':'y', 'area':'mass'}, inplace=True)
print(data.columns)
print(data.head())


t = tp.link_df(data, 20, memory=3)
print(t.head())

t1 = tp.filter_stubs(t, 20)

print(t1.columns)

fTrajHeader = t1.columns

# Save results
trajectoriesFileName = os.path.join(pathName, (os.path.splitext(fileName)[0]).replace('blobs', 'trajectories') + '.csv')
print("trajectories file name: ", trajectoriesFileName)

t1.to_csv(trajectoriesFileName)

# plt.figure()
# tp.plot_traj(t1);

# im = tp.imsd(tm, 752.23/1673.43, 30)  # microns per pixel = verify! , frames per second = 30



# fig, ax = plt.subplots()
# ax.plot(im.index, im, 'k-', alpha=0.1)  # black lines, semitransparent
# ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
       # xlabel='lag time $t$')
# ax.set_xscale('log')
# ax.set_yscale('log')




# em = tp.emsd(tm, 752.23/1673.43, 30)


# fig, ax = plt.subplots()
# ax.plot(em.index, em, 'o')
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
       # xlabel='lag time $t$')
# ax.set(ylim=(1e-2, 10));

# plt.figure()
# plt.ylabel(r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]')
# plt.xlabel('lag time $t$');
# tp.utils.fit_powerlaw(em)  # performs linear best fit in log space, plots]
