

# import the necessary packages
import argparse
import sys # this is for system paths
import numpy as np
import os
import pandas as pd
import math # for cos and sin

 
# Here I specify the path to cv2.
# On my computer the actual name is cv2.cpython-35m-darwin.so
# and I create a link that I name cv2.so
sys.path.append('/usr/local/lib/python3.7/site-packages/')
sys.path.append('/home/aperna/.local/lib/python3.7/dist-packages/')

# Import cv2
import cv2

# from matplotlib import cm
# myColourMap = cm.prism(np.arange(256))
# print(myColourMap.shape)
# print(myColourMap)

myColourMap =  255 * np.random.random((100,3))
myColourMap = myColourMap.astype(int)

# OPENCV is BGR!!!!

# If reddish colourmap:
myColourMap[:,0] = myColourMap[:,0]*0.5
myColourMap[:,1] = myColourMap[:,1]*0.5
myColourMap[:,2] = myColourMap[:,2]*0.5 + 127

# If yellowish colourmap:
# myColourMap[:,0] = myColourMap[:,0]*0.5
# myColourMap[:,1] = myColourMap[:,1]*0.5 + 127
# myColourMap[:,2] = myColourMap[:,2]*0.5 + 127

# If blueish colourmap:
# myColourMap[:,0] = myColourMap[:,0]*0.5 + 127
# myColourMap[:,1] = myColourMap[:,1]
# myColourMap[:,2] = myColourMap[:,2]*0.5


print(myColourMap)

# an empty coordinate
refPt = [float('nan'), float('nan')]

memory = 120 # number of frames for which we plot past coordinates
showImage = 1 # whether to show the image while processing
# allData = [("time", "frame", "lastFrameDetected", "lastFrameNotDetected", "entity", "assigned", "zone", "size", "x", "y")]



ap = argparse.ArgumentParser()
ap.add_argument("-v", "--video", type = str, default = '', help="video file name (including path)")
ap.add_argument("-t", "--trajectories", type = str, default = '', help="trajectory file name (including path)")
ap.add_argument("-m", "--memory", default = 120, help="Flag indicating how many frames backwards in time are plotted")
ap.add_argument("-i", "--imview", default = 1, help="Flag indicating whether to show tracked video frames while processing them")
args = vars(ap.parse_args())

if args.get("memory", True):
	memory = int(args["memory"])

if args.get("imview", True):
	showImage = int(args["imview"])
	
# video file
file_path = args["video"]
print("video file is ", file_path)

if os.path.isfile(file_path) == False:

	# file dialog 
	import tkinter as tk
	from tkinter import filedialog
	
	root = tk.Tk()
	root.withdraw()
	root.update()
	
	file_path = filedialog.askopenfilename()
	print(file_path)
	root.update()

drive, path_and_file = os.path.splitdrive(file_path)
pathName, fileName = os.path.split(path_and_file)


# trajectory file
traj_file_path = args["trajectories"]
print("trajectory file is ", traj_file_path)

if os.path.isfile(traj_file_path) == False:

	# file dialog 
	import tkinter as tk
	from tkinter import filedialog
	
	root = tk.Tk()
	root.withdraw()
	root.update()
	
	traj_file_path = filedialog.askopenfilename()
	print(traj_file_path)
	root.update()

traj_drive, traj_path_and_file = os.path.splitdrive(traj_file_path)
traj_pathName, traj_fileName = os.path.split(traj_path_and_file)


# read all the trajectory
data = pd.read_csv(traj_file_path, skipinitialspace = True)#, dtype = {"revenues" : "float64"})
print(data.columns)
# print(data.head())

dataByFrame = data.sort_values(by=['frame', 'particle']) # there should be no need for this!
# print("sorted array")
print(dataByFrame.head(20))

# dataByParticle = data.sort_values(by=['particle', 'frame'])
# print("sorted array")
# print(dataByParticle.head())

# # create an index, with the first and last element of each frame
# frameNumbers = dataByFrame['frame'].values
# print(frameNumbers[0:30])
# newFrames = np.diff(frameNumbers, n=1, axis = 0)
# newFrames = np.insert(newFrames, 0, 1, axis=0) # the first element is "new" and also shifts every other element
# newFrames = np.insert(newFrames, -1, 1, axis=0) # add an element at the end
# print(newFrames[0:30])

# newFrameIndices = np.where(newFrames == 1)
# print(newFrameIndices[0][0:10])
# print(newFrameIndices[0][-1])
# print(dataByFrame.shape)


camera = cv2.VideoCapture(file_path)

# get information about the original video (fps, width, height etc.)
frameWidth = int(camera.get(cv2.CAP_PROP_FRAME_WIDTH))   # float
frameHeight = int(camera.get(cv2.CAP_PROP_FRAME_HEIGHT)) # float
# videoNFrames = int(camera.get(cv2.CAP_PROP_FRAME_COUNT))
fps = int(camera.get(cv2.CAP_PROP_FPS))

startFrame = 0# min(frameNumColumn)
# print(startFrame)
camera.set(cv2.CAP_PROP_POS_FRAMES, startFrame);


outFileName = os.path.splitext(file_path)[0]+'_with_track.mp4'
print(outFileName)

# Define the codec and create VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'avc1') # other options avc1 (for h264), 'mp4v', 'raw '
out = cv2.VideoWriter(outFileName,fourcc, fps, (frameWidth,frameHeight))


# Read the entire video
while(camera.isOpened()):
	# read one frame; ret tells if it was successful
	ret, frame = camera.read()
	if ret == True:

		currFrameNumber = camera.get(cv2.CAP_PROP_POS_FRAMES) # get the current frame
		print("currentFrameNumber: ", currFrameNumber)
		#currTrajSoFar = np.zeros((1,2), dtype=int)
		#for i in range(len(allData)):
			#row = allData[i]
			# print(row[2])
			#if row[2] <= currFrame and np.isfinite(row[-2]) and np.isfinite(row[-1]):
			#	currTrajSoFar = np.append(currTrajSoFar, np.array([[int(row[-2]), int(frameHeight) - int(row[-1])]]), axis=0)
				
				# print(row[-2], row[-1])
				
		# for i in range(len(currTrajSoFar)):
			# print(currTrajSoFar[i])
		# cv2.polylines(frame,[currTrajSoFar[1:]],False,(255,200,0), thickness=4)


		currentFrameIndices = np.where(dataByFrame['frame'].values == currFrameNumber)
		# print("current frame indices")
		# print(currentFrameIndices[0][:])
		#print("here")
		for jj in currentFrameIndices[0][:]:
			# print("current jj:" + str(jj))
			currentParticle = int(dataByFrame.loc[jj, 'particle'])
			# print("processing particle: " + str(currentParticle))
			myInd = currentParticle % 100
			cv2.ellipse(frame, (int(dataByFrame.loc[jj, 'x']), int(dataByFrame.loc[jj, 'y'])),
			(int(dataByFrame.loc[jj, 'bEllipse']), int(dataByFrame.loc[jj, 'aEllipse'])), dataByFrame.loc[jj, 'orientationFromMovement'],
			0, 360, (int(myColourMap[myInd,0]),int(myColourMap[myInd,1]),int(myColourMap[myInd,2])), thickness=5, lineType=8, shift=0)
			#cv2.arrowedLine(frame, (int(dataByFrame.loc[jj, 'x']), int(dataByFrame.loc[jj, 'y'])), (int(dataByFrame.loc[jj, 'x'] + 40*math.sin(dataByFrame.loc[jj, 'orientation'] /180.0*math.pi)), int(dataByFrame.loc[jj, 'y'] - 40*math.cos(dataByFrame.loc[jj, 'orientation']/180.0*math.pi))), (255, 255, 127), 4, 8, 0, 0.3)
			cv2.putText(frame, str(dataByFrame.loc[jj, 'particle']), (int(dataByFrame.loc[jj, 'x']), int(dataByFrame.loc[jj, 'y'])), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0,0,0), 2)
            # cv2.circle(frame, (int(dataByFrame.loc[jj, 'x']), (int(dataByFrame.loc[jj, 'y']), int(5), (0, 255, 0), 2)

			# Here I get the full trajectory of each particle to plot the previous path
			dataCurrentParticle = dataByFrame[dataByFrame['particle'] == currentParticle]
			# print(dataCurrentParticle['frame'])
			#input('Before plot')
			particleTrajectory = dataCurrentParticle[dataCurrentParticle['frame'].isin(np.arange(currFrameNumber - memory, currFrameNumber + 1, 1))]
			#print('Ora')
			particleCoordinates = np.transpose(np.array([particleTrajectory['x'].values, particleTrajectory['y'].values], ndmin=2))
			# print('Dopo')
			# print(particleCoordinates)
			cv2.polylines(frame, [particleCoordinates], False,(int(myColourMap[myInd,0]),int(myColourMap[myInd,1]),int(myColourMap[myInd,2])), thickness=3)
			# input('Waiting for you to press a key')
        
		#cv2.polylines(frame,[staticObjectCoords[1:]],True,(128,255,255), thickness=3)


		# frame = cv2.flip(frame,0)

		# write the flipped frame
		out.write(frame)
		
		if showImage == 1:
			cv2.imshow('frame',frame)
    
		# exit the while loop if we press 'q' when on the video window
		if cv2.waitKey(1) & 0xFF == ord('q'):
			break

# clean everything
camera.release()
out.release()
cv2.destroyAllWindows()



# fileNameOut = ntpath.join(ntpath.dirname(fileName), "compressed_" + ntpath.basename(fileName))
# print(fileNameOut)
# fout = open(fileNameOut,"w")

# # use readline() to read the first line 
# line = f.readline()
# fout.write(line)

# # use the read line to read further.
# # If the file is not empty keep reading one line
# # at a time, till the file is empty
# lineCount = 1
# while line:
	# # in python 3 print is a builtin function, so
	# print(line)
	# # use realine() to read next line
	# line = f.readline()
	# lineElems = line.rsplit("\t")
	# if (len(lineElems) >= 5): # the sixth element is assigned (5 counting from zero)
		# print(lineElems[5])
		# if int(lineElems[5]) == 1:
			# print("ok")
			# fout.write(line)
	# print(len(lineElems))
	# print(lineElems[5])
	# lineCount += 1
	# if lineCount % 1000 == 0:
		# print(lineCount)
# f.close()
# fout.close()
