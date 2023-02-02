
# I run the script with >> python3 test_read_video.py


# import the necessary packages
import argparse
import sys # this is for system paths
import numpy as np
import os
 
# Here I specify the path to cv2.
# On my computer the actual name is cv2.cpython-35m-darwin.so
# and I create a link that I name cv2.so
sys.path.append('/usr/local/lib/python3.7/site-packages/')

# Import cv2
import cv2



# initialize the list of reference points and boolean indicating
# whether cropping is being performed or not
# allData = [(float('nan') for i in range (10))]
allData = [("time", "frame", "lastFrameDetected", "lastFrameNotDetected", "entity", "assigned", "zone", "size", "x", "y")]
frameCount = 0
entity = 0
currTraj = np.zeros((1,2), dtype=int)
startFrame = 1
refPt = [float('nan'), float('nan')]
yPos = 0 # I define the position of the smaller image 
xPos = 0
winWidth = 800
winHeight = 600
toggleSmallerWindow = False # when selecting a ROI use a smaller window



def add_coordinates(event, x, y, flags, param):
	# grab references to the global variables
	global refPt, currTraj

			
		
	# if the left mouse button was clicked, record the starting
	# (x, y) coordinates and indicate that cropping is being
	# performed
	# if event == cv2.EVENT_LBUTTONDOWN:
	#refPt = [(x, y)]

	# check to see if the left mouse button was released
	if event == cv2.EVENT_LBUTTONUP:
		# record the ending (x, y) coordinates and indicate that
		# the cropping operation is finished
		# refPt.append((frameCount, x, y))
		refPt = [x, y]
		currTraj = np.append(currTraj, np.array([[x + xPos, y + yPos]]), axis=0)
		# tempTraj = currTraj[1:]
		print(currTraj)
		# draw a line around the region of interest
		#cv2.line(frame,refPt[0],refPt[1],(0,255,255),2,8,0)
		cv2.polylines(frame,[currTraj],False,(0,255,255))
		cv2.imshow("frame", smallerRegion)



# file dialog 
import tkinter as tk
from tkinter import filedialog
from tkinter.simpledialog import askinteger

root = tk.Tk()
root.withdraw()
root.update()

file_path = filedialog.askopenfilename()
print(file_path)

root.update()

drive, path_and_file = os.path.splitdrive(file_path)
pathName, fileName = os.path.split(path_and_file)

# Open a videocatpure object
# This can be used to read a video from a camera or from a file
camera = cv2.VideoCapture(path_and_file)
#camera = cv2.VideoCapture('Manip_250_18-02-2016B.mp4')

fps = camera.get(cv2.CAP_PROP_FPS)
print("frames per second: ", fps)

cv2.namedWindow("frame", cv2.WINDOW_NORMAL)
cv2.setMouseCallback("frame", add_coordinates)

frameHeight = camera.get(cv2.CAP_PROP_FRAME_HEIGHT)
frameWidth = camera.get(cv2.CAP_PROP_FRAME_WIDTH)
videoNFrames = int(camera.get(cv2.CAP_PROP_FRAME_COUNT))

camera.set(cv2.CAP_PROP_POS_FRAMES, startFrame)

stopVideo = 0
# Read frame by frame
while (camera.isOpened()):
	print("camera is open")
	# Get an image
	ret, frame = camera.read()
	print(ret)
	if ret == False:
		break
	
	# I skip frames that were not read properly!
	
	# fgmask = cv2.morphologyEx(fgmask, cv2.MORPH_OPEN, kernel)

	
	# keep looping until the 'q' key is pressed
	while True:
		# isolate smaller region so that large videos can play on small monitors
		# print("yPos: " + str(yPos) + "; yPos + winHeight: " + str(yPos + winHeight) + "; xPos: " + str(xPos) +  "; xPos + winWidth: " + str(xPos + winWidth))
		if toggleSmallerWindow == True:
			smallerRegion = frame[yPos:yPos+winHeight, xPos:xPos+winWidth]
			cv2.imshow("frame", smallerRegion)
		else:
			# display the image and wait for a keypress
			cv2.imshow('frame',frame)		

		# display the image and wait for a keypress
		key = cv2.waitKey(1) & 0xFF
		

		# if the 'c' key is pressed, move to next frame
		if key == ord("c"):
			currFrame = camera.get(cv2.CAP_PROP_POS_FRAMES)
			currTime = camera.get(cv2.CAP_PROP_POS_MSEC)
			allData.append([round(currTime)/1000, frameCount, int(currFrame), float('nan'), entity, 1, 1, float('nan'), refPt[0]+xPos, int(frameHeight) - refPt[1] - yPos])
			refPt = [float('nan'), float('nan')]
			break
			
		# if the 'b' key is pressed, move backward
		elif key == ord("b"):
			currFrame = camera.get(cv2.CAP_PROP_POS_FRAMES)
			prevFrame = max(currFrame - 2, 1)
			camera.set(cv2.CAP_PROP_POS_FRAMES, prevFrame)
			break
			
		# if the 'e' key is pressed, move backward and erase last coordinates
		elif key == ord("e"):
			currFrame = camera.get(cv2.CAP_PROP_POS_FRAMES)
			prevFrame = max(currFrame - 2, 1)
			camera.set(cv2.CAP_PROP_POS_FRAMES, prevFrame)
			allData = allData[:-1]
			currTraj = currTraj[:-1]
			refPt = [float('nan'), float('nan')]
			break
			
		# if the 'w' key is pressed, move fast backward
		elif key == ord("w"):
			currFrame = camera.get(cv2.CAP_PROP_POS_FRAMES)
			prevFrame = max(currFrame - 21, 1)
			camera.set(cv2.CAP_PROP_POS_FRAMES, prevFrame)
			break
			
		# if the 'f' key is pressed, move fast forward
		elif key == ord("f"):
			currFrame = camera.get(cv2.CAP_PROP_POS_FRAMES)
			newFrame = currFrame + 19
			camera.set(cv2.CAP_PROP_POS_FRAMES, newFrame)
			# print("currFrame: " + str(currFrame) + "newFrame: " + str(newFrame)) 
			break
		
		elif key == ord("g"):
			frameNum = askinteger('jump to frame',  'go to frame number:', parent = root)
			frameNum = max(1, min(frameNum, videoNFrames - 2))
			camera.set(cv2.CAP_PROP_POS_FRAMES, frameNum)
			break
			
		# add new entity	
		elif key == ord("n"):
			entity += 1
			currTraj = np.zeros((1,2), dtype=int)
			#camera.set(cv2.CAP_PROP_POS_FRAMES, startFrame)
			frameCount = 0
			

			outFileName = os.path.join(pathName, os.path.splitext(fileName)[0] + '_coordinates_entity' + str(entity) + '.txt')
			print("output file name: ", outFileName)

			f = open(outFileName, 'a')

			for t in allData:
				line = ', '.join(str(x) for x in t)
				f.write(line + '\n')
			f.close()			
			allData = [("time", "frame", "lastFrameDetected", "lastFrameNotDetected", "entity", "assigned", "zone", "size", "x", "y")]

			break

		# if the 'q' key is pressed, exit completely from video
		elif key == ord("q"):
			stopVideo = 1
			break;
		
		elif key == ord("t"): # toggle small large window
			if toggleSmallerWindow == False:
				toggleSmallerWindow = True
			else:
				# reset yPos and xPos
				newYPos = 0
				yStep = newYPos - yPos
				refPt[1] = refPt[1] - yStep
				yPos = newYPos
				
				newXPos = 0
				xStep = newXPos - xPos
				refPt[0] = refPt[0] - xStep
				xPos = newXPos
				
				# toggle
				toggleSmallerWindow = False
	
		# if the 'd' key is pressed, move frame down
		elif key == ord("d") and toggleSmallerWindow == True:
			newYPos = min(yPos + 100, int(frameHeight) - winHeight)
			yStep = newYPos - yPos
			refPt[1] = refPt[1] - yStep
			yPos = newYPos
			# break
			
		# if the 'u' key is pressed, move frame up
		elif key == ord("u") and toggleSmallerWindow == True:
			newYPos = max(0, yPos - 100)
			yStep = newYPos - yPos
			refPt[1] = refPt[1] - yStep
			yPos = newYPos
			# break
			
		# if the 'l' key is pressed, move frame left
		elif key == ord("l") and toggleSmallerWindow == True:
			newXPos = max(0, xPos - 100)
			xStep = newXPos - xPos
			refPt[0] = refPt[0] - xStep
			xPos = newXPos
			# break
			
		# if the 'r' key is pressed, move frame right
		elif key == ord("r") and toggleSmallerWindow == True:
			newXPos = min(xPos + 100, int(frameWidth) - winWidth)
			xStep = newXPos - xPos
			refPt[0] = refPt[0] - xStep
			xPos = newXPos
			# break
		
		# if the 's' key is pressed, save one snapshot
		elif key == ord("s"):
			currFrame = camera.get(cv2.CAP_PROP_POS_FRAMES)
			outFileName = os.path.join(pathName, os.path.splitext(fileName)[0] + '_frame_' + str(currFrame) + '.png')
			print("output file name: ", outFileName)
			cv2.imwrite(outFileName, frame)     # save frame as JPEG file	
			


	frameCount += 1

	
	if stopVideo == 1:
		break



# print coordinates
for i in range(len(allData)):
	print(allData[i])


outFileName = os.path.join(pathName, os.path.splitext(fileName)[0] + '_coordinates.txt')
print("output file name: ", outFileName)

f = open(outFileName, 'a')

for t in allData:
    line = ', '.join(str(x) for x in t)
    f.write(line + '\n')
f.close()

# clean everything
camera.release()
cv2.destroyAllWindows()
