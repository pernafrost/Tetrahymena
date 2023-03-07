
# select the roi with the mouse, then press s, then press q
# tracking starts
# at the resolution of 40x, 779441.5 pixels are 4*10^-4 ml

# select the ROI, then press s and q


# import the necessary packages
import argparse
import sys # this is for system paths
import numpy as np
import os
import argparse
# from pathlib import Path
 
# Here I specify the path to cv2.
# On my computer the actual name is cv2.cpython-35m-darwin.so
# and I create a link that I name cv2.so
#sys.path.append('/usr/local/lib/python3.5/site-packages/')
# sys.path.append('/home/aperna/.local/lib/python3.7/dist-packages/')
# testsys.path.append('/Users/perna/Library/Python/3.7/lib/python/site-packages/')
sys.path.append('/usr/local/lib/python3.5/site-packages/')


# Import cv2
import cv2

from skimage.feature import peak_local_max
from skimage.morphology import watershed
from skimage.morphology import medial_axis
from skimage.morphology import skeletonize
from scipy import ndimage

# initialize the list of reference points and boolean indicating
# whether cropping is being performed or not
# allData = [(float('nan') for i in range (10))]
fBlobsHeader = "time, frame, xCircle, yCircle, xCentroid, yCentroid, orientation, area, xEllipse, yEllipse, aEllipse, bEllipse\n"
allDensities = []
refPt = []
frameCount = 0
entity = 0
currROI = np.zeros((1,2), dtype=int)
startFrame = 0
maxNumFrames = -1 # max number of frames after which it stops
trackBlobs = True
frameJump = 1 # how many frames it skips every time
refPt = [float('nan'), float('nan')]
minimumArea = 15 # minimum area of components in pixels used 15 for all those adapted at 15 degrees
bgThreshold = 4 # used 4 for all those adapted at 15 degrees
currentFrame = 0
showImage = 1 # I use 0 for false and 1 for true so that it is easier to parse the input argument

yPos = 0 # I define the position of the smaller image 
xPos = 0
winWidth = 800
winHeight = 600
toggleSmallerWindow = False # when selecting a ROI use a smaller window


# def moment_elongation(M):
    # x = M['m20'] + M['m02']
    # y = 4 * M['m11']**2 + (M['m20'] - M['m02'])**2
    # return (x + y**0.5) / (x - y**0.5)



def click_ROI(event, x, y, flags, param):
	# grab references to the global variables
	global refPt, currROI

			
		
	# if the left mouse button was clicked, record the starting
	# (x, y) coordinates and indicate that cropping is being
	# performed
	# if event == cv2.EVENT_LBUTTONDOWN:
	#refPt = [(x, y)]

	# check to see if the left mouse button was released
	if event == cv2.EVENT_LBUTTONDOWN:
		# record the ending (x, y) coordinates and indicate that
		# the cropping operation is finished
		# refPt.append((frameCount, x, y))
		refPt = [x, y]
		currROI = np.append(currROI, np.array([[x + xPos, y + yPos]]), axis=0)
		# print(currROI)
		# draw a line around the region of interest
		#cv2.line(frame,refPt[0],refPt[1],(0,255,255),2,8,0)
		cv2.polylines(frame,[currROI],False,(0,255,255))
		#cv2.imshow("frame", frame)


ap = argparse.ArgumentParser()
ap.add_argument("-m", "--mode", nargs='?',  default=1, type = int, help="running mode: 0 = do everything, 1 = only select ROI, 2 = only execute without selecting ROI")
ap.add_argument("-v", "--video", type = str, default = '', help="video file name (including path)")
ap.add_argument("-s", "--sample", default = 0, help="Flag indicating how many frames to sample use negative numbers for all frames")
ap.add_argument("-i", "--imview", default = 1, help="Flag indicating whether to show tracked video frames during tracking")

args = vars(ap.parse_args())

if args["mode"] != None:
	runningMode = args["mode"]
else:
	runningMode = 0

if args.get("sample", True):
	maxNumFrames = int(args["sample"])

if args.get("imview", True):
	showImage = int(args["imview"])

if runningMode == 1: # It must show the image in order to select the ROI
	showImage = 1
	

print("running mode is ", runningMode)

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

drive, path_and_file = os.path.splitdrive(file_path)
pathName, fileName = os.path.split(path_and_file)


roiFileName = os.path.join(drive, pathName, os.path.splitext(fileName)[0] + '_ROI.txt')
print("roi file name: ", roiFileName)

# Open a videocatpure object
# This can be used to read a video from a camera or from a file
camera = cv2.VideoCapture(file_path)
#camera = cv2.VideoCapture('Manip_250_18-02-2016B.mp4')


cv2.namedWindow("frame")
cv2.setMouseCallback("frame", click_ROI)

camera.set(cv2.CAP_PROP_POS_FRAMES, startFrame);

frameWidth = camera.get(cv2.CAP_PROP_FRAME_WIDTH)   # float
frameHeight = camera.get(cv2.CAP_PROP_FRAME_HEIGHT) # float
videoNFrames = int(camera.get(cv2.CAP_PROP_FRAME_COUNT))

if maxNumFrames > 0:
	frameJump = int(videoNFrames / maxNumFrames)

#currROI = np.append(currROI, np.array([[0, 0]]), axis=0)
#currROI = np.append(currROI, np.array([[0, frameWidth]]), axis=0)
#currROI = np.append(currROI, np.array([[frameHeight, frameWidth]]), axis=0)
#currROI = np.append(currROI, np.array([[frameHeight, 0]]), axis=0)



if runningMode == 2:
	f = open(roiFileName, 'r')
	myTempROI = []
	for line in f: # read all the lines
		myTempROI.append([int(x) for x in line.split(',')])
	tempROI = np.array(myTempROI)
	currROI = np.append(currROI, tempROI, axis = 0)

else:
	settingUp = 1
	# Read frame by frame
	while (camera.isOpened()):
		# Get an image
		ret, frame = camera.read();
		# fgmask = cv2.morphologyEx(fgmask, cv2.MORPH_OPEN, kernel)
		# keep looping until the 'q' key is pressed
		while True:
			
			if toggleSmallerWindow == True:
				smallerRegion = frame[yPos:yPos+winHeight, xPos:xPos+winWidth]
				cv2.imshow("frame", smallerRegion)
			else:
				# display the image and wait for a keypress
				cv2.imshow('frame',frame)
				
			key = cv2.waitKey(1) & 0xFF
	
			# if the 'c' key is pressed, move to next frame
			if key == ord("c"):
				currFrame = camera.get(cv2.CAP_PROP_POS_FRAMES)
				camera.set(cv2.CAP_PROP_POS_FRAMES, currFrame + 1);
				refPt = [float('nan'), float('nan')]
				break
				
			# start new ROI	
			elif key == ord("n"):
				currROI = np.zeros((1,2), dtype=int)
				camera.set(cv2.CAP_PROP_POS_FRAMES, startFrame);
				frameCount = 0
				break
					
			elif key == ord("g"): # go to frame
				frameNum = askinteger('jump to frame',  'go to frame number:', parent = root)
				frameNum = max(1, min(frameNum, videoNFrames - 2))
				camera.set(cv2.CAP_PROP_POS_FRAMES, frameNum)
				break
				
			elif key == ord("s"): # show the roi
				tempROI = currROI[1:]
				cv2.polylines(frame,[tempROI],True,(255,127,127))
				
			elif key == ord("a"):
				currROI = np.zeros((1,2), dtype=int)
				currROI = np.append(currROI, np.array([[0, 0]]), axis=0)
				currROI = np.append(currROI, np.array([[1, 1]]), axis=0)
				currROI = np.append(currROI, np.array([[1, int(frameHeight)-1]]), axis=0)
				currROI = np.append(currROI, np.array([[int(frameWidth)-1, int(frameHeight)-1]]), axis=0)
				currROI = np.append(currROI, np.array([[int(frameWidth)-1, 1]]), axis=0)
				cv2.polylines(frame,[currROI],True,(255,127,127))
				
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
				
			# go to frame
			#elif key == ord("g"):
			#	cap.set(cv2.CAP_PROP_POS_FRAMES, 500);
	
			# if the 'q' key is pressed, exit completely from video
			elif key == ord("q"):
				settingUp = 0;
				break;
				
	
			frameCount += 1
			
		if settingUp == 0:
			break
			
	
	cv2.destroyAllWindows()
	camera.set(cv2.CAP_PROP_POS_FRAMES, startFrame);
	
	tempROI = currROI[1:]
	
	
	f = open(roiFileName, 'w')
	
	for t in tempROI:
	    line = ', '.join(str(x) for x in t)
	    f.write(line + '\n')
	
	f.close()


roiArea = cv2.contourArea(tempROI)
print("roi area: ", roiArea)



if runningMode == 1:
	# clean everything
	camera.release()
	cv2.destroyAllWindows()
	sys.exit()

if trackBlobs == True:
	blobsFileName = os.path.join(drive, pathName, os.path.splitext(fileName)[0] + '_blobs.csv')
	print("blobs file name: %s", blobsFileName)
	fblobs = open(blobsFileName, 'w')
	fblobs.write(fBlobsHeader)
	

history = 350 # number of frames on which the background is calculated used 350 for the 15 degrees
threshold = 12
bs = cv2.createBackgroundSubtractorMOG2(history, threshold)
#bs = cv2.bgsegm.createBackgroundSubtractorMOG()

print("Learning background")
# Learn initial background model for a while
usedHistory = min(history, videoNFrames)
for t in range (1, usedHistory):
	# Get an image
	ret, im = camera.read();
	im = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY);
	# im = cv2.morphologyEx(im, cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(5,5)))
	fgmask = bs.apply(im, 0);
	#if t%10 == 0:
	#	print(t)

# rewind the video
camera.set(cv2.CAP_PROP_POS_FRAMES, startFrame);

frameCount = 0
while(camera.isOpened()):
	particleCounter = 0
	# Get an image
	ret, im = camera.read();
	
	if ret == False:
		break
	
	imGrey = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY);
	fgmask = bs.apply(imGrey, 0.001); # it was 0.001 for 15 degrees
	
	##blur the image
	blurred = cv2.GaussianBlur(fgmask, (3, 3), 0)
	# blurred = fgmask

	# threshold the image to reveal light regions in the
	# blurred image
	thresh = cv2.threshold(blurred, bgThreshold, 255, cv2.THRESH_BINARY)[1] # sometimes 6
	
	# perform a series of erosions and dilations to remove
	# any small blobs of noise from the thresholded image
	thresh = cv2.erode(thresh, None, iterations=2)
	thresh = cv2.dilate(thresh, None, iterations=1)



	# find contours
	## check to see if we are using OpenCV 2.X
	#if imutils.is_cv2():
		#(cnts, _) = cv2.findContours(thresh.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE) 
	## check to see if we are using OpenCV 3
	#elif imutils.is_cv3():
	#(_, cnts, _) = cv2.findContours(thresh.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

	
	if False: # if the density is low, then it is good to keep this as False; at high density the watershed algorithm can instead be useful
		# compute the exact Euclidean distance from every binary
		# pixel to the nearest zero pixel, then find peaks in this
		# distance map
		D = ndimage.distance_transform_edt(thresh)
		localMax = peak_local_max(D, indices=False, min_distance=10, labels=thresh)
		
		# perform a connected component analysis on the local peaks,
		# using 8-connectivity, then appy the Watershed algorithm
		markers = ndimage.label(localMax, structure=np.ones((3, 3)))[0]
		labels = watershed(-D, markers, mask=thresh)
		# print("[INFO] {} unique segments found".format(len(np.unique(labels)) - 1))
	
	else:
		labels = ndimage.label(thresh, structure=np.ones((3, 3)))[0]

	# tempImg = thresh.copy()
	# tempImg[tempImg > 1] = 1
	# #skel, distance = medial_axis(tempImg, return_distance=True)
	# #dist_on_skel = distance * skel
	
	# skel = skeletonize(tempImg)
	# skel.dtype='uint8'
	# skel[skel == 1] = 255
	# alpha = 0.5
	# # cv2.addWeighted(imGrey, alpha, skel, 1 - alpha, 0.0, im);
	# skelRGB = cv2.cvtColor(skel,cv2.COLOR_GRAY2RGB)
	# im[skelRGB >= 1] = 255
	
	# loop over the unique labels returned by the Watershed
	# algorithm
	
	# get current frame and time in the video
	currFrame = camera.get(cv2.CAP_PROP_POS_FRAMES)
	currTime = camera.get(cv2.CAP_PROP_POS_MSEC)
	
	for label in np.unique(labels):
		# if the label is zero, we are examining the 'background'
		# so simply ignore it
		if label == 0:
			continue
 
		# otherwise, allocate memory for the label region and draw
		# it on the mask
		mask = np.zeros(blurred.shape, dtype="uint8")		
		mask[labels == label] = 255
  
		# detect contours in the mask and grab the largest one
		cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)[-2]
		
		c = max(cnts, key=cv2.contourArea) # this is the best contour, the one with largest area
		
		if cv2.contourArea(c) > minimumArea: # only run analysis on large blobs
			cv2.drawContours(im, cnts, -1, (255, 200, 40), 2)

			tempImg = mask.copy()
			tempImg[tempImg > 1] = 1
			#skel, distance = medial_axis(tempImg, return_distance=True)
			#dist_on_skel = distance * skel
			
			# skel = skeletonize(tempImg)
			# skel.dtype='uint8'
			# skel[skel == 1] = 255
			# alpha = 0.5
			# # cv2.addWeighted(imGrey, alpha, skel, 1 - alpha, 0.0, im);
			# skelRGB = cv2.cvtColor(skel,cv2.COLOR_GRAY2RGB)
			# im[skelRGB >= 1] = 255
			
			
			# draw a circle enclosing the object
			((x, y), r) = cv2.minEnclosingCircle(c)
			
			
			
			if(cv2.pointPolygonTest(tempROI, (int(x), int(y)), False) > 0):
				# cv2.circle(im, (int(x), int(y)), int(r), (0, 255, 0), 2)			
				cv2.putText(im, "+", (int(x), int(y)), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255, 127), 2)
				particleCounter += 1
				
				if trackBlobs == True:
					M = cv2.moments(c)
					x1,y1 = int(M['m10']/M['m00']), int(M['m01']/M['m00'])
					area = M['m00']
					# elongation = moment_elongation(M)
					
					if len(c)>5:
						(center,axes,orientation) = cv2.fitEllipse(c)
						# angle,Qinit=tracking.GetOrientation(orientation,Qinit)
						# angle=angle*math.pi/180
						# print(center)
						ellipsex = center[0]
						ellipsey = center[1]
						ellipsea = axes[0]
						ellipseb = axes[1]
					else:
						# print("could not determine orientation")
						orientation = float('nan')
						ellipsex = float('nan')
						ellipsey = float('nan')
						ellipsea = float('nan')
						ellipseb = float('nan')
				
					# fblobs.write("{:.3f}".format(currTime) + ", " + str(currFrame) + "," + "{:.3f}".format(x) + ", " + "{:.3f}".format(y) + "," + "{:.3f}".format(x1) + ", " + "{:.3f}".format(y1) +  ", " + "{:.3f}".format(orientation) +  ", " + "{:.3f}".format(area) +  ", " + "{:.3f}".format(ellipsex) +  ", " + "{:.3f}".format(ellipsey) +  ", " + "{:.3f}".format(ellipsea) +  ", " + "{:.3f}".format(ellipseb) + ", " + "{:.3f}".format(elongation) + '\n')
					fblobs.write("{:.3f}".format(currTime) + ", " + str(currFrame) + "," + "{:.3f}".format(x) + ", " + "{:.3f}".format(y) + "," + "{:.3f}".format(x1) + ", " + "{:.3f}".format(y1) +  ", " + "{:.3f}".format(orientation) +  ", " + "{:.3f}".format(area) +  ", " + "{:.3f}".format(ellipsex) +  ", " + "{:.3f}".format(ellipsey) +  ", " + "{:.3f}".format(ellipsea) +  ", " + "{:.3f}".format(ellipseb) + '\n')
	
	
 
	# show the output image
	if showImage > 0:
		cv2.imshow("Output", im)
	#	print(particleCounter)
	allDensities.append([int(currFrame), particleCounter, int(roiArea)])
	camera.set(cv2.CAP_PROP_POS_FRAMES, currFrame + frameJump - 1)

	
	
	##ret, labels = cv2.connectedComponents(thresh)
	## Map component labels to hue val
	#label_hue = np.uint8(179*labels/np.max(labels))
	#blank_ch = 255*np.ones_like(label_hue)
	#labeled_img = cv2.merge([label_hue, blank_ch, blank_ch])

	## cvt to BGR for display
	#labeled_img = cv2.cvtColor(labeled_img, cv2.COLOR_HSV2BGR)

	## set bg label to black
	#labeled_img[label_hue==0] = 0

	# cv2.imshow('labeled.png', labeled_img)
	# cv2.resizeWindow('labeled.png', 600,600)
	
	#cv2.drawContours(imGrey, cnts, -1, (240, 0, 159), 3)

	# if the 'q' key is pressed, exit completely from video
	key = cv2.waitKey(1) & 0xFF
	if key == ord("q"):
		stopVideo = 1;
		break
		
	frameCount += 1
	if frameCount%100 == 0:
		print("frame: ", frameCount, " of ", videoNFrames)
		snapshotFileName = os.path.join(drive, pathName, os.path.splitext(fileName)[0] + '_frame_' + str(int(currFrame)) + '_particles_' + str(particleCounter) + '_ROI_area_' + str(roiArea) + '_snapshot.png')
		cv2.polylines(im,[tempROI],True,(255,127,127))
		#cv2.imshow("Output", im)
		cv2.imwrite(snapshotFileName, im)     # save frame as PNG file	
	
	# The lines below stop collecting densities after maxNumFrames frames used
	if maxNumFrames > 0:
		if frameCount >=maxNumFrames:
			stopVideo = 1;
			break
	if frameCount >= videoNFrames: # this is introduced so that it can still save a snapshot after the end of video (if it exits because camera does not read a frame it won't be able to save anything)
		stopVdieo = 1
		break
	
	
fblobs.close()
	
	
	
	

meanCount = np.mean(allDensities, axis=0)
print(meanCount)

# print the count for each frame 
#for t in allDensities:
#    line = ', '.join(str(x) for x in t)
#    print(line + '\n')
    



f = open('particle_count.csv', 'a+')

f.write(pathName + ', ' + fileName + ', ' + "{:.3f}".format(meanCount[1]) + ", " + "{:.3f}".format(meanCount[2]) + '\n')
f.close()


# export the frames	
snapshotFileName = os.path.join(drive, pathName, os.path.splitext(fileName)[0] + '_frame_' + str(int(currFrame)) + '_particles_' + str(particleCounter) + '_ROI_area_' + str(roiArea) + '_snapshot.png')
cv2.polylines(im,[tempROI],True,(255,127,127))
#cv2.imshow("Output", im)
cv2.imwrite(snapshotFileName, im)     # save frame as PNG file	
#cv2.imwrite("frame%05d.png" % frameCount, im)     # save frame as PNG file	
		
			
# clean everything
camera.release()
cv2.destroyAllWindows()

sys.exit
