
Run the various scripts in the following order:
1) Run run_video_tracking_in_all_subfolders.py, make sure that the first time you run with select (to only select the roi) and the second time you run with analyse (parameter -m)
2) Run run_link_blobs_in_all_subfolders.py
3) run extract_particle_information_from_trajectories (check FPS etc.)
4) If you want, run overlay_trajectory_on_video.py to check the result
5) If necessary, merge_files.sh to produce a single file with information on all the particles
6) fit_schoolfield_equation_to_long_term_speed_responses.R

List of files:
run_video_tracking_in_all_subfolders.py
This python script simply navigates all the subfolders and then every time that it finds a video with extension .avi it runs a system command which consists in running the video_tracking_tetrahymena_with_shape.py on that video file. The parameter -m  is important because it changes the behaviour of the video_tracking_tetrahymena_with_shape.py script:
If m is 1, the script asks to select a ROI and then quits
if m is 2, then the script expects a ROI already and runs the tracking
if m is 0, the script asks first to select a ROI and then runs the tracking
In general, if we work on multiple videos we probably want first to select a ROI in each video, that is, run with -m 1 and then we want to run the script again with -m 2 and we leave it to run unattended



video_tracking_tetrahymena_with_shape.py
This script runs the video tracking on one video and during the image processing steps it also collects some measurements of the cell shape of tetrahymena which are stored among the other pieces of information. In particular, the script reads a video and saves a blobs file that has the following columns:
time, frame, xCircle, yCircle, xCentroid, yCentroid, orientation, area, xEllipse, yEllipse, aEllipse, bEllipse
0.000, 1.0,1635.000, 536.000,1634.000, 536.000, 164.280, 43.000, 1634.951, 536.740, 5.605, 12.167
33.315, 2.0,1634.000, 536.500,1633.000, 536.000, 169.516, 43.500, 1633.158, 537.030, 5.317, 12.906
66.630, 3.0,1635.500, 536.000,1634.000, 536.000, 19.624, 66.500, 1635.256, 536.001, 9.756, 14.439

time is the time in milliseconds
frame is the frame number
xCircle and yCircle are the coordinates of a circle fitted around the shape of each tetrahymena (with the opencv function cv2.minEnclosingCircle)
xCentroid and yCentroid are the coordinates of the centroid of each tetrahymena (with the opencv function 
M = cv2.moments(c)
x1,y1 = int(M['m10']/M['m00']), int(M['m01']/M['m00']))
orientation is the orientation of the tetrahymena axis
(center,axes,orientation) = cv2.fitEllipse(c)
xEllipse and yEllipse are the coordinates of the centre of the ellipse
aEllispe and bEllipse are the length of the two axes of the ellipse
the code is:
(center,axes,orientation) = cv2.fitEllipse(c)
ellipsex = center[0]
ellipsey = center[1]
ellipsea = axes[0]
ellipseb = axes[1]



run_link_blobs_in_all_subfolders.py
This is another convenience script, to run the link_blobs_into_tracks.py on all the blobs files



link_blobs_into_tracks.py
This script simply puts together the blobs into trajectories based on proximity of positions between one frame and the next. The script doesn't do much, but it relies on trackpy

It produces a table with trajectories with the following columns:
frame	time	frame	xCircle	yCircle	x	y	orientation	mass	xEllipse	yEllipse	aEllipse	bEllipse	particle
1	0	1	679.5	1046	679	1046	169.468	185.5	679.252	1047.82	11.442	23.924	1
1	0	1	1698.5	1023	1698	1023	117.683	146	1697.786	1023.687	11.249	18.021	2
1	0	1	1465.5	979.5	1464	978	146.152	231.5	1464.964	977.962	12.718	24.307	3

The main difference between this file and the blobs is that now there is a particle column that has the same name for everything that is considered to be the same particle (the same tetrahymena)

Some columns have changed name: x instead of xCentroid, y instead of yCentroid and mass instead of area.

Inside the code of link_blobs_into_tracks.py we find some hidden parameters that determine how long a particle can disappear from the image and still be kept in memory, how many pixels it can move from one frame and the next, and what is the minimum duration of a trajectory in order to be considered valid. These parameters become important if the code is reused on videos with very different settings (e.g. different frame rate, etc.)

t = tp.link_df(data, 20, memory=3)
t1 = tp.filter_stubs(t, 20)





extract_particle_information_from_trajectories.R
This R script does a bit of extra analysis of trajectories, such as for instance calculating the speed of individual tetrahymenas. In doing so, it has to convert from the current units of trajectories (pixels and frames) into more meaningful units of micrometres and seconds. For this reason it is important to check at the beginning of the file that the conversions are correct for the frame rate of the video and for the magnification of the microscope.


The script produces another modified trajectories_with_aggregate_info file, which is similar to the original trajectory file but now the column of "mass", "aEllipse" and "bEllipse" are in units of micrometres, and are calculated as the median value over the entire trajectory rather than on each frame individually. In addition there are two columns with "orientation" and "autocorrelationTimeFrame". Orientation is the orientation calculated from the direction of movement (rather than from the shape), while autocorrelationTimeFrame is the typical time after which tetrahymena has completed one cycle of body oscillation. This is measured by first looking at the difference of orientation (cosine of the angle) between time T and time T plus or minus a lag. When tetrahymena has completed one cycle, it should have roughly the same orientation and so the cos(angle) should be close to 1. Now we want to find all the lags for which the autocorrelation is maximum again and see how distant they are from each other, or better we find the lag that corresponds to the first local maximum.
