macro "test" {
                run("Close All");
input = "/Users/perna/Tetrahymena_tracking/temp_experiment_files/body_size_25_degrees/";
output = "/Users/perna/Tetrahymena_tracking/temp_experiment_files/body_size_25_degrees/output_of_analysis/";
 
//setBatchMode(true);
list = getFileList(input);
for (i = 0; i < list.length; i++)
{
	print("current file: "+list[i]);
	if (endsWith(list[i], ".png"))
	{
        print("process file "+list[i]);
        action(input, output, list[i]);
	}
    else
        print("skip file "+list[i]);
}
saveAs("Results.csv");
run("Close All");
//setBatchMode(false);



 
function action(input, output, filename) {
	print("in action with file "+input+filename);
        open(input+filename);
        print(input + File.nameWithoutExtension);

// imTitle = getTitle();
//run("Duplicate...", " ");
//run("Set Scale...", "distance=1673.4276 known=752.23 pixel=1 unit=µm global");
//run("Find Edges");
//run("8-bit");
//run("Morphological Filters", "operation=Closing element=Octagon radius=8");
//run("Gaussian Blur...", "sigma=2");
////run("Threshold...");
//run("Threshold...");
//setThreshold(170, 255);
////setOption("BlackBackground", true);
//run("Convert to Mask");
//run("Fill Holes");
//selectWindow("Threshold");
//run("Close");
////run("Watershed");

//run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's median skewness kurtosis display add redirect=None decimal=3");
////run("Analyze Particles...", "size=300-Infinity circularity=0.30-1.00 show=Nothing display exclude");

//run("ROI Manager...");
//roiManager("Add");
//selectWindow(filename);
//roiManager("Measure");


nameOut = File.nameWithoutExtension;

run("8-bit");
imTitle = getTitle();
run("Duplicate...", " ");
run("Set Scale...", "distance=1673.4276 known=752.23 pixel=1 unit=µm global");
run("Find Edges");
run("8-bit");
run("Morphological Filters", "operation=Closing element=Octagon radius=8");
run("Gaussian Blur...", "sigma=2");
//run("Threshold...");
run("Threshold...");
setThreshold(170, 255);
//setOption("BlackBackground", true);
run("Convert to Mask");
run("Fill Holes");
selectWindow("Threshold");
run("Close");
//run("Watershed");
// run("Create Selection");
//selectWindow(imTitle);
//run("Restore Selection");

run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's median skewness kurtosis display add redirect=None decimal=3");
//run("Analyze Particles...", "size=300-Infinity circularity=0.30-1.00 show=Nothing display exclude");


run("ROI Manager...");
setTool("wand");

waitForUser("q"); 
selectWindow(imTitle);
if (roiManager("Count") > 0){

roiManager("Measure");
// run("Elliptic Fourier D.", "number=10 results reconstruction");
run("RGB Color");
setForegroundColor(255, 0, 0);
roiManager("Show All");
roiManager("Draw");
saveAs("Jpeg",  output+filename);
roiManager("Delete");

selectWindow("Results");
saveAs("Results", output + nameOut + ".csv");
}


//waitForUser("Select good ROI")
//run("ROI Manager...");
//roiManager("Add");

imageList = getList("image.titles");
for (i = 0; i < imageList.length; i++)
	{
		print("closing: "+imageList[i]);
	if (startsWith(imageList[i], File.nameWithoutExtension))
	{
		selectWindow(imageList[i]);
		print("closing: "+imageList[i]);
		run("Close");
	}
	}
}
}
