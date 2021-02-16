w = getWidth();
h = getHeight();
selectWindow("ChangCardioTemplate_bin4.tif");
path = getDirectory("image");
if (path==""){
    path = "not available";
}    
else{
	myDir = path + "roiExports";
	File.makeDirectory(myDir);
}
	
print("Path: " + myDir);

nROIs = roiManager("count");
for(i=0;i<nROIs;i++){
	newImage("Untitled", "8-bit black", w, h, 1);
	selectWindow("Untitled");
	roiManager("Select",i);
	roiName = Roi.getName();
	roiManager("Fill");
	oName = myDir + File.separator + roiName;
	saveAs("Tiff", oName);
	close();
}

