//BROWSES FOLDERS ORGANIZED BY '01 ORGANIZE IMAGES BY WELL', 
//counts the number of cycles, and allows the user to select neurons at a specific cycle from a desired range of wells
//creates a folder of text data 'neuronPos'

// Define selection variables (box size)
h = 8;
w = 8;
//Prompt User for Organized Data Folder - Input
sourceDir = getDirectory("Select the Organized Well Folder");
folderList = getFileList(sourceDir);
//Sort folders from low to high - necessary?
Array.sort(folderList);
//Create Folder for Neuron Position Data - Output
newDir= sourceDir + "\\.." + "/neuronPos/";
File.makeDirectory(newDir);
saveDir = newDir;
//Ask User what well to start with
startFile = getNumber("found "+folderList.length+" wells, start selecting at well?", 1);
startFile = startFile - 1;
//Ask User what well to end with
endFile = getNumber("found "+folderList.length+" wells, end selecting at well?", folderList.length);
endFile = endFile - 1;
//Ask user which cycle they want to see [ADD cycle# check and prompt]
clickFile = getNumber("Which video do you want to use for clicking?", 1);
clickFile = clickFile - 1;
//Loop from one well folder to the next, displaying an image for manual Neuron Selection
for (i=startFile; i<=endFile; i++) {
	//set path to well specific folder
	path = sourceDir+folderList[i];
	imageList = getFileList(path);
	//UNNECESSARY?
	posFile = path+"neuronPos.txt";
	//Sort cycles form high to low - necessary?
  	Array.sort(imageList);
	//Open Specific Cycle image from Well Folder 
	open(path+imageList[clickFile],100);
	//Improve neuron visibility
        	run("Enhance Contrast", "saturated=0.35");
	//DISABLED SCALING
	//run("Set Scale...", "distance=1 known=1 pixel=1 unit=um");
	//Neuron Selection Variables
	//dummy variables for capturing mouse data
	leftButton=16;
	rightButton=4;
	x2=-1;
	y2=-1;
	z2=-1;
	flags2=-1;
	getCursorLoc(x, y, z, flags);
	wasLeftPressed = false;
	//Animal Count
	a=0;


	//WHILE loop to detect position of left clicks when single image is displayed by FOR loop
	while (flags&rightButton==0){
		getCursorLoc(x, y, z, flags);

		//Detect left Mouse click
		if (flags&leftButton!=0) {
		// Wait for it to be released
		wasLeftPressed = true;
		} else if (wasLeftPressed) {

			wasLeftPressed = false;
			if (x!=x2 || y!=y2 || z!=z2 || flags!=flags2) {
			//Increment animal counter
			a=a+1;

			//Draw Rectangle Around clicked Neuron - size defined at start of file
			makeRectangle(x - w/2, y - h/2, w, h);

			run("Add Selection...", "stroke=yellow width=1 fill=0");

			//Display Coordinates clicked in Log
			print("well " + (i+1) + " x " + x + " y " + y + " a " + a + " v " + clickFile);
			}
		//Displayed Rectangle around pixel left clicked by mouse
		}
//While Loop interrupted by clicking right button (!=0).
//Note – doesn’t advance until subsequent left click – glitch with nested IF?
	}
	//focus on Log window
	selectWindow("Log");
	//Save Log of clicked XY coordinates and animal count as text file
	saveAs("Text", saveDir+"neuronPos"+(i+1));
	//Clear Log and close windows
	run("Clear Results");
	print("\\Clear");
	run("Close All");
//go to next well image
}
