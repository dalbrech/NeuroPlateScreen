//ORGANIZES IMAGES in WELL NUMBERED FOLDERS AFTER COMPLETION OF IMAGING CONTROLLED BY SCRIPT '02 WELLPLATE IMAGING'
//THIS IS THE FIRST MACRO TO RUN TO BEGIN ANALYSIS

//Prompt the user to select the input folder
dir = getDirectory("Select the Folder of Unsorted Well Images"); 
//prompt the user to select an output folder � must be done manually with the �New Folder� command�
SaveDir = getDirectory("Choose a Folder to save Sorted Well Images "); 
//Create file list 
list = getFileList(dir); 
print("number of files in the Directory ="+ list.length); 
//initialize variable for checking names in loop
PREVIOUSWELL="X";
//Loop through each file of the input directory
for (i=0; i<list.length; i++ ) { 
	//unnecessary statement? - TESTING FOR TXT FILES - FEATURE TO BE FIXED
	//before running this plugin the user must manually move the Experiment Info .txt generated by GUIRL
	if (endsWith(list[i], ".txt")){ 	
	}
	//Find the Well Number in the file name and move it to the right folder
	else{ 
		//Split the Name by "." and look at the first Entry 	
		FILENAME = split(list[i],"."); 	
		// find the 001, 002, 003 designation for the file and isolate it
		WELLID=substring(FILENAME[0],28,lengthOf(FILENAME[0])-14);			 
//EL note � convert WELLID to  INT in an array of wells to keep track, more efficient folder creation?
		print(WELLID);	
			//Check well number against check string
			if(WELLID!=PREVIOUSWELL){ 
				//Establish directory for a well
				// Make directory name 
				DirName = WELLID;
				//Create that directory 
//EL note � this seems redundant/inefficient if the list of files is by cycle � it would remake the directory every time it goes to another cycle -see how the file list is sorted
				File.makeDirectory(SaveDir+File.separator+DirName); 
				// Move the file to that Dir 
				//THIS IS CUTTING NOT COPYING
				File.rename(dir + list[i],SaveDir + DirName + File.separator + list[i]); 
//Update check string to well that has been given a folder and had a file moved into it
				PREVIOUSWELL=WELLID; 
			}
//if this file has a well number that already has a folder it can just be put there
else{ 
				File.rename(dir + list[i],SaveDir + DirName + File.separator + list[i]); 
			} 
	} 
//End of file loop � go to next file
}
// http://imagej.1557.x6.nabble.com/Select-every-other-image-and-move-into-specific-folder-td5005416.html