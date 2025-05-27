# E_faecalis_erm_treatment


Code is divided by the figure it produced.

Figure 2: Original image processing (use 3D data in treated/untreated tiff folder)
Use data from Untreated/Treated files also provided 
This code slices the images along the x-axis and computes the protective structure volume of each slice
It also registers how many slices have a protective structure within them
There are 2 main editable parameters:
1. Grey-scale threshold: tune the grey-scale threshold so that the density plot starts at approximately the 
	observed base density of the microscopy images (larger structures need lower thresholds as the 
	dye has trouble penetrating)
2. Epsilon: this is a parameter used in Matlab’s function dbscan
	epsilon is usually 15 (stress distance) unless structures seem to the user to we well separated
	but close together. In this case lower epsilon to 10.

Figure 5: testing BC’s (use substacks in Figure 5 folder)
This code can run tiff slices on the three antibiotic boundary conditions tested in figure 5 of the paper
The output is a movie and a binary matrix
The input needs to be a binary matrix constructed using Figure5_ICcreation.m from a 2-D tiff z-slice

Figure 6: Validation (use substacks in Figure 6 folder)
Use provided 2D tiff z-slices (or your own!) to run a simulation on them. The volume of protective structures before and after antibiotic simulated treatment are computed. Be sure to modify the grey-scale value and 
density cut off as directed in the file. This file does not produce movies, but will tell you the average beginning and end volumes of the protective structures. 
