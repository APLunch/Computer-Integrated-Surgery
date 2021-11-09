/*******************************Source Files**************************************/
cismath.py: 
	The math utility package that contains custom data types Vec3D, Rot3D, Frame, .etc) and useful functions to solve cis problems
assignment1_utilities.py: 
	Utility functions specified to programming assignment 1. Including functions that read and compare data from files.
assignment2_utilities.py:
	Utility functions specified to programming assignment 2. Including functions that read and compare data from files, distortion corrections, tooltip position calculations using EM markers data frames.
registration.py: 
	Contains 1 function, registration(), that computes 3D Points to 3D Points registration between 2 lists of Vec3D objects. The function returns a Frame object.
pivot_calibration.py: 
	Contains 2 function that calculates the positions of tool tips using pivot calibration method. One is for the EM probe, the other one is for the optical probe.
main.py: 
	The main file that calls functions to complete tasks of programming assignment 1.
main_PA2.py:
	The main file that calls functions to complete tasks of programming assignment 2.
plotter.py:
	The plotting functions used to visualize points/vectors in space for debugging purposes. Not involved in actual demonstrations of PA assignments.

/********************************How to run***************************************/
Requirement: 
	Please install NumPy in the enviornment if not already.
	Plotter function requires mpl_toolkits and matplotlib. You can ignore it if not using plotters in the source code.

Using commandline prompt:
	For PA 1:
		python3 main.py
	For PA 2:
		python3 main_PA2.py

Using IDE:
	For PA1:
		Open main.py in an IDE and run the script
	For PA2:
		Open main_PA2.py in an IDE and run the script

Then the program will automatically read datasets a-j(k, if PA 1) with appropriate filenames in "Input Data" folder, print out results of each step, and finally store output in the"OUTPUT" folder in the main directore.

/************************************Credits*****************************************/
Yunxin Chen:
	assignment1_utilities.py
	assignment2_utilities.py
	registration.py
	main.py
	main_PA2.py

HongYi Fan:
	assignment2_utilities.py
	cismath.py
	pivot_calibration.py
	main.py
	main_PA2.py

External resources:
	NumPy (numpy.org) : numpy linear algebra and least sqare root function
	matplotlib: Used for DEBUG - purposed plotting

