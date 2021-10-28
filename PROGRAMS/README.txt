/*******************************Source Files**************************************/
cismath.py: 
	The math utility package that contains custom data types Vec3D, Rot3D, Frame, .etc) and useful functions to solve cis problems
assignment1_utilities.py: 
	Utility functions specified to programming assignment 1. Including functions that read and compare data from files.
registration.py: 
	Contains 1 function, registration(), that computes 3D Points to 3D Points registration between 2 lists of Vec3D objects. The function returns a Frame object.
pivot_calibration.py: 
	Contains 2 function that calculates the positions of tool tips using pivot calibration method. One is for the EM probe, the other one is for the optical probe.
main.py: 
	The main file that calls functions to complete tasks of programming assignmen 1.


/********************************How to run***************************************/
Requirement: 
	Please install NumPy in the enviornment if not already.

Using commandline prompt:
	python3 main.py

Using IDE:
	Open main.py in an IDE and run the script

Then the program will automatically read datasets a-k with appropriate filenames in "Input Data" folder, print out results of each step, and finally store output in the"OUTPUT" folder in the main directore.

/************************************Credits*****************************************/
Yunxin Chen:
	assignment1_utilities.py
	registration.py
	main.py

HongYi Fan:
	cismath.py
	pivot_calibration.py
	main.py

External resources:
	NumPy (numpy.org) : numpy linear algebra and least sqare root function.

