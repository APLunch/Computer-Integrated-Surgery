/*******************************Source Files**************************************/
cismath.py: 
	The math utility package that contains custom data types Vec3D, Rot3D, Frame, Mesh, Triangle .etc) and useful functions to solve cis problems
assignment3_utilities.py: 
	Utility functions specified to programming assignment 3 and 4. Including functions that read and compare data from files, and mesh loading.
registration.py: 
	Contains 1 function, registration(), that computes 3D Points to 3D Points registration between 2 lists of Vec3D objects. The function returns a Frame object.
main_PA3.py:
	The main file that calls functions to complete tasks of programming assignment 3.
main_PA4.py:
	The main file that calls functions to complete tasks of programming assignment 4.
plotter.py:
	The plotting functions used to visualize points/vectors in space for debugging purposes. Not involved in actual demonstrations of PA assignments.

/********************************How to run***************************************/
Requirement: 
	Please install NumPy in the enviornment if not already.
	Plotter function requires mpl_toolkits and matplotlib. You can ignore it if not using plotters in the source code.

Using commandline prompt:
		python3 main_PA4.py

Using IDE:
		Open main_PA4.py in an IDE and run the script

Then the program will automatically read datasets A-H and J with appropriate filenames in "Input Data" folder, print out results of each step, and finally store output in the"OUTPUT" folder in the main directore.

/************************************Credits*****************************************/
Yunxin Chen:
	assignment3_utilities.py
	registration.py
	main_PA3.py
	main_PA4.py
	cismath.py
	
HongYi Fan:
	cismath.py
	main_PA3.py
	main_PA4.py
	plotter.py

External resources:
	NumPy (numpy.org) : numpy linear algebra and least sqare root function
	matplotlib: Used for DEBUG plotting
	kDTree Diagram: 
    	Figure 1 in PA3_Report.pdf. source = [GPL, https://commons.wikimedia.org/w/index.php?curid=589909]

