-----------------------------------------------
The FROGS package contains the following files: 
-----------------------------------------------

frogs.cpp: The core of FROGS written in C++

FROGS_utils.py: Contains some functions used by FROGS_validation_example.py. These are

- get_plane: Function that creates a plane built from scattering elements. 

- get_pipe: Function that creates a segment of a pipe or a halfpipe built from scattering elements. Originally created to model channels. 

- get_corners: Function that calculates the position of the corners of the square scattering elements. This is needed to plot the geometry. 

- mat2ascii: Writes a matrix to the disc. It is used by the python script to write the geometry of the problem to the disc to be read by the C++ FROGS. 

FROGS_validation_example.py: Runs the validation example of the paper "Fast 3D ground penetrating radar simulations for glaciers" by J. Hunziker, E.C. Slob and J. Irving, Computers & Geosciences, 2022 (submitted). Run this file. 

LICENSE: The license under which FROGS can be used, modified and distributed. 

parameter_class.py: A python class to hold all kinds of parameters used by FROGS. 

readme.txt: This file that you are currently looking at. 

trace_ES.txt: GPR-trace created with the semi-analytical code used to validate FROGS. 


-----------------------------------------------
How to run FROGS:
-----------------------------------------------

1.) Download the Fast Fourier Transform library FFTW from https://www.fftw.org 
2.) Install the FFTW library following the instructions given in the README that you find in the FFTW package.
3.) Run FROGS_validation_example.py with the docompile flag set to 1 in order to compile the C++ code. You can do this either in an interactive developer environment such as for example Spyder or directly in the terminal by typing "python FROGS_validation_example.py"


For details about the code, we refer to the paper "Fast 3D ground penetrating radar simulations for glaciers" by J. Hunziker, E.C. Slob and J. Irving, Computers & Geosciences, 2022 (submitted).
