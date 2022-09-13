# FROGS - Fast Radar-On-Glacier Simulations
# Simulate Ground Penetrating Radar data on glaciers in a fast and memory-efficient way.  
#
# Validation example of the paper "Fast 3D ground penetrating radar simulations for glaciers" 
# by J. Hunziker, E.C. Slob and J. Irving, Computers & Geosciences, 2022 (submitted). 
# 
# This file is part of FROGS. FROGS is free software: you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the Free Software Foundation, 
# either version 3 of the License, or (at your option) any later version. FROGS is distributed 
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
# more details. You should have received a copy of the GNU General Public License along with FROGS. 
# If not, see <https://www.gnu.org/licenses/>. 
# 
# When refering to FROGS in any publication, please cite the paper "Fast 3D ground penetrating 
# radar simulations for glaciers" by J. Hunziker, E.C. Slob and J. Irving, Computers & Geosciences, 
# 2022 (submitted).

import os
import numpy as np
import numpy.matlib
import math
import matplotlib.pyplot as plt
from struct import unpack
from FROGS_utils import get_plane, get_corners, mat2ascii
from parameter_class import Parameters


doplotgeo = 0 # Show geometry (1) or not (0)
doplotres = 1 # Show results (1) or not (0)
docompile = 1 # (Re)Compile C++ code (1) or not (0)
foldername = "." # Folder containing the FROGS source code
ellength = 0.5 # Side length of isquare scattering elements [m]
numcore = 16 # Amount of CPUs used for C++ FROGS


# Set FROGS parameters
params = Parameters()
params.epsr1 = 1.0 # Relative electric permittivity of the upper halfspace [-] (1.0 = air)
params.epsr2 = 3.2 # Relative electric permittivity of the lower halfspace [-] (3.2 = ice)
params.I = 1.0 # Current Amplitude [A]
params.dz = 0.5 # Infinitesimal dipole element length [m]
params.freq_c = 1e8 # Center frequency of antenna [Hz]
params.wavelet_shift = 0.0 # Time shift of wavelet [s]
params.noa = 1 # amount of source positions [-]
params.srcstartx = 0.0 # start position of source in x-direction [m]
params.srcstarty = 0.0 # start position of source in y-direction [m]
params.srcstartz = 0.0 # start position of source in z-direction [m]
params.srcdx = 0.0 # increment of source position in x-direction [m]
params.srcdy = 0.5 # increment of source position in y-direction [m]
params.srcdz = 0.0 # increment of source position in z-direction [m]
params.recstartx = 0.0 # start position of receiver in x-direction [m]
params.recstarty = 0.0 # start position of receiver in y-direction [m]
params.recstartz = 0.0 # start position of receiver in z-direction [m]
params.recdx = 0.0 # increment of receiver position in x-direction [m]
params.recdy = 0.5 # increment of receiver position in y-direction [m]
params.recdz = 0.0 # increment of receiver position in z-direction [m]
params.nf = 8192 # number of frequencies [-]
params.maxfreqin = 2e9 # highest frequency in data [Hz]
params.verbose = 2 # show no info (0), show just timing info (1), or show all info (2)
params.writemode = 2 # write output in ASCII-format (0), binary-format (1), or both (2)
params.doradpat = 1 # Circular radiation pattern (0) or actual dipole radiation pattern (1)
params.critdist = 20.0 # Critical distance for scattering element to be taken into account [m]
params.tapw = 10.0 # Width of taper inside which the weight of the scattering element is reduced [m]
params.phi_ant_src = 165.0 # Source antenna orientation [°]
params.phi_ant_rec = 165.0 # Reciever antenna orientation [°]
params.write2disc(foldername)

# Plane
xlen = 60.0 # Length of plane in x-direction [m]
ylen = 96.0 # Length of plane in y-direction [m]
xshift = 0.0 # Shift the plane in positive x-direction [m]
yshift = 0.0 # Shift the plane in positive y-direction [m]
depth = 50.0 # Depth of plane in negative z-direction [m]
dip = 0.0/180.0*math.pi # Dip rotation of plane [rad]
azi = 0.0/180.0*math.pi # Azimuthal rotation of plane [rad]
thick = 0.5 # Thickness of plane [m] 
            # (If the thickness is zero, the standard reflection 
            # coefficient is chosen by FROGS, meaning that the 
            # layer is infinitely thick.)
epsr_plane = 25.0 # Relative electric permittivity of plane [-] (81.0 = water)
epsr_below =  7.0 # Relative electric permittivity of material below plane [-] (3.2 = ice)
                  # If the plane has no thickness, it becomes infinitely
                  # thick. In that case, there is no material below the
                  # plane. This value is then disregarded. 

# The first object initializes the elempos matrix.
elempos =  get_plane(xlen,ylen,azi,dip,xshift,yshift,depth,ellength,thick,epsr_plane,epsr_below)


# # Pipe or Halfpipe
# dohalfpipe = 0 # Do halfpipe open upwards (1), open downwards (-1) or fullpipe (0)
# thick = 0.0 # Thickness of pipe [m] 
#             # (If the thickness is zero, the standard reflection 
#             # coefficient is chosen automatically, meaning that the 
#             # layer is infinitely thick.)
# epsr_pipe = 81.0 # Relative electric permittivity of pipe [-] (81.0 = water)
# epsr_below =  3.2 # Relative electric permittivity of material below pipe [-] (3.2 = ice)
#                   # If the pipe has no thickness, it becomes infinitely
#                   # thick. In that case, there is no material below the
#                   # pipe. This value is then disregarded. 
# xshift = 0.0 # Shift the pipe in positive x-direction [m]
# yshift = 1.0 # Shift the pipe in positive y-direction [m]
# depth = 2.8 # Depth of pipe in negative z-direction [m]
# rad = 0.45 # Radius of pipe [m]
# length = 2 # Length of pipe [m]
# azi = 90.0/180.0*math.pi # Azimuthal rotation of pipe [rad]

# # The second and further objects are stored in a separate matrix
# # and are then concatenated with the existing elempos matrix.
# temp =  get_pipe(rad,length,azi,xshift,yshift,depth,ellength,thick,epsr_pipe,epsr_below,dohalfpipe)
# elempos = np.concatenate((elempos,temp),axis=1)


# Write input_geometry.txt file
mat2ascii(elempos,foldername+'/input_geometry.txt')

if doplotgeo==1:
    # Ploting all the objects
    plotx, ploty, plotz =  get_corners(elempos,ellength)

    fig = plt.figure()
    ax = plt.axes(projection="3d")
    for il in range(elempos.shape[1]):
        ax.plot(plotx[il*5:il*5+5],ploty[il*5:il*5+5],plotz[il*5:il*5+5],color='blue')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')
    plt.show()


# Compile C++ FROGS
if docompile==1:
    os.system('g++ frogs.cpp -fopenmp -lfftw3 -O3 -o frogs')


# Write script to run C++ FROGS
filename_run = "run_FROGS.scr"
frun = open(filename_run,"w")
frun.write('#!/bin/bash\n')
frun.write('export OMP_NUM_THREADS=' + str(numcore) + '\n')
frun.write('./frogs\n')
frun.close()


# Make the script file executable
os.system('chmod u+x run_FROGS.scr')


# Run the C++ code
os.system('./run_FROGS.scr')


# Read the binary output file
# Open file
filename = "frogs_time_output.bin"
f = open(foldername+"/"+filename,"rb")
# Read amount of time samples andamount of source positions
ntime = unpack('=i', f.read(4))[0]
nsrc = unpack('=i', f.read(4))[0]
# Read time vector
timevec = np.zeros(ntime)
for it in range(ntime):
    timevec[it] = unpack('=d', f.read(8))[0]
# Read data
output_real = np.zeros(ntime*nsrc)
output_imag = np.zeros(ntime*nsrc)
for it in range(ntime*nsrc):
    output_real[it] = unpack('=d', f.read(8))[0]
    output_imag[it] = unpack('=d', f.read(8))[0]
output_real = output_real.reshape(nsrc,ntime).T
output_imag = output_imag.reshape(nsrc,ntime).T
# Close file
f.close()


# Load data from semi-analytical code of Evert Slob, TU Delft
traceES = np.loadtxt('trace_ES.txt', dtype='float')
ntES = 8192 # number of samples [-]
dtES = 0.3125 # dt [ns]
tvecES = np.linspace(0,ntES-1,num=ntES)*dtES;


if doplotres==1:
    # Plot the real part of the data
    fig, ax = plt.subplots()
    ax.plot(timevec*1e9, output_real[:,0]/max(output_real[:,0]), linewidth=2.0, label='FROGS')
    ax.plot(tvecES, traceES/max(traceES), linewidth=2.0, linestyle='dashed', label='semi-analytical code')
    ax.set(xlim=(580,680), xticks=np.arange(580,681,20))
    ax.legend()
    plt.xlabel('Time [ns]')
    plt.show()
