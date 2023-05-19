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
# Vol. 173, 2023: https://doi.org/10.1016/j.cageo.2023.105320


# Define the Parameters class to hold all kinds of parameters needed to run FROGS. 
class Parameters:
    
    epsr1 = 0.0 # Relative electric permittivity of the upper halfspace [-] (1.0 = air)
    epsr2 = 0.0 # Relative electric permittivity of the lower halfspace [-] (3.2 = ice)
    I = 0.0 # Current Amplitude [A]
    dz = 0.0 # Infinitesimal dipole element length [m]
    freq_c = 0.0 # Center frequency of antenna [Hz]
    wavelet_shift = 0.0 # Time shift of wavelet [s]
    noa = 1 # amount of source positions [-]
    srcstartx = 0.0 # start position of source in x-direction [m]
    srcstarty = 0.0 # start position of source in y-direction [m]
    srcstartz = 0.0 # start position of source in z-direction [m]
    srcdx = 0.0 # increment of source position in x-direction [m]
    srcdy = 0.0 # increment of source position in y-direction [m]
    srcdz = 0.0 # increment of source position in z-direction [m]
    recstartx = 0.0 # start position of receiver in x-direction [m]
    recstarty = 0.0 # start position of receiver in y-direction [m]
    recstartz = 0.0 # start position of receiver in z-direction [m]
    recdx = 0.0 # increment of receiver position in x-direction [m]
    recdy = 0.0 # increment of receiver position in y-direction [m]
    recdz = 0.0 # increment of receiver position in z-direction [m]
    nf = 1 # number of frequencies [-]
    maxfreqin = 0.0 # highest frequency in data [Hz]
    verbose = 1 # show no info (0), show just timing info (1), or show all info (2)
    writemode = 1 # write output in ASCII-format (0), binary-format (1), or both (2)
    doradpat = 1 # Circular radiation pattern (0) or actual dipole radiation pattern (1)
    critdist = 0.0 # Critical distance for scattering element to be taken into account [m]
    tapw = 0.0 # Width of taper inside which the weight of the scattering element is reduced [m]
    phi_ant_src = 0.0 # Source antenna orientation [°]
    phi_ant_rec = 0.0 # Reciever antenna orientation [°]
    
    
    # Function to write parameter file to disc to be read by C++ FROGS
    def write2disc(self, folder):
        f = open(folder + '/input_parameters.txt',"w")
        f.write(str(self.epsr1) + ' // Relative electric permittivity of the upper halfspace [-] (1.0 = air)\n')
        f.write(str(self.epsr2) + ' // Relative electric permittivity of the lower halfspace [-] (3.2 = ice)\n')
        f.write(str(self.I) + ' // Current Amplitude [A]\n')
        f.write(str(self.dz) + ' // Infinitesimal dipole element length [m]\n')
        f.write(str(self.freq_c) + ' // Center frequency of antenna [Hz]\n')
        f.write(str(self.wavelet_shift) + ' // Time shift of wavelet [s]\n')
        f.write(str(self.noa) + ' // amount of source positions [-]\n')
        f.write(str(self.srcstartx) + ' // start position of source in x-direction [m]\n')
        f.write(str(self.srcstarty) + ' // start position of source in y-direction [m]\n')
        f.write(str(self.srcstartz) + ' // start position of source in z-direction [m]\n')
        f.write(str(self.srcdx) + ' // increment of source position in x-direction [m]\n')
        f.write(str(self.srcdy) + ' // increment of source position in y-direction [m]\n')
        f.write(str(self.srcdz) + ' // increment of source position in z-direction [m]\n')
        f.write(str(self.recstartx) + ' // start position of receiver in x-direction [m]\n')
        f.write(str(self.recstarty) + ' // start position of receiver in y-direction [m]\n')
        f.write(str(self.recstartz) + ' // start position of receiver in z-direction [m]\n')
        f.write(str(self.recdx) + ' // increment of receiver position in x-direction [m]\n')
        f.write(str(self.recdy) + ' // increment of receiver position in y-direction [m]\n')
        f.write(str(self.recdz) + ' // increment of receiver position in z-direction [m]\n')
        f.write(str(self.nf) + ' // number of frequencies [-]\n')
        f.write(str(self.maxfreqin) + ' // highest frequency in data [Hz]\n')
        f.write(str(self.verbose) + ' // show no info (0), show just timing info (1), or show all info (2)\n')
        f.write(str(self.writemode) + ' // write output in ASCII-format (0), binary-format (1), or both (2)\n')
        f.write(str(self.doradpat) + ' // Circular radiation pattern (0) or actual dipole radiation pattern (1)\n')
        f.write(str(self.critdist) + ' // Critical distance for scattering element to be taken into account [m]\n')
        f.write(str(self.tapw) + ' // Width of taper inside which the weight of the scattering element is reduced [m]\n')
        f.write(str(self.phi_ant_src) + ' // Source antenna orientation [°]\n')
        f.write(str(self.phi_ant_rec) + ' // Reciever antenna orientation [°]\n')
        f.write('\n')
        f.write('// The order of the values has to be maintained and no lines can be inserted before the parameters!\n')
        f.close()
