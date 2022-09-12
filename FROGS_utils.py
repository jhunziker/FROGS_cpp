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

import numpy as np
import numpy.matlib
import math

def get_plane(xlen,ylen,azi,dip,xshift,yshift,depth,ellength,thick,epsr_plane,epsr_below):
    '''
    Create a plane constructed of square planar elements
    
    Input arguments: 
    - xlen: Length of plane in x-direction [m]
    - ylen: Length of plane in y-direction [m]
    - azi: Azimuthal rotation of plane [rad]
    - dip: Dip of plane [rad]
    - xshift: Shift of plane in positive x-direction [m]
    - yshift: Shift of plane in positive y-direction [m]
    - depth: Depth of plane in negative z-direction [m]
    - ellength: Length of scattering element [m]
    - thick: Thickness of plane [m] 
             (If the thicness is zero, the standard reflection 
             coefficient is chosen automatically, meaning that the 
             layer is infinitely thick.)
    - epsr_plane: Relative electric permittivity of plane [-]
    - epsr_below: Relative electric permittivity of material below plane [-]
                  If the plane has no thickness, it becomes infinitely thick. 
                  In that case, there is no material below the plane. This 
                  value is then disregarded, but still needs to be provided. 
    
    Output: 
    - elempos: 9 x nel matrix, where nel is the number of square scattering 
               elements. The nine rows are the x-, y-, and z-coordinate of 
               the center of the squares as well as the azimuth phi, the dip
               theta, the surface area and the thickness of the scattering 
               element. The last two entries are the relative electric
               permittivity of the plane and of the material below. 
    '''

    xel = int(xlen/ellength)
    yel = int(ylen/ellength)
    startx = (xlen-ellength)/2
    starty = (ylen-ellength)/2
    xvec=np.linspace(-startx,startx,xel)
    yvec=np.linspace(-starty,starty,yel)

    elempos = np.zeros((9,xel*yel))
    elempos[0,:] = np.matlib.repmat(xvec,1,yel) # x-pos of center of element [m]
    elempos[1,:] = (np.matlib.repmat(yvec.T,xel,1)).T.reshape(1,xel*yel) # y-pos of center of element [m]
    elempos[2,:] = np.zeros((1,xel*yel)) # z-pos of center of element [m]
    elempos[3,:] =  azi # azimuth of normal vector of element from x-axis [rad]
    elempos[4,:] =  dip # dip of normal vector of element from vertical [rad]
    elempos[5,:] = pow(ellength,2.0) # area of element [m^2]
    elempos[6,:] = thick # thickness of layer [m]
    elempos[7,:] = epsr_plane # Relative electric permittivity of plane [-]
    elempos[8,:] = epsr_below # Relative electric permittivity of material below plane [-]

    diprot = np.array([[math.cos(dip),0.0,math.sin(dip)],
             [0.0,1.0,0.0],
             [-math.sin(dip),0.0,math.cos(dip)]])
    elempos[0:3,:] = diprot.dot(elempos[0:3,:])

    azirot = np.array([[math.cos(azi),-math.sin(azi),0.0],
             [math.sin(azi),math.cos(azi),0.0],
             [0.0,0.0,1.0]])
    elempos[0:3,:] = azirot.dot(elempos[0:3,:])

    elempos[0,:] = elempos[0,:]+xshift
    elempos[1,:] = elempos[1,:]+yshift
    elempos[2,:] = elempos[2,:]-depth

    return elempos


def get_pipe(rad,length,azi,xshift,yshift,depth,ellength,thick,epsr_pipe,epsr_below,dohalfpipe):
    '''
    Create a pipe or a halfpipe constructed of square planar elements
    
    Input arguments: 
    - rad: Radius of pipe [m]
    - length: Length of pipe [m]
    - azi: Azimuthal rotation of pipe [rad]
    - xshift: Shift of pipe in positive x-direction [m]
    - yshift: Shift of pipe in positive y-direction [m]
    - depth: Depth of pipe in negative z-direction [m]
    - ellength: Length of scattering element [m]
    - thick: Thickness of plane [m] 
             (If the thicness is zero, the standard reflection 
             coefficient is chosen automatically, meaning that the 
             layer is infinitely thick.)
    epsr_pipe: Relative electric permittivity of pipe [-]
    epsr_below: Relative electric permittivity of material below pipe [-]
                If the pipe has no thickness, it becomes infinitely thick. 
                In that case, there is no material below the pipe. This value 
                is then disregarded, but still needs to be provided. 
    - dohalfpipe: Do halfpipe open upwards (1), open downwards (-1) or fullpipe (0)
    
    Output: 
    - elempos: 9 x nel matrix, where nel is the number of square scattering 
               elements. The nine rows are the x-, y-, and z-coordinate of 
               the center of the squares as well as the azimuth phi, the dip
               theta, the surface area and the thickness of the scattering 
               element. The last two entries are the relative electric
               permittivity of the pipe and of the material below. 
    '''

    circ = 2.0*math.pi*rad
    num_of_el = int(round(circ/ellength))
    dalpha = 2.0*math.pi/num_of_el # angular interval to have an element [rad]
    
    yel = int(length/ellength)
    start = (length-ellength)/2
    yvec = np.linspace(-start,start,yel)
    
    elempos = np.zeros((9,num_of_el*yel))
    for iy in range(yel):
        for iel in range(num_of_el):
            elempos[0,iel+iy*num_of_el] = -rad*math.cos(iel*dalpha+math.pi/2.0)
            elempos[1,iel+iy*num_of_el] = yvec[iy]
            elempos[2,iel+iy*num_of_el] = rad*math.sin(iel*dalpha+math.pi/2.0)
            elempos[3,iel+iy*num_of_el] = azi
            elempos[4,iel+iy*num_of_el] = iel*dalpha
            elempos[5,iel+iy*num_of_el] = pow(ellength,2.0)
            elempos[6,iel+iy*num_of_el] = thick
            elempos[7,iel+iy*num_of_el] = epsr_pipe
            elempos[8,iel+iy*num_of_el] = epsr_below
    
    if dohalfpipe==1:
        indices = []
        for iel in range(elempos.shape[1]):
            if elempos[2,iel]<=-depth:
                indices.append(iel)
        elempos = elempos[:,indices]
    elif dohalfpipe==-1:
        indices = []
        for iel in range(elempos.shape[1]):
            if elempos[2,iel]>=-depth:
                indices.append(iel)
        elempos = elempos[:,indices]
    
    azirot = np.array([[math.cos(azi),-math.sin(azi),0.0],
             [math.sin(azi),math.cos(azi),0.0],
             [0.0,0.0,1.0]])
    elempos[0:3,:] = azirot.dot(elempos[0:3,:])

    elempos[0,:] = elempos[0,:]+xshift
    elempos[1,:] = elempos[1,:]+yshift
    elempos[2,:] = elempos[2,:]-depth

    return elempos


def get_corners(elempos,ellength):
    '''
    For a set of square scattering elements contained in the 9 x nel matrix
    elempos, where nel is the number of elements, the location of the four 
    corners is calculated and returned in three matrices plotx, ploty and plotz. 
    They contain the x-, y- and z-coordinates of the corners, respectively. 
    Note, for plotting reasons, the first point is stored twice, so that the 
    plotting routine draws also a line from the last corner back to the first. 
    '''
    
    c2m = ellength/2.0 # Distance from center point to middle of side
    
    plotx = np.zeros(5*elempos.shape[1])
    ploty = np.zeros(5*elempos.shape[1])
    plotz = np.zeros(5*elempos.shape[1])
    for il in range(elempos.shape[1]):
        thetatemp = elempos[4,il]+math.pi/2.0
        phitemp = elempos[3,il];
    
        # Calculate the position of the midpoint of all sides of the square
        midlow = np.zeros(3)
        midlow[0] = c2m*math.sin(thetatemp)*math.cos(phitemp)
        midlow[1] = c2m*math.sin(thetatemp)*math.sin(phitemp)
        midlow[2] = c2m*math.cos(thetatemp)
        midtop = np.zeros(3)
        midtop[0] = c2m*math.sin(thetatemp-math.pi)*math.cos(phitemp)
        midtop[1] = c2m*math.sin(thetatemp-math.pi)*math.sin(phitemp)
        midtop[2] = c2m*math.cos(thetatemp-math.pi)
        # Before the shift, the midright and midleft points are in the
        # xy-plane (i.e., z = 0)
        midright = np.zeros(3)
        midright[0] = c2m*math.cos(phitemp+math.pi/2.0)
        midright[1] = c2m*math.sin(phitemp+math.pi/2.0)
        midright[2] = 0.0
        midleft = np.zeros(3)
        midleft[0] = c2m*math.cos(phitemp-math.pi/2.0)
        midleft[1] = c2m*math.sin(phitemp-math.pi/2.0)
        midleft[2] = 0.0
    
        # Calculate the corners
        botright = midlow+midright
        botleft = midlow+midleft
        topright = midtop+midright
        topleft = midtop+midleft
    
        plotx[0+il*5] = elempos[0,il]+botright[0]
        plotx[1+il*5] = elempos[0,il]+botleft[0]
        plotx[2+il*5] = elempos[0,il]+topleft[0]
        plotx[3+il*5] = elempos[0,il]+topright[0]
        plotx[4+il*5] = elempos[0,il]+botright[0]
        ploty[0+il*5] = elempos[1,il]+botright[1]
        ploty[1+il*5] = elempos[1,il]+botleft[1]
        ploty[2+il*5] = elempos[1,il]+topleft[1]
        ploty[3+il*5] = elempos[1,il]+topright[1]
        ploty[4+il*5] = elempos[1,il]+botright[1]
        plotz[0+il*5] = elempos[2,il]+botright[2]
        plotz[1+il*5] = elempos[2,il]+botleft[2]
        plotz[2+il*5] = elempos[2,il]+topleft[2]
        plotz[3+il*5] = elempos[2,il]+topright[2]
        plotz[4+il*5] = elempos[2,il]+botright[2]
    
    return plotx, ploty, plotz


def mat2ascii(mat,filename):
    '''
    Write a matrix mat to file whose name is filename. 
    '''
    
    f = open(filename,'wb')
    for istring in range(mat.shape[1]):
        outstring = ""
        for iel in range(mat.shape[0]):
            if iel==mat.shape[0]-1:
                outstring = outstring + str(mat[iel,istring]) + "\n"
            else:
                outstring = outstring + str(mat[iel,istring]) + ", "
        writestring = outstring.encode('ascii')
        f.write(writestring)
    f.close()

