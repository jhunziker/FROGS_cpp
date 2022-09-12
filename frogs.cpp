// FROGS - Fast Radar-On-Glacier Simulations
// Simulate Ground Penetrating Radar data on glaciers in a fast and memory-efficient way.  
//
// This file is part of FROGS. FROGS is free software: you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation, 
// either version 3 of the License, or (at your option) any later version. FROGS is distributed 
// in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
// more details. You should have received a copy of the GNU General Public License along with FROGS. 
// If not, see <https://www.gnu.org/licenses/>. 
// 
// When refering to FROGS in any publication, please cite the paper "Fast 3D ground penetrating 
// radar simulations for glaciers" by J. Hunziker, E.C. Slob and J. Irving, Computers & Geosciences, 
// 2022 (submitted).

#include <omp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <math.h>
#include <fftw3.h>
#include <chrono>
using namespace std;
using namespace std::chrono;

#define PI 3.141592653589793

// Electric permittivity in vaccum [F/m]
#define EPS0 8.854187817620389e-12

// Magnetic permeability in vaccum [H/m] 
#define MU0 4.0*PI*1e-7

// Speed of light in vacuum 
#define LIGHTSPEED 299792458.0

// Impedance of free space [Ohm]
#define ETA MU0*LIGHTSPEED

class parameters {
    public:
    // Data members:
    int noa, nf; // Number of antenna positions [-], number of frequencies [-]
    int verbose; // Flag for printing info to screen [-]
    int writemode; // Write output in ASCII-format (0), binary format (1), or both (2)
    int doradpat; // Circular radiation pattern (0) or actual dipole radiation pattern (1)
    double epsr1, epsr2; // Relative electric permittivity of the upper and lower halfspace 
    double I, dz; // Current Amplitude [A] and infinitesimal dipole element length [m]
    double freq_c, maxfreqin; // Center frequency of antenna [Hz], higest frequency in data [Hz]
    double wavelet_shift; // Time shift of wavelet [s] 
    double srcstartx, srcstarty, srcstartz; // Starting position of antenna [m]
    double srcdx, srcdy, srcdz; // Increments for antenna position [m]
    double recstartx, recstarty, recstartz; // Starting position of antenna [m]
    double recdx, recdy, recdz; // Increments for antenna position [m]
    double critdist; // Critical distance for scattering element to be taken into account [m]
    double tapw; // Width of taper inside which the weight of the scattering element is reduced [m]
    double phi_ant_src, phi_ant_rec; // Orientation of antennas [°]

    // Constructors:
    parameters(double val1, double val2, double val3, double val4, double val5, double val6, int val7, double val8, double val9, double val10, double val11, double val12, double val13, double val14, double val15, double val16, double val17, double val18, double val19, int val20, double val21, int val22, int val23, int val24, double val25, double val26, double val27, double val28){
        epsr1 = val1; 
        epsr2 = val2; 
        I = val3; 
        dz = val4;
        freq_c = val5;
        wavelet_shift = val6;
        noa = val7;
        srcstartx = val8;
        srcstarty = val9;
        srcstartz = val10;
        srcdx = val11;
        srcdy = val12;
        srcdz = val13;
        recstartx = val14;
        recstarty = val15;
        recstartz = val16;
        recdx = val17;
        recdy = val18;
        recdz = val19;
        nf = val20;
        maxfreqin = val21;
        verbose = val22;
        writemode = val23;
        doradpat = val24;
        critdist = val25;
        tapw = val26;
        phi_ant_src = val27;
        phi_ant_rec = val28;
    }
    
    parameters(string file_param){
        size_t found;
        string line;
        ifstream filein(file_param);

        // read epsr1
        getline(filein,line);
        found = line.find("//");
        epsr1=stof(line.substr(0,found-1));
        // read epsr2
        getline(filein,line);
        found = line.find("//");
        epsr2=stof(line.substr(0,found-1));
        // read I 
        getline(filein,line);
        found = line.find("//");
        I=stof(line.substr(0,found-1));
        // read dz
        getline(filein,line);
        found = line.find("//");
        dz=stof(line.substr(0,found-1));
        // read freq_c
        getline(filein,line);
        found = line.find("//");
        freq_c=stof(line.substr(0,found-1));
        // read wavelet_shift 
        getline(filein,line);
        found = line.find("//");
        wavelet_shift=stof(line.substr(0,found-1));
        // read noa 
        getline(filein,line);
        found = line.find("//");
        noa=stoi(line.substr(0,found-1));
        // read srcstartx 
        getline(filein,line);
        found = line.find("//");
        srcstartx=stof(line.substr(0,found-1));
        // read srcstarty 
        getline(filein,line);
        found = line.find("//");
        srcstarty=stof(line.substr(0,found-1));
        // read srcstartz 
        getline(filein,line);
        found = line.find("//");
        srcstartz=stof(line.substr(0,found-1));
        // read srcdx 
        getline(filein,line);
        found = line.find("//");
        srcdx=stof(line.substr(0,found-1));
        // read srcdy 
        getline(filein,line);
        found = line.find("//");
        srcdy=stof(line.substr(0,found-1));
        // read srcdz 
        getline(filein,line);
        found = line.find("//");
        srcdz=stof(line.substr(0,found-1));
        // read recstartx 
        getline(filein,line);
        found = line.find("//");
        recstartx=stof(line.substr(0,found-1));
        // read recstarty 
        getline(filein,line);
        found = line.find("//");
        recstarty=stof(line.substr(0,found-1));
        // read recstartz 
        getline(filein,line);
        found = line.find("//");
        recstartz=stof(line.substr(0,found-1));
        // read recdx 
        getline(filein,line);
        found = line.find("//");
        recdx=stof(line.substr(0,found-1));
        // read recdy 
        getline(filein,line);
        found = line.find("//");
        recdy=stof(line.substr(0,found-1));
        // read recdz 
        getline(filein,line);
        found = line.find("//");
        recdz=stof(line.substr(0,found-1));
        // read nf 
        getline(filein,line);
        found = line.find("//");
        nf=stoi(line.substr(0,found-1));
        // read maxfreqin
        getline(filein,line);
        found = line.find("//");
        maxfreqin=stof(line.substr(0,found-1));
        // read verbose
        getline(filein,line);
        found = line.find("//");
        verbose=stoi(line.substr(0,found-1));
        // read writemode
        getline(filein,line);
        found = line.find("//");
        writemode=stoi(line.substr(0,found-1));
        // read doradpat 
        getline(filein,line);
        found = line.find("//");
        doradpat=stoi(line.substr(0,found-1));
        // read critdist
        getline(filein,line);
        found = line.find("//");
        critdist=stof(line.substr(0,found-1));
        // read tapw
        getline(filein,line);
        found = line.find("//");
        tapw=stof(line.substr(0,found-1));
        // read phi_ant_src
        getline(filein,line);
        found = line.find("//");
        phi_ant_src=stof(line.substr(0,found-1));
        // read phi_ant_rec
        getline(filein,line);
        found = line.find("//");
        phi_ant_rec=stof(line.substr(0,found-1));
  
        filein.close();
    }

    // Member functions:
    parameters& print_info(){
        cout << " " << endl;
        cout << "The following parameters have been read:" << endl;
        cout << "(The order is the same as in the file input_parameters.txt)" << endl;
        cout << "epsr1 = " << epsr1 << endl;
        cout << "epsr2 = " << epsr2 << endl;
        cout << "I = " << I << endl;
        cout << "dz = " << dz << endl;
        cout << "freq_c = " << freq_c << endl;
        cout << "wavelet_shift =  " << wavelet_shift << endl;
        cout << "noa = " << noa << endl;
        cout << "srcstartx = " << srcstartx << endl;
        cout << "srcstarty = " << srcstarty << endl;
        cout << "srcstartz = " << srcstartz << endl;
        cout << "srcdx = " << srcdx << endl;
        cout << "srcdy = " << srcdy << endl;
        cout << "srcdz = " << srcdz << endl;
        cout << "recstartx = " << recstartx << endl;
        cout << "recstarty = " << recstarty << endl;
        cout << "recstartz = " << recstartz << endl;
        cout << "recdx = " << recdx << endl;
        cout << "recdy = " << recdy << endl;
        cout << "recdz = " << recdz << endl;
        cout << "nf = " << nf << endl;
        cout << "maxfreqin = " << maxfreqin << endl;
        cout << "verbose = " << verbose << endl;
        cout << "writemode = " << writemode << endl;
        cout << "doradpat = " << doradpat << endl;
        cout << "critdist = " << critdist << endl;
        cout << "tapw = " << tapw << endl;
        cout << "phi_ant_src = " << phi_ant_src << endl;
        cout << "phi_ant_rec = " << phi_ant_rec << endl;
        cout << " " << endl;
        return *this;
    }
};

class geometry_data {
    public:
    // Data members: 
    vector<double> midx;
    vector<double> midy;
    vector<double> midz;
    vector<double> azimuth;
    vector<double> dip;
    vector<double> elarea;
    vector<double> thick;
    vector<double> epsr_el;
    vector<double> epsr_below;

    // Constructors:
    geometry_data(double val1, double val2, double val3, double val4, double val5, double val6, double val7, double val8, double val9){
        midx.emplace_back(val1);
        midy.emplace_back(val2);
        midz.emplace_back(val3);
        azimuth.emplace_back(val4);
        dip.emplace_back(val5);
        elarea.emplace_back(val6);
        thick.emplace_back(val7);
        epsr_el.emplace_back(val8);
        epsr_below.emplace_back(val9);
    }

    geometry_data(string file_geo){
        ifstream filein(file_geo);
        streampos end;
        filein.seekg(0, ios::end);
        end = filein.tellg();
        filein.seekg(0, ios::beg);

        int nov = 9;
        string line;
        size_t found1, found2;
        vector<double> temp(nov);
        while (filein.tellg()<end){
            getline(filein,line);
            found1 = 0;
            for (int il=0; il<nov; il++){
                if (il<nov-1){
                    found2 = found1+1+line.substr(found1+1,line.size()).find(",");
                    temp[il]=stof(line.substr(found1,found2-found1));
                    found1 = found2+1;
                } else {
                    temp[il]=stof(line.substr(found1,line.size()));
                }
            }
            midx.emplace_back(temp[0]);
            midy.emplace_back(temp[1]);
            midz.emplace_back(temp[2]);
            azimuth.emplace_back(temp[3]);
            dip.emplace_back(temp[4]);
            elarea.emplace_back(temp[5]);
            thick.emplace_back(temp[6]);
            epsr_el.emplace_back(temp[7]);
            epsr_below.emplace_back(temp[8]);
        }
        filein.close();
    }

    // Member functions:
    void get_element(int isel, vector<double> &curel){
        curel[0] = midx[isel];
        curel[1] = midy[isel];
        curel[2] = midz[isel];
        curel[3] = azimuth[isel];
        curel[4] = dip[isel];
        curel[5] = elarea[isel];
        curel[6] = thick[isel];
        curel[7] = epsr_el[isel];
        curel[8] = epsr_below[isel];
    }

    // Calculate the horizontal distance between a scattering element 
    // and an antenna position
    void calc_horz_dist(int isel, const vector<double> &ant, double &horz_dist){
        horz_dist = sqrt(pow(ant[0]-midx[isel],2)+pow(ant[1]-midy[isel],2));
    }

    // Calculate taper value for scattering element
    void calc_tapval(int isel, double horz_dist, const parameters params, double &tapval){
        if (horz_dist<=params.critdist-params.tapw) {// Element is taken fully into account
            tapval = 1.0;
        } else if (horz_dist>=params.critdist) {// Element is not taken into account
            tapval = 0.0;
        } else {// Inside taper
            double tapx = (params.critdist-horz_dist)/params.tapw*PI;
            tapval = pow(1.0-((cos(tapx)+1.0)/2.0),2);
        }
    }

    void check_consistency(){
        if (midx.size()!=midy.size()) {
            cerr << "Inconsistent amount of y-coordinates received." << endl;;
            throw "size inconsistency";
        } else if (midx.size()!=midz.size()) {
            cerr << "Inconsistent amount of z-coordinates received." << endl;
            throw "size inconsistency";
        } else if (midx.size()!=azimuth.size()) {
            cerr << "Inconsistent amount of azimuth angles received."<< endl;
            throw "size inconsistency";
        } else if (midx.size()!=dip.size()) {
	    cerr << "Inconsistent amount of dip angles received." << endl;
            throw "size inconsistency";
        } else if (midx.size()!=elarea.size()) {
	    cerr << "Inconsistent amount of side lengths received." << endl;
            throw "size inconsistency";
        } else if (midx.size()!=thick.size()) {
	    cerr << "Inconsistent amount of thicknesses received." << endl;
            throw "size inconsistency";
        } else if (midx.size()!=epsr_el.size()) {
	    cerr << "Inconsistent amount of layer permittivities received." << endl;
            throw "size inconsistency";
        } else if (midx.size()!=epsr_below.size()) {
	    cerr << "Inconsistent amount of below-layer permittivities received." << endl;
            throw "size inconsistency";
        }
    }
};

class radar_data {
    public: 
    // Data members: 
    vector<complex<double>> E1; // Ex or Etheta for current scattering element
    vector<complex<double>> E2; // Ey or Ephi for current scattering element
    vector<complex<double>> E3; // Ez or Er for current scattering element
    vector<complex<double>> E1src; // Source radiation for current scattering element
    vector<complex<double>> E2src; // Source radiation for current scattering element
    vector<complex<double>> E3src; // Source radiation for current scattering element
    vector<complex<double>> E1rec; // Receiver radiation for current scattering element
    vector<complex<double>> E2rec; // Receiver radiation for current scattering element
    vector<complex<double>> E3rec; // Receiver radiation for current scattering element
    vector<double> freq; // frequency vector
    vector<double> omega; // angular frequency vector
    int nf; // number of frequencies
    double maxfreq; // maximum frequency

    // Constructors: 

    radar_data(int nfin, double maxfreqin){
        nf = nfin;
        maxfreq = maxfreqin;
        double df = maxfreq/(nf-1);
	for (int il=0; il<nf-1; il++) {// zero frequency is not considered
            E1.emplace_back(1.0,0.0);
            E2.emplace_back(1.0,0.0);
            E3.emplace_back(1.0,0.0);
            E1src.emplace_back(0.0,0.0);
            E2src.emplace_back(0.0,0.0);
            E3src.emplace_back(0.0,0.0);
            E1rec.emplace_back(0.0,0.0);
            E2rec.emplace_back(0.0,0.0);
            E3rec.emplace_back(0.0,0.0);
            freq.emplace_back((il+1)*df);
            omega.emplace_back(2.0*PI*(il+1)*df);
        }
    }

    // Member functions:
    radar_data& reinitialize(){
        // Reinitialize for loop over scattering elements
        for (int il=0; il<nf-1; il++) {
            E1[il] = {1.0,0.0};
            E2[il] = {1.0,0.0};
            E3[il] = {1.0,0.0};
            E1src[il] = {0.0,0.0};
            E2src[il] = {0.0,0.0};
            E3src[il] = {0.0,0.0};
            E1rec[il] = {0.0,0.0};
            E2rec[il] = {0.0,0.0};
            E3rec[il] = {0.0,0.0};
        }
        return *this;
    }

    radar_data& radiation_pattern(const vector<double> orientvec, const parameters params, int dosrc){
        // Equations 1 to 7 from Arcone, S. A., Numerical studies of the radiation patterns of resistively loaded dipoles, 1995, Journal of Applied Geophysics, 33, p. 39-52. 
        complex<double> imag_unit(0.0,1.0);
        double k1, k2, r, phi, theta, theta_c, v, n, frac1, frac2, temp1, temp2;
        complex<double> factor, K1, K2, cfrac1, cfrac2, E1tmp, E2tmp, E3tmp;

        r = orientvec[2];
        if (dosrc==1){
            phi = orientvec[1]+params.phi_ant_src/180.0*PI;
        } else {
            phi = orientvec[1]+params.phi_ant_rec/180.0*PI;
        }
        theta = orientvec[0]; 
        v = 1.0/sqrt(params.epsr2*EPS0*MU0);
        n = LIGHTSPEED/v;
        theta_c = asin(1.0/n);
        factor = imag_unit*params.I*params.dz*ETA;
        for (int il=0; il<nf-1; il++){
            k1 = omega[il]*sqrt(params.epsr1*EPS0*MU0);
            k2 = omega[il]*sqrt(params.epsr2*EPS0*MU0);
            K1 = factor*k1*exp(imag_unit*(k1*r))/(2.0*PI*r);
            K2 = factor*k2*exp(imag_unit*(k2*r))/(2.0*PI*r);
            if (params.doradpat==1) {// Actual radiation pattern
                if (theta < PI/2) {// The halfspace of air
                    frac1 = pow(cos(theta),2)/(cos(theta)+sqrt(n*n-pow(sin(theta),2)));
                    temp1 = pow(sin(theta),2)*cos(theta)*(cos(theta)-sqrt(n*n-pow(sin(theta),2)));
                    temp2 = n*n*cos(theta)+sqrt(n*n-pow(sin(theta),2));
                    frac2 = temp1/temp2;
                    E1tmp = K1*cos(phi)/n*(frac1-frac2);
                    E2tmp = (-K1)/n*cos(theta)*sin(phi)/(cos(theta)+sqrt(n*n-pow(sin(theta),2)));
                } else {// The lower halfspace
                    if (theta >= PI-theta_c){
                       temp1 = sqrt(1-n*n*pow(sin(theta),2)) + n*cos(theta);
                       temp2 = n*sqrt(1-n*n*pow(sin(theta),2)) - cos(theta);
                       frac1 = temp1/temp2;
                       frac2 = pow(cos(theta),2)/(sqrt(1-n*n*pow(sin(theta),2)) - n*cos(theta));
                       E1tmp = K2*cos(phi)*(pow(sin(theta),2)*cos(theta)*frac1-frac2);
                       E2tmp = K2*cos(theta)*sin(phi)/(sqrt(1-n*n*pow(sin(theta),2))-n*cos(theta));
                    } else { // theta < pi-theta_c
                        cfrac1 = sqrt(n*n*pow(sin(theta),2)-1)-imag_unit*n*cos(theta);
                        cfrac1 = cfrac1/(n*sqrt(n*n*pow(sin(theta),2)-1)+imag_unit*cos(theta));
                        cfrac2 = pow(cos(theta),2)/(sqrt(n*n*pow(sin(theta),2)-1)+imag_unit*n*cos(theta));
                        E1tmp = K2*cos(phi)*(pow(sin(theta),2)*cos(theta)*cfrac1+imag_unit*cfrac2);
                        E2tmp = -imag_unit*K2*cos(theta)*sin(phi)/(sqrt(n*n*pow(sin(theta),2)-1)+imag_unit*n*cos(theta));
                    }
                }
            } else { // Equal radiation in all directions
                if (theta < PI/2) {
                    E1tmp = K1;
                    E2tmp = K1;
                } else {
                    E1tmp = K2;
                    E2tmp = K2;
                } 
            }
            if (dosrc==1){
                E1src[il] = E1tmp;
                E2src[il] = E2tmp;
                E3src[il] = E3tmp;
            } else {
                E1rec[il] = E1tmp;
                E2rec[il] = E2tmp;
                E3rec[il] = E3tmp;
            }
        }
        return *this;
    }

    radar_data& rotate(double angle, int axis, int field){
        double rotmat[3][3];
        if (axis==0) {// rotation around x-axis
            rotmat[0][0] = 1.0;
            rotmat[0][1] = 0.0;
            rotmat[0][2] = 0.0;
            rotmat[1][0] = 0.0;
            rotmat[1][1] =  cos(angle);
            rotmat[1][2] = -sin(angle);
            rotmat[2][0] = 0.0;
            rotmat[2][1] =  sin(angle);
            rotmat[2][2] =  cos(angle);
        } else if (axis==1) {// rotation around y-axis
            rotmat[0][0] =  cos(angle); 
            rotmat[0][1] = 0.0;
            rotmat[0][2] =  sin(angle);
            rotmat[1][0] = 0.0;
            rotmat[1][1] = 1.0;
            rotmat[1][2] = 0.0;
            rotmat[2][0] = -sin(angle);
            rotmat[2][1] = 0.0;
            rotmat[2][2] =  cos(angle);
        } else if (axis==2) {// rotation around z-axis
            rotmat[0][0] =  cos(angle); 
            rotmat[0][1] = -sin(angle);
            rotmat[0][2] = 0.0;
            rotmat[1][0] =  sin(angle);
            rotmat[1][1] =  cos(angle);
            rotmat[1][2] = 0.0;
            rotmat[2][0] = 0.0;
            rotmat[2][1] = 0.0;
            rotmat[2][2] = 1.0;
        } else {
            cerr << "Rotation axis has to be 0, 1 or 2." << endl;
            throw "rotation axis error";
        }
        if (field==1) {//rotate E1src, E2src, E3src
            vector<complex<double>> tmpE(3);
	    for (int il=0; il<nf-1; il++) {
                tmpE[0] = E1src[il];
                tmpE[1] = E2src[il];
                tmpE[2] = E3src[il];
                E1src[il] = -1.0*(rotmat[0][0]*tmpE[0] + rotmat[0][1]*tmpE[1] + rotmat[0][2]*tmpE[2]);
                E2src[il] = -1.0*(rotmat[1][0]*tmpE[0] + rotmat[1][1]*tmpE[1] + rotmat[1][2]*tmpE[2]);
                E3src[il] =  1.0*(rotmat[2][0]*tmpE[0] + rotmat[2][1]*tmpE[1] + rotmat[2][2]*tmpE[2]);
            }
        } else if (field==2) {//rotate E1rec, E2rec, E3rec
            vector<complex<double>> tmpE(3);
	    for (int il=0; il<nf-1; il++) {
                tmpE[0] = E1rec[il];
                tmpE[1] = E2rec[il];
                tmpE[2] = E3rec[il];
                E1rec[il] = -1.0*(rotmat[0][0]*tmpE[0] + rotmat[0][1]*tmpE[1] + rotmat[0][2]*tmpE[2]);
                E2rec[il] = -1.0*(rotmat[1][0]*tmpE[0] + rotmat[1][1]*tmpE[1] + rotmat[1][2]*tmpE[2]);
                E3rec[il] =  1.0*(rotmat[2][0]*tmpE[0] + rotmat[2][1]*tmpE[1] + rotmat[2][2]*tmpE[2]);
            }
        } else {//rotate E1, E2, E3
            vector<complex<double>> tmpE(3);
	    for (int il=0; il<nf-1; il++) {
                tmpE[0] = E1[il];
                tmpE[1] = E2[il];
                tmpE[2] = E3[il];
                E1[il] = rotmat[0][0]*tmpE[0] + rotmat[0][1]*tmpE[1] + rotmat[0][2]*tmpE[2];
                E2[il] = rotmat[1][0]*tmpE[0] + rotmat[1][1]*tmpE[1] + rotmat[1][2]*tmpE[2];
                E3[il] = rotmat[2][0]*tmpE[0] + rotmat[2][1]*tmpE[1] + rotmat[2][2]*tmpE[2];
            }
        }
        return *this;
    }

    radar_data& combine_radpat(){
	for (int il=0; il<nf-1; il++) {
            E1[il] = E1src[il]*E1rec[il];
            E2[il] = E2src[il]*E2rec[il];
            E3[il] = E3src[il]*E3rec[il];
        }
        return *this;
    }

    radar_data& reflect(const parameters params, const double elarea, const double thick, const double epsr3, const double epsr4){
        double gamma2, gamma3, gamma4;

        if (fabs(thick)<=1e-6) {// Standard reflection coefficient
            double rte, rtm;
            for (int il=0; il<nf-1; il++){
                gamma2 = omega[il]*sqrt(params.epsr2*EPS0*MU0); // vertical wavenumber of ice
                gamma3 = omega[il]*sqrt(epsr3*EPS0*MU0); // vertical wavenumber of scatterer
                rte = (gamma2-gamma3)/(gamma2+gamma3);
                rtm = (gamma2-gamma3)/(gamma2+gamma3);
                E1[il] = E1[il]*rte*elarea;
                E2[il] = E2[il]*rte*elarea;
                E3[il] = E3[il]*rtm*elarea;
            }
        } else {// Bradford and Deeds (2006) reflection coefficient
            complex<double> rte, rtm, temp1, temp2;
            complex<double> imag_unit(0.0,1.0);
            for (int il=0; il<nf-1; il++){
                gamma2 = omega[il]*sqrt(params.epsr2*EPS0*MU0); // vertical wavenumber of ice
                gamma3 = omega[il]*sqrt(epsr3*EPS0*MU0); // vertical wavenumber of scatterer
                gamma4 = omega[il]*sqrt(epsr4*EPS0*MU0); // vertical wavenumber below scatterer
                temp1 =  gamma2-gamma4-imag_unit*(gamma2*gamma4/gamma3-gamma3)*tan(gamma3*thick);
                temp2 = (gamma2+gamma4-imag_unit*(gamma2*gamma4/gamma3+gamma3)*tan(gamma3*thick));
                rte = temp1/temp2;
                //cout << il << " " << temp1 << " " << temp2 << " " << rte << endl;
                temp1 = gamma2*epsr4*EPS0-gamma4*params.epsr2*EPS0-imag_unit*(gamma2*gamma4*epsr3*EPS0/gamma3-gamma3*params.epsr2*EPS0*epsr4/epsr3)*tan(gamma3*thick);
                temp2 = gamma2*epsr4*EPS0+gamma4*params.epsr2*EPS0-imag_unit*(gamma2*gamma4*epsr3*EPS0/gamma3+gamma3*params.epsr2*EPS0*epsr4/epsr3)*tan(gamma3*thick);
                rtm = temp1/temp2;
                E1[il] = E1[il]*rte*elarea;
                E2[il] = E2[il]*rte*elarea;
                E3[il] = E3[il]*rtm*elarea;
            }
        }
        return *this;
    }

    radar_data& sum(double tapval, vector<complex<double>> &fEtot){
        for (int il=0; il<nf-1; il++){
            fEtot[il] = fEtot[il] + tapval*(E1[il] + E2[il] + E3[il]);
        }
        return *this;
    }

};

class output_data {
    public: 
    // Data members: 
    vector<double> fEtotreal, fEtotimag;
    vector<complex<double>> fEtot; // total electric field in frequency domain 
    vector<complex<double>> tE; // total electric field in time domain 
    vector<double> freq; // frequency vector
    vector<double> omega; // angular frequency vector
    int nf; // number of frequencies
    double maxfreq; // maximum frequency

    // Constructors:
    output_data(int nfin, double maxfreqin){
        nf = nfin;
        maxfreq = maxfreqin;
        double df = maxfreq/(nf-1);
        for (int il=0; il<nf-1; il++) {// zero frequency is not considered
            fEtotreal.emplace_back(0.0);
            fEtotimag.emplace_back(0.0);
            fEtot.emplace_back(0.0,0.0);
            freq.emplace_back((il+1)*df);
            omega.emplace_back(2.0*PI*(il+1)*df);
        }
        for (int il=0; il<2*nf; il++) {
            tE.emplace_back(0.0,0.0);
        }
    }

    // Member functions:
    output_data& reinitialize(){
        // Reinitialize for loop over antenna positions
        for (int il=0; il<nf-1; il++) {
            fEtotreal.emplace_back(0.0);
            fEtotimag.emplace_back(0.0);
            fEtot[il] = {0.0,0.0};
        }
        return *this;
    }

    output_data& merge_real_imag(){
        complex<double> imag_unit(0.0,1.0);
        for (int il=0; il<nf-1; il++) {
            fEtot[il] = fEtotreal[il]+imag_unit*fEtotimag[il];
        }
        return *this;
    }

    output_data& convolve_wavelet(const parameters params){
        complex<double> imag_unit(0.0,1.0);
        complex<double> wavelet;
        double freqrat;
        for (int il=0; il<nf-1; il++){
            freqrat = -1.0*pow(freq[il]/params.freq_c,2.0);
            wavelet = freqrat*exp(freqrat)*exp(imag_unit*omega[il]*params.wavelet_shift);
            fEtot[il] = fEtot[il]*wavelet;
        }
        return *this;
    }

    output_data& fft(){// wrapper for FFTW
        double df = maxfreq/(nf-1);
        fftw_complex fft_in[2*nf], fft_out[2*nf];
        fftw_plan plan;
        plan = fftw_plan_dft_1d(2*nf, fft_in, fft_out, FFTW_BACKWARD, FFTW_ESTIMATE);
        for (int il=0; il<2*nf; il++) {
            if (il==0){
                fft_in[il][0] = 0.0;
                fft_in[il][1] = 0.0;
            } else if (il<nf){
                fft_in[il][0] = real(fEtot[il-1]);
                fft_in[il][1] = imag(fEtot[il-1]);
            } else if (il==nf){
                fft_in[il][0] = 0.0;
                fft_in[il][1] = 0.0;
            } else {
                fft_in[il][0] = real(conj(fEtot[nf-(il-nf+1)]));
                fft_in[il][1] = imag(conj(fEtot[nf-(il-nf+1)]));
            }
        }
        fftw_execute(plan);
        for (int il=0; il<2*nf; il++) {// normalize
            tE[il] = {fft_out[il][0]*df,fft_out[il][1]*df};
        }
        fftw_destroy_plan(plan);
        fftw_cleanup();
        return *this;
    }

    output_data& fftshift(){
        vector<complex<double>> temp(2*nf);
        for (int il=0; il<2*nf; il++) {
            if (il<nf){
                temp[il] = tE[nf+il]; // nf:2*nf-1
            } else {
                temp[il] = tE[il-nf]; // 0:nf-1
            }
        }
        for (int il=0; il<2*nf; il++) {
            tE[il] = temp[il];
        }
        return *this;
    }

    output_data& flip(){
        vector<complex<double>> temp(2*nf);
        for (int il=0; il<2*nf; il++) {
            temp[il] = tE[2*nf-il-1];
        }
        for (int il=0; il<2*nf; il++) {
            tE[il] = temp[il];
        }
        return *this;
    }

};

void antenna_scatterer_orientation(vector<double> &orientvec, const vector<double> &ant, const vector<double> &curel){
    // r
    orientvec[2] = sqrt(pow(ant[0]-curel[0],2.0)+pow(ant[1]-curel[1],2.0)+pow(ant[2]-curel[2],2.0));
    // phi
    orientvec[1] = PI-(atan2(curel[1]-ant[1],curel[0]-ant[0]));
    // theta
    orientvec[0] = acos((curel[2]-ant[3])/orientvec[2]);
}

void write_output(vector<complex<double>> &outputE, parameters params){
    double df = params.maxfreqin/(params.nf-1);
    double dt = 1.0/(2.0*params.nf*df);
    double time;
    if (params.writemode==0 || params.writemode==2){
        ofstream fileout;
        fileout.open ("frogs_time_output.txt");
        fileout << "src-number time real(E1) imag(E1)" << endl;
        for (int isrc=0; isrc<params.noa; isrc++) {
            time = -1.0*params.nf*dt;
            for (int itime=0; itime<2*params.nf; itime++) {
                fileout << isrc << " " << time << " ";
                fileout << real(outputE[isrc*2*params.nf+itime]) << " ";
                fileout << imag(outputE[isrc*2*params.nf+itime]) << " ";
                time += dt;
            }
        }
        fileout.close();
    } 
    if (params.writemode==1 || params.writemode==2){
        int ntime = 2*params.nf;
        ofstream fileout;
        fileout.open("frogs_time_output.bin", ios::out | ios::binary);
        fileout.write(reinterpret_cast<char*>(&ntime), sizeof(int));
        fileout.write(reinterpret_cast<char*>(&params.noa), sizeof(int));
        time = -1.0*params.nf*dt;
        for (int itime=0; itime<ntime; itime++) {
            fileout.write(reinterpret_cast<char*>(&time), sizeof(double));
            time += dt;
        }
        fileout.write(reinterpret_cast<char*>(&outputE[0]), outputE.size()*sizeof(complex<double>));
        fileout.close();
    }
}

int main(){
    // Measure time at program start
    auto start = high_resolution_clock::now(); 

    // Read input parameters
    parameters params("input_parameters.txt");

    // Print startup info
    if (params.verbose >= 1){
        cout << "Fast Radar-On-Glacier Simulation" << endl;
        cout << "FROGS - last updated on 10.9.2022" << endl;
    }

    // Provide some information
    if (params.verbose == 2) params.print_info();

    // Read geometry of scattering elements
    geometry_data elempos("input_geometry.txt"); 
    elempos.check_consistency();

    // Set up vector of all source positions
    vector<double> srcxvec(params.noa);
    vector<double> srcyvec(params.noa);
    vector<double> srczvec(params.noa);
    srcxvec[0] = params.srcstartx;
    srcyvec[0] = params.srcstarty;
    srczvec[0] = params.srcstartz;
    for (int isrc=1; isrc<params.noa; isrc++) {
        srcxvec[isrc] = srcxvec[isrc-1]+params.srcdx;
        srcyvec[isrc] = srcyvec[isrc-1]+params.srcdy;
        srczvec[isrc] = srczvec[isrc-1]+params.srcdz;
    }

    // Set up vector of all receiver positions
    vector<double> recxvec(params.noa);
    vector<double> recyvec(params.noa);
    vector<double> reczvec(params.noa);
    recxvec[0] = params.recstartx;
    recyvec[0] = params.recstarty;
    reczvec[0] = params.recstartz;
    for (int irec=1; irec<params.noa; irec++) {
        recxvec[irec] = recxvec[irec-1]+params.recdx;
        recyvec[irec] = recyvec[irec-1]+params.recdy;
        reczvec[irec] = reczvec[irec-1]+params.recdz;
    }

    // Print antenna information
    if (params.verbose == 2) {
        cout << "The following source-receiver positions are considered:" << endl;
        for (int isrc=0; isrc<params.noa; isrc++) {
            cout<<"("<<srcxvec[isrc]<<","<<srcyvec[isrc]<<","<<srczvec[isrc]<<") | ";
            cout<<"("<<recxvec[isrc]<<","<<recyvec[isrc]<<","<<reczvec[isrc]<<")"<<endl;
        }
        cout << " " << endl;
    }

    // Initialize matrices to contain output electromagnetic fields
    vector<complex<double>> outputE(2*params.nf*params.noa);

    // Initialize radar data and frequency vector
    radar_data data(params.nf,params.maxfreqin);

    // Initialize output data
    output_data data_out(params.nf,params.maxfreqin);

    // Initialize source and receiver vectors 
    vector<double> src(3);
    vector<double> rec(3);

    // Initialize vector to hold the geometry of the current scattering element
    vector<double> curel(9);

    // Initialize scalars to hold the horizontal distance between the current scattering element
    // and the current antenna positions as well as the corresponding taper value. 
    double srcdist, recdist, tapval;

    // Initialize vectors for geometry relation between the antennas and the scattering element
    vector<double> src_orientvec(3);
    vector<double> rec_orientvec(3);

    for (int isrc=0; isrc<params.noa; isrc++) {// Loop over source positions
        // Reinitialize the data for the new source position
        data.reinitialize();
        data_out.reinitialize();

        // Update source and receiver position
        src[0] = srcxvec[isrc]; 
        src[1] = srcyvec[isrc]; 
        src[2] = srczvec[isrc]; 
        rec[0] = recxvec[isrc]; 
        rec[1] = recyvec[isrc]; 
        rec[2] = reczvec[isrc]; 

        #pragma omp parallel for firstprivate(data,curel,srcdist,recdist,src_orientvec,rec_orientvec) shared(data_out,params,elempos,src,rec)
        for (int il=0; il<elempos.midx.size(); il++) {// Loop over scattering elements
            // Check if current scattering element is inside critical distance from antennas
            elempos.calc_horz_dist(il, src, srcdist);
            elempos.calc_horz_dist(il, rec, recdist);
            if (min(srcdist,recdist)<params.critdist) {
                // Run if current scattering element is inside critical distance from antennas

                // Reinitialize the data for the new element
                data.reinitialize();

                // Get current scattering element
                elempos.get_element(il, curel);

                // Calculate distance and angles between source and scattering element
                antenna_scatterer_orientation(src_orientvec, src, curel);

                // Calculate the source-antenna radiation pattern
                data.radiation_pattern(src_orientvec, params, 1);

                // Calculate distance and angles between receiver and scattering element
                antenna_scatterer_orientation(rec_orientvec, rec, curel);

                // Take the receiver-antenna radiation pattern into account
                data.radiation_pattern(rec_orientvec, params, 2);

                // Rotate the source radiation pattern
                // from the spherical coordinates into cartesian coordinates
                data.rotate(src_orientvec[0],1,1).rotate(-PI+src_orientvec[1],2,1);

                // Rotate the receiver radiation pattern
                // from the spherical coordinates into cartesian coordinates
                data.rotate(src_orientvec[0],1,2).rotate(-PI+src_orientvec[1],2,2);

                // Combine the source and receiver radiation patterns
                data.combine_radpat();

                // Rotate the electric field into the coordinate system of the plane of the reflector
                data.rotate(elempos.dip[il],1,0).rotate(elempos.azimuth[il],2,0);

                // Reflect the electric field at the reflector
                data.reflect(params, curel[5], curel[6], curel[7], curel[8]);
        
                // Rotate the scattered electric field back into the original cartesian coordinates
                data.rotate(-elempos.azimuth[il],2,0).rotate(-elempos.dip[il],1,0);

                // Calculate taper value for scattering element
                elempos.calc_tapval(il, min(srcdist,recdist), params, tapval);

                // Sum up the fields from the various scattering elements
                for (int il=0; il<data.nf-1; il++){
                    #pragma omp atomic update
                    data_out.fEtotreal[il] = data_out.fEtotreal[il] + tapval*(data.E1[il].real() + data.E2[il].real() + data.E3[il].real());
                    #pragma omp atomic update
                    data_out.fEtotimag[il] = data_out.fEtotimag[il] + tapval*(data.E1[il].imag() + data.E2[il].imag() + data.E3[il].imag());
                    // Atomic update apparently does not work with complex variables
                    // data_out.fEtot[il] = data_out.fEtot[il] + tapval*(data.E1[il] + data.E2[il] + data.E3[il]);
                }
            }
        } 

        // Merge real and imaginary parts of fEtot
        data_out.merge_real_imag();

        // Convolve the data in the frequency domain with a wavelet
        data_out.convolve_wavelet(params);

        // Fourier transform to time domain
        data_out.fft();

        // Rearrange data
        data_out.flip().fftshift();

        // Copy data to output matrix
        for (int itime=0; itime<2*params.nf; itime++) {
            outputE[isrc*2*params.nf+itime] = data_out.tE[itime];
        }
    } 

    // Write output matrix to disc
    write_output(outputE, params);

    // Measure time at program completion
    if (params.verbose >= 1){
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start); 
        cout << "FROGS runtime: " << duration.count()*1e-6 << " s" << endl; 
    }

    return 1;
}
