/*
 *  Copyright (c) 2011       Marius Cautun
 *
 *                           Kapteyn Astronomical Institute
 *                           University of Groningen, the Netherlands
 *
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <vector>

#include "DTFE.h"
#include "input_output.cc"

using namespace std;


int main( int argc, char *argv[] )
{
    User_options userOptions;		// stores the program options supplied by the user plus additional program wide constants
    userOptions.readOptions( argc, argv );  // read the user supplied options
    
    
    
    // read the input data
    vector<Particle_data> particles;	// vector that keeps track of particle positions and their data
    vector<Sample_point> samplingCoordinates; // keep track of the sampling coordinates if the grid interpolation is done to a non-regular grid
    readInputData( &particles, &samplingCoordinates, &userOptions ); // reads particle positions, non-regular sampling coordinates (if any). Can also set members of the 'User_options' class like the box dimensions or grid size.
    
    
    
    // compute the DTFE
    Quantities uQuantities;	// structure that will store the output quantities - i.e. the fields at the sampling points (NOTE: "uQuantities" stands for "unaveraged quantities" meaning that this fields were NOT averaged over the sampling cell)
    Quantities aQuantities;     // structure that will store the volume averaged output quantities - i.e. the fields volume averaged over the sampling cells (NOTE: "aQuantities" stands for "averaged quantities")
    DTFE( &particles, samplingCoordinates, userOptions, &uQuantities, &aQuantities );  // this function deletes the information in 'particles'
    
    
    
    // output the desired quantities to file/files
    writeOutputData( uQuantities, aQuantities, userOptions );
}
