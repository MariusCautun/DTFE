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


/** This file contain an example of how one can access the Delaunay tessellation of the point set using the DTFE function.

NOTE: For this code to work you must compile with the Makefile option: TRIANGULATION
*/

#include <vector>

#include "../../src/DTFE.h"
#include "../../src/input_output.cc"

using namespace std;



/* This function computes the density at the given sample point using the Delaunay cell as input. */
Real computeDensity(Cell_handle &cell,
                    Real* samplePoint)
{
    // get the vertex position difference matrix (= position(vertices!=base) - position(base); where "base" = vertex 0 of the Delaunay cell)
    Point base = cell->vertex(0)->point();  // not base stores the position of the 0-th vertex of the Delaunay cell
    Real A[NO_DIM][NO_DIM];
    vertexPositionMatrix( cell, A );
    
    // compute the inverse of A
    Real AInverse[NO_DIM][NO_DIM];
    matrixInverse( A, AInverse );
    
    // compute the density gradient
    Real dens[NO_DIM];  //matrix to store the density differences
    for (int i=0; i<NO_DIM; ++i)
        dens[i] = cell->vertex(i+1)->info().density() - cell->vertex(0)->info().density();
    Real densGrad[NO_DIM];
    matrixMultiplication( AInverse, dens, densGrad );    //computes the density gradient
    
    // get the density at the sample point
    Real result = cell->vertex(0)->info().density();
    for (int i=0; i<NO_DIM; ++i)
        result += (samplePoint[i]-base[i]) * densGrad[i];
    
    return result;
}



int main( int argc, char *argv[] )
{
    User_options userOptions;		// stores the program options supplied by the user plus additional program wide constants
    userOptions.field.triangulation = true;
    userOptions.readOptions( argc, argv );  // read the user supplied options
//     userOptions.field.density = false;
    
    
    // read the input data
    vector<Particle_data> particles;	// vector that keeps track of particle positions and their data
    vector<Sample_point> samplingCoordinates; // keep track of the sampling coordinates if the grid interpolation is done to a non-regular grid
    readInputData( &particles, &samplingCoordinates, &userOptions ); // reads particle positions, non-regular sampling coordinates (if any). Can also set members of the 'User_options' class like the box dimensions or grid size.
    
    
    // before computing the DTFE, store the particle positions in a vector since the DTFE function will delete the data in "particles"
    vector< Pvector<Real,NO_DIM> > positions;
    positions.reserve( particles.size() );
    for (size_t i=0; i<particles.size(); ++i)
        positions.push_back( particles[i].position() );
    /* // Or alternative way:
    for (size_t i=0; i<particles.size(); ++i)
    {
        Pvector<Real,NO_DIM> temp;
        for (int j=0; j<NO_DIM; ++j)
            temp[j] = particles[i].position(j);
        positions.push_back( temp );
    } */
    
    
    // compute the DTFE
    DT delaunayTriangulation;   // the Delaunay Triangulation of the point set
    Quantities quantities;	// structure that keeps vectors for all allowed field computations - only for grid quantities, not for the Delaunay triangulation itself
    DTFE( &particles, samplingCoordinates, userOptions, &quantities, delaunayTriangulation );  // this function deletes the information in 'particles'
    
    
    // find the density at the positions stored in vector "positions"
    std::cout << "Computing the density at particles' positions ... " << std::flush;
    vector<Real> particleDensity;   // variable to store the density values at each particle position
    particleDensity.reserve( positions.size() );
    for (size_t i=0; i<positions.size(); ++i)
    {
        // variables used by CGAL to locate the Delaunay cell where the sample point lies
        Locate_type lt;
        int li, lj;
        // locate the Delaunay cell where the sample point lies
#if NO_DIM==2
        Point samplePoint = Point( positions[i][0], positions[i][1] );
        Cell_handle cell = delaunayTriangulation.locate( samplePoint , lt, li );
#elif NO_DIM==3
        Point samplePoint = Point( positions[i][0], positions[i][1], positions[i][2] );
        Cell_handle cell = delaunayTriangulation.locate( samplePoint , lt, li, lj );
#endif
        
        //compute the density at the sample point
        particleDensity.push_back( computeDensity( cell, &(positions[i][0]) ) );
    }
    std::cout << "Done.\n";
    
    
    // output the desired quantities to file/files
    writeOutputData( quantities, userOptions );
    //now writting the particle positions to a file
    std::fstream outputFile;
    string filename = userOptions.outputFilename+".particleDensity";
    openOutputTextFile( outputFile, filename );
    outputFile << "# particle id\t pos_X\t pos_Y\t pos_Z\t density\n";
    for (size_t i=0; i<positions.size(); ++i)
        outputFile << i << "\t" << positions[i][0] << "\t" << positions[i][1] << "\t" << positions[i][2] << "\t" << particleDensity[i] << "\n";
    outputFile.close();
}


