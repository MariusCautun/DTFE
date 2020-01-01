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


/*!
This header contains the functions for reading the user options from the command line and also the class that keeps track of this options.
*/


#ifndef USER_OPTIONS_HEADER
#define USER_OPTIONS_HEADER

#include <string>
#include <vector>
#include <cstdlib>
#include <ctime> 

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "define.h"
#include "miscellaneous.h"
#include "box.h"
#include "message.h"



// structure to keep track what fields to compute
struct Field
{
    bool triangulation;
    bool density;
    bool velocity;
    bool velocity_gradient;
    bool velocity_divergence;
    bool velocity_shear;
    bool velocity_vorticity;
    bool velocity_std;
    bool scalar;
    bool scalar_gradient;
    
    Field()
    { triangulation = false; density = false; velocity = false; velocity_gradient = false; velocity_divergence = false; velocity_shear = false; velocity_vorticity = false; velocity_std = false; scalar = false; scalar_gradient = false; }
    
    bool updateChoices(std::string choice,
                       std::string str_triang, std::string str_den, std::string str_vel, std::string str_grad, std::string str_div, std::string str_shear, std::string str_vort, std::string str_velstd, std::string str_scalar, std::string str_scalarGrad)
    {
        if ( choice.compare(str_triang)==0 ) triangulation = true;
        else if ( choice.compare(str_den)==0 ) density = true;
        else if ( choice.compare(str_vel)==0 ) velocity = true;
        else if ( choice.compare(str_grad)==0 ) velocity_gradient = true;
        else if ( choice.compare(str_div)==0 ) velocity_divergence = true;
        else if ( choice.compare(str_shear)==0 ) velocity_shear = true;
        else if ( choice.compare(str_vort)==0 ) velocity_vorticity = true;
        else if ( choice.compare(str_velstd)==0 ) velocity_std = true;
        else if ( choice.compare(str_scalar)==0 ) scalar = true;
        else if ( choice.compare(str_scalarGrad)==0 ) scalar_gradient = true;
        else return false;
        return true;
    }
    
    bool selected()
    { return ( density or velocity or velocity_gradient or velocity_divergence or velocity_shear or velocity_vorticity or velocity_std or scalar or scalar_gradient ); }
    
    bool selectedVelocityDerivatives()
    { return ( velocity_divergence or velocity_shear or velocity_vorticity ); }
    void deselectVelocityDerivatives()
    { velocity_divergence = false; velocity_shear = false; velocity_vorticity = false; }
    bool selectedVelocity()
    { return ( velocity or velocity_gradient or velocity_divergence or velocity_shear or velocity_vorticity or velocity_std ); }
    void deselectVelocity()
    { velocity = false; velocity_gradient = false; velocity_divergence = false; velocity_shear = false; velocity_vorticity = false;  velocity_std = false; }
    
    bool selectedScalar()
    { return ( scalar or scalar_gradient ); }
    void deselectScalar()
    { scalar = false; scalar_gradient = false; }
};


// structure to hold user specified options
struct User_options
{
    std::string inputFilename; // name of input particle positions
    std::string outputFilename;// name/root name of output file/files
    
    
    std::vector<size_t> gridSize;// grid size along the 3 directions
    bool   periodic;      // true if the particle data is in a periodic boxBox
    Box    boxCoordinates;// keeps track of the coordinates of the particle data box (of the full box)
    bool   userGivenBoxCoordinates;// true if the user inserted the coordinates of the box as an option to the program
    int    inputFileType;   // stores the type of the input file
    std::vector<int>   readParticleData;       // stores which particle properties to read (true = read those properties, false = don't read)
    std::vector<int>   readParticleSpecies;    // stores which particle species to read (true = read the data for those species, false = don't read)
    int    outputFileType;  // stores the type of the output file
    
    
    Field  uField;        // "unaveraged field" - specifies which field to interpolate to grid at the position of the sampling points
    Field  aField;        // "averaged field" - specifies which field to interpolate to grid by taking the field volume average over the full grid cell
    
    
    bool   partitionOn;   // true if the user selects the partition option
    std::vector<size_t> partition;// gives the grid along which to split the particle data if the Delaunay triangulation is to be applied on parts of the data (due to time and memory problems)  
    int    partNo;        // which file should be outputed by the program (index in 'fileGrid')
    
    
    bool   paddingOn;     // true if there is set a padding value
    bool   paddingMpcOn;  // true if the padding value is set in Mpc and not as relative values with respect to the box length
    Real   paddingParticles;// gives the number of particles to be used for padding when partitioning the particle data
    Box    paddingLength; // this gives the the padding length (in Mpc or relative box lengths) along each direction (both for left and right edges)
    bool   testPaddedBoundaries; // if true, it inserts dummy points on the boundaries to find out if the padded region is enough to give a full Delaunay triangulation inside the region of interest
    
    
    bool   regionOn;      // true if to compute the quantities only in a user defined region of the full box
    bool   regionMpcOn;   // true if the region coordinates where inserted in Mpc and not as a fraction of absolute box length
    Box    region;        // the coordinates of the user defined region where to compute the quantities of interest
    
    
    // variables for the averaging computation
    int    method;        // method to volume average the fields to grid
    int    noPoints;      // how many points to use to sample the average of the field in each grid cell
    bool   noPointsOn;    // true if the user specified the number of sampling points
    Real   averageDensity;// normalization value for the density in the simulation box (if not given by the user, computed by the program as the average density)
    size_t randomSeed;    // random seed used for the Monte Carlo interpolation
    
    
    
    // variables for light cone computations
    std::vector<Real> originPosition;   // the position of the light cone origin
    bool              redshiftConeOn;   // true if to interpolate the fields to a light cone grid
    Box               redshiftCone;     // gives the boundaries of the light cone grid (r_min, r_max, theta_min, theta_max, psi_min, psi_max)
    
    
    // additional options
    std::string configFilename;// the program options can be supplied in a file too. This keeps track of the file name where the options where supplied.
    bool   CIC;           // true if to use CIC grid interpolation instead of DTFE
    bool   TSC;           // true if to use TSC grid interpolation instead of DTFE
    bool   SPH;           // true if to use SPH grid interpolation instead of DTFE
    int    SPH_neighbors; // keep track of the number of neighbors used for the SPH grid interpolation
    Real   MpcValue;      // the value of 1Mpc in input units
    bool   extensive;     // variable used for scalar field computations using the TSC or SPH method - true if the scalar fields are extensive (if false than the scalar fields are intensive)
    int    verboseLevel;  // keep track of the verbose level of the program (see header 'message.h' for additional details)
    Real   randomSample;  // the size of the random subsample if the users wishes to use only a subset of all the data for computations
    size_t poisson;       // if !=0 - generate a random sample of particles instead of reading the positions from a file
    std::vector<std::string> additionalOptions; // variable that can be used by the user to supply additional options to the code in a very simple way
    bool   transformToRedshiftSpaceOn;          // true if to transform from position to redshift space 
    std::vector<Real> transformToRedshiftSpace; // which coordinates of the velocity to use to transform to redshift space (it uses H=100 h km/s /Mpc)
    
    
    
    //! this options cannot be set from the command line - they are set during runtime
    Box    paddedBox;      // this gives the dimensions of the padded box
    std::vector<Real> fullBoxOffset; // stores the positions of the full box left side along each axis (redundant information with 'boxCoordinates')
    std::vector<Real> fullBoxLength; // stores the full box length along each axis (redundant information with 'boxCoordinates')
    
    bool   DTFE;           // true if DTFE interpolation method
    bool   userDefinedSampling; // true if the user inserted his on sampling coordinates
    int    noProcessors;   // the number of threads executing the program in parallel
    int    threadId;       // the thread id for parallel processes
    Real   totalTime;      // constant to keep track of the total CPU time used by the program
    std::string programOptions; // string that stores all program options supplied by the user
    
    
    
    
    //! functions - functions' body is found in the file "user_options.cc"
    User_options();
    void readOptions(int argc, char *argv[], bool getFileNames=true, bool showOptions=true );   // reads the user supplied program options
    void addOptions(po::options_description &allOptions,
                    po::options_description &visibleOptions,
                    po::positional_options_description &p); // defines what are the available program options
    void helpInformation( po::options_description &visibleOptions, char *progName ); //print help information with what are the available user options
    void shortHelp( char *progName );
    void printOptions(); // prints what values the program will use for the program options
    void updateFullBox(Box &newBox); // updates the value of the full box: 'boxCoordinates', 'fullBoxOffset' and 'fullBoxLength'
    void updateEntries(size_t const noTotalParticles,
                       bool userSampling);   // used to update, after reading the input data file, some of the members of this class: paddedBox, fullBoxOffset, fullBoxLength; Also does error checking of member values in the class.
    void updatePadding(size_t const noParticles); // computes the value of the padding length
};





#endif


