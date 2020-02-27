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


#include <fstream>
#include <cstdio>
#include "user_options.h"


#ifndef INPUT_FILE_DEFAULT
#define INPUT_FILE_DEFAULT 101
#endif
#ifndef OUTPUT_FILE_DEFAULT
#define OUTPUT_FILE_DEFAULT 101
#endif
#ifndef MPC_UNIT
#define MPC_UNIT 1.
#endif




// class constructor
User_options::User_options()
{
    periodic = false;
    boxCoordinates.assign(0.);
    userGivenBoxCoordinates = false;
    inputFileType = INPUT_FILE_DEFAULT;     // the default input file - see makefile for definition
    outputFileType = OUTPUT_FILE_DEFAULT;   // the default output file - see makefile for definition
    readParticleData.assign( 10, 0 );
    readParticleSpecies.assign( 10, 0 );
    readParticleData[0] = 1; readParticleData[1] = 1; readParticleData[2] = 1;
    readParticleSpecies[1] = 1;
    
    regionOn = false;
    regionMpcOn = false;
    
    partitionOn = false;
    partNo = -1;
    
    paddingOn = false;
    paddingMpcOn = false;
    paddingParticles = Real(5.);
    testPaddedBoundaries = true;
    
    // density options
    method = 1;
    noPoints = 100;
    noPointsOn = false;
    averageDensity = Real(-1.);
    
    // variables for light cone computation
    redshiftConeOn = false;
    
    // additional options
    DTFE = true;
    NGP = false;
    CIC = false;
    TSC = false;
    SPH = false;
    MpcValue = Real(MPC_UNIT);
    extensive = false;
    verboseLevel = 3;
    randomSample = Real(-1.);
    poisson = 0;
    
    
    // internal variables
    paddedBox.assign(0.);
    userDefinedSampling = false;
    noProcessors = 1;
    totalTime = Real(0.);
    programOptions = "";
}



/* Use the boost program options library to store the available user options.
The options are stored in two different structures:
    allOptions - stores all available options
    visibleOptions - stores the options visible to the user - this will be printed via the help message
*/
void User_options::addOptions(po::options_description &allOptions,
                              po::options_description &visibleOptions,
                              po::positional_options_description &p)
{
    // divide the options in several categories
    po::options_description mainOptions("Main options");
    mainOptions.add_options()
            ("help,h", "produce summary of help message. For more detailed help use the '--full_help' option.")
            ("full_help", "produce detailed help message. For more detailed help information consult the documentation.")
            ("grid,g", po::value< std::vector<size_t> >(&(this->gridSize))->multitoken(), "choose grid size along each direction (e.g. '-g 256' for a 256^3 grid; '-g 256 128 512' for different gridsize along each direction).")
            ("box", po::value< std::vector<Real> >()->multitoken(), "the coordinates of the box encompasing all the particles. It needs 6 arguments for 3D (4 for 2D) which give 'x_left', 'x_right', 'y_left', 'y_right', etc... (where 'x_left' is the left box coordinates along x-direction). For example '--box 0 1 0.5 1.5 0 10'.")
            ("input,i", po::value< std::vector<int> >()->multitoken(), "give the type of the input file (101=gadget multiple file, 102=gadget single file, 105=gadget HDF5 file, see documenation for more options). If present, a 2nd argument gives the data to be read from file (1=positions, 2=weights, 4=masses, ..., 2^n=the n+1 data) - e.g. to read positions, masses and velocities insert 1+2+4=7. If present, a 3rd argument gives the particle species to be read from file (1=1st species, 2=2nd species, ..., 2^n=the n+1 species) - e.g. to read the data of species 2,3 and 4 insert 2+4+8=14.")
            ("output,o", po::value< int >(&(this->outputFileType)), "give the type of the output file (101=binary file, 111=text file, see documenation for more options).")
            ("periodic,p", "particle data is in a periodic box with box coordinates given by option '--box' or read from input file.")
            ;
    
    
    po::options_description fieldOptions("Field choices");
    fieldOptions.add_options()
            ("field,f", po::value<std::vector<std::string> >()->multitoken(), "specify which field is to be interpolated to grid - e.g. '-f density'. You can specify multiple fields at a time. You can choose to output the field value at the sampling point or the volume averaged field value inside the grid cell associated to the sampling point (field names end with '_a' with 'a' standing for averaged). Available options are:\n"
                    "  density = \tcompute the density at the sampling point position.\n"
                    "  density_a = \tcompute the volume averaged density inside the sampling cell associated to the sampling point [DEFAULT choice if none provided]. Each of the options below have also a '*_a' version which will be left out to minimize the help messages.\n"
#ifdef VELOCITY
                    "  velocity = \tcompute the velocity at the sampling point position (use 'velocity_a' to get the averaged velocity inside the sampling cell).\n"
                    "  gradient = \tcompute the velocity gradient at the sampling point position (use 'gradient_a' to get the averaged velocity gradient inside the sampling cell).\n"
                    "  divergence = \tcompute velocity divergence at the sampling point position (use 'divergence_a' to get the averaged velocity divergence inside the sampling cell).\n"
                    "  shear = \tcompute velocity shear at the sampling point position (use 'shear_a' to get the averaged velocity shear inside the sampling cell).\n"
                    "  vorticity = \tcompute velocity vorticity at the sampling point position (use 'vorticity_a' to get the averaged velocity vortivity inside the sampling cell).\n"
                    "  velocityStd_a = \tcompute velocity standard deviation inside the sampling cell (NOTE: there is no 'velocityStd' of this option and this option works only with averaging method 2 '--method 2').\n"
#endif
#ifdef SCALAR
                    "  scalar = \tcompute scalar quantities at the sampling point position (use 'scalar_a' to get the averaged field components inside the sampling cell).\n"
                    "  scalarGradient = \tcompute the gradient of the scalar quantities at the sampling point position (use 'scalarGradient_a' to get the averaged field gradient inside the sampling cell).\n"
#endif
            );
    
    
    po::options_description regionOptions("Region options");
    regionOptions.add_options()
            ("region", po::value< std::vector<Real> >()->multitoken(), "choose this option to compute the field interpolation only in a given part of the full box. Use this option to specify the region of interest in terms of fractions of box length (i.e. '--region 0.4 0.6 0.3 0.7 0.45 0.55' computes the density for the box that extends from 0.4 to 0.6 of the box length along direction x, 0.3 to 0.7 along direction y and 0.45 to 0.55 along direction z).")
            ("regionMpc", po::value< std::vector<Real> >()->multitoken(), "choose this option to compute the field interpolation only in a given part of the box. Specify the limits of that region in Mpc units; see previous option for additional information.")
            ;
    
    
    po::options_description partitionOptions("Partition options");
    partitionOptions.add_options()
            ("partition", po::value< std::vector<size_t> >(&(this->partition))->multitoken(), "choose this option if the particle data is too large to compute the Delaunay triangulation for the full data at once. Specify here in how many parts to split the box along each direction (e.g. '--partition 3 3 3' splits the data in 27 chuncks).")
            ("partNo", po::value<int>(&(this->partNo)), "choose to compute the density only for this partition number (from 0 to 'maximum partitions'-1). This options is usefull if you would like to compute the density for a large particle number on several different machines at the same time - need to run the program on each machine independently.")
            ;
    
    
    po::options_description paddingOptions("Padding options");
    paddingOptions.add_options()
            ("padding", po::value< std::vector<Real> >()->multitoken(), "give the size of the padding need to make sure that the Delaunay triangulation fully covers the region of interest. There are two ways to give the padding size:\n"
                    "  1) \t by giving one value which is the average number of particles that will be copied along each face of the region of interest. The actual computation uses all the particle that are within 'padding number' * 'particle grid spacing' distance from the region of interest. For example '--padding 5' will add an average of 5 particles on both the left and right sides for each dimension.\n"
                    "  2) \t by giving the size of the padding for each face of the box of interest. This size is given with respect to the box length along each coordinate (i.e. '-padding 0.1 0.2 0.5 0.5 0.1 0.1' means that box will be padded with '0.1*x box length' on the left of the x-coordiante and by '0.2*x box length' on the right of the x-coordinate, similar for the y and z dimensions).")
            ("paddingMpc", po::value< std::vector<Real> >()->multitoken(), "give the size of the padding need to make sure that the Delaunay triangulation fully covers the region of interest. Similar to option 'padding' choice '2)' with the difference that the padding size is given in Mpc and not box lengths.")
#ifdef TEST_PADDING
            ("noTest", "do not test for the efficiency of the padding when computing the Delaunay tesselation. The default is to use dummy particles positioned at the boundary of the extended padded box to test if the Delaunay tesselation fully covers the region of interest (i.e. the unpadded box).")
#endif
            ;
    
    
    po::options_description averagingOptions("Averaging options");
    averagingOptions.add_options()
            ("method,m", po::value<int>(&(this->method))->default_value(this->method), "choose volume averaging method (only for fields inserted using the '_a' ending):\n"
                    "  1 = \tvolume average the fields using a Monte Carlo method with pseudo-random numbers inside the Delaunay cell.\n"
                    "  2 = \tvolume average the fields using the Monte Carlo method inside the grid cell.\n"
                    "  3 = \tvolume average the fields using uniformly distributed volume sampling points in each grid cell (recommended only for testing purposes).")
            ("samples,s", po::value<int>(&(this->noPoints)), "specify the number of sampling points when volume averaging the fields inside each grid cell (e.g. '-s 20'). DEFAULT values if none specified:\n"
                    "  1st method: \tan average of 100 sample points per grid cell.\n"
                    "  2nd method: \t20 random points per grid cell.\n"
                    "  3rd method: \t27 random points per grid cell.")
            ("density0", po::value<Real>(&averageDensity), "supply a value to be used to scale the density. If none is supplied, the average density will be used for this task.")
            ("seed", po::value<size_t>(&(this->randomSeed)), "integer value to be used for the random seed generator when interpolating to the grid using Monte Carlo methods. Generated randomly if not supplied by the user.")
            ;
    
    
    po::options_description redshiftConeOptions("Redshift cone options");
    redshiftConeOptions.add_options()
            ("redshiftCone", po::value< std::vector<Real> >()->multitoken(), "specify to interpolate the fields on a redshift cone grid (i.e. on spherical coordinates). Must give 6 arguments in 3D (4 in 2D) which give 'r_min', 'r_max', 'theta_min', 'theta_max', 'psi_min' and 'psi_max' (where 'r' is the distance in Mpc and 'theta' and 'psi' the two angles expressed in degrees). For example '--redshiftCone 1 10 0 90 0 360' gives you a half a sphere shell." )
            ("origin", po::value< std::vector<Real> >( &(this->originPosition) )->multitoken(), "specify the origin of the spherical coordinate system used in option '--redshiftCone'. It gives the x, y and z values of the origin point.")
            ;
    
    
    po::options_description additionalOptions("Additional options");
    additionalOptions.add_options()
            ("config,c", po::value<std::string>(&(this->configFilename)), "supply all/part of the program options in a configuration file. The option syntax is the same as at the command line. Can insert comments using the '#' symbol (everything after this symbol until the end of the line will be considered a comment)." )
            ("NGP", "choose the NGP (nearest grid point) as the grid interpolation method instead of DTFE. This method is available only for: density and velocity fields.")
            ("CIC", "choose the CIC (Cloud In Cell) as the grid interpolation method instead of DTFE. This method is available only for: density and velocity fields.")
            ("TSC", "choose the TSC (Triangular Shape Cloud) as the grid interpolation method instead of DTFE. This method is available only for: density and velocity fields.")
            ("SPH", po::value<int>(&(this->SPH_neighbors)), "choose the SPH (Smoothed Particle Hydrodynamics) as the grid interpolation method instead of DTFE. This method is available only for: density, velocity and scalar fields.")
            ("MpcUnit", po::value<Real>(&(this->MpcValue)), "specify the value of 1Mpc in units of the input particle position data. [DEFAULT value is the one given in the 'DMPC_UNIT' variable in the Makefile.]")
            ("extensive", "specify that all the fields under 'scalar fields' are extensive quantities. If this option is missing than the code treats the variables as intensive fields. This option is important only when using the TSC or SPH interpolation methods applied to the scalar variable.")
            ("verbose,v", po::value<int>(&(this->verboseLevel))->default_value(verboseLevel), "choose the verbosity level of the program (a value from 0 to 3). See the documentation for additional help.")
            ("randomSample", po::value<Real>(&(this->randomSample)), "generates a random subsample of the input data. The size of the subsample is given by value supplied to the option (with values from 0. to 1.). Only this random subsample of the full data set will be used in any further computations. For example '--randomSample 0.1' will keep only 10\% of the data set for further computations.")
            ("poisson", po::value<size_t>(&(this->poisson)), "generate the particle positions randomly. The argument gives the root 3 in 3D (and root 2 in 2D) of the random number of particles. The particles have the same weight and are in a box of size unity. For example '-poisson 256' will generate 256^3 particles in 3D, while only 256^2 particles in 2D.")
#ifdef REDSHIFT_SPACE
            ("redshiftSpace", po::value< std::vector<Real> >()->multitoken(), "specify this option to transform the particle positions from position-space to redhsift-space. This option takes 3 arguments that specifies the direction (d1,d2,d3) along which to tranform to redshift-space. For example '--redshiftSpace d1 d2 d3' specifies to transform to 'redshift-space = position-space + (d1,d2,d3)*velocity / H', with H=100 h km/s /Mpc and (d1,d2,d3) normalized to a unit vector." )
#endif
            ("options", po::value< std::vector<std::string> >( &(this->additionalOptions) )->multitoken(), "variable used to supply additional options to the program in a very simple way. Each additional option will be stored as a string in 'User_options.additionalOptions' (this variable is a vector of strings).")
            ;
    
    
    po::options_description hidden("Hidden options");
    hidden.add_options()
            ("inputFile", po::value<std::string>(&(this->inputFilename)), "name of the input position file")
            ("outputFile",po::value<std::string>(&(this->outputFilename)), "root name of the output file/files")
            ;
    
    
    // add the options to 'visibleOptions' and to 'allOptions'
    visibleOptions.add(mainOptions);
#ifdef FIELD_OPTIONS
    visibleOptions.add(fieldOptions);
#endif
#ifdef REGION_OPTIONS
    visibleOptions.add(regionOptions);
#endif
#ifdef PARTITION_OPTIONS
    visibleOptions.add(partitionOptions);
#endif
#ifdef PADDING_OPTIONS
    visibleOptions.add(paddingOptions);
#endif
#ifdef AVERAGING_OPTIONS
    visibleOptions.add(averagingOptions);
#endif
#ifdef REDSHIFT_CONE_OPTIONS
    visibleOptions.add(redshiftConeOptions);
#endif
#ifdef ADDITIONAL_OPTIONS
    visibleOptions.add(additionalOptions);
#endif
    
    
    allOptions.add(mainOptions);
    allOptions.add(fieldOptions);
    allOptions.add(regionOptions);
    allOptions.add(partitionOptions);
    allOptions.add(paddingOptions);
    allOptions.add(averagingOptions);
    allOptions.add(redshiftConeOptions);
    allOptions.add(additionalOptions);
    allOptions.add(hidden);
    
    // now add the hidden options to 'positional_options_description'
    p.add( "inputFile", 1 );
    p.add( "outputFile", 1 );
}


//print help information to the user
void User_options::helpInformation( po::options_description &visibleOptions, char *progName )
{
    MESSAGE::Message message( 3 );
    message << "Use this program to interpolate fields on a grid using the DTFE method - there are multiple options which allow for increased flexibility of the program. The computations are done using the Delaunay Triangulation module from the CGAL library.\n";
    message << "Usage:    " << progName << "  name_position_file  output(root)_file  'options - see below' \n";
    message << "On top of the above, the user can add any of the following options:\n";
    message << visibleOptions << "\n" << MESSAGE::Flush;
    exit( EXIT_SUCCESS );
}
void User_options::shortHelp( char *progName )
{
    Real temp;
    po::options_description mainOptions("Main options");
    mainOptions.add_options()
            ("help,h", "produce this help message.")
            ("full_help", "produce detailed help message.")
            ("grid,g", po::value< Real >(&temp), "specify grid size along each direction.")
            ("box", po::value< Real >(&temp), "specify the coordinates of the box encompasing all the particles.")
            ("input,i", po::value< std::vector<int> >()->multitoken(), "give the type of the input file, which data to read and for which particle species. See full help for details.")
            ("output,o", po::value< std::vector<int> >(), "give the type of the output file. See full help for details.")
            ("periodic,p", "specify the data is in a periodic box.")
            ;
    po::options_description fieldOptions("Field choices");
    fieldOptions.add_options()
            ("field,f", po::value< Real >(&temp), "specify which field to interpolate to grid. Available options are:\n"
                    "  density = \tdensity at the sampling point position.\n"
                    "  density_a = \tvolume averaged density inside the sampling cell[DEFAULT]. Each of the options below also have a '*_a' version which is left out to minimize the help messages.\n"
#ifdef VELOCITY
                    "  velocity = \tnon-averaged velocity (use 'velocity_a' to get the volume averaged value).\n"
                    "  gradient = \tnon-averaged velocity gradient (use 'gradient_a' to get the volume averaged value).\n"
                    "  divergence = \tnon-averaged velocity divergence (use 'divergence_a' to get the volume averaged value).\n"
                    "  shear = \tnon-averaged velocity shear (use 'shear_a' to get the volume averaged value).\n"
                    "  vorticity = \tnon-averaged velocity vorticity (use 'vorticity_a' to get the volume averaged value).\n"
#endif
#ifdef SCALAR
                    "  scalar = \tnon-averaged scalar quantities (use 'scalar_a' to get the volume averaged value).\n"
                    "  scalarGradient = \tnon-averaged gradient of the scalar quantities (use 'scalarGradient_a' to get the volume averaged value).\n"
#endif
            );
    po::options_description regionOptions("Region options");
    regionOptions.add_options()
            ("region", po::value< Real >(&temp), "specify to interpolate the fields only in a region of the box given in terms of fractions of box length.")
            ("regionMpc", po::value< Real >(&temp), "specify to interpolate the fields only in a region given in Mpc coordinates.")
            ;
    po::options_description partitionOptions("Partition options");
    partitionOptions.add_options()
            ("partition", po::value< Real >(&temp), "specify in how many parts to split the box along each direction.")
            ("partNo", po::value< Real >(&temp), "choose to compute the interpolation only for this partition number.")
            ;
    po::options_description paddingOptions("Padding options");
    paddingOptions.add_options()
            ("padding", po::value< Real >(&temp), "give the size of the buffer zone to make sure that the triangulation fully covers the region of interest.")
            ("paddingMpc", po::value< Real >(&temp), "give the size of the buffer zone in Mpc units.")
#ifdef TEST_PADDING
            ("noTest", "do NOT test for the efficiency of the padding.")
#endif
            ;
    po::options_description averagingOptions("Averaging options");
    averagingOptions.add_options()
            ("method,m", po::value< Real >(&temp), "choose the MC averaging method: 1 = volume average inside the Delaunay cell OR 2 = volume average inside the grid cell.")
            ("samples,s", po::value< Real >(&temp), "specify the number of MC sampling points for volume averaged interpolation.")
            ("density0", po::value< Real >(&temp), "value to scale the density [DEFAULT: use average density].")
            ("seed", po::value< Real >(&temp), "integer value used as seed for the random generator.")
            ;
    po::options_description redshiftConeOptions("Redshift cone options");
    redshiftConeOptions.add_options()
            ("redshiftCone", po::value< Real >(&temp), "specify to interpolate the fields on a redshift cone grid and give the coordinates  ('r_min', 'r_max', 'theta_min', 'theta_max', 'psi_min' and 'psi_max')." )
            ("origin", po::value< Real >(&temp), "specify the origin of the spherical coordinate system.")
            ;
    po::options_description additionalOptions("Additional options");
    additionalOptions.add_options()
            ("config,c", po::value< Real >(&temp), "name the configuration file from which to read the program options." )
            ("NGP", "choose NGP (Nearest Grid Point) for grid interpolation.")
            ("CIC", "choose CIC (Cloud In Cell) for grid interpolation.")
            ("TSC", "choose TSC (Triangular Shape Cloud) for grid interpolation.")
            ("SPH", po::value< Real >(&temp), "choose SPH for grid interpolation (argument = number nearest neighbors).")
            ("MpcUnit", po::value< Real >(&temp), "value of 1Mpc in units of the input position data.")
            ("extensive", "the fields under 'scalar fields' are extensive quantities [DEFAULT: intensive variables].")
            ("verbose,v", po::value< Real >(&temp), "choose the verbosity level (from 0 to 3).")
            ("randomSample", po::value< Real >(&temp), "generates a random subsample of the input data (argument = from 0 to 1, gives fraction of particles).")
            ("poisson", po::value< Real >(&temp), "generate the particle positions randomly.")
#ifdef REDSHIFT_SPACE
            ("redshiftSpace", po::value< std::vector<Real> >()->multitoken(), "transform the particle positions from position-space to redhsift-space -- need to give 3 values that give the direction for the shift." )
#endif
            ("options", po::value< Real >(&temp), "variable used to supply additional options.")
            ;
    
    po::options_description visibleOptions;
    visibleOptions.add(mainOptions);
#ifdef FIELD_OPTIONS
    visibleOptions.add(fieldOptions);
#endif
#ifdef REGION_OPTIONS
    visibleOptions.add(regionOptions);
#endif
#ifdef PARTITION_OPTIONS
    visibleOptions.add(partitionOptions);
#endif
#ifdef PADDING_OPTIONS
    visibleOptions.add(paddingOptions);
#endif
#ifdef AVERAGING_OPTIONS
    visibleOptions.add(averagingOptions);
#endif
#ifdef REDSHIFT_CONE_OPTIONS
    visibleOptions.add(redshiftConeOptions);
#endif
#ifdef ADDITIONAL_OPTIONS
    visibleOptions.add(additionalOptions);
#endif
    
    
    MESSAGE::Message message( 3 );
    message << "Use this program to interpolate fields on a grid using the DTFE method - there are multiple options which allow for increased flexibility of the program. The computations are done using the Delaunay Triangulation module from the CGAL library.\n";
    message << "Usage:    " << progName << "  name_position_file  output(root)_file  'options - see below' \n";
    message << "On top of the above, the user can add any of the following options:\n";
    message << visibleOptions << "\n";
    message << "\nUse '--full_help' for a more detailed help information.\n" << MESSAGE::Flush;
    exit( EXIT_SUCCESS );
}




/* Print to the user what options the program will use. */
void User_options::printOptions()
{
    std::string uField;
    if ( this->uField.density ) uField = " density,";
    if ( this->uField.velocity ) uField += " velocity,";
    if ( this->uField.velocity_gradient ) uField += " velocity gradient,";
    if ( this->uField.velocity_divergence ) uField += " velocity divergence,";
    if ( this->uField.velocity_shear ) uField += " velocity shear,";
    if ( this->uField.velocity_vorticity ) uField += " velocity vorticity,";
    if ( this->uField.scalar ) uField += " scalar,";
    if ( this->uField.scalar_gradient ) uField += " scalar gradient,";
    if ( this->uField.triangulation ) uField += " triangulation,";
    if ( not this->uField.selected() and this->NGP ) uField = "none since NGP interpolation";
    else if ( not this->uField.selected() and this->CIC ) uField = "none since CIC interpolation";
    else if ( not this->uField.selected() and this->TSC ) uField = "none since TSC interpolation";
    else if ( not this->uField.selected() and this->SPH ) uField = "none since SPH interpolation";
    else if ( not this->uField.selected() ) uField = "none";
    
    std::string aField;
    if ( this->aField.density ) aField = " density,";
    if ( this->aField.velocity ) aField += " velocity,";
    if ( this->aField.velocity_gradient ) aField += " velocity gradient,";
    if ( this->aField.velocity_divergence ) aField += " velocity divergence,";
    if ( this->aField.velocity_shear ) aField += " velocity shear,";
    if ( this->aField.velocity_vorticity ) aField += " velocity vorticity,";
    if ( this->aField.velocity_std ) aField += " velocity standard deviation,";
    if ( this->aField.scalar ) aField += " scalar,";
    if ( this->aField.scalar_gradient ) aField += " scalar gradient,";
    if ( not this->aField.selected() ) aField = "none";
    
    std::string interpolationMethod = "DTFE";
    if ( this->CIC ) interpolationMethod = "CIC";
    else if ( this->NGP ) interpolationMethod = "NGP";
    else if ( this->TSC ) interpolationMethod = "TSC";
    else if ( this->SPH ) 
    {
        char temp[100];
        sprintf( temp, "%d neighbors", this->SPH_neighbors );
        interpolationMethod = "SPH ";
        interpolationMethod += temp;
    }
    
    std::string inFileType = "unknown";
    if ( this->inputFileType==101 ) inFileType = "Gadget multiple files";
    else if ( this->inputFileType==101 ) inFileType = "Gadget single file";
    else if ( this->inputFileType==105 ) inFileType = "Gadget HDF5 file/files";
    else if ( this->inputFileType==111 ) inFileType = "text file with positions (first 3 columns and weights in 4th column)";
    
    std::string outFileType = "unknown";
    if ( this->outputFileType==101 ) outFileType = "binary file";
    else if ( this->outputFileType==100 ) outFileType = "density binary file";
    else if ( this->outputFileType==110 ) outFileType = "text file";
    
    
    MESSAGE::Message message( verboseLevel );
    message << "RUNNING: " << this->programOptions << "\n\n";
    message << "The program will interpolate to grid the chosen field using the DTFE method with the following input parameters:\n"
            << "\t unaveraged field(s)    : " << uField << "\n"
            << "\t averaged field(s)      : " << aField << "\n"
            << "\t interpolation method   : " << interpolationMethod << "\n"
            << "\t input data file        : " << this->inputFilename << "\n"
            << "\t input data file type   : " << this->inputFileType << " - " << inFileType << "\n"
            << "\t input data blocks      : " << MESSAGE::printElements( this->readParticleData, "  " ) << "\n"
            << "\t input particle species : " << MESSAGE::printElements( this->readParticleSpecies, "  " ) << "\n"
            << "\t output file            : " << this->outputFilename << "\n"
            << "\t output file type       : " << this->outputFileType << " - " << outFileType << "\n";
    if ( not this->gridSize.empty() )
        message << "\t grid size              : " << MESSAGE::printElements( this->gridSize, "  " ) << (this->regionOn ? "   for the box region selected by the user\n" : "   for the full particle box\n" );
    else
        message << "\t grid size              : none specifed at the moment";
    if ( not boxCoordinates.isNullBox() )
        message << "\t box coordinates        : [" << MESSAGE::printElements( boxCoordinates, ", " ) << "]\n";
    if ( this->periodic )
        message << "\t computing the grid interpolation in a PERIODIC box\n";
    
    
    if ( this->regionOn )
    {
        message << "\t computing the grid interpolation only for the user specified region of coordinates:\n"
            << "\t\t x extension : " << region[0] << "   " << region[1] << (regionMpcOn?"  Mpc":"  box length") << "\n"
            << "\t\t y extension : " << region[2] << "   " << region[3] << (regionMpcOn?"  Mpc":"  box length") << "\n";
#if NO_DIM==3
        message << "\t\t z extension : " << region[4] << "   " << region[5] << (regionMpcOn?"  Mpc":"  box length") << "\n";
#endif
    }
    
    
    if ( this->partitionOn )
    {
        message << "\t splitting the data in  : " << MESSAGE::printElements( partition, "  ") << "   separate data sets on which to apply the Delaunay triangulation.\n";
        if ( this->partNo>=0 )
            message << "\t computing grid interpolation ONLY for data set " << partNo << " of the partitioned data\n";
    }
    
    
    if ( not this->paddingLength.isNullBox() )
    {
        message << "\t padding length:\n"
            << "\t\t x axis : " << paddingLength[0] << "   " << paddingLength[1] << (paddingMpcOn?"  Mpc":"  box length") << "\n"
            << "\t\t y axis : " << paddingLength[2] << "   " << paddingLength[3] << (paddingMpcOn?"  Mpc":"  box length") << "\n";
#if NO_DIM==3
        message << "\t\t z axis : " << paddingLength[4] << "   " << paddingLength[5] << (paddingMpcOn?"  Mpc":"  box length") << "\n";
#endif
    }
    else
        message << "\t number padding particles: " << this->paddingParticles << "\n";
    if ( this->DTFE and this->testPaddedBoundaries )
        message << "\t the DTFE computation will add DUMMY TEST particles to test the efficiency of the padding\n";
    
    
    if ( this->aField.selected() and this->DTFE )
    {
        std::string averagingMethod;
        if ( this->method==1 ) averagingMethod = "Monte Carlo sampling using quasi-random points inside the Delaunay cells";
        else if ( this->method==2 ) averagingMethod = "Monte Carlo sampling using random points inside the sampling cells";
        else if ( this->method==3 ) averagingMethod = "equidistant sampling points inside the sampling cells";
        else averagingMethod = "unknown";
        message << "\t volume averaging method: " << this->method << " - " << averagingMethod << "\n"
                << "\t number sampling points : " << this->noPoints << "\n";
        if ( this->method==2 )
            message << "\t random generator seed  : " << this->randomSeed << "\n";
        if ( averageDensity>Real(0.) )
            message << "\t density scaling value  : " << averageDensity << "\n";
    }
    
    
    if ( redshiftConeOn )
        message << "\t Computing the grid interpolation to a redshift cone grid of coordinates " << (NO_DIM==2? "[r_min, r_max, psi_min, psi_max]" : "[r_min, r_max, theta_min, theta_max, psi_min, psi_max]") << " = [" << MESSAGE::printElements( redshiftCone, ", " ) << "]\n"
            << "\t The origin of the light cone is (x,y,z) : (" << MESSAGE::printElements( originPosition, ", " ) << ")\n";
    
    message << "\t 1Mpc = " << MpcValue << " in units of input data.\n";
    if ( poisson>0 )
        message << "\t Particle positions will be generated randomly. Particle number = " << (NO_DIM==2? poisson*poisson : poisson*poisson*poisson) << ".\n";

#ifdef REDSHIFT_SPACE
    if ( transformToRedshiftSpaceOn )
        message << "\t Transforming from position-space to redshift space using the velocity along the direction :  ( " << MESSAGE::printElements( transformToRedshiftSpace, ", " ) << " )\n";
#endif
    
    
    if ( this->gridSize.empty() )
        message << "\n\n~~~WARNING~~~ No grid size specified. Unless the input data file specifies the grid size for the output result, the program will end with and error message!\n\n";
    
    message << "\n" << MESSAGE::Flush;
}




/* Read the user supplied options and check that they satify some restrictions. */
void User_options::readOptions(int argc, char *argv[], bool getFileNames, bool showOptions)
{
    po::options_description visibleOptions("Allowed options"), allOptions("All options");
    po::positional_options_description p;
    
    this->addOptions( allOptions, visibleOptions, p );  //read the options available to the program
    
    po::variables_map vm;
    MESSAGE::Message showMessage(3);
    try
    {
        po::store( po::command_line_parser(argc, argv).options(allOptions).positional(p).run(), vm );
        po::notify(vm);
    }
    catch (exception& e)
    {
        throwError( "When reading the command line program options:\n\t\"", e.what(), "\"\n" );
    }
    catch (...)
    {
        throwError( "Unknown error when reading the command line program options!" );
    }
    
    
    
    if ( vm.count("config") )
    {
        // read the options from the file
        ifstream ifs( configFilename.c_str() );
        if (!ifs)
            showMessage << "~~~ERROR~~~ Can not open the configuration file '" << configFilename << "'.\n";
        else
        {
            try
            {
                store( parse_config_file(ifs, allOptions), vm );
                notify(vm);
            }
            catch (exception& e)
            {
                throwError( "When reading the configuration file program options:\n\t\"", e.what(), "\"\n" );
            }
            catch (...)
            {
                throwError( "Unknown error when reading the configuration file program options!" );
            }
        }
    }
    
    
    if ( not (vm.count("help") or vm.count("full_help")) and not vm.count("inputFile") and getFileNames )
        showMessage << "~~~ERROR~~~ No input file detected.\n";
    if ( not (vm.count("help") or vm.count("full_help")) and not vm.count("outputFile") and getFileNames )
        showMessage << "~~~ERROR~~~ No output file detected.\n";
    if ( not (vm.count("help") or vm.count("full_help")) and (not vm.count("inputFile") or not vm.count("outputFile")) and getFileNames )
        exit( EXIT_SUCCESS );
    
    if ( vm.count("help") )
        this->shortHelp( argv[0] ); // print available options and exit
    else if ( vm.count("full_help") )
        this->helpInformation( visibleOptions, argv[0] ); // print available options and exit
    
    
    conflicting_options(vm, "SPH", "TSC");
    conflicting_options(vm, "SPH", "CIC");
    conflicting_options(vm, "CIC", "TSC");
    conflicting_options(vm, "NGP", "CIC");
    conflicting_options(vm, "NGP", "TSC");
    conflicting_options(vm, "NGP", "SPH");
    conflicting_options(vm, "region", "regionMpc");
    conflicting_options(vm, "padding", "paddingMpc");
    
    option_dependency(vm, "partNo", "partition");
    option_dependency(vm, "redshiftCone", "origin");
    
    
    
    // Read the options supplied to the program
    if ( vm.count("grid") )
    {
        int const temp = this->gridSize.size();
        if ( temp==1 ) this->gridSize.assign( NO_DIM, this->gridSize.at(0) );
        else if ( temp!=NO_DIM ) throwError( "You can only insert 1 or ", NO_DIM, " values for the '-g [ --grid ]' option (e.g. '-g 256' or '-g 128 256 256')." );
        
        for (size_t i=0; i<this->gridSize.size(); ++i)
            lowerBoundCheck( this->gridSize[i], size_t(1), "the values of option '-g [ --grid ]'" );
    }
    if ( vm.count("box") )      // read the box coordinates
    {
        int const temp = vm["box"].as< std::vector<Real> >().size();
        if ( temp!=2*NO_DIM ) throwError( "You need to specify ", 2*NO_DIM, " arguments with the option '--box'." );
        
        boxCoordinates.coords = vm["box"].as< std::vector<Real> >();
        for (int i=0; i<NO_DIM; ++i)
            if( boxCoordinates[2*i+1]<boxCoordinates[2*i]) throwError( "The right coordinate of the box encompasing the data along axis " , i+1, " (1=x, 2=y, 3=z) must be larger than the left coordinate of the box. This is not the case in the arguments supplied to '--box' option." );
        userGivenBoxCoordinates = true;
    }
    if ( vm.count("input") )    // read the input file type, the particle data to read and the particle species to read
    {
        std::vector<int> temp = vm["input"].as< std::vector<int> >();
        inputFileType = temp[0];    // the file type
        if ( temp.size()>=2 )   // set what data blocks to read
            for (size_t i=0; i<readParticleData.size(); ++i)
            {
                readParticleData[i] = temp[1] % 2;
                temp[1] /= 2;
            }
        if ( temp.size()>=3 )   // set what particle species to read
            for (size_t i=0; i<readParticleSpecies.size(); ++i)
            {
                readParticleSpecies[i] = temp[2] % 2;
                temp[2] /= 2;
            }
    }
    if ( vm.count("periodic") )
        this->periodic = true;
    
    
    // read which field to compute using the DTFE method
    if ( vm.count("field") )
    {
        for (size_t i=0; i<vm["field"].as< std::vector<std::string> >().size(); ++i)
        {
            std::string field = vm["field"].as< std::vector<std::string> >().at(i);
            bool uFieldOption = this->uField.updateChoices( field, "triangulation", "density", "velocity", "gradient", "divergence", "shear", "vorticity", "", "scalar", "scalarGradient" );
            bool aFieldOption = this->aField.updateChoices( field, "", "density_a", "velocity_a", "gradient_a", "divergence_a", "shear_a", "vorticity_a",  "velocityStd_a", "scalar_a", "scalarGradient_a" );
            if ( not(uFieldOption or aFieldOption) ) throwError( "Unknown value '" + field + "' for the option '--field'." );
        }
    }
    else
        this->aField.density = true;
    
    
    // read the user defined region options
    if ( vm.count("regionMpc") )
        this->regionMpcOn = true;
    if ( vm.count("regionMpc") or vm.count("region") )
    {
        this->regionOn = true;
        size_t temp = 0;
        if (this->regionMpcOn)
        {
            temp = vm["regionMpc"].as< std::vector<Real> >().size();
            this->region.coords = vm["regionMpc"].as< std::vector<Real> >();
        }
        else
        {
            temp = vm["region"].as< std::vector<Real> >().size();
            this->region.coords = vm["region"].as< std::vector<Real> >();
        }
        if ( temp!=2*NO_DIM ) throwError( "You have to insert ", 2*NO_DIM, " values for the '--region' or '--regionMpc' option (e.g. '--region 0.4 0.6 0.3 0.7 0.45 0.55')." );
        for (int i=0; i<NO_DIM; ++i)
            if( region[2*i+1]<region[2*i]) throwError( "The right coordinate of the region of interest along axis ", i+1, " (1=x, 2=y, 3=z) must be larger than the left coordinate. This is not the case in the arguments supplied to '--region' option." );
    }
    
    
    // read the partition options
    if ( vm.count("partition") )
    {
        this->partitionOn = true;
        int const temp = this->partition.size();
        if ( temp==1 ) this->partition.assign( NO_DIM, this->partition.at(0) );
        else if ( temp!=3 ) throwError( "You can only insert 1 or ", NO_DIM, " values for the '--partition' option (e.g. '--partition 3' or '--partition 2 3 3')." );
        
        for (size_t i=0; i<this->partition.size(); ++i)
            lowerBoundCheck( this->partition[i], size_t(1), "the values for the option '--partition'" );
        
        size_t temp2 = 1;
        for (size_t i=0; i<this->partition.size(); ++i)
            temp2 *= this->partition[i];
        if ( vm.count("partNo") ) intervalCheck( this->partNo, 0, int(temp2-1), "'--partNo' program option" );
        
        if ( temp2==1 ) //if no actual partition is done
            this->partitionOn = false;
    }
    
    
    // read the padding options
    if ( vm.count("paddingMpc") )
        this->paddingMpcOn = true;
    if ( vm.count("paddingMpc") or vm.count("padding") )
    {
        this->paddingOn = true;
        size_t temp = 0;
        if (this->paddingMpcOn)
        {
            temp = vm["paddingMpc"].as< std::vector<Real> >().size();
            this->paddingLength.coords = vm["paddingMpc"].as< std::vector<Real> >();
        }
        else
        {
            temp = vm["padding"].as< std::vector<Real> >().size();
            if ( temp==1 )
            {
                paddingParticles = vm["padding"].as< std::vector<Real> >().at(0);
                lowerBoundCheck( paddingParticles, Real(0.), "'--padding' program option" );
                temp = 2*NO_DIM;
            }
            else
                paddingLength.coords = vm["padding"].as< std::vector<Real> >();
        }
        if ( temp!=2*NO_DIM )
            throwError( "You have to insert 1 or ", 2*NO_DIM, " values for the '--padding' or '--paddingMpc' option (e.g. '--padding 0.1 0.2 0.1 0.1 0.3 0.5')." );
    }
    if ( vm.count("noTest") ) this->testPaddedBoundaries = false;
    
    
    // read the family averaging options
    intervalCheck( this->method, 1, 3, "'--method' can have only 3 values (from 1 to 3) since there are implemented only 3 methods for field averaging inside the sampling cell (see '--help' for additional details)" );
    if ( vm.count("samples") )
    {
        this->noPointsOn = true;
        lowerBoundCheck( this->noPoints, 1, "value of the '-s' ['--samples'] option in the program command line options" );
    }
    else if ( this->method==2 )
        this->noPoints = 20;  //default values for method 2
    else if ( this->method==3 )
        this->noPoints = 27;  //default values for method 3
    if ( vm.count("density0") )
        lowerBoundCheck( averageDensity, Real(0.), "the value of '--density0'" );
    if ( vm.count("seed") )
        lowerBoundCheck( this->randomSeed, size_t(0), "value suplied with program option '--seed'" );
    else
    {
        std::srand( (unsigned)time(0) );
        this->randomSeed = std::rand();
    }
    
    
    // read the redshift cone options
    if ( vm.count("redshiftCone") )
    {
        redshiftConeOn = true;
        size_t temp = vm["redshiftCone"].as< std::vector<Real> >().size();
        if ( temp!=2*NO_DIM ) throwError( "The '--redshiftCone' option must be followed by ", 2*NO_DIM, " values which give the extension ", (NO_DIM==2 ? "(r_min, r_max, psi_min, psi_max)" : "(r_min, r_max, theta_min, theta_max, psi_min, psi_max)" ), " for the spherical coordinates region used to interpolate to grid."  );
        redshiftCone.coords = vm["redshiftCone"].as< std::vector<Real> >();
        
        for (int i=0; i<NO_DIM; ++i)
            if ( redshiftCone[2*i]>=redshiftCone[2*i+1] ) throwError( "When inserting the lower and upper values for option '--redshiftCone'. A lower value is higher than an upper value. If error was due to an angular interval value, increase the upper bound by 360 degrees." );
        if (NO_DIM==3)
        {
            intervalCheck( redshiftCone[2], Real(0.), Real(180.), "3rd value of option '--redshiftCone'" );
            intervalCheck( redshiftCone[3], Real(0.), Real(180.), "4th value of option '--redshiftCone'" );
        }
        Real tempRes = redshiftCone[2*(NO_DIM-1)+1] - redshiftCone[2*(NO_DIM-1)]; //(psi_max - psi_min)
        if ( tempRes>=Real(360.) ) throwError( "The 'psi' angle interval can strech at most 360 degrees." );
    }
    if ( vm.count("origin") )
        if ( originPosition.size()!=NO_DIM ) throwError( "The option '--origin' must be followed by ", NO_DIM, " values." );
    
    
    // read the additional options
    if ( vm.count("NGP") )
    {
        this->NGP = true;
        this->DTFE = false;
    }
    if ( vm.count("CIC") )
    {
        this->CIC = true;
        this->DTFE = false;
    }
    if ( vm.count("TSC") )
    {
        this->TSC = true;
        this->DTFE = false;
    }
    if ( vm.count("SPH") )
    {
        this->SPH = true;
        this->DTFE = false;
    }
    if ( vm.count("MpcUnit") )
        lowerBoundCheck( MpcValue, Real(0.), "value of option '--MpcUnit'" );
    if ( vm.count("extensive") )
        this->extensive = true;
    if ( vm.count("verbose") )
        intervalCheck( verboseLevel, 0, 3, "value of option '--verbose'" );
    if ( vm.count("randomSample") )
        intervalCheck( randomSample, Real(0.), Real(1.), "value of option '--randomSample'" );
    if ( vm.count("poisson") )
        lowerBoundCheck( poisson, size_t(1), "value of option '--poisson'" );
#ifdef REDSHIFT_SPACE
    if ( vm.count("redshiftSpace") )
    {
        transformToRedshiftSpaceOn = true;
        transformToRedshiftSpace = vm["redshiftSpace"].as< std::vector<Real> >();
        size_t temp = vm["redshiftSpace"].as< std::vector<Real> >().size();
        if ( temp!=NO_DIM ) throwError( "The '--redshiftSpace' option must be followed by ", NO_DIM, " values which give the direction ", (NO_DIM==2 ? "(d1,d2)" : "(d1,d2,d3)" ), " of the vector used to transform from position-space to redshift-space."  );
        
        //normalize the direction to a unitless vector
        Real length = transformToRedshiftSpace[0]*transformToRedshiftSpace[0] + transformToRedshiftSpace[1]*transformToRedshiftSpace[1];
        length += ( NO_DIM==2 ? 0 : transformToRedshiftSpace[2]*transformToRedshiftSpace[2] );
        length = std::sqrt( length );
        for (int i=0; i<NO_DIM; ++i)
            transformToRedshiftSpace[i] /= length;
    }
    else
#endif
    {
        transformToRedshiftSpaceOn = false;
        transformToRedshiftSpace.assign( NO_DIM, Real(0.) );
    }
    
    
    // set values to userOptions->programOptions
    for (int i=0; i<argc; ++i)
        this->programOptions += string( argv[i] ) + " ";
    
    
    // check that if some options are disabled during compilation, that the user did not insert them by mistake
#ifndef TEST_PADDING
    this->testPaddedBoundaries = false;
#endif
#ifndef VELOCITY
    this->uField.deselectVelocity();
    this->aField.deselectVelocity();
#endif
#ifndef SCALAR
    this->uField.deselectScalar();
    this->aField.deselectScalar();
#endif
    
    // some special settings only for the TSC or SPH methods
    if ( this->NGP or this->CIC or this->TSC or this->SPH )
    {
        // for CIC, TSC and SPH only one type of fields are available - the averaged ones
        if ( this->uField.density ) this->aField.density = true;
        if ( this->uField.velocity ) this->aField.velocity = true;
        if ( this->uField.velocity_gradient ) this->aField.velocity_gradient = true;
        if ( this->uField.velocity_divergence ) this->aField.velocity_divergence = true;
        if ( this->uField.velocity_shear ) this->aField.velocity_shear = true;
        if ( this->uField.velocity_vorticity ) this->aField.velocity_vorticity = true;
        if ( this->uField.scalar ) this->aField.scalar = true;
        if ( this->uField.scalar_gradient ) this->aField.scalar_gradient = true;
        this->uField = Field();
        
        // switch gradient computations off for the TSC and SPH methods
        if ( this->aField.velocity_gradient or this->aField.selectedVelocityDerivatives() )
        {
            MESSAGE::Warning warning( verboseLevel );
            warning << "The NGP, CIC, TSC or SPH grid interpolation methods do not have implemented a method for computing the velocity gradient. The velocity gradient will not be computed!" << MESSAGE::EndWarning;
            this->aField.velocity_gradient = false;
            this->aField.deselectVelocityDerivatives(); // switches off the divergence, shear and vorticity computations
        }
        if ( this->aField.scalar_gradient )
        {
            MESSAGE::Warning warning( verboseLevel );
            warning << "The NGP, CIC, TSC or SPH grid interpolation methods do not have implemented a method for computing the velocity gradient. The velocity gradient will not be computed!" << MESSAGE::EndWarning;
            this->aField.scalar_gradient = false;
        }
        if ( (this->NGP or this->CIC or this->TSC) and this->aField.scalar )
        {
            MESSAGE::Warning warning( verboseLevel );
            warning << "The CIC and TSC grid interpolation method does not have implemented a method for computing the scalar fields on the grid. The scalar field cannot be computed!" << MESSAGE::EndWarning;
            this->aField.scalar = false;
        }
    }
    
    
    // check to see that there is at least one task selected
    if ( not (this->uField.selected() or this->uField.triangulation or this->aField.selected() or this->aField.triangulation) )
        throwError( "No valid field interpolation quantities were selected for the computation. The program cannot continue since it does not compute anything." );
    
    
    if (showOptions)
        this->printOptions();
}


/* Updates the values of the full box. */
void User_options::updateFullBox(Box &newBox)
{
    boxCoordinates = newBox;
    fullBoxOffset.clear();
    fullBoxOffset.reserve(NO_DIM);
    fullBoxLength.clear();
    fullBoxLength.reserve(NO_DIM);
    for (size_t i=0; i<NO_DIM; ++i)
    {
        fullBoxOffset.push_back( boxCoordinates[2*i] );
        fullBoxLength.push_back( boxCoordinates[2*i+1] - boxCoordinates[2*i] );
    }
}




/* Updates some of the members of the class after reading the input data file.
NOTE: This function does additional error checks on the values of the 'User_options' class members. */
void User_options::updateEntries(size_t const noTotalParticles,
                                 bool userSampling)
{
    // check full box size and grid values
    if ( this->boxCoordinates.size()!=2*NO_DIM )
        throwError( "Failed a consistency check. The box encompasing the data should have ", 2*NO_DIM, " coordinates, but it has ", this->boxCoordinates.size(), " coordinates. Check again the values supplied as the coordinates of the full data box." );
    if ( this->boxCoordinates.volume()==0. )
        throwError( "Failed a consistency check. The box encompasing the data has 0. volume. Probably you forget to initialize the coordinates of the full particle data box." );
    
    if ( (this->gridSize.empty() and not userSampling) or this->gridSize.size()!=NO_DIM )
        throwError( "Failed a consistency check. The array storing the interpolation grid should have ", NO_DIM, " values, but it has ", this->gridSize.empty()?0:this->gridSize.size(), " values. Check again the values supplied as the size of the interpolation grid." );
    if ( not gridSize.empty() )
        for (size_t i=0; i<this->gridSize.size(); ++i)
            lowerBoundCheck( this->gridSize[i], size_t(1), "the values of option '--grid'" );
    
    
    // update some member values
    // update 'paddedBox' and also 'fullBoxOffset' and 'fullBoxLength'
    paddedBox = boxCoordinates;
    this->updateFullBox( boxCoordinates );  //computes 'fullBoxOffset' and 'fullBoxLength'
    
    // update the 'paddingLength' values
    updatePadding( noTotalParticles );
    
    
    // Set the 'region' variable
    if ( regionOn and not regionMpcOn )
    {
        for (size_t i=0; i<region.size(); ++i )
            region[i] = fullBoxOffset[i%2] + region[i] * fullBoxLength[i%2];
        regionMpcOn = true;
    }
    else if ( not regionOn )
        region = boxCoordinates;
    
    
    // update the user sampling points
    userDefinedSampling = userSampling;
    if ( userSampling )
    {
        redshiftConeOn = false;
        if ( method!=2 and aField.selected() )
        {
            method = 2;
            if ( not noPointsOn ) noPoints = 20;
            MESSAGE::Warning warning( verboseLevel );
            warning << "When computing the volume average of the fields on a user defined grid only volume averaging method 2 is available. The program will use volume averaging method '2 - Monte Carlo sampling using random points inside the grid cells' using '"<< noPoints <<"' random samples in each user given grid cell." << MESSAGE::EndWarning;
        }
    }
    else if ( redshiftConeOn )
    {
        if ( method!=2 and aField.selected() )
        {
            method = 2;
            if ( not noPointsOn ) noPoints = 20;
            MESSAGE::Warning warning( verboseLevel );
            warning << "When computing the volume average of the fields on a redshift cone grid only volume averaging method 2 is available. The program will use volume averaging method '2 - Monte Carlo sampling using random points inside the grid cells' using '"<< noPoints <<"' random samples in each grid cell." << MESSAGE::EndWarning;
        }
    }
    
    
    // check that the option "--partition" is disabled when dealing with redshift cone or used defined sampling coordinates
    if ( partitionOn and (userSampling or redshiftConeOn) )
    {
        MESSAGE::Warning warning( verboseLevel );
        warning << "The option '--partition' is not available when interpolating the fields to a redshift cone grid or to user defined sampling points. This '--partition' option will be disabled in the rest of the program. If the computation is too large for the available RAm we advise that you mannually split the computation in manageable data portions." << MESSAGE::EndWarning;
        partitionOn = false;
        for (int i=0; i<NO_DIM; ++i)
            partition[i] = 1;
        partNo = -1;
    }
    
    
    // additional error checking
    if ( redshiftConeOn )
    {
        if (NO_DIM==3)
        {
            intervalCheck( redshiftCone[2], Real(0.), Real(180.), "3rd value of option '--redshiftCone'" );
            intervalCheck( redshiftCone[3], Real(0.), Real(180.), "4th value of option '--redshiftCone'" );
        }
        Real tempRes = redshiftCone[2*(NO_DIM-1)+1] - redshiftCone[2*(NO_DIM-1)]; //(psi_max - psi_min)
        if ( tempRes>=Real(360.) ) throwError( "The 'psi' angle interval can strech at most 360 degrees." );
        if ( originPosition.empty() or originPosition.size()!=NO_DIM ) throwError( "The vector 'User_options::originPosition' must have ", NO_DIM, " entries." );
    }
    if ( randomSample>=Real(0.) )
        intervalCheck( randomSample, Real(0.), Real(1.), "value of '--randomSample'" );
    
    
    //check that for the equidistant points, the number of sample points is a perfect square/cube
    if ( method==3 )
        rootN( noPoints, NO_DIM );
    
    // check the values of the program constants defined in "define.h"
    if ( noVelComp!=NO_DIM )
        throwError( "The program constant 'noVelComp' must have the same value as the number of spatial dimensions. Please check the value of the constant in file 'define.h'." );
    if ( noGradComp != NO_DIM*NO_DIM )
        throwError( "The number of components of the velocity gradient must have the value (number of spatial dimensions)^2. Please check the value of the constant in file 'define.h'." );
    if ( noShearComp != (NO_DIM*(NO_DIM+1))/2-1 )
        throwError( "The number of components of the velocity shear must be 2 in 2D and 5 in 3D. Please check the value of the constant in file 'define.h'." );
    if ( noVortComp != ((NO_DIM-1)*NO_DIM)/2 )
        throwError( "The number of components of the velocity vorticity must be 1 in 2D and 3 in 3D. Please check the value of the constant in file 'define.h'." );
    if ( noScalarGradComp != noScalarComp*NO_DIM )
        throwError( "The number of components of the scalar values gradient must be the number of dimensions times the number of scalar components (which is ", noScalarComp*NO_DIM, " for the given parameters). Please check the value of the constant in file 'define.h'." );
    
    
    // switch gradient computations off for the TSC and SPH methods
    if ( NGP or CIC or TSC or SPH )
    {
        if ( aField.velocity_gradient or aField.selectedVelocityDerivatives() )
            throwError( "The NGP, CIC, TSC and SPH grid interpolation methods do not have implemented a method for computing the velocity gradient. The velocity gradient cannot be computed!");
        if ( aField.scalar_gradient )
            throwError( "The NGP, CIC, TSC and SPH grid interpolation methods do not have implemented a method for computing the scalar fields gradients. The scalar gradient cannot be computed!");
        if ( (CIC or TSC) and aField.scalar )
            throwError( "The NGP, CIC and TSC grid interpolation method does not have implemented a method for computing the scalar fields on the grid. The scalar field cannot be computed!");
    }
    
    
    // check that the options are disabled in the absence of the VELOCITY and SCALAR compiler directives
#ifndef VELOCITY
    if ( uField.selectedVelocity() or aField.selectedVelocity() )
        throwError( "Compiler directive 'VELOCITY' not detected. You cannot interpolate the velocity and/or velocity related quantities if the compiler directive 'VELOCITY' is not activated." );
#endif
#ifndef SCALAR
    if ( uField.selectedScalar() or aField.selectedScalar() )
        throwError( "Compiler directive 'SCALAR' not detected. You cannot interpolate the scalar data and/or scalar gradient if the compiler directive 'SCALAR' is not activated." );
#endif
#ifndef TEST_PADDING
    if ( testPaddedBoundaries )
        throwError( "Compiler directive 'TEST_PADDING' not detected. You cannot test the padding efficiency if the compiler directive 'TEST_PADDING' is not activated." );
#endif
}



/* Computes the values of the padding length. */
void User_options::updatePadding(size_t const noTotalParticles)
{
    if ( paddingLength.isNullBox() )
    {
        if ( paddingParticles<=Real(0.) )   // do not add padding
            paddingLength.assign( Real(0.) );
        else        // add padding using the number of particles
        {
            if ( this->DTFE or this->SPH )    //padding for the DTFE and SPH methods
                for (int i=0; i<2*NO_DIM; ++i)
                    paddingLength[i] = paddingParticles / Real( pow(noTotalParticles,1./NO_DIM) ) * fullBoxLength[i%2];
            else if ( this->TSC or this->NGP )
                for (int i=0; i<2*NO_DIM; ++i)
                    paddingLength[i] = fullBoxLength[i%2] / gridSize[i%2];
            else if ( this->CIC )
                for (int i=0; i<2*NO_DIM; ++i)
                    paddingLength[i] = 0.;
        }
        paddingOn = true;
        paddingMpcOn = true;
    }
    else if ( paddingOn and not paddingMpcOn )
    {
        for (int i=0; i<2*NO_DIM; ++i)
            paddingLength[i] *= fullBoxLength[i%2];
        paddingMpcOn = true;
    }
}


