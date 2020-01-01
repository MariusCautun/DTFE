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
Here are the functions used for data input and output.
To easily find how to read the input data go to the function:
    


*/

#ifndef INPUT_OUTPUT_HEADER
#define INPUT_OUTPUT_HEADER

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <boost/filesystem.hpp>
namespace bfs=boost::filesystem;

#include "define.h"
#include "particle_data.h"
#include "user_options.h"
#include "message.h"


// contains the definitions of some classes and fucntions used only for input and output purposes
#include "io/input_output.h"

// different data format readers and writters
#include "io/gadget_reader.cc"  // reader for Gadget snapshots
#include "io/binary_io.cc"      // writer for binary file format
#include "io/text_io.cc"        // reader/writer for text files
#include "io/my_io.cc"          // reader/writer for text files
#include "io/density_file_io.cc"    // own output file format



void openInputBinaryFile(std::fstream & inputFile,
                         std::string & fileName);
void openOutputBinaryFile(std::fstream & outputFile,
                          std::string & fileName);
void openInputTextFile(std::fstream & inputFile,
                       std::string & fileName);
void openOutputTextFile(std::fstream & outputFile,
                        std::string & fileName);





//! Functions for reading the input data

/* Read only positions from a text file and on top of that also reads user defined sampling points from an additional file (USE ONLY WHEN YOU NEED CUSTOM SAMPLING POINTS) - see the "src/io/text_io.cc" for the definition of this function */
//#define DATA_INPUT_FUNCTION readTextFile_userDefinedSampling

/* Read the input data from a binary file - see the "src/io/binary_io.cc" for the definition of this function */
//#define DATA_INPUT_FUNCTION readBinaryFile


typedef void (*FunctionReadInputData)(std::string, Read_data<float> *, User_options *);

/*! This function uses the 'inputFileType' entry in the userOptions strcuture to decide on the function that will be used to read the input data. */
FunctionReadInputData chooseInputDataReadFunction(int const inputFileType)
{
    if ( inputFileType==101 )
        return &readGadgetFile;    // Read the input data from a single/multiple Gadget snapshot file (works only for snapshots type 1 or 2). See the "src/io/gadget_reader.cc" for the definition of this function.
#ifdef HDF5
    else if ( inputFileType==105 )
        return &HDF5_readGadgetFile;        // Read the input data from a HDF5 Gadget snapshot file. See the "src/io/gadget_reader.cc" for the definition of this function.
#endif
    else if ( inputFileType==111 )
        return &readTextFile;               // Read the input data from a text file. See the "src/io/text_io.cc" for the definition of this function.
    else if ( inputFileType==112 )
        return &readTextFile_positions;     // Read the input data from a text file. The text file has only particle positions. See the "src/io/text_io.cc" for the definition of this function.
    else if ( inputFileType==121 )
        return &readBinaryFile;             // Define your custom function to read the input data. Do this in the "src/io/binary_io.cc" file.
    else if ( inputFileType==131 )
        return &readMyFile;                 // Define your custom function to read the input data. Do this in the "src/io/my_io.cc" file.
    else
        throwError( "Unknow value for the 'inputFileType' argument in function 'chooseInputDataReadFunction'. The program could not recognize the input data file type." );
}




/*! This function reads the particle data used for the DTFE computation. You can modify this function as you please. */
void readInputData(std::vector<Particle_data> *p,
                   std::vector<Sample_point> *samplingCoordinates,
                   User_options *userOptions)
{
    std::string filename = userOptions->inputFilename;  // the name of the input file
    Read_data<Real> readData; // will store the data read from the input file
    
    
    // Read the data from the input file - see the function 'chooseInputDataReadFunction' that selects the input function used to read in the input data file
    FunctionReadInputData functionReadInputData = chooseInputDataReadFunction( userOptions->inputFileType );
    (*functionReadInputData)( filename, &readData, userOptions );
    
    
    // 'userOptions->MpcValue' is the conversion factor from the units in the input data file to Mpc units - do the next computation only if userOptions->MpcValue!=1
    if ( userOptions->MpcValue!=Real(1.) )
    {
        for (size_t i=0; i<userOptions->boxCoordinates.size(); ++i)
            userOptions->boxCoordinates[i] /= userOptions->MpcValue;
        
        size_t noParticles = readData.noParticles();  // returns the number of particles
        float *positions = readData.position();  //returns pointer to array storing particle positions
        for (size_t i=0; i<noParticles*NO_DIM; ++i)
            positions[i] /= userOptions->MpcValue;
        
        // if the user inserted user defined sampling points, divide the coordinates of those points by the normalization factor
        if ( readData.noSamples()!=size_t(0) )
        {
            float *sampling = readData.sampling();  // returns pointer to user defined sampling coordinates
            float *delta = readData.delta();        // returns pointer to user defined sampling cell sizes
            for (size_t i=0; i<readData.noSamples()*NO_DIM; ++i)
            {
                sampling[i] /= userOptions->MpcValue;
                delta[i] /= userOptions->MpcValue;
            }
        }
    }
    
    
    // now store the data in the 'Particle_data list'. It also copies the user given sampling coordinates, if any - none in this case.
    redshiftSpaceOn = userOptions->transformToRedshiftSpaceOn;
    redshiftSpaceVector = Pvector<Real,NO_DIM>( &(userOptions->transformToRedshiftSpace[0]) );

//    size_t noParticles = readData.noParticles();
//    readData.scalar( noParticles );
//    for (size_t i=0; i<noParticles; ++i)
//        readData.scalar()[i] = 1.;
    
    readData.transferData( p, samplingCoordinates );
    
//    for (int i=0; i<10; ++i) std::cout << "\nscalar " << (*p).at(i).scalar(0);  
}










//! Functions for writing the output data


//! short class to easily change between writing the output to different file types while the program is running. The user can choose the output file type using the '--output' option.
class OutputData
{
    public:
    int outputFileType;
    
    // class constructor - initializes the type of the output file
    OutputData(int fileType)
    {
        this->outputFileType = fileType;
    }
    
    // the function that calls the function doing the actual writing
    template <typename T>
    void write(T &dataToWrite,
                std::string filename,
                std::string variableName,
                User_options const &userOptions)
    {
        if ( this->outputFileType==101 )
            writeBinaryFile( dataToWrite, filename, variableName, userOptions );        // writes the data to a binary file (see "binary_io.cc" for function definition)
        else if ( this->outputFileType==111 )
            writeTextFile( dataToWrite, filename, variableName, userOptions );          // writes the data to a text file (see "text_io.cc" for function definition)
        else if ( this->outputFileType==112 )
            writeTextFile_gridIndex( dataToWrite, filename, variableName, userOptions );// writes the data to a text file, but on each line writes also the coordinates of the grid cell corresponding to the result being written (see "text_io.cc" for function definition)
        else if ( this->outputFileType==113 )
            writeTextFile_samplingPosition( dataToWrite, filename, variableName, userOptions ); // writes the data to a text file, but on each line writes also the sampling point coordinates corresponding to the result being written  (see "text_io.cc" for function definition)
        else if ( this->outputFileType==114 )
            writeTextFile_redshiftConePosition( dataToWrite, filename, variableName, userOptions ); // writes the data to a text file, but on each line writes also the sampling point coordinates for a redshift cone grid  (see "text_io.cc" for function definition)
        else if ( this->outputFileType==121 )
            writeMyFile( dataToWrite, filename, variableName, userOptions );   // the format that I use where it write to a binary file with a header that describes the data stored in it
        else if ( this->outputFileType==100 )
            writeSpecialFile( dataToWrite, filename, variableName, userOptions );   // the format that I use where it write to a binary file with a header that describes the data stored in it
    }
};


/*! This function writes the output data to a file, each different 'variable' being written to a separate file. 
You can modify this function as you please. */
void writeOutputData(Quantities &uQuantities,
                     Quantities &aQuantities,
                     User_options const &userOptions)
{
    // select the file type for output
    OutputData output( userOptions.outputFileType );
    
    // output the desired quantities to file/files
    // outputs the density
    if ( userOptions.uField.density )
        output.write( uQuantities.density, userOptions.outputFilename + ".den", "density", userOptions );
    if ( userOptions.aField.density )
        output.write( aQuantities.density, userOptions.outputFilename + ".a_den", "volume averaged density", userOptions );
    
    // outputs the velocity
    if ( userOptions.uField.velocity )
        output.write( uQuantities.velocity, userOptions.outputFilename + ".vel", "velocity", userOptions );
    if ( userOptions.aField.velocity )
        output.write( aQuantities.velocity, userOptions.outputFilename + ".a_vel", "volume averaged velocity", userOptions );
    
    // outputs the velocity gradient
    if ( userOptions.uField.velocity_gradient )
        output.write( uQuantities.velocity_gradient, userOptions.outputFilename + ".velGrad", "velocity gradient", userOptions );
    if ( userOptions.aField.velocity_gradient )
        output.write( aQuantities.velocity_gradient, userOptions.outputFilename + ".a_velGrad", "volume averaged velocity gradient", userOptions );
    
    // outputs the velocity divergence
    if ( userOptions.uField.velocity_divergence )
        output.write( uQuantities.velocity_divergence, userOptions.outputFilename + ".velDiv", "velocity divergence", userOptions );
    if ( userOptions.aField.velocity_divergence )
        output.write( aQuantities.velocity_divergence, userOptions.outputFilename + ".a_velDiv", "volume averaged velocity divergence", userOptions );
    
    // outputs the velocity shear
    if ( userOptions.uField.velocity_shear )
        output.write( uQuantities.velocity_shear, userOptions.outputFilename + ".velShear", "velocity shear", userOptions );
    if ( userOptions.aField.velocity_shear )
        output.write( aQuantities.velocity_shear, userOptions.outputFilename + ".a_velShear", "volume averaged velocity shear", userOptions );
    
    // outputs the velocity vorticity
    if ( userOptions.uField.velocity_vorticity )
        output.write( uQuantities.velocity_vorticity, userOptions.outputFilename + ".velVort", "velocity vorticity", userOptions );
    if ( userOptions.aField.velocity_vorticity )
        output.write( aQuantities.velocity_vorticity, userOptions.outputFilename + ".a_velVort", "volume averaged velocity vorticity", userOptions );
    
    // outputs the velocity standard deviation
    if ( userOptions.aField.velocity_std )
        output.write( aQuantities.velocity_std, userOptions.outputFilename + ".a_velStd", "volume averaged velocity standard deviation", userOptions );
    
    // outputs the scalar fields
    if ( userOptions.uField.scalar )
        output.write( uQuantities.scalar, userOptions.outputFilename + ".scalar", "scalar", userOptions );
    if ( userOptions.aField.scalar )
        output.write( aQuantities.scalar, userOptions.outputFilename + ".a_scalar", "volume averaged scalar", userOptions );
    
    // outputs the scalar fields gradient
    if ( userOptions.uField.scalar_gradient )
        output.write( uQuantities.scalar_gradient, userOptions.outputFilename + ".scalarGrad", "scalar gradient", userOptions );
    if ( userOptions.aField.scalar_gradient )
        output.write( aQuantities.scalar_gradient, userOptions.outputFilename + ".a_scalarGrad", "volume averaged scalar gradient", userOptions );
}























#endif
