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



/* You can use this file to easily customize the DTFE software data input/output. The functions below were left intentionally empty such that you can put your read/write to file code there. Use the rest of files in this directory for examples on how to read from / write to different types of files (I would recommend as starting point the "text_io.cc" file). */



/* Your custom functions for reading the input data. */
void readMyFile(std::string filename,
                Read_data<float> *readData,
                User_options *userOptions)
{
    MESSAGE::Message message( userOptions->verboseLevel );
    message << "Reading the input data from my custom type file '" << filename << "' ... " << MESSAGE::Flush;
    
    // open the file
    std::fstream inputFile;
    openInputBinaryFile( inputFile, filename );
//     openInputTextFile( inputFile, filename );
    
    
    // Read the number of particles and the coordinates of the box encompassing the data
    int noParticles;    // the number of particles
    float boxCoordinates[2*NO_DIM];
    // read 'noParticles' and 'boxCoordinates'
    for (size_t i=0; i<2*NO_DIM; ++i) userOptions->boxCoordinates[i] = boxCoordinates[i];
    
    
    // reserve memory for the input data - select only the ones you use
    Real *positions = readData->position(noParticles);  //particle positions
//     Real *weights = readData->weight(noParticles);      //particle weights (e.g. weights = particle masses)
//     Real *velocities = readData->velocity(noParticles); //particle velocities
//     Real *scalars = readData->scalar(noParticles);      //scalar component for each particle
    
    
    // read the rest of particle data from the file
    
    
    checkFileOperations( inputFile, "read from" );   // check that the data reading was succesful
    inputFile.close();
    
    message << "Done.\n";
}


/* Your custom function for writting the output data to a file. */
void writeMyFile(std::vector<Real> &dataToWrite,
                 std::string filename,
                 std::string variableName,
                 User_options const &userOptions)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "Writing the " << variableName << " to my custom type file '" << filename << "' ...  " << MESSAGE::Flush;
    
    // open the file
    std::fstream outputFile;
    openOutputBinaryFile( outputFile, filename );
//     openOutputTextFile( outputFile, filename );
    
    
    // write the data to file
    
    
    checkFileOperations( outputFile, "write to" );   // check that the data writing was succesful
    outputFile.close();
    message << "Done.\n";
}
template <size_t N>
void writeMyFile(std::vector< Pvector<Real,N> > &dataToWrite,
                 std::string filename,
                 std::string variableName,
                 User_options const &userOptions)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "Writing the " << variableName << " to my custom type file '" << filename << "' ...  " << MESSAGE::Flush;
    
    // open the file
    std::fstream outputFile;
    openOutputBinaryFile( outputFile, filename );
//     openOutputTextFile( outputFile, filename );
    
    
    // write the data to file
    
    
    checkFileOperations( outputFile, "write to" );   // check that the data writing was succesful
    outputFile.close();
    message << "Done.\n";
}


