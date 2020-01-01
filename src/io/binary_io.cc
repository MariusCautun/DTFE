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



/* This file contains functions for reading and writing the data to a binary file. */


/* This function reds the input data from a binary file.
This function works when the binary file has the following format:
    1) file begins with an int value that gives the number of particles
    2) continues with 6 float values that give the coordinates of the box encompassing the data (i.e. xMin, xMax, yMin, yMax, zMin, zMax)
    3) the next part of the file gives the positions of the particles which are saved as x1, y1, z1, x2, y2, z2, x3, ... (with x1 = x coordinate for particle 1, y1 = y coordinate for particle 1, etc ...)
    4) after the particle positions, come the particle weights array saved as w1, w2, w3, ... (with w1 = weight of particle 1, etc ...)
    5) after the weights array, comes the velocity array which are is saved as vx1, vy1, vz1, vx2, vy2, vz2, ... (with vx1 = the velocity along the x direction of particle 1, ...)

NOTE: In this example the data is saved in the binary file in single precision format (float 32 bytes).
*/
void readBinaryFile(std::string filename,
                    Read_data<float> *readData,
                    User_options *userOptions)
{
    MESSAGE::Message message( userOptions->verboseLevel );
    message << "Reading the input data from the binary file '" << filename << "' ... " << MESSAGE::Flush;
    
    // open the file
    std::fstream inputFile;
    openInputBinaryFile( inputFile, filename );
    
    
    // read the number of particles and the box coordinates from the binary file
    int noParticles;
    float boxCoordinates[2*NO_DIM];
    inputFile.read( reinterpret_cast<char *>(&noParticles), sizeof(noParticles) );
    inputFile.read( reinterpret_cast<char *>(boxCoordinates), sizeof(float) );
    for (size_t i=0; i<2*NO_DIM; ++i)
        userOptions->boxCoordinates[i] = boxCoordinates[i];
    
    
    // reserve memory for the input data
    Real *positions = readData->position(noParticles);  //particle positions
    Real *weights = readData->weight(noParticles);      //particle weights (e.g. weights = particle masses)
    Real *velocities = readData->velocity(noParticles); //particle velocities
    
    
    // read the rest of the input data: positions, weights and velocities
    size_t dataSize = noParticles * sizeof(float) * NO_DIM;    // number of data bytes that store the particle positions (3*4 bytes per particle)
    inputFile.read( reinterpret_cast<char *>(positions), dataSize );
    
    dataSize = noParticles * sizeof(float);    // number of data bytes that store the particle weights (1*4 bytes per particle)
    inputFile.read( reinterpret_cast<char *>(positions), dataSize );
    
    dataSize = noParticles * sizeof(float) * NO_DIM;    // number of data bytes that store the particle velocities (3*4 bytes per particle)
    inputFile.read( reinterpret_cast<char *>(positions), dataSize );
    
    checkFileOperations( inputFile, "read from" );   // check that the data reading was succesful
    inputFile.close();
    
    message << "Done.\n";
}



/* This function writes the data to a binary file.*/
template <typename T>
void writeBinaryFile(T const &dataToWrite,
                     std::string filename,
                     std::string variableName,
                     User_options const &userOptions)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "Writing the " << variableName << " to the binary file '" << filename << "' ...  " << MESSAGE::Flush;
    
    // open the file
    std::fstream outputFile;
    openOutputBinaryFile( outputFile, filename );
    
    
    // write the data to file
    // it is very simple to write data to a binary file, just do: outputFile.write( reinterpret_cast<char const *>(&(dataToWrite[0])), dataSize ). The following implements a recursive method that deals with possible issuess in writing very large data sets.
    size_t dataSize = dataToWrite.size()*sizeof(dataToWrite[0]);    //total number of data bytes to be written
    size_t maxSize = 256*256*256;   // write at most 256^3 elements at once, otherwise the write function might fail
    size_t noRepeats = size_t( dataToWrite.size() / maxSize ), currentPosition = 0;
    size_t tempBuffer = maxSize * sizeof(dataToWrite[0]);
    for (size_t i=0; i<noRepeats; ++i)  //write in blocks of 256^3 elements
    {
        outputFile.write( reinterpret_cast<char const *>(&(dataToWrite[currentPosition])), tempBuffer );
        currentPosition += maxSize;
    }
    tempBuffer = (dataToWrite.size() - currentPosition) * sizeof(dataToWrite[0]);    // write everything else that is left
    outputFile.write( reinterpret_cast<char const *>(&(dataToWrite[currentPosition])), tempBuffer );
    
    checkFileOperations( outputFile, "write to" );   // check that the data writing was succesful
    outputFile.close();
    
    message << "Done.\n";
}


