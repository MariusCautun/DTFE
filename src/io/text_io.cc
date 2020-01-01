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



/* This file implements functions for reading and writing the data from/to text files. */
#include <cmath>



//! Functions for reading in the data from a text file

/* This function shows an example of how to read the particle data from a text file. The text file has the following data:
    1 line: number of particles
    2 line: "xMin xMax yMin yMax zMin zMax" for the box fully encompassing the given particle data
    3rd line: "posX posY posZ velX velY velZ weight scalar" for the 1st particle
    rest of the lines: same as line 3, but for the other particles
    
To read the file one must open in. It can do so using functions:
    openInputTextFile( fstream file, string filename ) - to open a text file for reading
    openInputBinaryFile( fstream file, string filename ) - to open a binary file for reading

After reading the number of particle in the file, memory can be reserved for reading the data. One can do so using the function "Read_data::position(size_t noParticles)" which allocates memory for storing NO_DIM*noParticles 'Reals' that are the particle positions. The function returns a pointer to the allocated emory block, so a typical call of the function is:
    Real *positions = readData->position( noParticles );
Once the memory for storing particle positions was allocated, the pointer to the memory can also be called using:
    Real *positions2 = readData->position();    // this will give an error if no memory was allocated beforehand
One has similar functions for allocating memory for the rest of the particle properties quantities:
    position:
                Real* Read_data :: position( size_t noParticles );  //allocate memory and return pointer
                Real* Read_data :: position();  // return pointer (can be called only after calling the previous function)
    velocity:
                Real* Read_data :: velocity( size_t noParticles );  //allocate memory and return pointer
                Real* Read_data :: velocity();  // return pointer
    weight:
                Real* Read_data :: weight( size_t noParticles );  //allocate memory and return pointer
                Real* Read_data :: weight();  // return pointer
    scalar:
                Real* Read_data :: scalar( size_t noParticles );  //allocate memory and return pointer
                Real* Read_data :: scalar();  // return pointer
NOTE: When calling the above functions the number of particle has to be the same for position, velocity, weight and scalar.
*/
void readTextFile(std::string filename,
                  Read_data<float> *readData,
                  User_options *userOptions)
{
    MESSAGE::Message message( userOptions->verboseLevel );
    message << "Reading the input data from the text file '" << filename << "' ... " << MESSAGE::Flush;
    // open the text file for reading
    std::fstream inputFile;
    openInputTextFile( inputFile, filename );
    
    
    // read the first line that gives the number of particles and the second line that gives the coordinates of the box encompassing the data
    size_t noParticles;
    inputFile >> noParticles;
    for (int i=0; i<2*NO_DIM; ++i)
        inputFile >> userOptions->boxCoordinates[i];
    
    
    // assign memory to store the particle data being read from file
    // the following assumes that the text file has a line for each particle with: posX, posY, posZ, velX, velY, velZ, weight(=particle mass), scalar(1 component)
    Real *positions = readData->position(noParticles);  //particle positions
    Real *velocities = readData->velocity(noParticles); //particle velocities
    Real *weights = readData->weight(noParticles);      //particle weights (e.g. weights = particle/galaxy masses)
    Real *scalars = readData->scalar(noParticles);      //scalar component for each particle
    
    
    // now read the particle data
    for (int i=0; i<noParticles; ++i)
    {
        // read positions
        for (int j=0; j<NO_DIM; ++j)
            inputFile >> positions[NO_DIM*i+j];
        // read velocities
//         for (int j=0; j<noVelComp; ++j)
//             inputFile >> velocities[noVelComp*i+j];
        // read weight
        inputFile >> weights[i];
        // read scalar components
//         for (int j=0; j<noScalarComp; ++j)  // must have noScalarComp==1 otherwise there will be a runtime error or the data will be read wrongly
//             inputFile >> scalars[noScalarComp*i+j];
    }
    
    // close the input file
    checkFileOperations( inputFile, "read from" );   // check that the data reading was succesful
    inputFile.close();
    message << "Done.\n";
}


/* This function reads only the particle positions from a text file that contains only the positions. */
void readTextFile_positions(std::string filename,
                            Read_data<float> *readData,
                            User_options *userOptions)
{
    MESSAGE::Message message( userOptions->verboseLevel );
    message << "Reading the particle position data from the text file '" << filename << "' ... " << MESSAGE::Flush;
    // open the text file for reading
    std::fstream inputFile;
    openInputTextFile( inputFile, filename );
    
    // read the first line that gives the number of particles and the second line that gives the coordinates of the box encompassing the data
    size_t noParticles;
    inputFile >> noParticles;
    for (int i=0; i<2*NO_DIM; ++i)
        inputFile >> userOptions->boxCoordinates[i];
    
    // assign memory to store the particle data being read from file
    // the following assumes that the text file has a line for each particle with: posX, posY, posZ
    Real *positions = readData->position(noParticles);  //particle positions
    
    // now read the particle data
    for (int i=0; i<noParticles; ++i)
    {
        // read positions
        for (int j=0; j<NO_DIM; ++j)
            inputFile >> positions[NO_DIM*i+j];
    }
    
    // close the input file
    checkFileOperations( inputFile, "read from" );   // check that the data reading was succesful
    inputFile.close();
    message << "Done.\n";
}



/* Example of how to read the particle positions from a text file and also how to read in user defined sampling points for the grid where to interpolate the output fields. */
void readTextFile_userDefinedSampling(std::string filename,
                                      Read_data<float> *readData,
                                      User_options *userOptions)
{
    // 1st part of the functions: read the particle positions from the input file
    MESSAGE::Message message( userOptions->verboseLevel );
    message << "Reading the particle position data from the text file '" << filename << "' ... " << MESSAGE::Flush;
    // open the text file for reading
    std::fstream inputFile;
    openInputTextFile( inputFile, filename );
    
    // read the first line that gives the number of particles and the second line that gives the coordinates of the box encompassing the data
    size_t noParticles;
    inputFile >> noParticles;
    for (int i=0; i<2*NO_DIM; ++i)
        inputFile >> userOptions->boxCoordinates[i];
    
    // assign memory to store the particle data being read from file
    // the following assumes that the text file has a line for each particle with: posX, posY, posZ
    Real *positions = readData->position(noParticles);  //particle positions
    
    // now read the particle data
    for (int i=0; i<noParticles; ++i)
    {
        // read positions
        for (int j=0; j<NO_DIM; ++j)
            inputFile >> positions[NO_DIM*i+j];
    }
    
    // close the input file
    checkFileOperations( inputFile, "read from" );   // check that the data reading was succesful
    inputFile.close();
    message << "Done.\n";
    
    
    // 2nd part of the function: read the user defined sampling points from a different input text file
    // for the purpose of this example I suppose that you inserted the name of the file giving the user defined sampling points using the program option "--options fileName" (this option can be used to give to the program additional arguments on top of the default ones). Each additional value supplied to the "--options" option is stored in a vector of strings called "userOptions->additionalOptions".
     message << "Reading the user defined sampling points from the text file '" << filename << "' ... " << MESSAGE::Flush;
    std::string filename2 = userOptions->additionalOptions[0];  // the name of the file containing the user defined sampling points
    openInputTextFile( inputFile, filename2 );
    
    // first line of this file should contain the number of user defined sampling points
    // on each other lines one should give 6 values: X, Y, Z, dX, dY, dZ where (X,y,Z) is the coordinate of the center of the user defined sampling point and (dX,dY,dZ) is the size of the grid cell associated to each sampling point. The (dX,dY,dZ) values are needed only when computing volume averaged values for the output fields.
    // In the following we suppose that the file has the above format: 1st line = number of sampling points, all other lines have 6 columns giving: X, Y, Z, dX, dY and dZ
    // read the number of sampling points
    size_t noSamplingPoints;
    inputFile >> noSamplingPoints;
    
    // reserve memory for reading in the user defined sampling points
    float *samples = readData->sampling(noSamplingPoints);  // reserve memory for the coordinates of the sampling points
    float *delta = readData->delta(noSamplingPoints);       // reserve memory for the cell size associated to each sampling point
    
    // read the data from the file
    for (int i=0; i<noSamplingPoints; ++i)
    {
        // read sampling point positions
        for (int j=0; j<NO_DIM; ++j)
            inputFile >> samples[NO_DIM*i+j];
        // read sampling cell sizes
        for (int j=0; j<NO_DIM; ++j)
            inputFile >> delta[NO_DIM*i+j];
    }
    
    // close the input file
    checkFileOperations( inputFile, "read from" );   // check that the data reading was succesful
    inputFile.close();
    message << "Done.\n";
}




//! Functions for writing the data to a text file

/* This function writes the data to a text file. 
 The field components at each sampling point are written in a single line, in different columns (i.e. number of columns = number of field components at a given point).*/
void writeTextFile(std::vector<Real> &dataToWrite,
                    std::string filename,
                    std::string variableName,
                    User_options const &userOptions)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "Writing the " << variableName << " to the text file '" << filename << "' ...  " << MESSAGE::Flush;
    
    // open the file
    std::fstream outputFile;
    openOutputTextFile( outputFile, filename );
    
    // write the data to file
    for (size_t i=0; i<dataToWrite.size(); ++i)
        outputFile << dataToWrite[i] << "\n";
    
    checkFileOperations( outputFile, "write to" );   // check that the data writing was succesful
    outputFile.close();
    message << "Done.\n";
}
template <size_t N>
void writeTextFile(std::vector< Pvector<Real,N> > &dataToWrite,
                    std::string filename,
                    std::string variableName,
                    User_options const &userOptions)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "Writing the " << variableName << " to the text file '" << filename << "' ...  " << MESSAGE::Flush;
    
    // open the file
    std::fstream outputFile;
    openOutputTextFile( outputFile, filename );
    
    // write the data to file
    for (size_t i=0; i<dataToWrite.size(); ++i)
    {
        for (size_t j=0; j<N; ++j)
            outputFile << dataToWrite[i][j] << "\t";
        outputFile << "\n";
    }
    
    checkFileOperations( outputFile, "write to" );   // check that the data writing was succesful
    outputFile.close();
    message << "Done.\n";
}



/* This function writes the data to a text file. 
On each line it writes:
        1) the grid indices
        2) the field value at that grid position (writes all the field components in different columns)
*/
void writeTextFile_gridIndex(std::vector<Real> &dataToWrite,
                             std::string filename,
                             std::string variableName,
                             User_options const &userOptions)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "Writing the " << variableName << " to the text file '" << filename << "' ...  " << MESSAGE::Flush;
    
    if ( userOptions.userDefinedSampling )
        throwError( "You cannot use the function 'writeTextFile_gridIndex' to write the data to a text file when using user defined coordinates since there are no grid indices associated to each sampling point." );
    
    
    // open the file
    std::fstream outputFile;
    openOutputTextFile( outputFile, filename );
    
    
    // write the data to file
    size_t const *grid = &(userOptions.gridSize[0]);
#if NO_DIM==2   // for the 2D case
    for (size_t i=0; i<grid[0]; ++i)
        for (size_t j=0; j<grid[1]; ++j)
        {
            size_t index = i*grid[1] + j;
            outputFile << i << "\t" << j << "\t";
            outputFile << dataToWrite[index] << "\n";
        }
#elif NO_DIM==3 // for the 3D case
    for (size_t i=0; i<grid[0]; ++i)
        for (size_t j=0; j<grid[1]; ++j)
            for(size_t k=0; k<grid[2]; ++k)
            {
                size_t index = i*grid[1]*grid[2] + j*grid[2] + k;
                outputFile << i << "\t" << j << "\t" << k << "\t";
                outputFile << dataToWrite[index] << "\n";
            }
#endif
    
    checkFileOperations( outputFile, "write to" );   // check that the data writing was succesful
    outputFile.close();
    message << "Done.\n";
}
template <typename T, size_t N>
void writeTextFile_gridIndex(std::vector< Pvector<T,N> > &dataToWrite,
                             std::string filename,
                             std::string variableName,
                             User_options const &userOptions)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "Writing the " << variableName << " to the text file '" << filename << "' ...  " << MESSAGE::Flush;
    
    if ( userOptions.userDefinedSampling )
        throwError( "You cannot use the function 'writeTextFile_gridIndex' to write the data to a text file when using user defined coordinates since there are no grid indices associated to each sampling point." );
    
    
    // open the file
    std::fstream outputFile;
    openOutputTextFile( outputFile, filename );
    
    
    // write the data to file
    size_t const *grid = &(userOptions.gridSize[0]);
#if NO_DIM==2   // for the 2D case
    for (size_t i=0; i<grid[0]; ++i)
        for (size_t j=0; j<grid[1]; ++j)
        {
            size_t index = i*grid[1] + j;
            outputFile << i << "\t" << j << "\t";
            for (size_t i1=0; i1<N; ++i1)
                outputFile << dataToWrite[index][i1] << "\t";
            outputFile << "\n";
        }
#elif NO_DIM==3 // for the 3D case
    for (size_t i=0; i<grid[0]; ++i)
        for (size_t j=0; j<grid[1]; ++j)
            for(size_t k=0; k<grid[2]; ++k)
            {
                size_t index = i*grid[1]*grid[2] + j*grid[2] + k;
                outputFile << i << "\t" << j << "\t" << k << "\t";
                for (size_t i1=0; i1<N; ++i1)
                    outputFile << dataToWrite[index][i1] << "\t";
                outputFile << "\n";
            }
#endif
    
    checkFileOperations( outputFile, "write to" );   // check that the data writing was succesful
    outputFile.close();
    message << "Done.\n";
}



/* This function writes the data to a text file. This particular function works only for regular rectangular grids (you can write a similar version for redshift coordinate grids).
On each line it writes:
        1) the coordinates of the sampling point
        2) the field value at that sampling point (writes all the field components in different columns)
*/
void writeTextFile_samplingPosition(std::vector<Real> &dataToWrite,
                                    std::string filename,
                                    std::string variableName,
                                    User_options const &userOptions)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "Writing the " << variableName << " to the text file '" << filename << "' ...  " << MESSAGE::Flush;
    
    if ( userOptions.userDefinedSampling or userOptions.redshiftConeOn )
        throwError( "You cannot use the function 'writeTextFile_samplingPosition' to write the data to a text file when using redshift cone or user defined coordinates since teh sampling coordinates are incorrect in this case." );
    
    
    // open the file
    std::fstream outputFile;
    openOutputTextFile( outputFile, filename );
    
    
    // write the data to file
    size_t const *grid = &(userOptions.gridSize[0]);
    Box boxCoordinates = userOptions.region;    // 'userOptions.region' stores the boundaries of the box in which the fields where interpolated to a regular grid
    Real dx[NO_DIM];
    for (size_t i=0; i<NO_DIM; ++i) dx[i] = (boxCoordinates[2*i+1]-boxCoordinates[2*i]) / grid[i];
    
#if NO_DIM==2   // for the 2D case
    for (size_t i=0; i<grid[0]; ++i)
        for (size_t j=0; j<grid[1]; ++j)
        {
            size_t index = i*grid[1] + j;
            Real x = boxCoordinates[0] + dx[0] * (i+0.5);
            Real y = boxCoordinates[2] + dx[1] * (j+0.5);
            outputFile << x << "\t" << y << "\t";
            outputFile << dataToWrite[index] << "\n";
        }
#elif NO_DIM==3 // for the 3D case
    for (size_t i=0; i<grid[0]; ++i)
        for (size_t j=0; j<grid[1]; ++j)
            for(size_t k=0; k<grid[2]; ++k)
            {
                size_t index = i*grid[1]*grid[2] + j*grid[2] + k;
                Real x = boxCoordinates[0] + dx[0] * (i+0.5);
                Real y = boxCoordinates[2] + dx[1] * (j+0.5);
                Real z = boxCoordinates[4] + dx[2] * (k+0.5);
                outputFile << x << "\t" << y << "\t" << z << "\t";
                outputFile << dataToWrite[index] << "\n";
            }
#endif
    
    checkFileOperations( outputFile, "write to" );   // check that the data writing was succesful
    outputFile.close();
    message << "Done.\n";
}
template <typename T, size_t N>
void writeTextFile_samplingPosition(std::vector< Pvector<T,N> > &dataToWrite,
                                    std::string filename,
                                    std::string variableName,
                                    User_options const &userOptions)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "Writing the " << variableName << " to the text file '" << filename << "' ...  " << MESSAGE::Flush;
    
    if ( userOptions.userDefinedSampling or userOptions.redshiftConeOn )
        throwError( "You cannot use the function 'writeTextFile_samplingPosition' to write the data to a text file when using redshift cone or user defined coordinates since teh sampling coordinates are incorrect in this case." );
    
    
    // open the file
    std::fstream outputFile;
    openOutputTextFile( outputFile, filename );
    
    
    // write the data to file
    size_t const *grid = &(userOptions.gridSize[0]);
    Box boxCoordinates = userOptions.region;    // 'userOptions.region' stores the boundaries of the box in which the fields where interpolated to a regular grid
    Real dx[NO_DIM];
    for (size_t i=0; i<NO_DIM; ++i) dx[i] = (boxCoordinates[2*i+1]-boxCoordinates[2*i]) / grid[i];
    
#if NO_DIM==2   // for the 2D case
    for (size_t i=0; i<grid[0]; ++i)
        for (size_t j=0; j<grid[1]; ++j)
        {
            size_t index = i*grid[1] + j;
            Real x = boxCoordinates[0] + dx[0] * (i+0.5);
            Real y = boxCoordinates[2] + dx[1] * (j+0.5);
            outputFile << x << "\t" << y << "\t";
            for (size_t i1=0; i1<N; ++i1)
                outputFile << dataToWrite[index][i1] << "\t";
            outputFile << "\n";
        }
#elif NO_DIM==3 // for the 3D case
    for (size_t i=0; i<grid[0]; ++i)
        for (size_t j=0; j<grid[1]; ++j)
            for(size_t k=0; k<grid[2]; ++k)
            {
                size_t index = i*grid[1]*grid[2] + j*grid[2] + k;
                Real x = boxCoordinates[0] + dx[0] * (i+0.5);
                Real y = boxCoordinates[2] + dx[1] * (j+0.5);
                Real z = boxCoordinates[4] + dx[2] * (k+0.5);
                outputFile << x << "\t" << y << "\t" << z << "\t";
                for (size_t i1=0; i1<N; ++i1)
                    outputFile << dataToWrite[index][i1] << "\t";
                outputFile << "\n";
            }
#endif
    
    checkFileOperations( outputFile, "write to" );   // check that the data writing was succesful
    outputFile.close();
    message << "Done.\n";
}




/* This function writes the data to a text file. This particular function works only for regular rectangular grids (you can write a similar version for redshift coordinate grids).
On each line it writes:
        1) the coordinates of the sampling point (can choose between (r,theta,psi) or (x,y,z) )
        2) the field value at that sampling point (writes all the field components in different columns)
*/
void writeTextFile_redshiftConePosition(std::vector<Real> &dataToWrite,
                                        std::string filename,
                                        std::string variableName,
                                        User_options const &userOptions)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "Writing the " << variableName << " to the text file '" << filename << "' ...  " << MESSAGE::Flush;
    
    if ( userOptions.userDefinedSampling )
        throwError( "You cannot use the function 'writeTextFile_redshiftConePosition' to write the data to a text file when using user defined coordinates since the program doesn't know how to compute the redshift cone coordinates." );
    if ( not userOptions.redshiftConeOn )
        throwError( "The function 'writeTextFile_redshiftConePosition' can be used to write the data to a text file only when using redshift cone coordinates since it writes the redshift cone coordinates in the file too." );
    
    // open the file
    std::fstream outputFile;
    openOutputTextFile( outputFile, filename );
    
    
    // write the data to file
    size_t const *grid = &(userOptions.gridSize[0]);
    Box coneCoordinates = userOptions.redshiftCone;    // 'userOptions.redshiftCone' stores the redshift cone coordinates
    std::vector<Real> origin = userOptions.originPosition;// origin of the redshift cone
    Real dx[NO_DIM];
    for (size_t i=0; i<NO_DIM; ++i) dx[i] = (coneCoordinates[2*i+1]-coneCoordinates[2*i]) / grid[i];
    Real factor = 3.14/180.; // to transform from degrees to radians
    
#if NO_DIM==2   // for the 2D case
    for (size_t i=0; i<grid[0]; ++i)
        for (size_t j=0; j<grid[1]; ++j)
        {
            size_t index = i*grid[1] + j;
            Real r = coneCoordinates[0] + dx[0] * (i+0.5);
            Real theta = coneCoordinates[2] + dx[1] * (j+0.5);
            Real x = origin[0] + r*cos(theta*factor);
            Real y = origin[1] + r*sin(theta*factor);
            
            // to write the position r and angle coordinate theta
            outputFile << r << "\t" << theta << "\t";
            // to write the position x and y
//             outputFile << x << "\t" << y << "\t";
            outputFile << dataToWrite[index] << "\n";
        }
#elif NO_DIM==3 // for the 3D case
    for (size_t i=0; i<grid[0]; ++i)
        for (size_t j=0; j<grid[1]; ++j)
            for(size_t k=0; k<grid[2]; ++k)
            {
                size_t index = i*grid[1]*grid[2] + j*grid[2] + k;
                Real r = coneCoordinates[0] + dx[0] * (i+0.5);
                Real theta = coneCoordinates[2] + dx[1] * (j+0.5);
                Real phi = coneCoordinates[4] + dx[2] * (k+0.5);
                Real x = origin[0] + r*sin(theta*factor)*cos(phi*factor);
                Real y = origin[1] + r*sin(theta*factor)*sin(phi*factor);
                Real z = origin[2] + r*cos(theta*factor);
                
                // to write the position r and angle coordinate theta
//                 outputFile << r << "\t" << theta << "\t" << phi << "\t";
                // to write the position x and y
                outputFile << x << "\t" << y << "\t" << z << "\t";
                outputFile << dataToWrite[index] << "\n";
            }
#endif
    
    checkFileOperations( outputFile, "write to" );   // check that the data writing was succesful
    outputFile.close();
    message << "Done.\n";
}
template <typename T, size_t N>
void writeTextFile_redshiftConePosition(std::vector< Pvector<T,N> > &dataToWrite,
                                        std::string filename,
                                        std::string variableName,
                                        User_options const &userOptions)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "Writing the " << variableName << " to the text file '" << filename << "' ...  " << MESSAGE::Flush;
    
    if ( userOptions.userDefinedSampling )
        throwError( "You cannot use the function 'writeTextFile_redshiftConePosition' to write the data to a text file when using user defined coordinates since the program doesn't know how to compute the redshift cone coordinates." );
    if ( not userOptions.redshiftConeOn )
        throwError( "The function 'writeTextFile_redshiftConePosition' can be used to write the data to a text file only when using redshift cone coordinates since it writes the redshift cone coordinates in the file too." );
    
    // open the file
    std::fstream outputFile;
    openOutputTextFile( outputFile, filename );
    
    
    // write the data to file
    size_t const *grid = &(userOptions.gridSize[0]);
    Box coneCoordinates = userOptions.redshiftCone;    // 'userOptions.redshiftCone' stores the redshift cone coordinates
    std::vector<Real> origin = userOptions.originPosition;// origin of the redshift cone
    Real dx[NO_DIM];
    for (size_t i=0; i<NO_DIM; ++i) dx[i] = (coneCoordinates[2*i+1]-coneCoordinates[2*i]) / grid[i];
    Real factor = 3.14/180.; // to transform from degrees to radians
    
#if NO_DIM==2   // for the 2D case
    for (size_t i=0; i<grid[0]; ++i)
        for (size_t j=0; j<grid[1]; ++j)
        {
            size_t index = i*grid[1] + j;
            Real r = coneCoordinates[0] + dx[0] * (i+0.5);
            Real theta = coneCoordinates[2] + dx[1] * (j+0.5);
            Real x = origin[0] + r*cos(theta*factor);
            Real y = origin[1] + r*sin(theta*factor);
            
            // to write the position r and angle coordinate theta
//             outputFile << r << "\t" << theta << "\t";
            // to write the position x and y
            outputFile << x << "\t" << y << "\t";
            for (size_t i1=0; i1<N; ++i1)
                outputFile << dataToWrite[index][i1] << "\t";
            outputFile << "\n";
        }
#elif NO_DIM==3 // for the 3D case
    for (size_t i=0; i<grid[0]; ++i)
        for (size_t j=0; j<grid[1]; ++j)
            for(size_t k=0; k<grid[2]; ++k)
            {
                size_t index = i*grid[1]*grid[2] + j*grid[2] + k;
                Real r = coneCoordinates[0] + dx[0] * (i+0.5);
                Real theta = coneCoordinates[2] + dx[1] * (j+0.5);
                Real phi = coneCoordinates[4] + dx[2] * (k+0.5);
                Real x = origin[0] + r*sin(theta*factor)*cos(phi*factor);
                Real y = origin[1] + r*sin(theta*factor)*sin(phi*factor);
                Real z = origin[2] + r*cos(theta*factor);
                
                // to write the position r and angle coordinate theta
//                 outputFile << r << "\t" << theta << "\t" << phi << "\t";
                // to write the position x and y
                outputFile << x << "\t" << y << "\t" << z << "\t";
                for (size_t i1=0; i1<N; ++i1)
                    outputFile << dataToWrite[index][i1] << "\t";
                outputFile << "\n";
            }
#endif
    
    checkFileOperations( outputFile, "write to" );   // check that the data writing was succesful
    outputFile.close();
    message << "Done.\n";
}
