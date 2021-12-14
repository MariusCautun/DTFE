/*
 *  Copyright (c) 2013       Marius Cautun
 *
 *                           Institute for Computational Cosmology
 *                           Durham University, Durham, UK
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



// reads the modified gravity forces
void readForces_MOG(std::string fileName,
                    Read_data<Real> *readData,
                    User_options &userOptions,
                    size_t const noParticlesThisFile,
                    int const noScalars,
                    size_t *numberParticlesRead);



/*! This function reads a Gadget snapshots saved in a single or multiple files.
 *  It also reads the GR and MOG forces from a separate file.
 * */
void readGadgetFile_MOG(std::string filename,
                        Read_data<Real> *readData,
                        User_options *userOptions)
{
    int gadgetFileType;     // stores the type of the gadget file (1 or 2)
    int noBytesPos, noBytesVel;     // stores the data type used to write positions and velocities
    bool swapEndian;        // true if need to swap the endianess of the data
    size_t noParticles;     // total number of particles to be read from the file
    Gadget_header gadgetHeader; // the gadget header for the first file
    
    
    // get the input gadget data characteristics: file type, endianess, particle number and number of files. The following function also reserves memory for the positions, weights and velocity data (if requested so).
    initializeGadget( filename, readData, userOptions, &gadgetHeader, &gadgetFileType, &noBytesPos, &noBytesVel, &swapEndian, &noParticles );
    
    
    MESSAGE::Message message( userOptions->verboseLevel );
    std::string fileName = filename;
    std::string fileNameMOG;
    
    
    // reserve memory for the scalar particle data
    int noScalars = atoi( userOptions->additionalOptions[1].c_str() );
    if ( noScalars>0  and  noScalars<=noScalarComp )
        readData->scalar( readData->noParticles() );  // particle scalar quantity
    else if ( noScalars>noScalarComp )
    {
        std::cout << "~~~ ERROR ~~~ The file required number of scalar components = " << noScalars << " is larger than the one specified when the program has been compiled = " << noScalarComp << "! The program ended unsucessfully!\n";
        exit( EXIT_FAILURE );
    }
    
    
    // read the data in each file
    size_t numberParticlesRead  = 0;   // the total number of particles read after each file
    size_t numberParticlesRead2 = 0;   // the total number of particles read after each file
    
    // iterate over all the files and read the data
    for (int i=0; i<gadgetHeader.num_files; ++i )
    {
        fileName = gadgetHeader.filename( filename, i );
        message << "Reading GADGET snapshot file '" << fileName << "' which is file " << i+1 << " of " << gadgetHeader.num_files << " files...\n" << MESSAGE::Flush;

        // call the function that reads the data
        readGadgetData( fileName, readData, *userOptions, gadgetFileType, noBytesPos, noBytesVel, swapEndian, &numberParticlesRead );
        
        // read the forces from that file
        fileNameMOG = gadgetHeader.filename( userOptions->additionalOptions[0], i+1 );
        size_t noParticlesThisFile = numberParticlesRead - numberParticlesRead2;
        if ( noScalars > 0 )
            readForces_MOG( fileNameMOG, readData, *userOptions, noParticlesThisFile, noScalars, &numberParticlesRead2 );
    }
    
    
    // check if the data needs to be swapped
    if ( swapEndian )
    {
        if ( userOptions->readParticleData[0] )
            ByteSwapArray( readData->position(), NO_DIM*noParticles );
        if ( userOptions->readParticleData[1] )
            ByteSwapArray( readData->weight(), noParticles );
#ifdef VELOCITY
        if ( userOptions->readParticleData[2] )
            ByteSwapArray( readData->velocity(), NO_DIM*noParticles );
#endif
#ifdef SCALAR
        if ( noScalars>0 )
            ByteSwapArray( readData->scalar(), NO_DIM*noParticles );
#endif
    }
    
    return;
}








/*! This function reads the forces from a binary file. 
*/
void readForces_MOG(std::string fileName,
                    Read_data<Real> *readData,
                    User_options &userOptions,
                    size_t const noParticlesThisFile,
                    int const noScalars,
                    size_t *numberParticlesRead)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "\t reading force data from file '" << fileName << "' ... " << MESSAGE::Flush;
    
    // open the file
    std::fstream inputFile;
    openInputBinaryFile( inputFile, fileName );
    int buffer1, buffer2;
    
    
    // read the data block
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    
    double *data = new double[ noParticlesThisFile * noScalars ];
    size_t readBytes = noParticlesThisFile * sizeof(double) * noScalars;
    inputFile.read( reinterpret_cast<char *>( data ), readBytes );
    
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    if ( buffer1!=buffer2 )
        throwError( "The integers before and after the particle 'MOG forces' data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." );
    
    
    // copy the data to the scalar block
    size_t alreadyRead = (*numberParticlesRead) * noScalarComp;
    Real *scalar = &( readData->scalar()[ alreadyRead ] );  // pointer to the first place where to start writing
    for (size_t i=0; i<noParticlesThisFile; ++i)
        for (int j=0; j<noScalars; ++j)
            scalar[ i*noScalarComp + j ] = data[ i*noScalars + j ];
    delete[] data;
    message << "Done.";

    inputFile.close();
    message << "\n";
    (*numberParticlesRead) += noParticlesThisFile;
}








