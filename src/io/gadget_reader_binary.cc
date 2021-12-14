/*
 *  Copyright (c) 2021       Marius Cautun
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



/* This file defines a reader for binary Gadget snapshot files. The reader given here can read all the particles and properties stored in the Gadget snapshot file by giving the appropiate options to the program. See the documentation for the options that control this behavior. */




// function that counts the number of particles for snapshots saved as multiple files
void countGadgetParticleNumber(std::string filenameRoot,
                               int const noFiles,
                               int const gadgetFileType,
                               bool const swapEndian,
                               int const verboseLevel,
                               size_t numberTotalParticles[]);

// function that reads the data in a single gadget file - this function is called by the 'readGadgetFile_single' and 'readGadgetFile_multiple' to do the actual data reading
void readGadgetData(std::string fileName,
                    Read_data<Real> *readData,
                    User_options &userOptions,
                    int const gadgetFileType,
                    int const noBytesPos,
                    int const noBytesVel,
                    bool const swapEndian,
                    size_t *numberParticlesRead);




/*! This functions reads the gadget header for one of the input files and uses that to set all the properties need for reading the input data:
        - finding the number of files
        - finding the header type and endianess of the data
        - setting the box length in the header as dimensions of the box
        - finding the number of particles to be read from the input files
        - reserving memory for reading the data
 */
template <typename T>
void initializeGadget(std::string filename,
                      Read_data<T> *readData,
                      User_options *userOptions,
                      Gadget_header *gadgetHeader,
                      int *gadgetFileType,
                      int *noBytesPos,
                      int *noBytesVel,
                      bool *swapEndian,
                      size_t *noParticles)
{
    MESSAGE::Message message( userOptions->verboseLevel );
    std::string fileName = filename;
    bool singleFile = true;
    if ( not bfs::exists(fileName) ) //if this is true, than the input is in several files
    {
        fileName = gadgetHeader->filename( filename, 0 );
        singleFile = false;
    }


    // open the first binary file for reading and read some of the overall characteristics
    std::fstream inputFile;
    openInputBinaryFile( inputFile, fileName );


    // detect the Gadget file type -> gadget file type 1 or 2
    int buffer1, buffer2, buffer3, buffer4;       // variables to read the buffer before and after each gadget data block
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    bool validFile = gadgetHeader->detectSnapshotType( buffer1, gadgetFileType, swapEndian );
    if ( not validFile )
        throwError( "Unknown file type for the input Gadget snapshot. Tried Gadget snapshots type 1 and 2 as well as changing endianness, but none worked. Check that you inserted the correct input file." );
    if ( *swapEndian )
        message << "Detected that the input data file has a different endianness than the current system. The program will automatically change endianness for the data!" << MESSAGE::Flush;
    int offset = (*gadgetFileType)==2 ? 16 : 0;      // keep track if file type 2 to have a 16 bytes offset every time reading a new data block


    // now read the actual values of the gadget header
    inputFile.seekg( offset, std::ios::beg );
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    inputFile.read( reinterpret_cast<char *>(gadgetHeader), sizeof(*gadgetHeader) );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    SWAP_HEADER_ENDIANNESS( *swapEndian, buffer1, buffer2, (*gadgetHeader) ); //swap endianness if that is the case  
    
    // get the type (float/double) used to store position and velocity data
    inputFile.read( reinterpret_cast<char *>(&buffer3), sizeof(buffer3) );
    inputFile.seekg( buffer3, std::ios::cur );
    inputFile.read( reinterpret_cast<char *>(&buffer4), sizeof(buffer4) );
    inputFile.read( reinterpret_cast<char *>(&buffer4), sizeof(buffer4) );
    inputFile.close();
    
    SWAP_HEADER_ENDIANNESS( *swapEndian, buffer1, buffer2, (*gadgetHeader) ); //swap endianness if that is the case
    if ( buffer1!=buffer2 or buffer1!=256 )
        throwError( "The was an error while reading the header of the GADGET snapshot file. The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );
    
    // compute the number of bytes used to save each real value
    int thisNoParts = 0;
    for (int i=0; i<6; ++i)
        thisNoParts += gadgetHeader->npart[i];
    *noBytesPos = buffer3 / (3*thisNoParts);
    *noBytesVel = buffer4 / (3*thisNoParts);
    std::cout << "!!!! Number of bytes for position and velocity data: " << *noBytesPos << " and  " << *noBytesVel << ", respectively.\n" << std::flush;


    // check if to set the box coordinates from the information in the header
    if ( not userOptions->userGivenBoxCoordinates )
    {
        for (size_t i=0; i<NO_DIM; ++i)
        {
            userOptions->boxCoordinates[2*i] = 0.;                    // this is the left extension of the full box
            userOptions->boxCoordinates[2*i+1] = gadgetHeader->BoxSize;// right extension of the full box
        }
    }
    else
        message << "The box coordinates were set by the user using the program options. The program will keep this values and will NOT use the box length information from the Gadget file!" << MESSAGE::Flush;
#ifdef WOJTEK
    if ( userOptions->additionalOptions.size()!=0 ) //if inserted an option
        gadgetHeader->num_files = atoi( userOptions->additionalOptions[0].c_str() );
    gadgetHeader->print();
#endif


    // get the total number of particles to be read from the file/files
    size_t numberTotalParticles[6];
    if ( singleFile )    // if single file
    {
        gadgetHeader->num_files = 1;
        for (int i=0; i<6; ++i)
            numberTotalParticles[i] = gadgetHeader->npart[i];
    }
    else                  // for multiple files read the number of particles in each file and keep track of that
        countGadgetParticleNumber( filename, gadgetHeader->num_files, *gadgetFileType, *swapEndian, userOptions->verboseLevel, numberTotalParticles );

    // from that keep only the particles specified by the user (not all the particle species may be of interest)
    *noParticles = 0;   //total number of particles to be read from file
    for (int i=0; i<6; ++i)
    {
        if ( not userOptions->readParticleSpecies[i] )
            numberTotalParticles[i] = 0;
        *noParticles += numberTotalParticles[i];
    }
    message << "Reading " << *noParticles << " particle data from the input file. These particles are made from the particle species: " << numberTotalParticles[0] << " + "  << numberTotalParticles[1] << " + "  << numberTotalParticles[2] << " + "  << numberTotalParticles[3] << " + "  << numberTotalParticles[4] << " + "  << numberTotalParticles[5] << " .\n" << MESSAGE::Flush;



    // allocate memory for the particle data
    message << "Allocating memory for: positions... " << MESSAGE::Flush;
    if ( userOptions->readParticleData[0] )
        readData->position( *noParticles );  // particle positions
    message << "weights... " << MESSAGE::Flush;
    if ( userOptions->readParticleData[1] )
        readData->weight( *noParticles );    // particle weights (weights = particle masses from the snapshot file)
#ifdef VELOCITY
    message << "velocity... " << MESSAGE::Flush;
    if ( userOptions->readParticleData[2] )
        readData->velocity( *noParticles );  // particle velocities
#endif
#ifdef SCALAR
    message << "scalars... " << MESSAGE::Flush;
    int noScalars = 0;
    for (size_t i=3; i<userOptions->readParticleData.size(); ++i)
        if ( userOptions->readParticleData[i] )
            noScalars += 1;
    if ( noScalars>0 )
        readData->scalar( *noParticles );  // particle scalar quantity
#endif
    message << "Done.\n";

    // check that the 'readParticleData' and 'readParticleSpecies' options make sense
    if ( not userOptions->readParticleData[0] )
        throwError( "The program needs the particle position information to be able to interpolate the fields on a grid. Please add '1' to the integer number giving the data blocks to be read from the input Gadget snapshot." );
    if ( not userOptions->readParticleData[1] )
    {
        MESSAGE::Warning warning( userOptions->verboseLevel );
        warning << "You selected not to read the Gadget particle masses. This means that the program will treat all particles as having the same weight (mass)." << MESSAGE::EndWarning;
    }
    if ( (*noParticles)<=0 )
        throwError( "Please select again the particle species that you would like to read from the Gadget file. There are no particles in the current selection!" );
}



/*! This function reads a Gadget snapshots saved in a single or multiple files.

NOTE: This function is a not the best choice to understand how to read input data for users unfamiliar with the Gadget snapshot files. */
void readGadgetFile(std::string filename,
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


    // read the data in each file
    size_t numberParticlesRead = 0;   // the total number of particles read after each file

    // iterate over all the files and read the data
    for (int i=0; i<gadgetHeader.num_files; ++i )
    {
        fileName = gadgetHeader.filename( filename, i );
        message << "Reading GADGET snapshot file '" << fileName << "' which is file " << i+1 << " of " << gadgetHeader.num_files << " files...\n" << MESSAGE::Flush;

        // call the function that reads the data
        readGadgetData( fileName, readData, *userOptions, gadgetFileType, noBytesPos, noBytesVel, swapEndian, &numberParticlesRead );
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
    }

    return;
}




/*! Function that counts how many particles are in a Gadget file snapshot that was saved in multiple files.
*/
void countGadgetParticleNumber(std::string filenameRoot,
                               int const noFiles,
                               int const gadgetFileType,
                               bool const swapEndian,
                               int const verboseLevel,
                               size_t numberTotalParticles[])
{
    for (int i=0; i<6; ++i)
        numberTotalParticles[i] = 0;
    int offset = gadgetFileType==2 ? 16 : 0;

    // loop over all the files and count the number of particles
    for (int i=0; i<noFiles; ++i)
    {
        Gadget_header gadgetHeader;
        std::string fileName = gadgetHeader.filename( filenameRoot, i ); //get the filename for file i

        std::fstream inputFile;
        openInputBinaryFile( inputFile, fileName );

        // read the header
        int buffer1, buffer2;
        inputFile.seekg( offset, std::ios::beg );
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        inputFile.read( reinterpret_cast<char *>(&gadgetHeader), sizeof(gadgetHeader) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        inputFile.close();
        SWAP_HEADER_ENDIANNESS( swapEndian, buffer1, buffer2, gadgetHeader ); //swap endianness if that is the case
        if ( buffer1!=buffer2 or buffer1!=256 )
            throwError( "The was an error while reading the header of the GADGET snapshot file '" + fileName + "'. The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );

        // add the particle numbers
        for (int j=0; j<6; ++j)
            numberTotalParticles[j] += gadgetHeader.npart[j];
    }

    MESSAGE::Message message( verboseLevel );
    message << "The data is in " << noFiles << " files and contains the following number of particles: " << numberTotalParticles[0] << " + "  << numberTotalParticles[1] << " + "  << numberTotalParticles[2] << " + "  << numberTotalParticles[3] << " + "  << numberTotalParticles[4] << " + "  << numberTotalParticles[5] << " .\n" << MESSAGE::Flush;
}



/*! This function reads the gadget data from a single file. This function should be called from within a loop to read a gadget snapshot saved in multiple files.
*/
void readGadgetData(std::string fileName,
                    Read_data<Real> *readData,
                    User_options &userOptions,
                    int const gadgetFileType,
                    int const noBytesPos,
                    int const noBytesVel,
                    bool const swapEndian,
                    size_t *numberParticlesRead)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    int offset = gadgetFileType==2 ? 16 : 0;


    // open the file
    std::fstream inputFile;
    openInputBinaryFile( inputFile, fileName );


    // read the header
    int buffer1, buffer2;
    Gadget_header tempHeader;
    READ_DELIMETER;
    inputFile.read( reinterpret_cast<char *>(&tempHeader), sizeof(tempHeader) );
    DELIMETER_CONSISTANCY_CHECK("header");
    if ( buffer1!=256 )
        throwError( "The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );


    // read the position block
    READ_DELIMETER;
    if ( userOptions.readParticleData[0] )
    {
        Real *positions = readData->position();                 // returns a pointer to the particle positions array
        size_t dataOffset = (*numberParticlesRead) * NO_DIM;    // the offset in the positions array from where to start reading the new positions
        message << "\t reading positions of the particles... " << MESSAGE::Flush;

        // loop over each species and reads its data if the user requested so
        for (int i=0; i<6; ++i)
        {
            // check if to skip this data block
            if ( tempHeader.npart[i]==0 )       // skip since no particle of this species is present
                continue;
            else if ( not userOptions.readParticleSpecies[i] )    // skip this data block since don't need this data
            {
                size_t skipBytes = tempHeader.npart[i] * noBytesPos * NO_DIM;
                inputFile.seekg( skipBytes, std::ios::cur );
                continue;
            }

            // read this data block
            size_t readBytes = tempHeader.npart[i] * noBytesPos * NO_DIM;
            if ( sizeof(Real) == noBytesPos )
                inputFile.read( reinterpret_cast<char *>( &(positions[dataOffset]) ), readBytes );
            else if ( noBytesPos==4 )    // read float data
            {
                float *temp = new float[ tempHeader.npart[i] * NO_DIM ];
                inputFile.read( reinterpret_cast<char *>( temp ), readBytes );
                for (int u=0; u<tempHeader.npart[i] * NO_DIM; ++u)  positions[ dataOffset + u ] = temp[u];
                delete[] temp;
            }
            else if ( noBytesPos==8 )    // read double data
            {
                double *temp = new double[ tempHeader.npart[i] * NO_DIM ];
                inputFile.read( reinterpret_cast<char *>( temp ), readBytes );
                for (int u=0; u<tempHeader.npart[i] * NO_DIM; ++u)  positions[ dataOffset + u ] = temp[u];
                delete[] temp;
            }
            dataOffset += tempHeader.npart[i] * NO_DIM;
        }

        // check for consistency in the input file
        message << "Done.";
    }
    else    //skip the positions if not interested in them
    {
        message << "\n\t (skipping positions)" << MESSAGE::Flush;
        size_t skipBytes = 0;
        for (int i=0; i<6; ++i) skipBytes += tempHeader.npart[i] * noBytesPos * NO_DIM;
        inputFile.seekg( skipBytes, std::ios::cur );
    }
    DELIMETER_CONSISTANCY_CHECK("position");


    // read the velocities block
    READ_DELIMETER;
#ifdef VELOCITY
    if ( userOptions.readParticleData[2] )
    {
        Real *velocities = readData->velocity();              // returns a pointer to the particle velocity array
        size_t dataOffset = (*numberParticlesRead) * NO_DIM;    // the offset in the velocities array from where to start reading the new velocities
        message << "\n\t reading velocities of the particles... " << MESSAGE::Flush;

        // loop over each species and reads its data if the user requested so
        for (int i=0; i<6; ++i)
        {
            // check if to skip this data block
            if ( tempHeader.npart[i]==0 )       // skip since no particle of this species is present
                continue;
            else if ( not userOptions.readParticleSpecies[i] )    // skip this data block since don't need this data
            {
                size_t skipBytes = tempHeader.npart[i] * noBytesVel * NO_DIM;
                inputFile.seekg( skipBytes, std::ios::cur );
                continue;
            }

            // read this data block
            size_t readBytes = tempHeader.npart[i] * noBytesVel * NO_DIM;
            if ( sizeof(Real) == noBytesVel )
                inputFile.read( reinterpret_cast<char *>( &(velocities[dataOffset]) ), readBytes );
            else if ( noBytesVel==4 )    // read float data
            {
                float *temp = new float[ tempHeader.npart[i] * NO_DIM ];
                inputFile.read( reinterpret_cast<char *>( temp ), readBytes );
                for (int u=0; u<tempHeader.npart[i] * NO_DIM; ++u)  velocities[ dataOffset + u ] = temp[u];
                delete[] temp;
            }
            else if ( noBytesVel==8 )    // read double data
            {
                double *temp = new double[ tempHeader.npart[i] * NO_DIM ];
                inputFile.read( reinterpret_cast<char *>( temp ), readBytes );
                for (int u=0; u<tempHeader.npart[i] * NO_DIM; ++u)  velocities[ dataOffset + u ] = temp[u];
                delete[] temp;
            }
            dataOffset += tempHeader.npart[i] * NO_DIM;
        }
        message << "Done.";
    }
    else
    {   //skip the velocities if not interested in them
        message << "\n\t (skipping velocities)" << MESSAGE::Flush;
        size_t skipBytes = 0;
        for (int i=0; i<6; ++i) skipBytes += tempHeader.npart[i] * noBytesVel * NO_DIM;
        inputFile.seekg( skipBytes, std::ios::cur );
    }
#else
    message << "\n\t (skipping velocities)" << MESSAGE::Flush;
    size_t skipBytes = 0;
    for (int i=0; i<6; ++i)
        skipBytes += tempHeader.npart[i] * noBytesVel * NO_DIM;
    inputFile.seekg( skipBytes, std::ios::cur );
#endif
    DELIMETER_CONSISTANCY_CHECK("velocity");


    // skip the particle ID block
    message << "\n\t (skipping ids)" << MESSAGE::Flush;
    READ_DELIMETER;
    inputFile.seekg( buffer1, std::ios::cur );
    DELIMETER_CONSISTANCY_CHECK("id");


    // read the masses (or weights) if different
    bool massBlockPresent = false;
    for (int i=0; i<6; ++i)
        if ( tempHeader.mass[i]==0. and tempHeader.npart[i]!=0 )
            massBlockPresent = true;
    
    if ( massBlockPresent )
    {
        READ_DELIMETER;
    }
    if ( userOptions.readParticleData[1] )
    {
        Real *weights = readData->weight();          // returns a pointer to the particle weights array
        size_t dataOffset = (*numberParticlesRead);    // the offset in the masses array from where to start reading the new masses
        message << "\n\t reading masses of the particles... " << MESSAGE::Flush;

        // loop over each species and reads its data if the user requested so
        for (int i=0; i<6; ++i)
        {
            // check if to skip this data block
            if ( tempHeader.npart[i]==0 )       // skip since no particle of this species is present
                continue;
            else if ( tempHeader.mass[i]==0. and not userOptions.readParticleSpecies[i] )    // skip this data block if it exists since don't need this data
            {
                size_t skipBytes = tempHeader.npart[i] * sizeof(float);
                inputFile.seekg( skipBytes, std::ios::cur );
                continue;
            }

            // read the masses
            if ( tempHeader.mass[i]==0. )       // each particle has a different mass
            {
                size_t readBytes = tempHeader.npart[i] * noBytesPos;
                if ( sizeof(Real) == noBytesPos )
                    inputFile.read( reinterpret_cast<char *>( &(weights[dataOffset]) ), readBytes );
                else if ( noBytesPos==4 )    // read float data
                {
                    float *temp = new float[ tempHeader.npart[i] ];
                    inputFile.read( reinterpret_cast<char *>( temp ), readBytes );
                    for (int u=0; u<tempHeader.npart[i]; ++u)  weights[ dataOffset + u ] = temp[u];
                    delete[] temp;
                }
                else if ( noBytesPos==8 )    // read double data
                {
                    double *temp = new double[ tempHeader.npart[i] ];
                    inputFile.read( reinterpret_cast<char *>( temp ), readBytes );
                    for (int u=0; u<tempHeader.npart[i]; ++u)  weights[ dataOffset + u ] = temp[u];
                    delete[] temp;
                }
            }
            else if ( userOptions.readParticleSpecies[i] ) // all particles have the same mass
            {
                float mass = tempHeader.mass[i];
                if ( swapEndian ) BYTESWAP(mass);       // keep the same endianess for all the data
                for (size_t j=dataOffset; j<dataOffset+tempHeader.npart[i]; ++j)
                    weights[j] = mass;
            }

            dataOffset += tempHeader.npart[i];
        }
        message << "Done.\n";
    }
    if ( massBlockPresent )
    {
        DELIMETER_CONSISTANCY_CHECK("mass");
    }


    // read the internal energy
    size_t noScalarsRead = 0;
    if ( userOptions.readParticleData[3] and tempHeader.npart[0]>0) 
    { // only makes sense for gas particles
        READ_DELIMETER;
        Real *scalar = readData->scalar();          // returns a pointer to the particle scalar properties array
        Real *weights = readData->weight();         // returns a pointer to the particle weights array
        size_t dataOffset = (*numberParticlesRead);  // the offset in the array to store current files' particles' U
        message << "\t reading internal energy of the particles... " << MESSAGE::Flush;

        // read the masses
        size_t readBytes = tempHeader.npart[0] * sizeof(float);
        float *tempData = new float[ tempHeader.npart[0] ];
        inputFile.read( reinterpret_cast<char *>( tempData ), readBytes );

        float mean = 0;
        for (size_t i; i<size_t(tempHeader.npart[0]); ++i)
        {
            size_t index1 = dataOffset + i;
            size_t index2 = index1 * NO_SCALARS + noScalarsRead;
            scalar[index2] = weights[index1] * tempData[i]; // U is given per unit mass in Gadget
            mean += scalar[index2];
        }
        mean /= tempHeader.npart[0];
        // or read file on the fly (more memory efficient but slower):
        //for (size_t i=dataOffset; i<size_t(tempHeader.npart[0])+dataOffset; ++i ) {
        //    float tempData;
        //    inputFile.read( reinterpret_cast<char *>( &tempData ), sizeof(tempData) );
        //    scalar[i*NO_SCALARS+noScalarsRead] = weights[i]*tempData;
        //}

        delete[] tempData;
        noScalarsRead += 1;
        message << "\r\t reading internal energy of the particles (mean energy: " << mean << ")... Done.\n";
        DELIMETER_CONSISTANCY_CHECK("internal energy U");
    }
    //~ else
    //~ {
        //~ message << "\n\t (skipping internal energy U)" << MESSAGE::Flush;
        //~ size_t skipBytes = 0;
        //~ for (int i=0; i<6; ++i) skipBytes += tempHeader.npart[i] * sizeof(float) * NO_DIM;
        //~ inputFile.seekg( skipBytes, std::ios::cur );
        //~ DELIMETER_CONSISTANCY_CHECK("internal energy U");
    //~ }

    inputFile.close();
    message << "\n";
    for (int i=0; i<6; ++i)     // update the number of read particles
        if ( userOptions.readParticleSpecies[i] )
            (*numberParticlesRead) += tempHeader.npart[i];
}








