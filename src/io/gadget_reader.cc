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



/* This file defines a reader for binary Gadget snapshot files and also for HDF5 files (need to have the HDF5 library for that). The reader given here can read all the particles and properties stored in the Gadget snapshot file by giving the appropiate options to the program. See the documentation for the options that control this behavior. */



#define SWAP_HEADER_ENDIANNESS(x1,x2,x3,x4) { if( x1 ) {BYTESWAP( x2 ); BYTESWAP( x3 ); x4.swapBytes();} }
#define SWAP_ENDIANNESS(x1,x2,x3)           { if( x1 ) {BYTESWAP( x2 ); BYTESWAP( x3 );} }


// Header structure for reading Gadget snapshots
struct Gadget_header
{
    int      npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[6];
    int      flag_cooling;
    int      num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
    
    
    // return the file name for a Gadget snapshot saved in single or multiple files - note that the name must contain a '%i' or '%s' character
    std::string filename(std::string fileRoot, int const fileNumber, bool checkFileExists=true )
    {
        char buf[500];
        sprintf( buf, fileRoot.c_str(), fileNumber );
        std::string fileName( buf );
        if ( not bfs::exists(fileName) and checkFileExists )
            throwError( "The program could not open the input GADGET snapshot file/files: '" + fileName + "'. It cannot find the file/files." );
        return fileName;
    }
    
    // Function that prints the Gadget header.
    void print()
    {
        std::cout << "\nThe header of the Gadget file contains the following info:\n" 
            << "npart[6]     =  " << npart[0] << "  " << npart[1] << "  " << npart[2] << "  " << npart[3] << "  " <<  npart[4] << "  " <<  npart[5] << "\n"
            << "mass[6]      =  " << mass[0] << "  " << mass[1] << "  " << mass[2] << "  " << mass[3] << "  " << mass[4] << "  " << mass[5] << "\n"
            << "time         =  " << time << "\n"
            << "redshift     =  " << redshift << "\n"
            << "flag_sfr     =  " << flag_sfr << "\n"
            << "flag_feedback=  " << flag_feedback << "\n"
            << "npartTotal[6]=  " << npartTotal[0] << "  " << npartTotal[1] << "  " << npartTotal[2] << "  " << npartTotal[3] << "  " << npartTotal[4] << "  " << npartTotal[5] << "  " << "\n"
            << "flag_cooling =  " << flag_cooling << "\n"
            << "num_files    =  " << num_files << "\n"
            << "BoxSize      =  " << BoxSize << "\n"
            << "Omega0       =  " << Omega0 << "\n"
            << "OmegaLambda  =  " << OmegaLambda << "\n"
            << "h            =  " << HubbleParam << "\n\n";
    }
    
    // Swap endianness
    void swapBytes()
    {
        ByteSwapArray( npart, 6 );
        ByteSwapArray( mass, 6 );
        BYTESWAP( time );
        BYTESWAP( redshift );
        BYTESWAP( flag_sfr );
        BYTESWAP( flag_feedback );
        ByteSwapArray( npartTotal, 6 );
        BYTESWAP( flag_cooling );
        BYTESWAP( num_files );
        BYTESWAP( BoxSize );
        BYTESWAP( Omega0 );
        BYTESWAP( OmegaLambda );
        BYTESWAP( HubbleParam );
    }
    
    // Checks for the type of the Gadget file -> can detected Gadget file type 1 & 2. Returns true if it could identify the gadget file type.
    bool detectSnapshotType(int const bufferValue,
                            int *gadgetFileType,
                            bool *swapEndian)
    {
        int buffer1 = bufferValue;
        *swapEndian = false;
        
        if ( buffer1 == 8 )             // gadget file format 2
            *gadgetFileType = 2;
        else if ( buffer1 == 256 )      // gadget file format 1
            *gadgetFileType = 1;
        else                            // check for swapped endianness
        {
            BYTESWAP( buffer1 );
            *swapEndian = true;
            if ( buffer1 == 8 )	        // gadget file format 2
                *gadgetFileType = 2;
            else if ( buffer1 == 256 )  // gadget file format 1
                *gadgetFileType = 1;
            else                        // could not detect the file type
                return false;
        }
        return true;
    }
};




// function that counts the number of particles for snapshots saved as multiple files
void countGadgetParticleNumber(std::string filenameRoot,
                               int const noFiles,
                               int const gadgetFileType,
                               bool const swapEndian,
                               int const verboseLevel,
                               size_t numberTotalParticles[]);

// function that reads the data in a single gadget file - this function is called by the 'readGadgetFile_single' and 'readGadgetFile_multiple' to do the actual data reading
void readGadgetData(std::string fileName,
                    Read_data<float> *readData,
                    User_options &userOptions,
                    int const gadgetFileType,
                    bool const swapEndian,
                    size_t *numberParticlesRead);




/*! This functions reads the gadget header for one of the input files and uses that to set all the properties need for reading the input data:
        - finding the number of files
        - finding the header type and endianess of the data
        - setting the box length in the header as dimensions of the box
        - finding the number of particles to be read from the input files
        - reserving memory for reading the data
 */
void initializeGadget(std::string filename,
                      Read_data<float> *readData,
                      User_options *userOptions,
                      Gadget_header *gadgetHeader,
                      int *gadgetFileType,
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
    int buffer1, buffer2;       // variables to read the buffer before and after each gadget data block
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
    inputFile.close();
    SWAP_HEADER_ENDIANNESS( *swapEndian, buffer1, buffer2, (*gadgetHeader) ); //swap endianness if that is the case
    if ( buffer1!=buffer2 or buffer1!=256 )
        throwError( "The was an error while reading the header of the GADGET snapshot file. The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );
    
    
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
    if ( userOptions->readParticleData[0] )
        readData->position( *noParticles );  // particle positions
    if ( userOptions->readParticleData[1] )
        readData->weight( *noParticles );    // particle weights (weights = particle masses from the snapshot file)
#ifdef VELOCITY
    if ( userOptions->readParticleData[2] )
        readData->velocity( *noParticles );  // particle velocities
#endif
    
    
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
                    Read_data<float> *readData,
                    User_options *userOptions)
{
    int gadgetFileType;     // stores the type of the gadget file (1 or 2)
    bool swapEndian;        // true if need to swap the endianess of the data
    size_t noParticles;     // total number of particles to be read from the file
    Gadget_header gadgetHeader; // the gadget header for the first file
    
    
    // get the input gadget data characteristics: file type, endianess, particle number and number of files. The following function also reserves memory for the positions, weights and velocity data (if requested so).
    initializeGadget( filename, readData, userOptions, &gadgetHeader, &gadgetFileType, &swapEndian, &noParticles );
    
    
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
        readGadgetData( fileName, readData, *userOptions, gadgetFileType, swapEndian, &numberParticlesRead );
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
                    Read_data<float> *readData,
                    User_options &userOptions,
                    int const gadgetFileType,
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
    inputFile.seekg( offset, std::ios::cur );      // jump the block id data if type 2 Gadget file
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    inputFile.read( reinterpret_cast<char *>(&tempHeader), sizeof(tempHeader) );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    SWAP_HEADER_ENDIANNESS( swapEndian, buffer1, buffer2, tempHeader );  // swap endianness if that is the case
    if ( buffer1!=buffer2 or buffer1!=256 )
        throwError( "The was an error while reading the header of the GADGET file '" + fileName + "'. The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );
    
    
    // read the position block
    if ( userOptions.readParticleData[0] )
    {
        float *positions = readData->position();               // returns a pointer to the particle positions array
        size_t dataOffset = (*numberParticlesRead) * NO_DIM;    // the offset in the positions array from where to start reading the new positions
        message << "\t reading positions of the particles ... " << MESSAGE::Flush;
        inputFile.seekg( offset, std::ios::cur );      // jump the block id data if type 2 Gadget file
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        
        // loop over each species and reads its data if the user requested so
        for (int i=0; i<6; ++i)
        {
            // check if to skip this data block
            if ( tempHeader.npart[i]==0 )       // skip since no particle of this species is present
                continue;
            else if ( not userOptions.readParticleSpecies[i] )    // skip this data block since don't need this data
            {
                size_t skipBytes = tempHeader.npart[i] * sizeof(float) * NO_DIM;
                inputFile.seekg( skipBytes, std::ios::cur );
                continue;
            }
            
            // read this data block
            size_t readBytes = tempHeader.npart[i] * sizeof(float) * NO_DIM;
            inputFile.read( reinterpret_cast<char *>( &(positions[dataOffset]) ), readBytes );
            dataOffset += tempHeader.npart[i] * NO_DIM;
        }
        
        // check for consistency in the input file
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        SWAP_ENDIANNESS( swapEndian, buffer1, buffer2 );
        if ( buffer1!=buffer2 )
            throwError( "The integers before and after the particle position data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." );
        message << "Done\n";
    }
    else    //skip the positions if not interested in them
    {
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        SWAP_ENDIANNESS( swapEndian, buffer1, buffer2 );
        inputFile.seekg( buffer1+sizeof(int), std::ios::cur );
    }
    
    
    // read the velocities block
#ifdef VELOCITY
    if ( userOptions.readParticleData[2] )
    {
        float *velocities = readData->velocity();              // returns a pointer to the particle velocity array
        size_t dataOffset = (*numberParticlesRead) * NO_DIM;    // the offset in the velocities array from where to start reading the new velocities
        message << "\t reading velocities of the particles ... " << MESSAGE::Flush;
        inputFile.seekg( offset, std::ios::cur );      // jump the block id data if type 2 Gadget file
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        
        // loop over each species and reads its data if the user requested so
        for (int i=0; i<6; ++i)
        {
            // check if to skip this data block
            if ( tempHeader.npart[i]==0 )       // skip since no particle of this species is present
                continue;
            else if ( not userOptions.readParticleSpecies[i] )    // skip this data block since don't need this data
            {
                size_t skipBytes = tempHeader.npart[i] * sizeof(float) * NO_DIM;
                inputFile.seekg( skipBytes, std::ios::cur );
                continue;
            }
            
            // read this data block
            size_t readBytes = tempHeader.npart[i] * sizeof(float) * NO_DIM;
            inputFile.read( reinterpret_cast<char *>( &(velocities[dataOffset]) ), readBytes );
            dataOffset += tempHeader.npart[i] * NO_DIM;
        }
        
        // check for consistency in the input file
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        SWAP_ENDIANNESS( swapEndian, buffer1, buffer2 );
        if ( buffer1!=buffer2 )
            throwError( "The integers before and after the particle velocity data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." );
        message << "Done\n";
    }
    else    //skip the velocities if not interested in them
    {
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        SWAP_ENDIANNESS( swapEndian, buffer1, buffer2 );
        inputFile.seekg( buffer1+sizeof(int), std::ios::cur );
    }
#else
    // skip the velocity data block if code not compiled with enabled velocities
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    SWAP_ENDIANNESS( swapEndian, buffer1, buffer2 );
    inputFile.seekg( buffer1+sizeof(int), std::ios::cur );
#endif
    
    
    // skip the particle ID block
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    SWAP_ENDIANNESS( swapEndian, buffer1, buffer2 );
    inputFile.seekg( buffer1+sizeof(int), std::ios::cur );
    
    
    // read the masses (or weights) if different
    if ( userOptions.readParticleData[1] )
    {
        float *weights = readData->weight();          // returns a pointer to the particle weights array
        size_t dataOffset = (*numberParticlesRead);    // the offset in the masses array from where to start reading the new masses
        message << "\t reading masses of the particles ... " << MESSAGE::Flush;
        
        // find if there is a mass block in the file
        bool massBlockPresent = false;
        for (int i=0; i<6; ++i)
            if ( tempHeader.mass[i]==0. and tempHeader.npart[i]!=0 )
                massBlockPresent = true;
        
        if ( massBlockPresent )
        {
            inputFile.seekg( offset, std::ios::cur );      // jump the block id data if type 2 Gadget file
            inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        }
        
        
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
                size_t readBytes = tempHeader.npart[i] * sizeof(float);
                inputFile.read( reinterpret_cast<char *>( &(weights[dataOffset]) ), readBytes );
            }
            else                                // all particles have the same mass
            {
                float mass = tempHeader.mass[i];
                if ( swapEndian ) BYTESWAP(mass);       // keep the same endianess for all the data
                for (size_t j=dataOffset; j<dataOffset+tempHeader.npart[i]; ++j)
                    weights[j] = mass;
            }
                    
            dataOffset += tempHeader.npart[i];
        }
        message << "Done\n";
    }
    
    
    inputFile.close();
    for (int i=0; i<6; ++i)     // update the number of read particles
        if ( userOptions.readParticleSpecies[i] )
            (*numberParticlesRead) += tempHeader.npart[i];
}






// The following are the functions used to read an HDF5 gadget file
#ifdef HDF5
#include <H5Cpp.h>
using namespace H5;


/*! Function to check if attribute exists. */
extern "C"
{
    bool doesAttributeExist(hid_t obj_id, const char* name)
    {
        return( H5Aexists( obj_id, name ) > 0 ? true : false );
    }
}




/*! Reads some of the entries for the Gadget header from an HDF5 file.
NOTE: it does not read all the entries, it only reads the particle number in the file, the mass array, the box length and the number of files per snapshot. 
*/
void HDF5_readGadgetHeader(std::string filename,
                           Gadget_header *gadgetHeader)
{
    // the name of the HDF5 file
    const H5std_string FILE_NAME( filename );
    
    // open the HDF5 file and the header group
    H5File *file = new H5File( FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT );
    Group *group = new Group( file->openGroup("/Header") );
    
    
    // start reading one header attribute at a time
    std::string name( "NumPart_ThisFile" );
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_INT, gadgetHeader->npart );
    else throwError( "No '" + name + "' attribute found in the HDF5 file '" + filename + "'. Cannot continue with the program!" );
    
    name = "MassTable";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, gadgetHeader->mass );
    else throwError( "No '" + name + "' attribute found in the HDF5 file '" + filename + "'. Cannot continue with the program!" );
    
    name = "NumFilesPerSnapshot";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_INT, &(gadgetHeader->num_files) );
    else throwError( "No '" + name + "' attribute found in the HDF5 file '" + filename + "'. Cannot continue with the program!" );
    
    name = "BoxSize";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, &(gadgetHeader->BoxSize) );
    
    name = "NumPart_Total";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_UINT, gadgetHeader->npartTotal );
    
    name = "Redshift";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, &(gadgetHeader->redshift) );
    
    name = "Time";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, &(gadgetHeader->time) );
    name = "Time_GYR";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, &(gadgetHeader->time) );
    
    name = "Omega0";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, &(gadgetHeader->Omega0) );
    
    name = "OmegaLambda";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, &(gadgetHeader->OmegaLambda) );
    
    name = "HubbleParam";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, &(gadgetHeader->HubbleParam) );
    
    
    // close the group and file
    delete group;
    delete file;
}



/*! Reads the Gadget data from a single HDF5 file. This function should be called to do the actual data reading for snapshots saved in single or multiple files.
*/
void HDF5_readGadgetData(std::string filename,
                         Read_data<float> *readData,
                         User_options &userOptions,
                         int const fileIndex,
                         size_t *numberParticlesRead)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    
    // get the header associated to this file
    Gadget_header gadgetHeader;
    HDF5_readGadgetHeader( filename, &gadgetHeader );
    
    //select only the species of interest
    for (int i=0; i<6; ++i)
        if ( not userOptions.readParticleSpecies[i] )
            gadgetHeader.npart[i] = 0;
    
    
    // open the HDF5 file
    const H5std_string FILE_NAME( filename );
    H5File *file = new H5File( FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT );
    Group *group;
    
    
    // read the coordinates
    if ( userOptions.readParticleData[0] )
    {
        float *positions = readData->position();               // returns a pointer to the particle positions array
        size_t dataOffset = (*numberParticlesRead) * NO_DIM;   // the offset in the positions array from where to start reading the new positions
        message << "\t reading the particles positions ... " << MESSAGE::Flush;
        for(int type=0; type<6; type++)                        // loop over the particle species
        {
            if ( gadgetHeader.npart[type]<=0 ) continue;
            char buf[500];
            sprintf( buf, "/PartType%d", type );
            group = new Group( file->openGroup(buf) );
            
            // open the data set
            DataSet dataset = group->openDataSet("Coordinates");
            
            dataset.read( &(positions[dataOffset]), PredType::NATIVE_FLOAT );
            delete group;
            
            dataOffset += gadgetHeader.npart[type] * NO_DIM;
        }
        message << "Done\n";
    }
    
    
    // read the masses (or weights) if different
    if ( userOptions.readParticleData[1] )
    {
        float *weights = readData->weight();          // returns a pointer to the particle weights array
        size_t dataOffset = (*numberParticlesRead);    // the offset in the masses array from where to start reading the new masses
        message << "\t reading the particles masses ... " << MESSAGE::Flush;
        for(int type=0; type<6; type++)                        // loop over the particle species
        {
            if ( gadgetHeader.npart[type]<=0 ) continue;
            if ( gadgetHeader.mass[type]!=0. )          // particle masses given in header
            {
                float mass = gadgetHeader.mass[type];
                for (size_t j=dataOffset; j<dataOffset+gadgetHeader.npart[type]; ++j)
                    weights[j] = mass;
            }
            else                                        // particle masses are variable
            {
                char buf[500];
                sprintf( buf, "/PartType%d", type );
                group = new Group( file->openGroup(buf) );
                
                DataSet dataset = group->openDataSet("Mass");
                
                dataset.read( &(weights[dataOffset]), PredType::NATIVE_FLOAT );
                delete group;
            }
            dataOffset += gadgetHeader.npart[type];
        }
        message << "Done\n";
    }
    
    
    // read the velocities
    if ( userOptions.readParticleData[2] )
    {
        float *velocities = readData->velocity();              // returns a pointer to the particle velocity array
        size_t dataOffset = (*numberParticlesRead) * NO_DIM;   // the offset in the velocity array from where to start reading the new data
        message << "\t reading the particles velocities ... " << MESSAGE::Flush;
        for(int type=0; type<6; type++)                        // loop over the particle species
        {
            if ( gadgetHeader.npart[type]<=0 ) continue;
            char buf[500];
            sprintf( buf, "/PartType%d", type );
            group = new Group( file->openGroup(buf) );
            
            // open the data set
            DataSet dataset = group->openDataSet("Velocity");
            
            dataset.read( &(velocities[dataOffset]), PredType::NATIVE_FLOAT );
            delete group;
            
            dataOffset += gadgetHeader.npart[type] * NO_DIM;
        }
        message << "Done\n";
    }
    
    
    // read the HI content
    size_t noScalarsRead = 0;
    if ( userOptions.readParticleData[3] and gadgetHeader.npart[0]>0 )
    {
        // the element abundance makes sense only for gas particles 
        group = new Group( file->openGroup( "/PartType0/ElementAbundance" ) );
        DataSet dataset = group->openDataSet("Hydrogen");
        float *hydrogenMassFraction = new float[ gadgetHeader.npart[0] ];
        dataset.read( hydrogenMassFraction, PredType::NATIVE_FLOAT );
        delete group;
        
        
        // open the file that gives the HI fraction
        if ( userOptions.additionalOptions.empty() )
            throwError( "You need to give the name of the HDF5 files giving the Urchin HI fraction if you want to compute the HI mass distribution. This is inserted using the '--options' program option." );
        std::string h1FileName = gadgetHeader.filename( userOptions.additionalOptions[0], fileIndex );
        H5File *file2 = new H5File( h1FileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
        group = new Group( file2->openGroup( "/PartType0" ) );
        DataSet datasetMolecular = group->openDataSet("MolecularHydrogenMassFraction");
        DataSet datasetHydrogen = group->openDataSet("HydrogenOneFraction");
        
        float *molecularMassFraction = new float[ gadgetHeader.npart[0] ];
        datasetMolecular.read( molecularMassFraction, PredType::NATIVE_FLOAT );
        float *hydrogenOneFraction   = new float[ gadgetHeader.npart[0] ];
        datasetHydrogen.read( hydrogenOneFraction, PredType::NATIVE_FLOAT );
        delete group;
        delete file2;
        
        
        // get the HI mass
        float *scalar = readData->scalar();          // returns a pointer to the particle scalar properties array
        float *weights = readData->weight();         // returns a pointer to the particle weights array
        size_t dataOffset = (*numberParticlesRead);  // the offset in the scalar array from where to start reading the new values
        for (size_t i=0; i<size_t(gadgetHeader.npart[0]); ++i)
        {
            size_t index1 = dataOffset + i;
            size_t index2 = index1 * NO_SCALARS + noScalarsRead;
            scalar[index2] = weights[index1] * hydrogenMassFraction[i] * (float(1.)-molecularMassFraction[i]) * hydrogenOneFraction[i];
        }
        
        delete[] hydrogenMassFraction;
        delete[] molecularMassFraction;
        delete[] hydrogenOneFraction;
        noScalarsRead += 1;
    }
    
    
    // read the temperatures
    if ( userOptions.readParticleData[4] and gadgetHeader.npart[0]>0 )
    {
        // the tempartures make sense only for gas particles
        group = new Group( file->openGroup( "/PartType0" ) );
        DataSet dataset = group->openDataSet("Temperature");
        float *tempData = new float[ gadgetHeader.npart[0] ];
        dataset.read( tempData, PredType::NATIVE_FLOAT );
        delete group;
        
        // get the mass weighted temperature
        float *scalar = readData->scalar();          // returns a pointer to the particle scalar properties array
        float *weights = readData->weight();         // returns a pointer to the particle weights array
        size_t dataOffset = (*numberParticlesRead);  // the offset in the scalar array from where to start reading the new values
        for (size_t i=0; i<size_t(gadgetHeader.npart[0]); ++i)
        {
            size_t index1 = dataOffset + i;
            size_t index2 = index1 * NO_SCALARS + noScalarsRead;
            scalar[index2] = weights[index1] * tempData[i];
        }
        
        delete[] tempData;
        noScalarsRead += 1;
    }
    
    
    delete file;
    for (int i=0; i<6; ++i)     // update the number of read particles
        (*numberParticlesRead) += gadgetHeader.npart[i];
}



/*! Function that counts how many particles are in a Gadget file snapshot that was saved in multiple files. Works only for HDF5 files.
 */
void HDF5_countGadgetParticleNumber(std::string filenameRoot,
                                    int const noFiles,
                                    int const verboseLevel,
                                    size_t numberTotalParticles[])
{
    for (int i=0; i<6; ++i)
        numberTotalParticles[i] = 0;
    
    // loop over all the files and count the number of particles
    for (int i=0; i<noFiles; ++i)
    {
        Gadget_header gadgetHeader;
        std::string fileName = gadgetHeader.filename( filenameRoot, i );
        
        // read the header
        HDF5_readGadgetHeader( fileName, &gadgetHeader );
        
        // add the particle numbers
        for (int j=0; j<6; ++j)
            numberTotalParticles[j] += gadgetHeader.npart[j];
    }
    
    MESSAGE::Message message( verboseLevel );
    message << "The data is in " << noFiles << " files and contains the following number of particles: " << numberTotalParticles[0] << " + "  << numberTotalParticles[1] << " + "  << numberTotalParticles[2] << " + "  << numberTotalParticles[3] << " + "  << numberTotalParticles[4] << " + "  << numberTotalParticles[5] << " .\n" << MESSAGE::Flush;
}




/*! Read the gadget header for one snapshot and intialize some of the values in the userOptions class.
*/
void HDF5_initializeGadget(std::string filename,
                           Read_data<float> *readData,
                           User_options *userOptions,
                           Gadget_header *gadgetHeader,
                           size_t *noParticles)
{
    MESSAGE::Message message( userOptions->verboseLevel );
    std::string fileName = filename;
    bool singleFile = true;
    
    // check to see if there is only one input file or several
    if ( not bfs::exists(fileName) ) //if this is true, than the input is in several files
    {
        fileName = gadgetHeader->filename( filename, 0 );
        singleFile = false;
    }
    
    
    // now read the actual values of the gadget header
    HDF5_readGadgetHeader( fileName, gadgetHeader );
    
    
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
    else                                // for multiple files read the number of particles in each file and keep track of that
        HDF5_countGadgetParticleNumber( filename, gadgetHeader->num_files, userOptions->verboseLevel, numberTotalParticles );
    
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
    if ( userOptions->readParticleData[0] )
        readData->position( *noParticles );  // particle positions
    if ( userOptions->readParticleData[1] )
        readData->weight( *noParticles );    // particle weights (weights = particle masses from the snapshot file)
#ifdef VELOCITY
    if ( userOptions->readParticleData[2] )
        readData->velocity( *noParticles );  // particle velocities
#endif
    // allocate memory for scalar quantities
#ifdef SCALAR
    int noScalars = 0;
    for (size_t i=3; i<userOptions->readParticleData.size(); ++i)
        if ( userOptions->readParticleData[i] )
            noScalars += 1;
    if ( noScalars>0 )
        readData->scalar( *noParticles );  // particle scalar quantity
#endif
    
    
    
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
#ifdef SCALAR
    if ( noScalars>NO_SCALARS )
        throwError( "You asked for too many scalars to be read from the file using the input data option. Either ask for less scalar quantities or increase the number of scalar field components using the Makefile option '-DNO_SCALARS'. You asked for ", noScalars, " scalar components, but 'NO_SCALARS' is only ", NO_SCALARS, "." );
#endif
}



/*! The function that reads the data from single or multiple HDF5 files.
*/
void HDF5_readGadgetFile(std::string filename,
                         Read_data<float> *readData,
                         User_options *userOptions)
{
    MESSAGE::Message message( userOptions->verboseLevel );
    
    // get the general properties about the data
    Gadget_header gadgetHeader;
    size_t noParticles = 0;
    HDF5_initializeGadget( filename, readData, userOptions, &gadgetHeader, &noParticles );
    
    
    // read the data from the files
    size_t numberParticlesRead = 0;
    std::string fileName;
    for (int i=0; i<gadgetHeader.num_files; ++i)
    {
        fileName = gadgetHeader.filename( filename, i );
        message << "Reading GADGET snapshot file '" << fileName << "' which is file " << i+1 << " of " << gadgetHeader.num_files << " files...\n" << MESSAGE::Flush;
        
        HDF5_readGadgetData( fileName, readData, *userOptions, i, &numberParticlesRead );
    }
}






#endif





