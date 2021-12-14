/*
 *  Copyright (c) 2013       Marius Cautun
 *
 *                           ICC, 
 *                           Durham University, UK
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
    
    
    int noScalarsRead = 0;
    // read the temperatures
    if ( userOptions.readParticleData[3] and gadgetHeader.npart[0]>0 )
    {
        // the tempartures make sense only for gas particles
        message << "\t reading the gas temperature ..." << MESSAGE::Flush;
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
        message << "Done\n";
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







/*! Reads the HI Gadget data using 2 HDF5 files: the first gives the gas content while the second gives the HI fraction. This function should be called to do the actual data reading for snapshots saved in single or multiple files.
 */
void HDF5_readGadgetData_HI(std::string filename,
                            std::string h1FileName,
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
    
    
    // open the HDF5 file - this gives the gas data
    H5File *file = new H5File( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    Group *group;
    // open the file with the HI fraction
    H5File *file2 = new H5File( h1FileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    
    
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
    
    
    // get the actual HI mass using the gas data and the HI fraction
    if ( gadgetHeader.npart[0]>0 )
    {
        message << "\t reading the HI fraction ... " << MESSAGE::Flush;
        
        // the element abundance makes sense only for gas particles 
        group = new Group( file->openGroup( "/PartType0/ElementAbundance" ) );
        DataSet dataset = group->openDataSet("Hydrogen");
        float *hydrogenMassFraction = new float[ gadgetHeader.npart[0] ];
        dataset.read( hydrogenMassFraction, PredType::NATIVE_FLOAT );
        delete group;
        
        group = new Group( file2->openGroup( "/PartType0" ) );
        DataSet datasetMolecular = group->openDataSet("MolecularHydrogenMassFraction");
        DataSet datasetHydrogen = group->openDataSet("HydrogenOneFraction");
        
        float *molecularMassFraction = new float[ gadgetHeader.npart[0] ];
        datasetMolecular.read( molecularMassFraction, PredType::NATIVE_FLOAT );
        float *hydrogenOneFraction   = new float[ gadgetHeader.npart[0] ];
        datasetHydrogen.read( hydrogenOneFraction, PredType::NATIVE_FLOAT );
        delete group;
        
        
        // get the HI mass
        float *weights = readData->weight();         // returns a pointer to the particle weights array
        size_t dataOffset = (*numberParticlesRead);  // the offset in the weights array from where to start reading the new values
        for (size_t i=0; i<size_t(gadgetHeader.npart[0]); ++i)
        {
            size_t index = dataOffset + i;
            weights[index] *= hydrogenMassFraction[i] * (float(1.)-molecularMassFraction[i]) * hydrogenOneFraction[i];
            weights[index] = weights[index]<float(0.) ? float(0.) : weights[index];
        }
        
        delete[] hydrogenMassFraction;
        delete[] molecularMassFraction;
        delete[] hydrogenOneFraction;
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
    
    
    int noScalarsRead = 0;
    // read the temperatures
    if ( userOptions.readParticleData[3] and gadgetHeader.npart[0]>0 )
    {
        // the tempartures make sense only for gas particles
        message << "\t reading the gas temperature ..." << MESSAGE::Flush;
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
        message << "Done\n";
    }
    
    
    delete file;
    delete file2;
    for (int i=0; i<6; ++i)     // update the number of read particles
        (*numberParticlesRead) += gadgetHeader.npart[i];
}



/*! The function that reads the HI data from single or multiple HDF5 files.

NOTE_this function reads the HI mass and not the gas mass.
 */
void HDF5_readGadgetFile_HI(std::string filename,
                            Read_data<float> *readData,
                            User_options *userOptions)
{
    MESSAGE::Message message( userOptions->verboseLevel );
    message << "Reading the HI mass fraction from an HDF5 file ...\n" << MESSAGE::Flush;
    
    // get the general properties about the data
    Gadget_header gadgetHeader;
    size_t noParticles = 0;
    HDF5_initializeGadget( filename, readData, userOptions, &gadgetHeader, &noParticles );
    
    // check that the user reads only the gas particles and hence only HI content
    bool otherSpecies = false;
    for (int i=1; i<6; ++i)
        if ( userOptions->readParticleSpecies[1] )
            otherSpecies = true;
    if ( otherSpecies )
        throwError( "You asked for the function that reads in the HI data yet you also asked for reading other particle species. This does not make sense." );
    if ( userOptions->additionalOptions.empty() )
        throwError( "You need to give the name of the HDF5 files giving the Urchin HI fraction if you want to compute the HI mass distribution. This is inserted using the '--options' program option." );
    
    
    // read the data from the files
    size_t numberParticlesRead = 0;
    std::string fileName, h1FractionFile;
    for (int i=0; i<gadgetHeader.num_files; ++i)
    {
        fileName = gadgetHeader.filename( filename, i );
        h1FractionFile = gadgetHeader.filename( userOptions->additionalOptions[0], i );
        message << "Reading GADGET snapshot file '" << fileName << "' and the HI fraction file '" << h1FractionFile << "' which are files " << i+1 << " of " << gadgetHeader.num_files << " files ...\n" << MESSAGE::Flush;
        
        HDF5_readGadgetData_HI( fileName, h1FractionFile, readData, *userOptions, i, &numberParticlesRead );
    }
}






#endif





