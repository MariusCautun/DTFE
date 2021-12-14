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
#include <math.h>
using namespace H5;




/*! Reads the Gadget data from a single HDF5 file. This function should be called to do the actual data reading for snapshots saved in single or multiple files.
*/
void HDF5_readGadgetData_Cristian(std::string filename,
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
            DataSet dataset = group->openDataSet("IntegerCoordinates");
            
            // read the positions -- these are saved as integers
            unsigned int *int_positions = new unsigned int[ gadgetHeader.npart[type] * NO_DIM ];
            dataset.read( int_positions, PredType::NATIVE_UINT );
            delete group;
            
            // convert teh positions to floats
            float scaling = gadgetHeader.BoxSize / pow(2.,32.);
            for(int i=0; i<gadgetHeader.npart[type] * NO_DIM; i++)
                positions[dataOffset+i] = int_positions[i] * scaling;
            
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
            DataSet dataset = group->openDataSet("Velocities");
            
            dataset.read( &(velocities[dataOffset]), PredType::NATIVE_FLOAT );
            delete group;
            
            dataOffset += gadgetHeader.npart[type] * NO_DIM;
        }
        message << "Done\n";
    }
    
    
    
    delete file;
    for (int i=0; i<6; ++i)     // update the number of read particles
        (*numberParticlesRead) += gadgetHeader.npart[i];
}





/*! The function that reads the data from single or multiple HDF5 files.
*/
void HDF5_readGadgetFile_Cristian(std::string filename,
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
        
        HDF5_readGadgetData_Cristian( fileName, readData, *userOptions, i, &numberParticlesRead );
    }
}








#endif





