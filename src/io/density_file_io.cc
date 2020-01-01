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




//! Define a header for the output file.

// constants to keep track of how the density was computed
static size_t const DTFE_METHOD = 1;
static size_t const TSC_METHOD = 2;
static size_t const SPH_METHOD = 3;
static size_t const UNKNOW_METHOD = -1;
// constants to keep track what fields the file contains
static int const DENSITY_FILE = 1;
static int const VELOCITY_FILE = 11;
static int const VELOCITY_GRADIENT_FILE = 12;
static int const VELOCITY_DIVERGENCE_FILE = 13;
static int const VELOCITY_SHEAR_FILE = 14;
static int const VELOCITY_VORTICITY_FILE = 15;
static int const VELOCITY_STD_FILE = 16;
static int const SCALAR_FIELD_FILE = 20;
static int const SCALAR_FIELD_GRADIENT_FILE = 21;
static int const UNKNOW_FILE = -1;


static int const fillSize = 1024 - 13*8 - 8*18 - 8*2;   //fill the header up to 1024 bytes

// data structure that stores information in the header of a density file output
struct Density_header
{
    // Information about the density computations
    size_t  gridSize[3];    // the size of the density grid along all 3D directions
    size_t  totalGrid;      // the total size of the grid = gridSize[0]*gridSize[1]*gridSize[2]
    int     fileType;       // keeps track of the field type stored in the file: DENSITY_FILE, VELOCITY_FILE, etc...
    int     noDensityFiles; // the number of density files corresponding to this run
    int     densityFileGrid[3];  // if density is saved within multiple files, each file stores only a patch of the full density (the patches are given by a regular grid with this variable giving the dimensions of that grid)
    int     indexDensityFile;    // keeps track of this density file index compared to the rest for multiple files (this tells the program what region of the density is saved in this file)
    double  box[6];         // keep track of the box coordinates in which the density was computed (xMin, xMax, yMin, yMax, zMin, zMax)
    
    
    // the following are the same as for the gadget file header and are meant to store information about the input file used to compute the density
    size_t   npartTotal[6];// gives the total number of particles in the gadget simulation for which we compute the density profile
    double   mass[6];      // the mass of the particles in the N-body code
    double   time;         // the expansion parameter 'a' of the snapshot file from which we computed the density
    double   redshift;     // the corresponding redshift
    double   BoxSize;      // the box size in kpc
    double   Omega0;       // Omega_matter
    double   OmegaLambda;  // Omega_Lambda
    double   HubbleParam;  // Hubble parameter h (where H=100 h km/s /Mpc )
    
    
    // additional information about files
    size_t  method;        // the method used to compute the density DTFE_METHOD, TSC_METHOD or SPH_METHOD
    char    fill[fillSize];// fill to 1024 bytes - left 760 - used to keep track of information on how the file was obtained
    size_t  FILE_ID;       // keep a unique id for this type of file
    
    
    // constructor - initializes to 0 or to non-assigned value (=-1) 
    Density_header()
    {
        for (int i=0; i<3; ++i)
        {
            gridSize[i] = size_t(0);
            densityFileGrid[i] = 1;
        }
        totalGrid = size_t(0);
        fileType = UNKNOW_FILE;
        noDensityFiles = 1;
        indexDensityFile = -1;
        for (int i=0; i<6; ++i)
        {
            box[i] = 0.;
            npartTotal[i] = 0;
        }
    
        time = -1.; redshift = -1.;
        BoxSize = -1.; Omega0 = -1.; OmegaLambda = -1.; HubbleParam = -1.;
    
        method = UNKNOW_METHOD;
        FILE_ID = 1;
    }
    
    // print the content of the density header to the standard output
    void print()
    {
        std::string densityMethod = "unknown";
        if ( method==DTFE_METHOD ) densityMethod = "DTFE";
        else if ( method==TSC_METHOD ) densityMethod = "TSC";
        else if ( method==SPH_METHOD ) densityMethod = "SPH";
    
        std::string fileTypeName = "unknown file type";
        if ( fileType==DENSITY_FILE ) fileTypeName = "the file stores a density field";
        else if ( fileType==VELOCITY_FILE ) fileTypeName = "the file stores a velocity field";
        else if ( fileType==VELOCITY_GRADIENT_FILE ) fileTypeName = "the file stores the gradient of a velocity field";
        else if ( fileType==VELOCITY_DIVERGENCE_FILE ) fileTypeName = "the file stores a velocity divergence";
        else if ( fileType==VELOCITY_SHEAR_FILE ) fileTypeName = "the file stores a velocity shear";
        else if ( fileType==VELOCITY_VORTICITY_FILE ) fileTypeName = "the file stores a velocity vorticity";
        else if ( fileType==SCALAR_FIELD_FILE ) fileTypeName = "the file stores a scalar field";
        else if ( fileType==SCALAR_FIELD_GRADIENT_FILE ) fileTypeName = "the file stores the gradient of a scalar field";
    
    
        std::cout << "\nThe header of the density file contains the following info:\n" <<
            "1) Information about the actual density computations:\n"
            << "gridSize      = " << gridSize[0] << "  " << gridSize[1] << "  " << gridSize[2] << "\n"
            << "totalGrid     = " << totalGrid << "\n"
            << "file type     = " << fileTypeName << "\n"
            << "# density file= " << noDensityFiles << "\n";
        if ( noDensityFiles>1 )
            std::cout << "file grid size= " << densityFileGrid[0] << "  " << densityFileGrid[1] << "  " << densityFileGrid[3] << "\n"
                << "file index    = " << indexDensityFile << "\n";
        std::cout << "box coords    = " << box[0] << "  " << box[1] << "  " << box[2] << "  " << box[3] << "  " << box[4] << "  " << box[5] << "\n";
            
        std::cout << "\n2) Information about the snapshot file used to compute the density:\n"
            << "npartTotal[6] =  " << npartTotal[0] << "  " << npartTotal[1] << "  " << npartTotal[2] << "  " << npartTotal[3] << "  " << npartTotal[4] << "  " << npartTotal[5] << "\n"
            << "mass[6]       =  " << mass[0] << "  " << mass[1] << "  " << mass[2] << "  " << mass[3] << "  " << mass[4] << "  " << mass[5] << "\n"
            << "time          =  " << time << "\n"
            << "redshift      =  " << redshift << "\n"
            << "BoxSize       =  " << BoxSize << "\n"
            << "Omega0        =  " << Omega0 << "\n"
            << "OmegaLambda   =  " << OmegaLambda << "\n"
            << "HubbleParam   =  " << HubbleParam << "\n";
    
    
        std::cout << "\n3) Information about files and additional remarks:\n"
            << "method          = " << densityMethod << "\n"
            << "fill            = " << fill << "\n\n";
    }
    
    // update entries in the header using the userOptions class
    void updateDensityHeader(User_options userOptions, std::string variableName)
    {
        if ( userOptions.regionOn and not userOptions.regionMpcOn )
        {
            for (size_t i=0; i<userOptions.region.size(); ++i )
                userOptions.region[i] = userOptions.boxCoordinates[i%2] + userOptions.region[i] * (userOptions.boxCoordinates[i%2+1]-userOptions.boxCoordinates[i%2]);
            userOptions.regionMpcOn = true;
        }
        
        for (int i=0; i<NO_DIM; ++i)        // update the grid dimensions
            this->gridSize[i] = userOptions.gridSize[i];
        if ( userOptions.partNo>=0 )        // update the data coordinates if file was split
            for (int i=0; i<NO_DIM; ++i)
                this->densityFileGrid[i] = userOptions.partition[i]; 
        for (int i=0; i<2*NO_DIM; ++i)      // update the box coordinates
            this->box[i] = userOptions.regionOn ? userOptions.region[i] : userOptions.boxCoordinates[i];
        
#if NO_DIM==2
        this->gridSize[2] = 1;
        if ( userOptions.partNo>=0 )
            this->densityFileGrid[2] = 1;
        this->box[4] = 0.;
        this->box[5] = 1.;
#endif
        this->totalGrid = this->gridSize[0]*this->gridSize[1]*this->gridSize[2];
        
        // copy the program options used to get the data
        for (int i=0; i<userOptions.programOptions.length(); ++i)
            this->fill[i] = userOptions.programOptions[i];
        for (int i=userOptions.programOptions.length(); i<fillSize; ++i)
            this->fill[i] = '\0';   // initialize the other elements to empty
        
        // get details about the methods
        if ( userOptions.partNo>=0 )
            this->indexDensityFile = userOptions.partNo;
        if ( userOptions.DTFE ) this->method = DTFE_METHOD;
        else if ( userOptions.TSC ) this->method = TSC_METHOD;
        else if ( userOptions.SPH ) this->method = SPH_METHOD;
        
        // get details about the output fields in the file
        if ( variableName.find( "density" )!=std::string::npos ) this->fileType = DENSITY_FILE;
        else if ( variableName.find( "velocity" )!=std::string::npos ) this->fileType = VELOCITY_FILE;
        else if ( variableName.find( "velocity gradient" )!=std::string::npos ) this->fileType = VELOCITY_GRADIENT_FILE;
        else if ( variableName.find( "velocity divergence" )!=std::string::npos ) this->fileType = VELOCITY_DIVERGENCE_FILE;
        else if ( variableName.find( "velocity shear" )!=std::string::npos ) this->fileType = VELOCITY_SHEAR_FILE;
        else if ( variableName.find( "velocity vorticity" )!=std::string::npos ) this->fileType = VELOCITY_VORTICITY_FILE;
        else if ( variableName.find( "velocity standard deviation" )!=std::string::npos ) this->fileType = VELOCITY_STD_FILE;
        else if ( variableName.find( "scalar" )!=std::string::npos ) this->fileType = SCALAR_FIELD_FILE;
        else if ( variableName.find( "scalar gradient" )!=std::string::npos ) this->fileType = SCALAR_FIELD_GRADIENT_FILE;
    }
    
    // copy information from the gadget header
    void copyGadgetHeaderInfo(User_options userOptions)
    {
        Gadget_header gadgetHeader;
        std::string filename = gadgetHeader.filename( userOptions.inputFilename, 0, false );
        if ( not bfs::exists(filename) ) return;    //file could not be found
        
        if ( userOptions.inputFileType==101 or userOptions.inputFileType==102 )
        {
            // open the first binary file for reading and read some of the overall characteristics
            std::fstream inputFile;
            openInputBinaryFile( inputFile, filename );
            
            // read the header
            int buffer, gadgetFileType;
            bool swapEndian = false;
            inputFile.read( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
            bool validFile = gadgetHeader.detectSnapshotType( buffer, &gadgetFileType, &swapEndian );
            if ( not validFile )
                return;
            int offset = gadgetFileType==2 ? 16+sizeof(buffer) : 0+sizeof(buffer);
            inputFile.seekg( offset, std::ios::beg );
            inputFile.read( reinterpret_cast<char *>(&gadgetHeader), sizeof(gadgetHeader) );
            inputFile.close();
        }
#ifdef HDF5
        else if ( userOptions.inputFileType==105 )
            HDF5_readGadgetHeader( filename, &gadgetHeader );
#endif
        else
            return;
        
        // copy the info from the gadget header
        for (int i=0; i<6; ++i)
        {
            this->npartTotal[i] = gadgetHeader.npartTotal[i];
            this->mass[i] = gadgetHeader.mass[i];
        }
        this->time = gadgetHeader.time;
        this->redshift = gadgetHeader.redshift;
        this->BoxSize = gadgetHeader.BoxSize;
        this->Omega0 = gadgetHeader.Omega0;
        this->OmegaLambda = gadgetHeader.OmegaLambda;
        this->HubbleParam = gadgetHeader.HubbleParam;
    }
};





/*! Writes a binary density file - it has a custom header that keeps track of the data in the file. */
template <typename T>
void writeSpecialFile(T const &dataToWrite,
                        std::string filename,
                        std::string variableName,
                        User_options const &userOptions)
{
    // update the density header information
    Density_header densityHeader;
    densityHeader.updateDensityHeader( userOptions, variableName );
    densityHeader.copyGadgetHeaderInfo( userOptions );
    
    
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "Writing the " << variableName << " to the file '" << filename << "' ...  " << MESSAGE::Flush;
    
    // open the file
    std::fstream outputFile;
    openOutputBinaryFile( outputFile, filename );
    
    
    // write the header
    size_t buffer = sizeof( densityHeader );
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    outputFile.write( reinterpret_cast<char *>(&densityHeader), sizeof(densityHeader) );
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    
    // write the position of the particles
    buffer = dataToWrite.size()*sizeof(dataToWrite[0]); //total number of data bytes written
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    size_t maxSize = 256*256*256;   //write at most 256^3 elements at once, otherwise the write function might fail
    size_t noRepeats = size_t( dataToWrite.size() / maxSize ), currentPosition = 0;
    size_t tempBuffer = maxSize * sizeof(dataToWrite[0]);
    for (size_t i=0; i<noRepeats; ++i)  //write in blocks of 256^3 elements
    {
        outputFile.write( reinterpret_cast<char const *>(&(dataToWrite[currentPosition])), tempBuffer );
        currentPosition += maxSize;
    }
    tempBuffer = (dataToWrite.size() - currentPosition) * sizeof(dataToWrite[0]);    // write everything else that is left
    outputFile.write( reinterpret_cast<char const *>(&(dataToWrite[currentPosition])), tempBuffer );
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    
    outputFile.close();
    message << "Done.\n";
}





