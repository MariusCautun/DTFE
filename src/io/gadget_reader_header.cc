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
 
 
#define SWAP_HEADER_ENDIANNESS(x1,x2,x3,x4) { if( x1 ) {BYTESWAP( x2 ); BYTESWAP( x3 ); x4.swapBytes();} }
#define SWAP_ENDIANNESS(x1,x2,x3)           { if( x1 ) {BYTESWAP( x2 ); BYTESWAP( x3 );} }
#define READ_DELIMETER \
    inputFile.seekg( offset, std::ios::cur ); \
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) )
#define DELIMETER_CONSISTANCY_CHECK(field) \
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) ); \
    SWAP_ENDIANNESS( swapEndian, buffer1, buffer2 ); \
    if ( buffer1!=buffer2 ) \
        throwError( "The integers before and after the particle " field " data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." )


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
            if ( buffer1 == 8 )         // gadget file format 2
                *gadgetFileType = 2;
            else if ( buffer1 == 256 )  // gadget file format 1
                *gadgetFileType = 1;
            else                        // could not detect the file type
                return false;
        }
        return true;
    }
};

