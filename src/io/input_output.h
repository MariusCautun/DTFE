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



/* This files defines 2 classes used for easy reading of the data from an input file.
 
It also defines 4 functions to open binary/text files for input/output. */

#include <algorithm>


// structure for pairs of pointer and bool used in structure 'Read_data'
template <typename T>
struct PairPtrBool
{
    T *_ptr;
    bool _assigned;
    
    PairPtrBool()   //constructor - initializes everything to NULL and false
    {
        _ptr = NULL;
        _assigned = false;
    }
    ~PairPtrBool()  //destructor - frees memory if assigned
    {
        if ( _assigned ) delete[] _ptr;
    }
    
    // assigns memory to the pointer
    T* assignMemory(size_t const size,
                    size_t *expectedSize,
                    size_t dimensions,
                    std::string name)
    {
        // if memory was already allocated, throw error
        if ( _assigned ) throwError( "When allocating memory for the variable 'Read_data::", name, "'. The memory was already allocated for the given variable." );
        // if 'size'!='expected size'
        if ( size!=(*expectedSize) and (*expectedSize)!=size_t(0) )
        {
            MESSAGE::Error error;
            error << "When allocating memory for the variable 'Read_data::" << name << "'. The size of the variable (which is " << size << ") is different from the expected size of " << *expectedSize << " (all variables: position, velocity, weight and scalar must have the same size)." << MESSAGE::EndError;
        }
        *expectedSize = size;
        _ptr = new T[dimensions*size];
        _assigned = true;
        return _ptr;
    }
    
    // returns a pointer to the data
    T* returnPointer(std::string name)
    {
        // if memory wasn't already allocated, throw error
        if ( not _assigned ) throwError( "When returning pointer for the variable 'Read_data::", name, "'. There was no memory allocation for the given variable." );
        return _ptr;
    }
};



// structure to transfer the read data between functions
template <typename T>
struct Read_data
{
    size_t _noParticles; // stores the number of particles 
    PairPtrBool<T>  _position;   // pointer to array storing the position data
    PairPtrBool<T>  _velocity;   // pointer to array storing the velocity data
    PairPtrBool<T>  _weight;     // pointer to array storing the weight data
    PairPtrBool<T>  _scalar;     // pointer to array storing the scalar data
    
    size_t _noSamples;   // stores the number of user given sample points (if any)
    PairPtrBool<T>  _sampling;   // pointer to array storing the user given sampling grid points
    PairPtrBool<T>  _delta;      // pointer to array storing the size of the user given sampling grids 
    
    
    Read_data()
    {
        _noParticles = 0;
        _noSamples = 0;
    }
    
    
    size_t noParticles() { return _noParticles;}
    size_t noSamples() { return _noSamples;}
    
     // return pointer to position data (valid only if assigned memory)
    T * position()
    { return _position.returnPointer("position"); }
    // assign memory to store position data and return pointer to assigned memory
    T * position( size_t const noParticles)
    { return _position.assignMemory( noParticles, &_noParticles, NO_DIM, "position" ); }
    
    T * velocity()
    { return _velocity.returnPointer("velocity"); }
    T * velocity( size_t const noParticles)
    { return _velocity.assignMemory( noParticles, &_noParticles, noVelComp, "velocity" ); }
    
    T * weight()
    { return _weight.returnPointer("weight"); }
    T * weight( size_t const noParticles)
    { return _weight.assignMemory( noParticles, &_noParticles, 1, "weight" ); }
    
    T * scalar()
    { return _scalar.returnPointer("scalar"); }
    T * scalar( size_t const noParticles)
    { return _scalar.assignMemory( noParticles, &_noParticles, noScalarComp, "scalar" ); }
    
    
    T * sampling()
    { return _sampling.returnPointer("sampling"); }
    T * sampling( size_t const noSamples)
    { return _sampling.assignMemory( noSamples, &_noSamples, NO_DIM, "sampling" ); }
    
    T * delta()
    { return _delta.returnPointer("delta"); }
    T * delta( size_t const noSamples)  // assign memory to store delta for user given sampling points and return pointer to assigned memory
    { return _delta.assignMemory( noSamples, &_noSamples, NO_DIM, "delta" ); }
    
    
    // writes the data kept in pointers of this structure to a 'Particle_data' and 'Sample_point' vectors
    void transferData(std::vector<Particle_data> *p,
                      std::vector<Sample_point>  *s)
    {
        // first write the particle data
        if ( _noParticles>size_t(0) )
        {
            p->clear();
            p->reserve( _noParticles );
            for (size_t i=0; i<_noParticles; ++i)
            {
                Particle_data temp;
                if ( _position._assigned )   // copy positions if assigned
                    for (size_t j=0; j<NO_DIM; ++j)
                        temp.position(j) = _position._ptr[NO_DIM*i+j];
                if ( _velocity._assigned )   // copy velocities if assigned
                    for (size_t j=0; j<noVelComp; ++j)
                        temp.velocity(j) = _velocity._ptr[noVelComp*i+j];
                if ( _weight._assigned )   // copy weight if assigned
                    temp.weight() = _weight._ptr[i];
                if ( _scalar._assigned )   // copy scalar fields if assigned
                    for (size_t j=0; j<noScalarComp; ++j)
                        temp.scalar(j) = _scalar._ptr[noScalarComp*i+j];
                p->push_back( temp );
            }
        }
        if ( _noSamples>size_t(0) )
        {
            s->clear();
            s->reserve( _noSamples );
            for (size_t i=0; i<_noSamples; ++i)
            {
                Sample_point temp;
                if ( _sampling._assigned )   // copy sample positions if assigned
                    for (int j=0; j<NO_DIM; ++j)
                        temp.position(j) = _sampling._ptr[NO_DIM*i+j];
                if ( _delta._assigned )      // copy grid point sizes if assigned
                    for (int j=0; j<NO_DIM; ++j)
                        temp.delta(j) = _delta._ptr[NO_DIM*i+j];
                s->push_back( temp );
            }
        }
    }
};




/* This function opens a binary file for input (to read from the file). It checks that the operation was done succesfully. */
void openInputBinaryFile(std::fstream & inputFile,
                         std::string & fileName)
{
    char const *temp = fileName.c_str();
    inputFile.open( temp, std::ios::in | std::ios::binary );
    
    if ( not inputFile.is_open() )
    {
        std::cout << "~~~ ERROR ~~~ The file '" << fileName << "' could not be opened for reading! The program ended unsucessfully!\n";
        exit( EXIT_FAILURE );
    }
}
/* This function opens a binary file for output (to write data to the file). It checks that the operation was done succesfully. */
void openOutputBinaryFile(std::fstream & outputFile,
                          std::string & fileName)
{
    char const *temp = fileName.c_str();
    outputFile.open( temp, std::ios::out | std::ios::binary );
    
    if ( not outputFile.is_open() )
    {
        std::cout << "~~~ ERROR ~~~ The file '" << fileName << "' could not be opened for writing! The program ended unsucessfully!\n";
        exit( EXIT_FAILURE );
    }
}
/* This function opens a text file for input (to read from the file). It checks that the operation was done succesfully. */
void openInputTextFile(std::fstream & inputFile,
                       std::string & fileName)
{
    char const *temp = fileName.c_str();
    inputFile.open( temp, std::ios::in );
    
    if ( not inputFile.is_open() )
    {
        std::cout << "~~~ ERROR ~~~ The file '" << fileName << "' could not be opened for reading! The program ended unsucessfully!\n";
        exit( EXIT_FAILURE );
    }
}
/* This function opens a text file for output (to write data to the file). It checks that the operation was done succesfully. */
void openOutputTextFile(std::fstream & outputFile,
                        std::string & fileName)
{
    char const *temp = fileName.c_str();
    outputFile.open( temp, std::ios::out );
    
    if ( not outputFile.is_open() )
    {
        std::cout << "~~~ ERROR ~~~ The file '" << fileName << "' could not be opened for writing! The program ended unsucessfully!\n";
        exit( EXIT_FAILURE );
    }
}


/* This function check that the file operations (reading/writting to file) were execute correctly. */
void checkFileOperations(std::fstream & file,
                         std::string operationName)
{
    MESSAGE::Error error;
    if ( file.bad() )
        error << "The " << operationName << " file operation has failed. The function ios::bad() has return an error that can be due to the file not being open or if the device has no more writing space left." << MESSAGE::EndError;
    else if ( file.fail() )
        error << "The " << operationName << " file operation has failed. A format error has happened. There are two possible causes for the error: you were trying to read a number, but an alphabetical character was found or you got to the end of the file before being able to read all the data." << MESSAGE::EndError;
//     else if ( file.eof() )
//         error << "The " << operationName << " to/from file operation has failed. The file ended before all the data could be " << operationName << " to/from file." << MESSAGE::EndError;
}


// Function used to swap bytes between different endianness
inline void ByteSwap(unsigned char * b, int n)
{
	register int i = 0;
	register int j = n-1;
	while (i<j)
	{
		std::swap(b[i], b[j]);
		i++, j--;
	}
}
template <typename T>
void ByteSwapArray(T *x, size_t const elements)
{
	int size = sizeof(x[0]);
	for (size_t i=0; i<elements; ++i)
		ByteSwap( (unsigned char *) &(x[i]), size );
}

#define BYTESWAP(x) ByteSwap( (unsigned char *) &x, sizeof(x) )

