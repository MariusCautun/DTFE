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


#include "../math_functions.h"
#define REAL_PRECISSION 1.e-6


/* Print information to the user about the amount of time taken. */
inline void printComputationTime(Timer *t, User_options *userOptions,
                                 string computationQuantityName)
{
    t->stop();
    userOptions->totalTime += t->time();
    MESSAGE::Message message( userOptions->verboseLevel );
    message << "  >>> Time: " << t->time()/userOptions->noProcessors << " sec. (" << computationQuantityName << ")\n" << MESSAGE::Flush;
    t->reset();
}



/* Computes the minor determinant coresponding to entry (row,column) of a 3x3 matrix. */
inline double minorDeterminant(double matrix[][NO_DIM],
                             size_t const row,
                             size_t const column)
{
    size_t i1=(row+1)%3,
    i2=(row+2)%3,
    j1=(column+1)%3,
    j2=(column+2)%3;
    return matrix[j1][i1]*matrix[j2][i2] - matrix[j1][i2]*matrix[j2][i1];
}


/* Computes the inverse of a 2x2 or 3x3 matrix. */
void matrixInverse(double matrix[][NO_DIM],
                   Real result[][NO_DIM])
{
    double tempM[NO_DIM][NO_DIM];
#if NO_DIM==2
    double det = matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
    result[0][0] = matrix[1][1] / det;
    result[0][1] = -matrix[0][1] / det;
    result[1][0] = -matrix[1][0] / det;
    result[1][1] = matrix[0][0] / det;
    
    // check that the determinat is well defined
    double temp = fabs(matrix[0][0]) + fabs(matrix[0][1]) + fabs(matrix[1][0]) + fabs(matrix[1][1]);
    if ( not (std::fabs(det/temp)>REAL_PRECISSION) )
        for (size_t row=0; row<NO_DIM; ++row)
            for (size_t column=0; column<NO_DIM; ++column)
                result[row][column] = Real(0.);
    else
        for (size_t row=0; row<NO_DIM; ++row)
            for (size_t column=0; column<NO_DIM; ++column)
                result[row][column] = Real(tempM);
    
#elif NO_DIM==3
    double det = matrix[0][0] * minorDeterminant( matrix, 0, 0)
            + matrix[0][1] * minorDeterminant( matrix, 1, 0)
            + matrix[0][2] * minorDeterminant( matrix, 2, 0);
    for (size_t row=0; row<NO_DIM; ++row)
        for (size_t column=0; column<NO_DIM; ++column)
            tempM[row][column] = minorDeterminant( matrix, row, column) / det;
    
    // check that the determinat is well defined
    double temp = fabs(matrix[0][0]) + fabs(matrix[0][1]) + fabs(matrix[0][2]) + fabs(matrix[1][0]) + fabs(matrix[1][1]) + fabs(matrix[1][2]) + fabs(matrix[2][0]) + fabs(matrix[2][1]) + fabs(matrix[2][2]);
    if ( not (std::fabs(det/temp)>REAL_PRECISSION) )
        for (size_t row=0; row<NO_DIM; ++row)
            for (size_t column=0; column<NO_DIM; ++column)
                result[row][column] = Real(0.);
    else
        for (size_t row=0; row<NO_DIM; ++row)
            for (size_t column=0; column<NO_DIM; ++column)
                result[row][column] = Real(tempM[row][column]);
#endif
}

/* Computes the vertex position difference matrix for a given Delaunay cell. */
void vertexPositionMatrix(Cell_handle &cell,
                          Real vertexMatrix[][NO_DIM])
{
    Point base = cell->vertex(0)->point();
    // Now store in 'vertexMatrix' the vertices position with respect to the 'base' vertex
    for (int v = 0; v<NO_DIM; ++v)	//loop over vertices != base
        for (int i=0; i<NO_DIM; ++i)	//loop over spatial dimensions
            vertexMatrix[v][i] = cell->vertex(v+1)->point()[i] - base[i];
}



/* Checks the number of vertices of the Delaunay triangulation. The interpolation computation cannot continue if there are no cellls, so it gives a warning message to the user and initializes all the values to 0. */
template <typename T> inline void assingZeroValues(vector<T> *quant, size_t const noElements )
{
    quant->assign( noElements, T::zero() );
}
template <> inline void assingZeroValues<Real>(vector<Real> *quant, size_t const noElements )
{
    quant->assign( noElements, Real(0.) );
}
template <typename T>
inline bool isTriangulationIncomplete(DT &dt,
                                      size_t gridSize,
                                      vector<T> *quant,
                                      string quantityName,
                                      bool showMessage = true)
{
    if ( dt.number_of_vertices()<NO_DIM+1 )
    {
        if ( showMessage )
        {
            MESSAGE::Warning warning(1);
            warning << "Because there are less than " << NO_DIM+1 << " vertices in the Delaunay triangulation there is no cell and hence there is no information that can be used to interpolate the '" << quantityName << "'. All '" << quantityName << "' values will be initialized to 0.\n" << MESSAGE::EndWarning;
        }
        assingZeroValues( quant, gridSize );
        return true;
    }
    return false;
}



/* Reserves memory for an output field and initializes all the elements to 0 if the Delaunay traingulation is not defined. */
template <typename T>
inline void reserveMemory(T *fieldStorage,
                          size_t const gridSize,
                          bool fieldToBeComputed,
                          DT &dt,
                          string quantityName)
{
    if ( fieldToBeComputed )
    {
        fieldStorage->clear();
        fieldStorage->reserve( gridSize );
        isTriangulationIncomplete( dt, gridSize, fieldStorage, quantityName );
    }
}




/* Checks the values of the angles to give correct interval values and transforms the angles from degrees to radians. */
void checkAngles(Real *angles,
                 size_t const size,
                 Real *offset = NULL)
{
    if ( size!=2*(NO_DIM-1) ) throwError( "The function 'checkAngles' in 'density_interpolation.cc' must have the 2nd argument = ", 2*(NO_DIM-1), ", but now that argument is ", size, "." );
    
    // rescale the angles to the interval 0 - 360 (only psi) and theta in the interval -90 - 90
#if NO_DIM==2
    angles[0] -= floor( angles[0]/Real(360.) ) * Real(360.);
    angles[1] -= floor( angles[1]/Real(360.) ) * Real(360.);
#elif NO_DIM==3
    intervalCheck( angles[0], Real(0.), Real(180.), "lower limit of the 'theta' spherical coordinate angle" );
    intervalCheck( angles[1], Real(0.), Real(180.), "upper limit of the 'theta' spherical coordinate angle" );
    angles[2] -= floor( angles[2]/Real(360.) ) * Real(360.);
    angles[3] -= floor( angles[3]/Real(360.) ) * Real(360.);
#endif
    
    // if psi_min>psi_max => psi_max += 360
    if ( offset!=NULL )
{
    offset[0] = Real(0.);
    offset[1] = Real(0.);
}
    if ( angles[size-2]>angles[size-1] )
{
    angles[size-1] += Real(360.);
    if ( offset!=NULL )
    {
        offset[0] = angles[size-2] / Real(2.);
        offset[1] = Real(2.*PI);
    }
}
    
    
    // transform the degrees to radians
    for (size_t i=0; i<size; ++i)
        angles[i] *= Real(PI/180.);
}



/* Computes the area/volume of a Delaunay triangle/tetrahedron in 2D/3D. */
template <typename DTCell>
inline Real volume(DT & dt,
                    DTCell &cell)
{
#if NO_DIM==2
    return dt.triangle( cell ).area();
#elif NO_DIM==3
    return dt.tetrahedron( cell ).volume();
#endif
}



/* Returns the index of the specified grid cell. */
template <typename T>
inline size_t gridCellIndex(T index1, T index2, T index3,
                            T *totalGridSize)
{
#if NO_DIM==2
    return index1*totalGridSize[1] + index2;
#elif NO_DIM==3
    return index1*totalGridSize[1]*totalGridSize[2] + index2*totalGridSize[2] + index3;
#endif
}
template <typename T>
inline size_t gridCellIndex(T *index,
                            size_t *totalGridSize)
{
#if NO_DIM==2
    return index[0]*totalGridSize[1] + index[1];
#elif NO_DIM==3
    return index[0]*totalGridSize[1]*totalGridSize[2] + index[1]*totalGridSize[2] + index[2];
#endif
}



/* Returns the standard deviation of a set of numbers/vectors. */
template <typename T>
inline T standardDeviation(T *data, size_t const size)
{
    T mean = T(0.), result = T(0.);
    for (size_t i=0; i<size; ++i)
        mean += data[i];
    mean /= size;
    for (size_t i=0; i<size; ++i)
    {
        T temp = data[i] - mean;
        result += temp*temp;
    }
    return sqrt( result/size );
}
template <typename T, size_t N>
inline T standardDeviation(Pvector<T,N> *data, size_t const size)
{
    Pvector<T,N> mean = Pvector<T,N>::zero();
    T result = T(0.);
    for (size_t i=0; i<size; ++i)
        mean += data[i];
    mean /= size;
    for (size_t i=0; i<size; ++i)
    {
        Pvector<T,N> temp = data[i] - mean;
        for (size_t j=0; j<N; ++j)
            result += temp[j]*temp[j];
    }
    return sqrt( result/size );
}
