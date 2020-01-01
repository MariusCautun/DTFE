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


#include "quantities.h"
#include "message.h"




/* This function copies one field of 'Quantities' at a time from the subgrid results to the main grid ones. */
template<typename T>
void copyField(T const &subgridResults,
                T *mainGridResults,
                bool toCopy,
                std::vector<size_t> const &mainGrid,
                std::vector<size_t> const &subgrid,
                std::vector<size_t> const &subgridOffset)
{
    if ( not toCopy )
        return;
    
    size_t start[NO_DIM], end[NO_DIM];
    for (size_t i=0; i<NO_DIM; ++i)
    {
        start[i] = subgridOffset[i];
        end[i] = subgridOffset[i] + subgrid[i];
    }
    
    
    for (size_t i1=0, i2=start[0]; i1<subgrid[0]; ++i1)
    {
        for (size_t j1=0, j2=start[1]; j1<subgrid[1]; ++j1)
        {
#if NO_DIM==2
            (*mainGridResults)[i2*mainGrid[1]+j2] = subgridResults[i1*subgrid[1]+j1];
#elif NO_DIM==3
            for (size_t k1=0, k2=start[2]; k1<subgrid[2]; ++k1)
            {
                size_t indexMain = i2*mainGrid[1]*mainGrid[2] + j2*mainGrid[2] + k2;
                size_t indexSec = i1*subgrid[1]*subgrid[2] + j1*subgrid[2] + k1;
                (*mainGridResults)[indexMain] = subgridResults[indexSec];
                ++k2;
            }
#endif
            ++j2;
        }
        ++i2;
    }
}

/* Function used to copy the results obtain on a subgrid using the option 'partion' to the results for the full grid.
NOTE: For additional details on how to compute the 'subgrid' and 'subgridOffset' check the function 'subgridAxesSize' in the file 'subpartition.h'. */
void Quantities::copyFromSubgrid(Quantities const &subgridResults,
                                 Field const &field,
                                 std::vector<size_t> const &mainGrid,
                                 std::vector<size_t> const &subgrid,
                                 std::vector<size_t> const &subgridOffset)
{
    copyField( subgridResults.density, &(this->density), field.density,  mainGrid, subgrid, subgridOffset );
    copyField( subgridResults.velocity, &(this->velocity), field.velocity,  mainGrid, subgrid, subgridOffset );
    copyField( subgridResults.velocity_gradient, &(this->velocity_gradient), field.velocity_gradient,  mainGrid, subgrid, subgridOffset );
    copyField( subgridResults.velocity_divergence, &(this->velocity_divergence), field.velocity_divergence,  mainGrid, subgrid, subgridOffset );
    copyField( subgridResults.velocity_shear, &(this->velocity_shear), field.velocity_shear,  mainGrid, subgrid, subgridOffset );
    copyField( subgridResults.velocity_vorticity, &(this->velocity_vorticity), field.velocity_vorticity,  mainGrid, subgrid, subgridOffset );
    copyField( subgridResults.velocity_std, &(this->velocity_std), field.velocity_std,  mainGrid, subgrid, subgridOffset );
    copyField( subgridResults.scalar, &(this->scalar), field.scalar,  mainGrid, subgrid, subgridOffset );
    copyField( subgridResults.scalar_gradient, &(this->scalar_gradient), field.scalar_gradient,  mainGrid, subgrid, subgridOffset );
}



/* Compares the size of a vector with the expected size and outputs error if this is not the case. */
template< typename T>
void fieldSize(T const & object,
                size_t *expectedSize)
{
    if( not object.empty() )
    {
        if ( (*expectedSize)!=0 and (*expectedSize)!=object.size() )
            throwError( "Two or more objects of class 'Quantities' have different sizes. All objects in this class should be empty or have the same size." );
        else
            (*expectedSize) = object.size();
    }
}

/* Returns what is the size of any non-empty object.
NOTE: All non-empty objects should have the same size. */
size_t Quantities::size() const
{
    size_t temp = 0;
    fieldSize( this->density, &temp );
    fieldSize( this->velocity, &temp );
    fieldSize( this->velocity_gradient, &temp );
    fieldSize( this->velocity_divergence, &temp );
    fieldSize( this->velocity_shear, &temp );
    fieldSize( this->velocity_vorticity, &temp );
    fieldSize( this->velocity_std, &temp );
    fieldSize( this->scalar, &temp );
    fieldSize( this->scalar_gradient, &temp );
    return temp;
}


/* This function reserves memory for the main grid quantities when using the 'partition' option. */
void Quantities::reserveMemory(size_t *gridSize, Field &field)
{
    size_t totalSize = 1;
    for (size_t i=0; i<NO_DIM; ++i)
        totalSize *= gridSize[i];
    
    if ( field.density )
        this->density.resize( totalSize, Real(0.) );
    if ( field.velocity )
        this->velocity.resize( totalSize, Pvector<Real,noVelComp>::zero() );
    if ( field.velocity_gradient )
        this->velocity_gradient.resize( totalSize, Pvector<Real,noGradComp>::zero() );
    if ( field.velocity_divergence )
        this->velocity_divergence.resize( totalSize, Real(0.) );
    if ( field.velocity_shear )
        this->velocity_shear.resize( totalSize, Pvector<Real,noShearComp>::zero() );
    if ( field.velocity_vorticity )
        this->velocity_vorticity.resize( totalSize, Pvector<Real,noVortComp>::zero() );
    if ( field.velocity_std )
        this->velocity_std.resize( totalSize, Real(0.) );
    if ( field.scalar )
        this->scalar.resize( totalSize, Pvector<Real,noScalarComp>::zero() );
    if ( field.scalar_gradient )
        this->scalar_gradient.resize( totalSize, Pvector<Real,noScalarGradComp>::zero() );
}


