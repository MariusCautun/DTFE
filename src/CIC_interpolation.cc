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


#include <vector>
#include <list>
#include <string>
#include <cmath>

#include "define.h"
#include "particle_data.h"
#include "quantities.h"
#include "box.h"
#include "user_options.h"
#include "message.h"

using namespace std;
typedef vector<Particle_data>::iterator     vectorIterator;


void CIC_interpolation_regular_grid(vector<Particle_data> &particles,
                                    User_options &userOptions,
                                    Quantities *q);




/* This function interpolates the density and velocity to grid using the CIC (cloud in cell) method. */
void CIC_interpolation(vector<Particle_data> *particles,
                       vector<Sample_point> &samples,
                       User_options &userOptions,
                       Quantities *q)
{
    if ( samples.empty() and not userOptions.redshiftConeOn )  // interpolate using a regular cubic grid
        CIC_interpolation_regular_grid( *particles, userOptions, q );
    else
        throwError( "The CIC method can interpolate the fields only on a regular rectangular grid. No CIC interpolation methods are implemented for redshift cone coordinates or for user defined sample points." );
    particles->clear();
}





/* This function uses the CIC method to interpolate quantities to a grid. It interpolates only the density and the velocity.
NOTE: It does not interpolate the velocity to the grid, but in fact the momentum. The velocity is than obtained as the momentum in the cell divided by the mass in the grid cell.
*/
void CIC_interpolation_regular_grid(vector<Particle_data> &particles,
                                    User_options &userOptions,
                                    Quantities *q)
{
    size_t const *nGrid = &(userOptions.gridSize[0]);
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "\nInterpolating the fields to the grid using the CIC method. The interpolation takes place inside the box of coordinates " << userOptions.region.print()
            << " on a " << MESSAGE::printElements( nGrid, NO_DIM, "*" ) << " grid ... " << MESSAGE::Flush;
    boost::timer t;
    t.restart();
    
    
    // allocate memory for the results
    size_t const reserveSize = (NO_DIM==2) ? nGrid[0]*nGrid[1] : nGrid[0]*nGrid[1]*nGrid[2];
    q->density.assign( reserveSize, Real(0.) );
    if ( userOptions.aField.velocity )
		q->velocity.assign( reserveSize, Pvector<Real,noVelComp>::zero() );
    
    
    // get the grid spacing
    Box box = userOptions.region;
    Real dx[NO_DIM];
    for (int i=0; i<NO_DIM; ++i)
        dx[i] = (box[2*i+1]-box[2*i]) / nGrid[i];
    
    
    // find the particles in box 'box'
    for (vectorIterator it=particles.begin(); it!=particles.end(); ++it)
	{
		int cell[NO_DIM];
		bool validCell = true;
		for (int j=0; j<NO_DIM; ++j)
		{
			cell[j] = int( floor( (it->position(j)-box[2*j])/dx[j] ) );
			if ( cell[j]<0 or cell[j]>=nGrid[j] )
				validCell = false;
		}
		if ( not validCell ) continue;
#if NO_DIM==2
		int index = cell[0] * nGrid[1] + cell[1];
#elif NO_DIM==3
		int index = cell[0] * nGrid[1]*nGrid[2] + cell[1] * nGrid[2] + cell[2];
#endif
		q->density[index] += it->weight();
		if ( userOptions.aField.velocity )
			q->velocity[index] += it->velocity() * it->weight();
	}
    
    
    // divide the momentum by the mass in the cells
    if ( userOptions.aField.velocity )
        for (size_t i=0; i<reserveSize; ++i)
            if ( q->density[i]!=Real(0.) )
                q->velocity[i] /= q->density[i];
            else
                q->velocity[i] = Pvector<Real,noVelComp>::zero();
    
    // normalize the density to average background density
    if ( userOptions.aField.density )
    {
        Real factor = Real( q->density.size() ) / box.volume() / userOptions.averageDensity;
        for (vector<Real>::iterator it=q->density.begin(); it!=q->density.end(); ++it)
            (*it) *= factor;
    }
    else
        q->density.clear();
    
    message << "Done.\n" << MESSAGE::Flush;
    printElapsedTime( &t, &userOptions, "CIC interpolation" );
}


/* This function counts how many particles are in each cell of a grid using the CIC method. */
void CIC_particle_count(vector<Particle_data> &particles,
                        size_t const *nGrid,
                        Box box,
                        vector<int> *counts)
{
    // allocate memory for the results
    size_t const reserveSize = (NO_DIM==2) ? nGrid[0]*nGrid[1] : nGrid[0]*nGrid[1]*nGrid[2];
    counts->assign( reserveSize, int(0) );
	
    // get the grid spacing
    Real dx[NO_DIM];
    for (int i=0; i<NO_DIM; ++i)
        dx[i] = (box[2*i+1]-box[2*i]) / nGrid[i];
    
    
    // find the particles in box 'box'
    for (vectorIterator it=particles.begin(); it!=particles.end(); ++it)
	{
		int cell[NO_DIM];
		bool validCell = true;
		for (int j=0; j<NO_DIM; ++j)
		{
			cell[j] = int( floor( (it->position(j)-box[2*j])/dx[j] ) );
			if ( cell[j]<0 or cell[j]>=nGrid[j] )
				validCell = false;
		}
		if ( not validCell ) continue;
#if NO_DIM==2
		int index = cell[0] * nGrid[1] + cell[1];
#elif NO_DIM==3
		int index = cell[0] * nGrid[1]*nGrid[2] + cell[1] * nGrid[2] + cell[2];
#endif
		(*counts)[index] += 1;
	}
}



