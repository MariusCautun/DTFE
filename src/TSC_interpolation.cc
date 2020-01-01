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


void TSC_interpolation_regular_grid(vector<Particle_data> &particles,
                                    User_options &userOptions,
                                    Quantities *q);




/* This function interpolates the density and velocity to grid using the TSC (triangular shape cloud) method. */
void TSC_interpolation(vector<Particle_data> *particles,
                       vector<Sample_point> &samples,
                       User_options &userOptions,
                       Quantities *q)
{
    if ( samples.empty() and not userOptions.redshiftConeOn )  // interpolate using a regular cubic grid
        TSC_interpolation_regular_grid( *particles, userOptions, q );
    else
        throwError( "The TSC method can interpolate the fields only on a regular rectangular grid. No TSC interpolation methods are implemented for redshift cone coordinates or for user defined sample points." );
    particles->clear();
}





/* This function uses the TSC method to interpolate quantities to a grid. It interpolates only the density and the velocity.
NOTE: It does not interpolate the velocity to the grid, but in fact the momentum. The velocity is than obtained as the momentum in the cell divided by the mass in the grid cell.
*/
void TSC_interpolation_regular_grid(vector<Particle_data> &particles,
                                    User_options &userOptions,
                                    Quantities *q)
{
    size_t const *nGrid = &(userOptions.gridSize[0]);
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "\nInterpolating the fields to the grid using the TSC method. The interpolation takes place inside the box of coordinates " << userOptions.region.print()
            << " on a " << MESSAGE::printElements( nGrid, NO_DIM, "*" ) << " grid ... " << MESSAGE::Flush;
    boost::timer t;
    t.restart();
    
    
    // allocate memory for the results
    size_t const reserveSize = (NO_DIM==2) ? nGrid[0]*nGrid[1] : nGrid[0]*nGrid[1]*nGrid[2];
    q->density.assign( reserveSize, Real(0.) );
    q->velocity.assign( reserveSize, Pvector<Real,noVelComp>::zero() );
    
    
    // get the grid spacing and find the extended box of particles that give contribution in the region of interest
    Box box = userOptions.region,
        outerBox = userOptions.region,
        innerBox = userOptions.region;
    Real dx[NO_DIM], outerPadding[2*NO_DIM], innerPadding[2*NO_DIM];
    for (int i=0; i<NO_DIM; ++i)
    {
        dx[i] = (box[2*i+1]-box[2*i]) / nGrid[i];
        outerPadding[2*i] = dx[i];
        outerPadding[2*i+1] = dx[i];
    }
    for (int i=0; i<2*NO_DIM; ++i)
        innerPadding[i] = -1.1*outerPadding[i];    // slightly larger than 1 grid cell to account for numerical uncertanties
    outerBox.addPadding( outerPadding );  // only particles in this box give contributions in 'userOptions.region'
    innerBox.addPadding( innerPadding );  // the contribution of particles in this box is limited to the region of interest
    
    
    // find the particles in box 'box'
    list< vectorIterator > innerParticles, outerParticles;
    for (vectorIterator it=particles.begin(); it!=particles.end(); ++it)
        if ( innerBox.isParticleInBox(*it) )
            innerParticles.push_back( it );     // keep track of the particles that give contributions only inside the region of interest
        else if ( outerBox.isParticleInBox(*it) )
            outerParticles.push_back( it );     // these particles give contributions also outside the region of interest
    
    
    
    // loop over all the particles in the inner box
    for (list< vectorIterator >::iterator it=innerParticles.begin(); it!=innerParticles.end(); ++it)
    {
        int cell[NO_DIM][3];    // the grid coordinates of the cell that get a contribution from the particle
        Real weight[NO_DIM][3]; // the weight associated to the particle for each cell that it contributes to
        for (int j=0; j<NO_DIM; ++j)
        {
            Real temp = ( (*it)->position(j) - box[2*j] ) / dx[j];
            cell[j][0] = int(floor( temp )) - 1;
            cell[j][1] = cell[j][0] + 1;
            cell[j][2] = cell[j][1] + 1;
            
            temp -= cell[j][1] + 0.5; //this gives the distance of the particle with respect to the center of the density grid in which the particle is located
            
            // compute the weight associated to this direction
            weight[j][0] = .5 * (0.5-temp) * (0.5-temp);
            weight[j][1] = .75 - temp * temp;
            weight[j][2] = .5 * (0.5+temp) * (0.5+temp);
        }
        
        // get the density contribution of the particle to the neighboring cells
#if NO_DIM==2
        for (int i1=0; i1<3; ++i1)
            for (int i2=0; i2<3; ++i2)
            {
                int index = cell[0][i1] * nGrid[1] + cell[1][i2];
                Real result = (*it)->weight() * weight[0][i1] * weight[1][i2];
                q->density[index] += result;
                q->velocity[index] += (*it)->velocity() * result;
            }
#elif NO_DIM==3
        for (int i1=0; i1<3; ++i1)
            for (int i2=0; i2<3; ++i2)
                for (int i3=0; i3<3; ++i3)
                {
                    int index = cell[0][i1] * nGrid[1]*nGrid[2] + cell[1][i2] * nGrid[2] + cell[2][i3];
                    Real result = (*it)->weight() * weight[0][i1] * weight[1][i2] * weight[2][i3];
                    q->density[index] += result;
                    q->velocity[index] += ( (*it)->velocity() * result );
                }
#endif
    }
    
    
    // now loop over the particles on the boundary - must check that all neighbors are valid cells
    for (list< vectorIterator >::iterator it=outerParticles.begin(); it!=outerParticles.end(); ++it)
    {
        int cell[NO_DIM][3];    // the grid coordinates of the cell that get a contribution from the particle
        int cellCount[NO_DIM];  // counts how many valid cells are along each dimension that get contribution from the particle
        for (int i=0; i<NO_DIM; ++i) cellCount[i] = 0;
        Real weight[NO_DIM][3]; // the weight associated to the particle for each cell that it contributes to
        for (int j=0; j<NO_DIM; ++j)
        {
            Real temp = ( (*it)->position(j) - box[2*j] ) / dx[j];
            int tempInt = int(floor( temp ));
            temp -= tempInt + 0.5; //this gives the distance of the particle with respect to the center of the density grid in which the particle is located
            if ( temp<-0.5 or temp>0.5 ) cout << temp << "\t" << j << "\t" << tempInt << "\n";
            
            // check that all neighbors must be valid cells
            if ( tempInt==-1 )
            {
                cell[j][0] = 0; cellCount[j] = 1;
                weight[j][0] = .5 * (0.5+temp) * (0.5+temp);
            }
            else if ( tempInt==0 )
            {
                cell[j][0] = 0; cell[j][1] = 1; cellCount[j] = 2;
                weight[j][0] = .75 - temp * temp;
                weight[j][1] = .5 * (0.5+temp) * (0.5+temp);
            }
            else if ( tempInt==int(nGrid[j])-1 )
            {
                cell[j][0] = nGrid[j]-2; cell[j][1] = nGrid[j]-1; cellCount[j] = 2;
                weight[j][0] = .5 * (0.5-temp) * (0.5-temp);
                weight[j][1] = .75 - temp * temp;
            }
            else if ( tempInt==int(nGrid[j]) )
            {
                cell[j][0] = nGrid[j]-1; cellCount[j] = 1;
                weight[j][0] = .5 * (0.5-temp) * (0.5-temp);
            }
            else if ( not( tempInt<-1 or tempInt>int(nGrid[j]) ) )
            {
                cell[j][0] = tempInt-1;
                cell[j][1] = tempInt;
                cell[j][2] = tempInt+1;
                cellCount[j] = 3;
                
                weight[j][0] = .5 * (0.5-temp) * (0.5-temp);
                weight[j][1] = .75 - temp * temp;
                weight[j][2] = .5 * (0.5+temp) * (0.5+temp);
            }
        }
        
        // get the density contribution of the particle to the neighboring cells
#if NO_DIM==2
        for (int i1=0; i1<cellCount[0]; ++i1)
            for (int i2=0; i2<cellCount[1]; ++i2)
            {
                int index = cell[0][i1] * nGrid[1] + cell[1][i2];
                Real result = (*it)->weight() * weight[0][i1] * weight[1][i2];
                q->density[index] += result;
                q->velocity[index] += (*it)->velocity() * result;
            }
#elif NO_DIM==3
        for (int i1=0; i1<cellCount[0]; ++i1)
            for (int i2=0; i2<cellCount[1]; ++i2)
                for (int i3=0; i3<cellCount[2]; ++i3)
                {
                    int index = cell[0][i1] * nGrid[1]*nGrid[2] + cell[1][i2] * nGrid[2] + cell[2][i3];
                    Real result = (*it)->weight() * weight[0][i1] * weight[1][i2] * weight[2][i3];
                    q->density[index] += result;
                    q->velocity[index] += (*it)->velocity() * result;
                }
#endif
    }
    
    
    // divide the momentum by the mass in the cells
    if ( userOptions.aField.velocity )
    {
        for (size_t i=0; i<reserveSize; ++i)
            if ( q->density[i]!=Real(0.) )
                q->velocity[i] /= q->density[i];
            else
                q->velocity[i] = Pvector<Real,noVelComp>::zero();
    }
    else
        q->velocity.clear();
    
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
    printElapsedTime( &t, &userOptions, "TSC interpolation" );
}




