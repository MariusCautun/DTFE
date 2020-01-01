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



/*
  This header contains the class which keeps track of all the possible different quantities that the DTFE interpolation to grid can compute:
        density             - a 'Real' value
        velocity            - a 'Pvector<Real,noVelComp>' value
        velocity_gradient   - a 'Pvector<Real,noGradComp>' value
        velocity_divergence - a 'Real' value
        velocity_shear      - a 'Pvector<Real,noShearComp>' value
        velocity_vorticity  - a 'Pvector<Real,noShearComp>' value
        scalar              - a 'Pvector<Real,noScalarComp>' value
        scalar_gradient     - a 'Pvector<Real,noScalarGradComp>' value

NOTE: After computation the quantities will be stored in a vector, hence memory will be allocated only to the quantities that will be computed and not to all possible variables.
NOTE 2: 'Pvector' stands for 'physical vector' and can be access as a normal array using the '[]' operator which takes values from 0 to N-1 (with N the number of components). Check the "Pvector.h" file for additional details and to see the matematical operations involving 'Pvector's.
NOTE 3: The program constants 'noVelComp', 'noGradComp', etc ... are defined in  "define.h".
*/



#ifndef QUANTITIES_HEADER
#define QUANTITIES_HEADER

#include <vector>

#include "define.h"
#include "particle_data.h"
#include "user_options.h"




/* This class keeps track of all the possible outcomes of the grid interpolation using the DTFE method. The results are stored in vectors, so memory is allocated only to the quantities that will be computed. */
struct Quantities
{
    std::vector<Real>                         density;              // vector that stores the density interpolated to grid using the DTFE method
    std::vector< Pvector<Real,noVelComp> >    velocity;             // vector that stores the velocity interpolated to grid using the DTFE method
    std::vector< Pvector<Real,noGradComp> >   velocity_gradient;    // vector that stores the velocity gradient map using the DTFE method
    std::vector<Real>                         velocity_divergence;  // vector that stores the velocity divergence map using the DTFE method
    std::vector< Pvector<Real,noShearComp> >  velocity_shear;       // vector that stores the velocity shear map using the DTFE method
    std::vector< Pvector<Real,noVortComp> >   velocity_vorticity;   // vector that stores the velocity vorticity map using the DTFE method
    std::vector<Real>                         velocity_std;         // vector that stores the velocity standard deviation map using the DTFE method
    std::vector< Pvector<Real,noScalarComp> > scalar;               // vector that stores a scalar field interpolated to grid using the DTFE method
    std::vector< Pvector<Real,noScalarGradComp> > scalar_gradient;  // vector that stores the gradient of the scalar field map using the DTFE method
        
    //Functions - you need to modify the below function if you add aditional members to this class ( - this is the case to be able to use the 'partition' option) 
    void copyFromSubgrid(Quantities const &subgridResults,
                        Field const &field,
                        std::vector<size_t> const &mainGrid,
                        std::vector<size_t> const &subgrid,
                        std::vector<size_t> const &subgridOffset); // used to copy the results obtain on a subgrid using the option 'partion' to the results for the full grid
    size_t size() const;	// returns the size of any non-empty object
    void reserveMemory(size_t *gridSize, Field &field); // reserve memory for the main grid quantities when using the 'partition' option
};


#endif
