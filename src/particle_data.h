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
  This header contains 3 classes:
        Data_structure - contains all the point/particle data with the exception of the position
        Particle_data - contains the point/particle position plus the data fields given in "Data_structure"
        Sample_point - contains a class to store sample points in case the user wants to interpolate the field at his choice of grid points
See below before each class for additional information.

NOTE: to add additional point/particle properties just modify the variable "noScalarComp".
*/


#ifndef PARTICLE_DATA_HEADER
#define PARTICLE_DATA_HEADER

#include "define.h"
#include "Pvector.h"




/* This class contains all the point/particle properties with the exception of the position vector.
NOTE: The "particle data" was split in two such that the point properties (with the exception of the position) can be easily transferred to the vertices of the Delaunay triangulation. */
class Data_structure
{
    protected:
        Real                       _weight; // weight of each particle
        Real                       _density;// the value of the density will be computed later using the DTFE density interpolation at each vertex position
#ifdef VELOCITY
    Pvector<Real,noVelComp>    _velocity;// stores particle's velocity NOTE: "Pvector" is a 1D array which behaves like a physical vector
#else
    static Pvector<Real,noVelComp> _velocity;//this will be defined only once for all 'Data_structure' classes - so it takes no memory at all
#endif
#ifdef SCALAR
    Pvector<Real,noScalarComp> _scalar; // this variable will store additional point/particle data - the size of these data can be change via the parameter "noScalarComp"
#else
    static Pvector<Real,noScalarComp> _scalar;//this will be defined only once for all 'Data_structure' classes - so it takes no memory at all
#endif
    
    public:
    Data_structure()
    {
        _weight = Real(1.);
        _density = Real(0.);
#ifdef VELOCITY
        _velocity = Pvector<Real,noVelComp>::zero();
#endif
#ifdef SCALAR
        _scalar = Pvector<Real,noScalarComp>::zero();
#endif
    }

    // functions to access the weight
    inline Real& weight() { return _weight;}
    inline void setWeight(Real const w) { _weight = w;}
    // functions to access the density
    inline Real& density() { return _density;}
    inline void setDensity(Real const d) { _density = d;}
    // functions to access the velocity and its components
    inline Pvector<Real,noVelComp>& velocity() { return _velocity;}
    inline Real& velocity(int const i) { return _velocity[i];}
    inline void setVelocity(Real * vel) { for(size_t j=0; j<noVelComp; ++j) _velocity[j] = vel[j];}
    // functions to access the scalar data
    inline Pvector<Real,noScalarComp>& scalar() { return _scalar;}
    inline Real& scalar(int const i) { return _scalar[i];}
    inline void setScalar(Real * s) { for(size_t j=0; j<noScalarComp; ++j) _scalar[j] = s[j];}
};



/* This is the class that will be used to transfer the point/particle's position and data between the different functions (units) of the program.
NOTE: Do not modify the name of the position variable from 'pos' since the program will not compile. */
struct Particle_data : public Data_structure
{
    Pvector<Real,NO_DIM>   pos;    // stores particle's position
    
    inline Pvector<Real,NO_DIM>& position() { return pos;}
    inline Real& position(int const i) { return pos[i];}
    inline void setPosition(Real *p) { for(size_t j=0; j<NO_DIM; ++j) pos[j] = p[j];}
};



/* In case the user whishes to know the interpolated fields as specific locations that do not follow the patern available with the code (spatial grid and light cone grid), than he can supply the positions of the sample points that he is interested of. The two components are:
        pos - the coordinates of the user-specified sample points
        delta - the grid size along each dimension associated to the given sample point (this is needed only for density interpolation since one needs to average over the full grid cell to filter the Poisson noise)
NOTE: For the case of density interpolation the user needs to supply also values for 'delta', otherwise this is not required. */
struct Sample_point
{
    Pvector<Real,NO_DIM>    pos;   // stores the sampling point position for non-uniform grids
    Pvector<Real,NO_DIM>    _delta; // stores the size of the grid cell along each direction for the given sampling point - only for non-unifrom grids
    
    inline Pvector<Real,NO_DIM>& position() { return pos;}
    inline Real& position(int const i) { return pos[i];}
    inline void setPosition(Real *p) { for(size_t j=0; j<NO_DIM; ++j) pos[j] = p[j];}
    inline Pvector<Real,NO_DIM>& delta() { return _delta;}
    inline Real& delta(int const i) { return _delta[i];}
    inline void setDelta(Real *p) { for(size_t j=0; j<NO_DIM; ++j) _delta[j] = p[j];}
};


#endif
