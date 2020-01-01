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


#ifndef VERTEX_DATA_HEADER
#define VERTEX_DATA_HEADER

#include "../define.h"
#include "../particle_data.h"


/* This class stores the properties of each delaunay triangulation vertex. These are the properties coming from the point/particle properties that the vertices represent (see the "particle_data.h" file) plus two additional bool variables. */
struct vertexData : public Data_structure
{
    protected:
    bool   dummy;         // true if the vertex is dummy test point - used to test padding efficiency
    bool   dummyNeighbor; // true if the vertex has at least one dummy point as neighbor - used to test padding efficiency for density computations
    
    
    public:
    // constructor
    vertexData(){ dummy=false; dummyNeighbor=false; }
    
    // Customizable function definition for the scalar data
    inline Pvector<Real,noScalarComp> myScalar()
    {
        return scalar();
    }
    
    // copy the point/particle data from "Particle_data" to this class
    inline void setData(Particle_data &other)
    {
        weight() = other.weight();
        density() = other.density();
#ifdef VELOCITY
        velocity() = other.velocity();
#endif
#ifdef SCALAR
        scalar() = other.scalar();
#endif
    }
    //! see the "particle_data.h" file for the rest of functions that can access and modify the values in 'Data_structure'
    
    
    // do not modify the following
    inline void setDummy() { dummy=true; dummyNeighbor=true; setDensity(0.); }
    inline void setDummyNeighbor() { dummyNeighbor=true; }
    inline bool isDummy() { return dummy; }
    inline bool hasDummyNeighbor() { return dummyNeighbor; }
};

#endif
