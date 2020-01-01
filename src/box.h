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



/*!
Defines a class that keeps track of the coordinates of a box.
The box coordinates are stored in a vector using the sequance:
        xMin, xMax, yMin, yMax, ...
*/

#ifndef BOX_HEADER
#define BOX_HEADER

#include <vector>
#include <sstream>
#include <string>

#include "define.h"
#include "message.h"
#include "miscellaneous.h"


struct Box
{
    std::vector<Real> coords;    //the coordinates of the box as: xMin, xMax, yMin, yMax, ...
    
    //! Class functions
    Box() {coords.assign( 2*NO_DIM, Real(0.) );}   // constructor
    
    // returns true if the particle is inside the box
    template <typename Particle> bool isParticleInBox(Particle &p) const
    {
#if NO_DIM==2
        if ( p.pos[0]>=coords[0] and p.pos[0]<=coords[1] and p.pos[1]>=coords[2] and p.pos[1]<=coords[3] )
            return true;
#elif NO_DIM==3
        if ( p.pos[0]>=coords[0] and p.pos[0]<=coords[1] and p.pos[1]>=coords[2] and p.pos[1]<=coords[3] and p.pos[2]>=coords[4] and p.pos[2]<=coords[5] )
            return true;
#endif
    return false;
    }
    // returns true if the point is inside the box
    template <typename Point> bool isPointInBox(Point &p) const
    {
#if NO_DIM==2
        if ( p[0]>=coords[0] and p[0]<=coords[1] and p[1]>=coords[2] and p[1]<=coords[3] )
            return true;
#elif NO_DIM==3
        if ( p[0]>=coords[0] and p[0]<=coords[1] and p[1]>=coords[2] and p[1]<=coords[3] and p[2]>=coords[4] and p[2]<=coords[5] )
            return true;
#endif
    return false;
    }
    
    // checks if second box is overlaping with this one
    bool isBoxOverlaping(Box const &other) const
    {
        for (size_t i=0; i<NO_DIM; ++i)
        {
            Real len = this->coords[2*i+1] - this->coords[2*i];
            Real xMin = other.coords[2*i] - this->coords[2*i];
            Real xMax = other.coords[2*i+1] - this->coords[2*i];
            if ( Real(0.)>=xMin or xMax>=Real(0.) )
                return true;
            if ( len>=xMin or xMax>=len )
                return true;
            if ( Real(0.)<=xMin or xMax<=len )
                return true;
        }
        return false;
    }
    
    // translates the box according to the offset
    void translate(Real *offset)
    {
        for (size_t i=0; i<NO_DIM; ++i)
        {
            coords[2*i] += offset[i];
            coords[2*i+1] += offset[i];
        }
    }
    
    // increase box size by padding[i] along each edge
    template <typename T>
    void addPadding(T padding)
    {
        for (size_t i=0; i<NO_DIM; ++i)
        {
            coords[2*i] -= padding[2*i];
            coords[2*i+1] += padding[2*i+1];
        }
    }
    
    // checks that this is a valid sub-box of a larger box (i.e. see exact requierments at function definition)
    void validSubBox(Box &mainBox, bool periodic) const
    {
        if (periodic)
            for (int i=0; i<NO_DIM; ++i)
            {
                if ( coords[2*i]<(2*mainBox[2*i]-mainBox[2*i+1]) )
                {
                    MESSAGE::Error error;
                    error << "The padded box selected by you must have left boundaries at most '-boxLength' along each axis from the main data box. Error in function 'Box::validSubBox'." << MESSAGE::EndError;
                }
                if ( coords[2*i+1]>(2*mainBox[2*i+1]-mainBox[2*i]) )
                {
                    MESSAGE::Error error;
                    error <<  "The box selected by you must have right boundaries at most '+boxLength' along each axis from the main data box. Error in function 'Box::validSubBox'." << MESSAGE::EndError;
                }
            }
        else
            for (int i=0; i<NO_DIM; ++i)
                if ( coords[2*i]<mainBox[2*i] or coords[2*i+1]>mainBox[2*i+1] )
                {
                    MESSAGE::Error error;
                    error <<  "The padded box selected by you must fit inside the main data box since the particle positions are NOT PERIODIC. Choose a new box inside the main box or choose periodic computations. Error in function 'Box::validSubBox'." << MESSAGE::EndError;
                }
    }
    
    // returns box volume
    Real volume() const
    {
        Real temp = 1.;
        for (size_t i=0; i<NO_DIM; ++i)
            temp *= coords[2*i+1] - coords[2*i];
        return temp;
    }
    
    // outputs the box boundaries
    std::string print() const
    {
        std::ostringstream buffer;
        buffer << "{";
        for (size_t i=0; i<NO_DIM; ++i)
            buffer << "[" << coords[2*i] << "," << coords[2*i+1] << "],";
        buffer << char(8) << "}";
        return buffer.str();
    }
    
     // overload the [] operator - returns access to vector entries directly
    Real& operator [](size_t const i)
    { return coords[i]; }
    
     // assigns a constant value to all the entries
    void assign(Real const value)
    { coords.assign( 2*NO_DIM, value ); }
    
    // returns the number of coordinates
    size_t size() const
    { return coords.size();}
    
    // check if there are specified coordinates for the box (i.e. if box !=[0,0] )
    bool isNullBox()
    {
        for (size_t i=0; i<coords.size(); ++i)
            if ( coords[i]!=Real(0.) )
                return false;
        return true;
    }
};

#endif
