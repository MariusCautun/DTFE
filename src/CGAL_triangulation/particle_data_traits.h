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


#ifndef PARTICLE_DATA_TRAITS_HEADER
#define PARTICLE_DATA_TRAITS_HEADER

#include "../define.h"
#include "../particle_data.h"



struct Particle_data_less_x
{
    bool operator()(Particle_data p, Particle_data q) const
    { return (p.pos[0] < q.pos[0]); }
};
struct Particle_data_less_y
{
    bool operator()(Particle_data p, Particle_data q) const
    { return (p.pos[1] < q.pos[1]); }
};
#if NO_DIM==2
struct Particle_data_sort_traits
{
    typedef Particle_data Point_2;
    typedef Particle_data_less_x Less_x_2;
    typedef Particle_data_less_y Less_y_2;
    
    Less_x_2 less_x_2_object() const
    { return Less_x_2(); }
    Less_y_2 less_y_2_object() const
    { return Less_y_2(); }
    Particle_data_sort_traits() {}
};

#elif NO_DIM==3
struct Particle_data_less_z
{
    bool operator()(Particle_data p, Particle_data q) const
    { return p.pos[2]<q.pos[2]; }
};

struct Particle_data_sort_traits
{
    typedef Particle_data Point_3;
    typedef Particle_data_less_x Less_x_3;
    typedef Particle_data_less_y Less_y_3;
    typedef Particle_data_less_z Less_z_3;
    
    Less_x_3 less_x_3_object() const
    { return Less_x_3(); }
    Less_y_3 less_y_3_object() const
    { return Less_y_3(); }
    Less_z_3 less_z_3_object() const
    { return Less_z_3(); }
    Particle_data_sort_traits() {}
};
#endif


// Compares two particles to be able to sort them according to positions
template <typename T>
bool compareParticles(T p1, T p2)
{
#if NO_DIM==2
    if ( p1.pos[0]!=p2.pos[0] )
        return p1.pos[0]<p2.pos[0];
    else
        return p1.pos[1]<p2.pos[1];
#elif NO_DIM==3
    if ( p1.pos[0]!=p2.pos[0] )
        return p1.pos[0]<p2.pos[0];
    else if ( p1.pos[1]!=p2.pos[1] )
        return p1.pos[1]<p2.pos[1];
    else
        return p1.pos[2]<p2.pos[2];
#endif
}

template <typename T>
bool sameParticle(T p1, T p2)
{
#if NO_DIM==2
    return p1.pos[0]==p2.pos[0] and p1.pos[1]==p2.pos[1];
#elif NO_DIM==3
    return p1.pos[0]==p2.pos[0] and p1.pos[1]==p2.pos[1] and p1.pos[2]==p2.pos[2];
#endif
}

#endif
