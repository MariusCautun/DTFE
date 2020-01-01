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


#ifndef DTFE_HEADER
#define DTFE_HEADER

#include <vector>

#include "define.h"
#include "user_options.h"
#include "particle_data.h"
#include "quantities.h"


// Computes the grid interpolation
void DTFE(std::vector<Particle_data> *allParticles,
          std::vector<Sample_point> &samples,
          User_options &userOptions,
          Quantities *uQuantities,
          Quantities *aQuantities);



// section used if you need to access the Delaunay triangulations
#ifdef TRIANGULATION
// include the headers
#if NO_DIM==2
#include "CGAL_triangulation/CGAL_include_2D.h"
#elif NO_DIM==3
#include "CGAL_triangulation/CGAL_include_3D.h"
#endif
#include "math_functions.h"

// Computes the grid interpolation (no partition options available as well as no parallel threads)
void DTFE(std::vector<Particle_data> *allParticles,
          std::vector<Sample_point> &samples,
          User_options &userOptions,
          Quantities *uQuantities,
          Quantities *aQuantities
          DT &delaunay_triangulation);


#endif


#endif

