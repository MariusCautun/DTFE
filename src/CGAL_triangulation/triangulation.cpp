/*
 *  Copyright (c) 2011       Marius Cautun
 *                           Erwin Platen
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

#include "../define.h"
#if NO_DIM==2
    #include "CGAL_include_2D.h"
#elif NO_DIM==3
    #include "CGAL_include_3D.h"
#endif

#include "../user_options.h"
#include "../quantities.h"
#include "../particle_data.h"
#include "particle_data_traits.h"
#include "../box.h"
#include "../message.h"
#include "../miscellaneous.h"

using namespace std;
#include "triangulation_miscellaneous.cc"   // defines various function used to compute triangulation related quantities
#include "my_function.h"                    // function to define custom computations
#include "padding_test.cc"                  // defines functions used to test the triangulation completness in the region of interest
#include "unaveraged_interpolation.cc"      // defines function used to interpolate the fields at the sampling points
#include "averaged_interpolation_1.cc"      // defines function used to compute the volume average of the fields using averaging method 1
#include "averaged_interpolation_2.cc"      // defines function used to compute the volume average of the fields using averaging methods 2 and 3


void delaunayTriangulation(DT *dt,
                           vector<Particle_data> *p,
                           int const verboseLevel);

void vertexDensity(DT & dt,
                User_options &userOptions);



/* This function computes the grid density in a periodic cosmological box of size [0,1] on each axis. If a non-periodic box is used, than points can have coordinates outside the [0,1] box to give a full Delaunay triangulation cover of the box. 
The function arguments are:
        p - vector which stores the particle positions and properties
        samples - vector storing the user defined sampling points, if any
        userOptions - structure to keep track of user given options and parameters for the interpolation
        uQuantities - struture that will hold the output fields, fields interpolated at the sampling points
        aQuantities - struture that will hold the output volume averaged fields, fields volume averaged over the grid cell associated to each sampling point
        dt - the CGAL Delaunay triangulation structure that stores the full triangulation of the point set
*/
void DTFE_interpolation(vector<Particle_data> *p,
                        vector<Sample_point> &samples,
                        User_options &userOptions,
                        Quantities *uQuantities,
                        Quantities *aQuantities,
                        DT &dt)
{
    // check that some conditions are satified (to don't get errors later on after all the heavy work was done)
    intervalCheck( userOptions.method, 1, 3, "'userOptions.method' in function 'DTFE_interpolation' must have values from 1 to 3");
    if ( userOptions.method==3 )
        rootN( userOptions.noPoints, NO_DIM );
    
    
    // Compute the Delaunay triangulation
    Timer t;     // construct a timer
    t.start();
    delaunayTriangulation( &dt, p, userOptions.verboseLevel );   // compute the triangulation
    printComputationTime( &t, &userOptions, "triangulation" );
    p->clear();  // delete the vector storing the particle data - not needed anymore
    
    
    // insert dummy points to test the padding
#ifdef TEST_PADDING
    if ( userOptions.testPaddedBoundaries )
    {
        t.start();
        insertDummyTestParticles( dt, userOptions );
        printComputationTime( &t, &userOptions, "insertion of dummy points" );
    }
#endif
    
    
    // compute the density associated to each Delaunay triangulation vertex (to each particle)
    t.start();
    vertexDensity( dt, userOptions );
    printComputationTime( &t, &userOptions, "vertex density computation" );
    
    
    
    // interpolate the required fields
    if ( userOptions.uField.selected() )    // interpolate the fields at the sampling points
    {
        t.start();
        if ( not userOptions.redshiftConeOn and not userOptions.userDefinedSampling )   // interpolate on a regular grid
            interpolateGrid( dt, userOptions, uQuantities );
        else if ( userOptions.redshiftConeOn and not userOptions.userDefinedSampling )  // interpolate on a redshift cone grid
            interpolateRedshiftCone( dt, userOptions, uQuantities );
        else if ( userOptions.userDefinedSampling )
            interpolateUserSampling( dt, samples, userOptions, uQuantities );
        printComputationTime( &t, &userOptions, "interpolation to sampling points" );
    }
    
    if ( userOptions.aField.selected() )    // interpolate the fields volume averaged inside the grid cells
    {
        t.start();
        if ( not userOptions.redshiftConeOn and not userOptions.userDefinedSampling )   // interpolate on a regular grid
        {
            if ( userOptions.method==1 )
                interpolateGrid_averaged_1( dt, userOptions, aQuantities );
            else if ( userOptions.method==2 )
                interpolateGrid_averaged_2( dt, userOptions, aQuantities );
            else if ( userOptions.method==3 )
                interpolateGrid_averaged_3( dt, userOptions, aQuantities );
            else throwError( "Unknown averaging method '", userOptions.method, "' when volume averaging the fields on a regular grid." );
        }
        else if ( userOptions.redshiftConeOn and not userOptions.userDefinedSampling )  // interpolate on a redshift cone grid
        {
            if ( userOptions.method==2 )
                interpolateRedshiftCone_averaged_2( dt, userOptions, aQuantities );
            else throwError( "Unknown averaging method '", userOptions.method, "' when volume averaging the fields on a redshift cone grid. The only available redshift cone grid averaging method is '2'." );
        }
        else if ( userOptions.userDefinedSampling )
        {
            if ( userOptions.method==2 )
                interpolateUserSampling_averaged_2( dt, samples, userOptions, aQuantities );
            else throwError( "Unknown averaging method '", userOptions.method, "' when volume averaging the fields on a user defined grid. The only available user defined grid averaging method is '2'." );
        }
        printComputationTime( &t, &userOptions, "interpolation to sampling points" );
    }
}

// Accesses the above 'DTFE_interpolation', but without returning the Delaunay triangulation
void DTFE_interpolation(vector<Particle_data> *p,
                        vector<Sample_point> &samples,
                        User_options &userOptions,
                        Quantities *uQuantities,
                        Quantities *aQuantities)
{
    DT dt;                // the triangulation
    DTFE_interpolation( p, samples, userOptions, uQuantities, aQuantities, dt );
}




/* Constructs the Delaunay triangulation using a set of points. */
void delaunayTriangulation(DT *dt,
                           vector<Particle_data> *p,
                           int const verboseLevel)
{
    // sort the particles to be spatially close - increase speed of Delaunay triangulation 
    MESSAGE::Message message( verboseLevel );
    message << "\nSorting the points to be spatially close yet randomly distributed ... " << MESSAGE::Flush;
    CGAL::spatial_sort( p->begin(), p->end(), Particle_data_sort_traits() );
    message << "Done.\n";
    
    
    // construct the actual triangulation
    message << "Constructing the Delaunay triangulation.\n\t Done: " << MESSAGE::Flush;
    size_t prev = 0, amount100 = 0, count = 0; // variable to show the user about the progress of the computation
    size_t const noPoints = p->size();
    Vertex_handle vh;     // vertex handle - points to each inserted point
    for (vector<Particle_data>::iterator it=p->begin(); it!=p->end(); ++it)
    {
#if NO_DIM==2
        vh = dt->insert( Point(it->pos[0],it->pos[1]) );
#elif NO_DIM==3
        vh = dt->insert( Point(it->pos[0],it->pos[1],it->pos[2]), vh );
#endif
        vh->info().setData( *it );      // set vertex quantities
        
        // show the progress of the computation
        amount100 = (100 * count++)/ noPoints;
        if (prev < amount100)
            message.updateProgress( ++prev );
    }
    
    message << "100\%.\n"
        << "The triangulation has " << dt->number_of_vertices() << " points.\n" << MESSAGE::Flush;
}



/* The function computes the volume of all Delaunay tetrahedra for each vertex and than computes the density at each vertex as the inverse of the total tetrahedra volumes incident on the respective vertex. */
void vertexDensity(DT & dt,
                  User_options &userOptions)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "\nComputing the density at each particle position.\n\t Done: " << MESSAGE::Flush;
    if ( userOptions.averageDensity<=0. )
    {
        MESSAGE::Error error;
        error << "The member 'averageDensity' of class 'User_options' must be positive since it represents the average density. Error found in function 'vertexDensity'." << MESSAGE::EndError;
    }
    if ( dt.number_of_vertices()<NO_DIM+1 )
    {
        MESSAGE::Warning warning(1);
        warning << "Because there are less than " << NO_DIM+1 << " vertices in the Delaunay triangulation there is no cell and hence there is no information that can be used to compute the density associated to each vertex. All vertex density values will be initialized to 0.\n" << MESSAGE::EndWarning;
        return;
        for (DT::Finite_vertices_iterator vIT = dt.finite_vertices_begin(); vIT != dt.finite_vertices_end(); ++vIT )
            vIT->info().setDensity( 0. );
    }
    
    Real const factor = (NO_DIM+1.) / userOptions.averageDensity;  //factor to normalize the density to the average background density; NO_DIM+1 comes from the number of vertices of a Delaunay tetrahedra - needed for mass conservation
    size_t prev = 0, amount100 = 0, count = 0; // variable to show the user about the progress of the computation
    size_t const noVertices = dt.number_of_vertices();  // total number of finite vertices
    
    
    DT::Finite_vertices_iterator vIT;
    for (vIT = dt.finite_vertices_begin(); vIT != dt.finite_vertices_end(); ++vIT )
    {
#ifdef TEST_PADDING
        if ( vIT->info().isDummy() )    // this is a dummy test point
            continue;
#endif
        
        vector<Cell_handle> cells;
#if NO_DIM==2
        DT::Face_circulator fc = dt.incident_faces( vIT );
        cells.push_back( fc++ );
        for (; fc!=cells[0]; ++fc)
            cells.push_back( fc );
#elif NO_DIM==3
        dt.incident_cells( vIT, back_inserter(cells) ); //get all the cells which have 'vIT' as vertex
#endif
        
        Real vol = 0.;
        bool infinite_volume = false;
        
        // compute the volume of those cells
        for ( vector< Cell_handle >::const_iterator itC = cells.begin(); itC!=cells.end(); ++itC )
        {
            if ( !dt.is_infinite(*itC) )
            {
                vol += volume( dt, *itC );
#ifdef TEST_PADDING
                if( hasDummyVertex(*itC) )  // the cell has a dummy point
                    vIT->info().setDummyNeighbor();
#endif
            }
            else
            {
                vIT->info().setDummyNeighbor();
                infinite_volume = true;
                break;
            }
        }
        
        // compute the density
        if ( not infinite_volume )
            vIT->info().setDensity( (vIT->info().weight()) * factor /vol );
        else
            vIT->info().setDensity( 0. );
        
        
        // show the progress of the computation
        amount100 = (100 * count++)/ noVertices;
        if (prev < amount100)
            message.updateProgress( ++prev );
    }
    
    message << "100\%.\n" << MESSAGE::Flush;
}




