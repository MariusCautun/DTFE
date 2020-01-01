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


#include "volume_split.h"


/* Computes the contribution to the input variable for the grid cells given in 'indices'.  */
/* Each entry of the array 'contributions' contains the first 2/3 entries the numbers that need to be multiplied with the gradient to give you the field while the last entry (NO_DIM) contains the intersection volume between the simplex and the grid cell. This needs to be multiplied with the value of the field corresponding to the base vertex of the Delaunay cell (e.g. density_base). */
inline void distributeContributions(vector<size_t> &indices,
                                    vector< Pvector<Real,NO_DIM+1> > &contributions,
                                    std::vector<Real> *density,
                                    Real const base_density,
                                    Real const densityGrad[NO_DIM],
                                    Vertex_handle &base)
{
    for (size_t i=0; i<indices.size(); ++i)
    {
        size_t index = indices[i];
        if ( index==size_t(-1) )
            continue;
        Real volume = contributions[i][NO_DIM];
        Real result = base_density * volume;
        for (int j=0; j<NO_DIM; ++j)
            result += densityGrad[j] * (contributions[i][j] - volume*base->point()[j]);
        (*density)[index] += result;
    }
}
template <size_t N1, size_t N2>
inline void distributeContributions(vector<size_t> &indices,
                                    vector< Pvector<Real,NO_DIM+1> > &contributions,
                                    std::vector< Pvector<Real,N1> > *field,
                                    Pvector<Real,N1> field_base,
                                    Real fieldGradient[NO_DIM][N2],
                                    Vertex_handle &base)
{
    for (size_t i=0; i<indices.size(); ++i)
    {
        size_t index = indices[i];
        if ( index==size_t(-1) )
            continue;
        Real volume = contributions[i][NO_DIM];
        Pvector<Real,N1> result = field_base * volume;
        for (int j=0; j<NO_DIM; ++j)
            for (int j2=0; j2<N1; ++j2)
                result[j2] += fieldGradient[j][j2] * (contributions[i][j] - volume*base->point()[j]);
        (*field)[index] += result;
    }
}
// the following function computes the average gradient - this is contant inside the Delaunay cell
template <size_t N1>
inline void distributeContributions(vector<size_t> &indices,
                                    vector< Pvector<Real,NO_DIM+1> > &contributions,
                                    std::vector< Pvector<Real,N1> > *fieldGrad,
                                    Pvector<Real,N1> fieldGradient)
{
    for (size_t i=0; i<indices.size(); ++i)
    {
        size_t index = indices[i];
        if ( index==size_t(-1) )
            continue;
        Real volume = contributions[i][NO_DIM];
        (*fieldGrad)[index] += fieldGradient * volume;
    }
}




/* Returns true if the Delaunay cell is outside the region of interest. */
inline bool cellOutsideRegion(DT & dt,
                              Finite_cells_iterator &cell,
                              Bbox &fullBox)
{
#if NO_DIM==2
    return not do_overlap( dt.triangle( cell ).bbox(), fullBox );
#elif NO_DIM==3
    return not do_overlap( dt.tetrahedron( cell ).bbox(), fullBox );
#endif
}

bool simplexInSingleCell(Finite_cells_iterator &cell,
                         Box &regionBox,
                         Real *dx,
                         int *baseGridCell,
                         Real *basePosition)
{
    // check if Delaunay cell is fully contained in the same grid cell
    Point base = cell->vertex(0)->point();
    for (int i=0; i<NO_DIM; ++i)    // get base point positions in the grid cell it is contained
    {
        basePosition[i] = base[i] - regionBox[2*i]; //position of the base vertex
        baseGridCell[i] = int( floor( (basePosition[i])/dx[i] ) ); //coordinates of grid cell that contains the base point
    }
    Box gridCellBox;    // keep track of coordinates of the grid cell in which 'base' is
    for (int i=0; i<NO_DIM; ++i)
    {
        gridCellBox[2*i] = baseGridCell[i] * dx[i] + regionBox[2*i];
        gridCellBox[2*i+1] = gridCellBox[2*i] + dx[i];
    }
    for (int i=1; i<NO_DIM+1; ++i)
        if ( not gridCellBox.isPointInBox( cell->vertex(i)->point() ) )
            return false;
    return true;
}





/* This function computes the field interpolation averaged over the sampling cell for each sampling point using a MC method over the cells of the Delaunay triangulation.
The function takes the following parameters:
    'dt' - the periodic Delaunay triangulation; with the density at the vertex points already computed
    'userOptions' - structure which store the number of grid points along each axis, the number of random sample points in each tetrahedra. This structure also gives information about the size of the box along each axis. 
    'quantities' - structure storing the vectors for the output fields on the grid
*/
void interpolateGrid_averaged_0(DT &dt,
                                User_options &userOptions,
                                Quantities *quantities)
{
    size_t *nGrid = &(userOptions.gridSize[0]); // size of the density grid along each axis
    Box boxCoordinates = userOptions.region;    // the box coordinates in which the density is computed
    vector<Real> boxLength, startPos;           // the size of the particle box of interest along each direction as well as the origin of the box
    for (size_t i=0; i<NO_DIM; ++i)
    {
        boxLength.push_back( boxCoordinates[2*i+1]-boxCoordinates[2*i] );
        startPos.push_back( boxCoordinates[2*i] );
    }
#if NO_DIM==2
    Bbox fullBox( boxCoordinates[0], boxCoordinates[2], boxCoordinates[1], boxCoordinates[3] );
#elif NO_DIM==3
    Bbox fullBox( boxCoordinates[0], boxCoordinates[2], boxCoordinates[4], boxCoordinates[1], boxCoordinates[3], boxCoordinates[5] );
#endif
    
    
    // define pointers to the vectors storing the output fields
    Field field = userOptions.aField;
    std::vector<Real>                         *density = &(quantities->density);
    std::vector< Pvector<Real,noVelComp> >    *velocity = &(quantities->velocity);
    std::vector< Pvector<Real,noGradComp> >   *velocity_gradient = &(quantities->velocity_gradient);
    std::vector< Pvector<Real,noScalarComp> > *scalar = &(quantities->scalar);
    std::vector< Pvector<Real,noScalarGradComp> > *scalar_gradient = &(quantities->scalar_gradient);
    
    
    // show a message to the user
    MESSAGE::Message message( userOptions.verboseLevel );   // structure used to output messages to the user ( the higher the verboseLevel, the more messages will be showns (check 'message.h' for details) )
    message << "\nComputing interpolation of the fields volume averaged over the sampling cell on a regular " 
            << MESSAGE::printElements( nGrid, NO_DIM, "*" ) << " grid in the region " 
            << boxCoordinates.print() 
            << ". The volume average is computed using geometric intersections up to a precision of " << VOLUME_TOL << ".\n" 
            << "\t Done: " << MESSAGE::Flush;
    
    
    // Initialize the grid field values to 0
    size_t const gridSize = (NO_DIM==2) ? nGrid[0]*nGrid[1] : nGrid[0]*nGrid[1]*nGrid[2];
    if ( field.density ) assingZeroValues<Real>( density, gridSize );
    if ( field.velocity ) assingZeroValues< Pvector<Real,noVelComp> >( velocity, gridSize );
    if ( field.velocity_gradient ) assingZeroValues< Pvector<Real,noGradComp> >( velocity_gradient, gridSize );
    if ( field.scalar ) assingZeroValues< Pvector<Real,noScalarComp> >( scalar, gridSize );
    if ( field.scalar_gradient ) assingZeroValues< Pvector<Real,noScalarGradComp> >( scalar_gradient, gridSize );
    if ( field.velocity_std )
    {
        MESSAGE::Warning warning(userOptions.verboseLevel);
        warning << "You cannot use the averaging method 0 to compute the velocity standard deviation inside the grid cell. Please use the averaging method 2 '--method 2' (Monte Carlo sampling inside the grid cell) to compute the velocity standard deviation.\n" << MESSAGE::EndWarning;
    }
    //checks that there are cells in the Delaunay triangulation, since otherwise there is no interpolation
    if ( dt.number_of_vertices()<NO_DIM+1 )
    {
        MESSAGE::Warning warning(1);
        warning << "Because there are less than " << NO_DIM+1 << " vertices in the Delaunay triangulation there is no cell and hence there is no information that can be used to interpolate the fields to a grid. All field values will be initialized to 0.\n" << MESSAGE::EndWarning;
        return;
    }
    
    
    
    // temporary variables needed for the interpolation
    Real dx[NO_DIM];        // grid spacing along each axis
    for (size_t i=0; i<NO_DIM; ++i)
        dx[i] = boxLength[i] / nGrid[i];
    Real const gridCellVolume = (NO_DIM==2) ? (dx[0]*dx[1]) : (dx[0]*dx[1]*dx[2]);    // the volume of a grid cell
#ifdef TEST_PADDING
    vector<size_t> incompleteCells_d;// keep track of grid cells where there is an error in the density field computation (because one of the vertices is a dummy point)
    vector<size_t> incompleteCells;  // keep track of grid cells where there is an error in the field (all fields except density) computation (because one of the vertices is a dummy point)
#endif
    
    
    // initialize the class that does the geometrical intersection by splitting the volume into triangles/tetrahedra
    VOLUME_SPLIT::VolumeSplit volSplit( &(startPos[0]), nGrid, dx );
    VOLUME_SPLIT::Simplex     simplex;      // will store the simplex: triangle or tetrahedron
    vector<size_t>            indices;      // stores the grid cell indices of the large grid that intersect the simplex
    vector< Pvector<Real,NO_DIM+1> > contributions;    //stores the contribution of the simplex to each of the grid cells in indices
    
    
    // iterate over all the cells of the Delaunay triangulation
    size_t prev = 0, amount100 = 0, count = 0;    // variable to show the user about the progress of the computation
#if NO_DIM==2
    size_t noTotalCells = dt.number_of_faces();
    for(Finite_cells_iterator itC = dt.finite_faces_begin(); itC!= dt.finite_faces_end(); ++itC)
#elif NO_DIM==3
    size_t noTotalCells = dt.number_of_finite_cells();
    for(Finite_cells_iterator itC = dt.finite_cells_begin(); itC!= dt.finite_cells_end(); ++itC)
#endif
    {
        // show the progress of the computation
        amount100 = (100 * count++)/ noTotalCells;
        if (prev < amount100)
            message.updateProgress( ++prev );
        
        
        // if tetrahedron outside the box of interest, continue to next cell
        if ( cellOutsideRegion( dt, itC, fullBox) )
            continue;
        
        
        // compute the inverse position matrix of the Delaunay cell
        Real posMatrixInverse[NO_DIM][NO_DIM];  //stores the inverse of the vertex position difference matrix
        positionMatrix( itC, posMatrixInverse );    //returns 'posMatrixInverse'
        Vertex_handle base = itC->vertex(0);    // stores the "base" vertex of the Delaunay cell
        
        
        // using the volume splitting procedure to get the intersection of the tetrahedron with the grid
        simplex.assign( itC );
        volSplit.findIntersection( simplex, indices, contributions );
        
        // compute the gradients inside the simplex
        Real densGrad[NO_DIM];
        densityGrad( itC, posMatrixInverse, densGrad );
        if (field.density) distributeContributions( indices, contributions, density, base->info().density(), densGrad, base );
#ifdef VELOCITY
        Real velGrad[NO_DIM][noVelComp];
        velocityGrad( itC, posMatrixInverse, velGrad );
        if (field.velocity) distributeContributions( indices, contributions, velocity, base->info().velocity(), velGrad, base );
        if (field.velocity_gradient) distributeContributions( indices, contributions, velocity_gradient, velocityGradient(velGrad) );
#endif
#ifdef SCALAR
        Real sGrad[NO_DIM][noScalarComp];
        scalarGrad( itC, posMatrixInverse, sGrad );
        if (field.scalar) distributeContributions( indices, contributions, scalar, base->info().myScalar(), sGrad, base );
        if (field.scalar_gradient) distributeContributions( indices, contributions, scalar_gradient, scalarGradient(sGrad) );
#endif
        
        
        // check if any of the cell vertices have dummy points as neighbors
#ifdef TEST_PADDING
        bool dummyNeighbors = hasDummyNeighbor( itC );
        bool dummyVertices = hasDummyVertex( itC );
        if ( dummyNeighbors ) updateDummyGridCells( indices, &incompleteCells_d );
        if ( dummyVertices ) updateDummyGridCells( indices, &incompleteCells );
#endif
    }
    
    
    
    // divide the values in each cell by the volume of each sampling cell to get the volume averaged field value at the sample point
    if (field.density)
        for (vector<Real>::iterator it=density->begin(); it!=density->end(); ++it)
            *it /= gridCellVolume;
#ifdef VELOCITY
    if (field.velocity)
        for (vector< Pvector<Real,noVelComp> >::iterator it=velocity->begin(); it!=velocity->end(); ++it)
            *it /= gridCellVolume;
    if (field.velocity_gradient)
        for (vector< Pvector<Real,noGradComp> >::iterator it=velocity_gradient->begin(); it!=velocity_gradient->end(); ++it)
            *it /= gridCellVolume;
#endif
#ifdef SCALAR
    if (field.scalar)
        for (vector< Pvector<Real,noScalarComp> >::iterator it=scalar->begin(); it!=scalar->end(); ++it)
            *it /= gridCellVolume;
    if (field.scalar_gradient)
        for (vector< Pvector<Real,noScalarGradComp> >::iterator it=scalar_gradient->begin(); it!=scalar_gradient->end(); ++it)
            *it /= gridCellVolume;
#endif
    
    message << "100%\n" << MESSAGE::Flush;
    
    
#ifdef TEST_PADDING
    if (field.density)
        showCellsContainingDummyPoints( &incompleteCells_d, userOptions, userOptions.outputFilename+"_density", "density" ); // check if there are any cells which may have an error in thedensity estimation due to an incomplete Delaunay tesselation over the region of interest
    if ( field.velocity or field.velocity_gradient or field.scalar or field.scalar_gradient )
        showCellsContainingDummyPoints( &incompleteCells, userOptions, userOptions.outputFilename+"_fields", "fields" ); // check if there are any cells which may have an error in the fields estimation (all fields except density) due to an incomplete Delaunay tesselation over the region of interest
#endif
}



