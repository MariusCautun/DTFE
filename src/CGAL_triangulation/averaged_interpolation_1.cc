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


/* This file contains most of the functions necessary to interpolate to grid volume averaged fields inside the sampling cell using averaging method 1 (= choose sampling points inside the Delaunay cell). */
#include <typeinfo>
#include <gsl/gsl_qrng.h>



/* Stores a sequence of quasi-random numbers using the GSL quasi-number generator. */
void quasiRandomSequence(Real quasiRandomNumbers[][NO_DIM],
                         size_t const arraySize)
{
    gsl_qrng * q = gsl_qrng_alloc( gsl_qrng_sobol, NO_DIM );
//     gsl_qrng * q = gsl_qrng_alloc( gsl_qrng_niederreiter_2, NO_DIM );
    
    double randomRoot[NO_DIM];
    for (size_t i=0; i<arraySize; )
    {
        gsl_qrng_get (q, randomRoot);
        Real sum = NO_DIM==2 ? randomRoot[0]+randomRoot[1] : randomRoot[0]+randomRoot[1]+randomRoot[2];
        if ( sum<=1. )
        {
            for (int j=0; j<NO_DIM; ++j)
                quasiRandomNumbers[i][j] = Real( randomRoot[j] );
            ++i;
        }
    }
    gsl_qrng_free (q);
}
/* Generates quasi-random sample points inside the Delaunay cell: traingle in 2D and tetrahedron in 3D. */
void quasiRandomPointsInCell(double vM[][NO_DIM],
                             size_t noRandomPoints,
                             Real qR[][NO_DIM],
                             Point *randomPoints)
{
#if NO_DIM==2
    for (size_t i=0; i<noRandomPoints; ++i)
        randomPoints[i] = Point( qR[i][0]*vM[0][0]+qR[i][1]*vM[1][0], qR[i][0]*vM[0][1]+qR[i][1]*vM[1][1] );
#elif NO_DIM==3
    for (size_t i=0; i<noRandomPoints; ++i)
        randomPoints[i] = Point( qR[i][0]*vM[0][0]+qR[i][1]*vM[1][0]+qR[i][2]*vM[2][0], qR[i][0]*vM[0][1]+qR[i][1]*vM[1][1]+qR[i][2]*vM[2][1], qR[i][0]*vM[0][2]+qR[i][1]*vM[1][2]+qR[i][2]*vM[2][2] );
#endif
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

static int const SINGLE_GRID_CELL = 0;
static int const MULTIPLE_GRID_CELLS = 1;
/* Checks the position of a Delaunay Cell and computes the vertices position matrix 'vertexMatrix'. 
NOTE: The vertex matrix position is stored as a multiple of the grid spacing along the given coordinate. */
int checkCellPosition(Finite_cells_iterator &cell,
                      Box &regionBox,
                      Real *dx,
                      double vertexMatrix[][NO_DIM],
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
    bool insideOneCell = true;
    for (int i=1; i<NO_DIM+1; ++i)
        if ( not gridCellBox.isPointInBox( cell->vertex(i)->point() ) )
        {
            insideOneCell = false;
            break;
        }


    // Now store in 'vertexMatrix' the vertices position with respect to the 'base' vertex
    for (int v = 0; v<NO_DIM; ++v)  //loop over vertices != base
        for (int i=0; i<NO_DIM; ++i)    //loop over spatial dimensions
            vertexMatrix[v][i] = double(cell->vertex(v+1)->point()[i]) - double(base[i]);

    if ( insideOneCell )
        return SINGLE_GRID_CELL;
    return MULTIPLE_GRID_CELLS;
}



/* Checks if the given grid cell is valid. */
inline bool isGridCellIndex(int *gridCell,
                            size_t *totalGridSize,
                            size_t *cellIndex)
{
    // check if grid cell is valid
    for (int i=0; i<NO_DIM; ++i)
        if ( gridCell[i]<0 or gridCell[i]>=int(totalGridSize[i]) )
            return false;
    
    *cellIndex = gridCellIndex( gridCell, totalGridSize );
    return true;
}
/* Finds the grid cell where the sample point lies and checks if it is a valid grid cell. */
inline bool isGridCellIndex(size_t *totalGridSize,
                            Real *dx,
                            Point &samplePoint,
                            Real *basePosition,
                            size_t *cellIndex)
{
    // compute the grid cell associated to the sample point and check if valid
    int pos[NO_DIM];
    for (int i=0; i<NO_DIM; ++i)
    {
        pos[i] = int(floor( (basePosition[i] + samplePoint[i])/dx[i] ));
        if ( pos[i]<0 or pos[i]>=int(totalGridSize[i]) )
            return false;
    }
    
    *cellIndex = gridCellIndex( pos, totalGridSize );
    return true;
}



/* Computes the average density in a Delaunay triangle/tetrahedron in 2D/3D. */
inline Real averageDensity(Finite_cells_iterator &cell)
{
    Real temp = 0.;
    for (int i=0; i<NO_DIM+1; ++i)
        temp += cell->vertex(i)->info().density();
    return temp / Real(NO_DIM+1);
}
/* Computes the average density in a Delaunay triangle/tetrahedron in 2D/3D. */
inline Pvector<Real,noVelComp> averageVelocity(Finite_cells_iterator &cell)
{
    Pvector<Real,noVelComp> temp = Pvector<Real,noVelComp>::zero();
    for (int i=0; i<NO_DIM+1; ++i)
        temp += cell->vertex(i)->info().velocity();
    return temp / Real(NO_DIM+1);
}
/* Computes the average density in a Delaunay triangle/tetrahedron in 2D/3D. */
inline Pvector<Real,noScalarComp> averageScalar(Finite_cells_iterator &cell)
{
    Pvector<Real,noScalarComp> temp = Pvector<Real,noScalarComp>::zero();
    for (int i=0; i<NO_DIM+1; ++i)
        temp += cell->vertex(i)->info().myScalar();
    return temp / Real(NO_DIM+1);
}





/* This function computes the field interpolation averaged over the sampling cell for each sampling point using a MC method over the cells of the Delaunay triangulation.
The function takes the following parameters:
    'dt' - the periodic Delaunay triangulation; with the density at the vertex points already computed
    'userOptions' - structure which store the number of grid points along each axis, the number of random sample points in each tetrahedra. This structure also gives information about the size of the box along each axis. 
    'quantities' - structure storing the vectors for the output fields on the grid
*/
void interpolateGrid_averaged_1(DT &dt,
                                User_options &userOptions,
                                Quantities *quantities)
{
    size_t *nGrid = &(userOptions.gridSize[0]); // size of the density grid along each axis
    size_t const NN = userOptions.noPoints;     // average number of random points per grid cell used to compute the density
    size_t const minNN = NN/6 + 1;              // minimum number of sample points in any given Delaunay cell
    Real const minRatio = 1./6.;                // minimum ratio of volume(Delaunay cell) / volume(grid cell) when to switch from 'minNN' sample points per Delaunay cell to 'minNN * volume ratio' sample points per Delaunay cell
    size_t const maxNN = 100*NN;                // maximum number of sample points per Delaunay cell
    Box boxCoordinates = userOptions.region;    // the box coordinates in which the density is computed
    vector<Real> boxLength;                     // the size of the particle box of interest along each direction
    for (size_t i=0; i<NO_DIM; ++i)
        boxLength.push_back( boxCoordinates[2*i+1]-boxCoordinates[2*i] );
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
            << ". The volume average is done using Monte Carlo sampling in each of the triangles/tetrahedra of the Delaunay triangulation. " 
            << "There are " << NN << " random sampling points on average in each grid cell of the output fields (the random samples in the Delaunay triangulation's triangles/tetrahedra are proportional to the area/volume of the cell).\n" 
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
        warning << "You cannot use the averaging method 1 (Monte Carlo sampling inside the Delaunay cell) to compute the velocity standard deviation inside the grid cell. Please use the averaging method 2 '--method 2' (Monte Carlo sampling inside the grid cell) to compute the velocity standard deviation.\n" << MESSAGE::EndWarning;
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
    
    
    // quasi-random sequence of numbers - using the GSL quasi-number generator
    Real quasiRandomNumbers[maxNN][NO_DIM];
    quasiRandomSequence( quasiRandomNumbers, maxNN );
    
    
    
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
        
        
        // check if any of the cell vertices have dummy points as neighbors
#ifdef TEST_PADDING
        bool dummyNeighbors = hasDummyNeighbor( itC );
        bool dummyVertices = hasDummyVertex( itC );
#endif
        
        
        double vertexMatrix[NO_DIM][NO_DIM]; // keep track of the relative differences between the vertices of a Delaunay cell
        int baseGridCell[NO_DIM]; // will store the grid cell where is the base point of the tetrahedron
        Real basePosition[NO_DIM];// will store the position of the base point of the tetrahedron with respect to the lower left corner of the box
        int cellPosition = checkCellPosition( itC, boxCoordinates, dx, vertexMatrix, baseGridCell, basePosition );  // this computes if the Delaunay cell is fully contained in a grid cell or not
        
        Real cellVolume = volume( dt, itC );
        Real posMatrixInverse[NO_DIM][NO_DIM];  // stores the inverse of the vertex position difference matrix
        matrixInverse( vertexMatrix, posMatrixInverse );  // takes the inverse of the 'vertexMatrix' matrix
        size_t index;
        
        // the tetrahedron is inside one grid cell - no need for MC sampling, can just add the total contribution
#ifndef MY_SCALAR
        if ( cellPosition==SINGLE_GRID_CELL and isGridCellIndex( baseGridCell, nGrid, &index ) )
        {
            if (field.density) (*density)[index] += averageDensity(itC) * cellVolume;
#ifdef VELOCITY
            if (field.velocity) (*velocity)[index] += averageVelocity(itC) * cellVolume;
            if (field.velocity_gradient)
            {
                Real velGrad[NO_DIM][noVelComp];
                velocityGrad( itC, posMatrixInverse, velGrad );
                (*velocity_gradient)[index] += velocityGradient(velGrad) * cellVolume;
            }
#endif
#ifdef SCALAR
            if (field.scalar) (*scalar)[index] += averageScalar(itC) * cellVolume;
            if (field.scalar_gradient)
            {
                Real sGrad[NO_DIM][noScalarComp];
                scalarGrad( itC, posMatrixInverse, sGrad );
                (*scalar_gradient)[index] += scalarGradient(sGrad) * cellVolume;
            }
#endif
#ifdef TEST_PADDING
            if ( dummyNeighbors ) updateDummyGridCells( index, &incompleteCells_d );
            if ( dummyVertices ) updateDummyGridCells( index, &incompleteCells );
#endif
            continue;
        }
#endif
        
        
        // the tetrahedron crosses multiple interpolation grid cells
        // get the quasi-random points inside the Delaunay cell
        size_t const tempInt = size_t(NN*cellVolume/gridCellVolume) + 1;
        size_t const noRandomPoints = (cellVolume/gridCellVolume>minRatio) ? (tempInt>maxNN ? maxNN:tempInt) : minNN;// number of random points
        Point randomPoints[noRandomPoints];
        quasiRandomPointsInCell( vertexMatrix, noRandomPoints, quasiRandomNumbers, randomPoints );   // get quasi-random points inside the Delaunay cell
        Real factor = cellVolume / noRandomPoints;  // volume associated with each random sample point
        Vertex_handle base = itC->vertex(0);  // stores the "base" vertex of the Delaunay cell
        
        Real densGrad[NO_DIM];
        densityGrad( itC, posMatrixInverse, densGrad );
#ifdef VELOCITY
        Real velGrad[NO_DIM][noVelComp];
        velocityGrad( itC, posMatrixInverse, velGrad );
#endif
#ifdef SCALAR
        Real sGrad[NO_DIM][noScalarComp];
        scalarGrad( itC, posMatrixInverse, sGrad );
#endif
        for (size_t i=0; i<noRandomPoints; ++i) // now for each sample point find the grid cell where the sample point lies and add the new field value to the field stored in that grid cell
        {
            // compute the grid cell associated to each sample point - if point outside region of interest, continue to next random sample point
            if ( not isGridCellIndex( nGrid, dx, randomPoints[i], basePosition, &index ) )
                continue;
            
            if (field.density) (*density)[index] += densityValue(densGrad,base,randomPoints[i]) * factor;
#ifdef VELOCITY
            if (field.velocity) (*velocity)[index] += velocityValue(velGrad,base,randomPoints[i]) * factor;
            if (field.velocity_gradient) (*velocity_gradient)[index] += velocityGradient(velGrad) * factor;
#endif
#ifdef SCALAR
#ifdef MY_SCALAR
            (*scalar)[index] +=  customScalar(itC, posMatrixInverse, base, randomPoints[i]) * factor;
#else
            if (field.scalar) (*scalar)[index] += scalarValue(sGrad,base,randomPoints[i]) * factor;
            if (field.scalar_gradient) (*scalar_gradient)[index] += scalarGradient(sGrad) * factor;
#endif
#endif
#ifdef TEST_PADDING
            if ( dummyNeighbors ) updateDummyGridCells( index, &incompleteCells_d );
            if ( dummyVertices ) updateDummyGridCells( index, &incompleteCells );
#endif
        }
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



