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


/* This file contains most of the functions necessary to interpolate to grid volume averaged fields inside the sampling cell using averaging method 2 (= choose sampling points inside grid cell).
 NOTE: it also contains the function the interpolate using averaging method 3 - see end of the file.*/



/* This function interpolates the fields to a uniform grid - it computes the volume averaged value of the fields within the grid cell. To compute the volume average in each grid cell it uses a Monte Carlo algorithm to sample 'noRandPoints' in each grid cell and take the average.
The arguments are:
    'dt' - the periodic Delaunay triangulation; with the density at the vertex points already computed
    'userOptions' - structure which store the number of grid points along each axis and the number of random sample points in each grid cell. This structure also gives information about the size of the box along each axis.
    'quantities' - structure storing the vectors for the output fields on the grid
*/
void interpolateGrid_averaged_2(DT &dt,
                                User_options &userOptions,
                                Quantities *quantities)
{
    size_t *nGrid = &(userOptions.gridSize[0]); // size of the density grid along each axis
    int const NN = userOptions.noPoints;        // number of random points in each grid cell used to compute the volume average
    Box boxCoordinates = userOptions.region;    // the box coordinates in which the density is computed
    vector<Real> boxLength;                     // the size of the particle box of interest along each direction
    for (size_t i=0; i<NO_DIM; ++i)
        boxLength.push_back( boxCoordinates[2*i+1]-boxCoordinates[2*i] );
    
    
    // define pointers to the vectors storing the output fields
    Field field = userOptions.aField;
    std::vector<Real>                         *density = &(quantities->density);
    std::vector< Pvector<Real,noVelComp> >    *velocity = &(quantities->velocity);
    std::vector< Pvector<Real,noGradComp> >   *velocity_gradient = &(quantities->velocity_gradient);
    std::vector<Real>                         *velocity_std = &(quantities->velocity_std);
    std::vector< Pvector<Real,noScalarComp> > *scalar = &(quantities->scalar);
    std::vector< Pvector<Real,noScalarGradComp> > *scalar_gradient = &(quantities->scalar_gradient);
    
    
    // show a message to the user
    MESSAGE::Message message( userOptions.verboseLevel );   // structure used to output messages to the user ( the higher the verboseLevel, the more messages will be showns (check 'message.h' for details) )
    message << "\nComputing interpolation of the fields volume averaged over the sampling cell on a regular " 
            << MESSAGE::printElements( nGrid, NO_DIM, "*" ) << " grid in the region " 
            << boxCoordinates.print() 
            << ". The volume average is done using " << NN << " Monte Carlo samples in each grid cell.\n"
            << "\t Done: " << MESSAGE::Flush;
    
    
    //reserve memory for the output fields
    size_t reserveSize = (NO_DIM==2) ? nGrid[0]*nGrid[1] : nGrid[0]*nGrid[1]*nGrid[2];      // the number of cell in the grid
    reserveMemory( density, reserveSize, field.density, dt, "density" );    // reserves memory for the density computations (if any) and initializes the density values to 0 if the Delaunay triangulation is not well defined
    reserveMemory( velocity, reserveSize, field.velocity, dt, "velocity" );
    reserveMemory( velocity_gradient, reserveSize, field.velocity_gradient, dt, "velocity gradient" );
    reserveMemory( velocity_std, reserveSize, field.velocity_std, dt, "velocity standard deviation" );
    reserveMemory( scalar, reserveSize, field.scalar, dt, "scalar" );
    reserveMemory( scalar_gradient, reserveSize, field.scalar_gradient, dt, "scalar gradient" );
    if ( dt.number_of_vertices()<NO_DIM+1 ) return;     // stop the computation if the Delaunay triangulation is not well defined
    
    
    // temporary variables needed for the interpolation
    Point samplePoint;      // the coordinates of the MC sampling point
    Real dx[NO_DIM];        // grid spacing along each axis
    for (size_t i=0; i<NO_DIM; ++i)
        dx[i] = boxLength[i] / nGrid[i];
    
#ifdef TEST_PADDING
    vector<size_t> incompleteCells_d;// keep track of grid cells where there is an error in the density field computation (because one of the vertices is a dummy point)
    vector<size_t> incompleteCells;  // keep track of grid cells where there is an error in the field (all fields except density) computation (because one of the vertices is a dummy point)
#endif
    int prev = 0, amount100 = 0;    // variable to show the user about the progress of the computation
    
    
    // parameters for searching in the triangulation
    Locate_type lt;
    int li, lj;
    Cell_handle current, previous;
    
    
    // the random number generator
    size_t seed = userOptions.randomSeed;
    if ( userOptions.partNo>=0 )
        seed += 1000*userOptions.partNo; //change seed between different partition computations
    base_generator_type generator( seed );
    boost::uniform_real<> uni_dist(-0.5,0.5);
    boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
    
    
    
    Real x = boxCoordinates[0] + dx[0]/2.;
    for (size_t i=0; i<nGrid[0]; ++i)
    {
        // print information to the user about the ongoing computation
        amount100 = (100 * i)/ nGrid[0];
        if (prev < amount100)
            message.updateProgress( ++prev );
        
        Real y = boxCoordinates[2] + dx[1]/2.;
        for (size_t j=0; j<nGrid[1]; ++j)
        {
#if NO_DIM==2
            previous = dt.locate( Point(x,y) , lt, li, previous );
#elif NO_DIM==3     // add additional z-loop
            Real z = boxCoordinates[4] + dx[2]/2.;
            for (size_t k=0; k<nGrid[2]; ++k)
            {
                previous = dt.locate( Point(x,y,z) , lt, li, lj, previous );
#endif
                
                // define some temporary variables used inside the MC loop
                bool dummyNeighbors = false, dummyVertices = false;
                Real tempD = Real(0.);
#ifdef VELOCITY
                Pvector<Real,noVelComp> tempV = Pvector<Real,noVelComp>::zero();
                Pvector<Real,noGradComp> tempVG = Pvector<Real,noGradComp>::zero();
                Pvector<Real,noVelComp> tempVelocities[NN];
#endif
#ifdef SCALAR
                Pvector<Real,noScalarComp> tempS = Pvector<Real,noScalarComp>::zero();
                Pvector<Real,noScalarGradComp> tempSG = Pvector<Real,noScalarGradComp>::zero();
#endif
                for (int w=0; w<NN; ++w)        // the Monte Carlo loop over the random samples inside each grid cell
                {
#if NO_DIM==2
                    samplePoint = Point( uni()*dx[0] + x, uni()*dx[1] + y );
                    current = dt.locate( samplePoint,lt,li,previous );
#elif NO_DIM==3
                    samplePoint = Point( uni()*dx[0] + x, uni()*dx[1] + y, uni()*dx[2] + z );
                    current = dt.locate( samplePoint,lt,li,lj, previous );
#endif
                    
#ifdef TEST_PADDING
                    if ( hasDummyNeighbor( current ) ) dummyNeighbors = true; // if one of the vertices is or has dummy neighbors, keep track of this cell
                    if ( hasDummyVertex( current ) ) dummyVertices = true; // if one of the vertices is a dummy test point, keep track of this cell
#endif
                    
                    // check if cell is infinite (one of the vertices is virtual) - don't add anything
                    if ( dt.is_infinite(current) )
                        continue;
                    
                    Real posMatrixInverse[NO_DIM][NO_DIM];  //stores the inverse of the vertex position difference matrix
                    positionMatrix( current, posMatrixInverse );
                    Vertex_handle base = current->vertex(0);  // stores the "base" vertex of the Delaunay cell
                    samplePoint = relativeSamplePoint( base, samplePoint ); // now sample point stores the relative position of the sample point with respect to the base vertex
                    
                    if (field.density)
                    {
                        Real densGrad[NO_DIM];
                        densityGrad(current, posMatrixInverse, densGrad);
                        tempD += densityValue(densGrad,base,samplePoint);
                    }
#ifdef VELOCITY
                    if (field.velocity or field.velocity_gradient or field.velocity_std)
                    {
                        Real velGrad[NO_DIM][noVelComp];
                        velocityGrad( current, posMatrixInverse, velGrad );
                        if (field.velocity) tempV += velocityValue(velGrad,base,samplePoint);
                        if (field.velocity_gradient) tempVG += velocityGradient(velGrad);
                        if (field.velocity_std) tempVelocities[w] = velocityValue(velGrad,base,samplePoint);
                    }
#endif
#ifdef SCALAR
                    if (field.scalar or field.scalar_gradient)
                    {
#ifdef MY_SCALAR
                        tempS += customScalar(current, posMatrixInverse, base, samplePoint);
#else
                        Real sGrad[NO_DIM][noScalarComp];
                        scalarGrad( current, posMatrixInverse, sGrad );
                        if (field.scalar) tempS += scalarValue(sGrad,base,samplePoint);
                        if (field.scalar_gradient) tempSG += scalarGradient(sGrad);
#endif
                    }
#endif
                }
                
                
                // update the field values
                if (field.density) density->push_back( tempD/NN );
#ifdef VELOCITY
                if (field.velocity) velocity->push_back( tempV/NN );
                if (field.velocity_gradient) velocity_gradient->push_back( tempVG/NN );
                if (field.velocity_std) velocity_std->push_back( standardDeviation(tempVelocities,NN) );
#endif
#ifdef SCALAR
                if (field.scalar) scalar->push_back( tempS/NN );
                if (field.scalar_gradient) scalar_gradient->push_back( tempSG/NN );
#endif
                
                // keep track of the cell index if the computation is incomplete
#ifdef TEST_PADDING
#if NO_DIM==2
                size_t k = 0;
#endif
                Real index = gridCellIndex( i,j,k, nGrid );
                if ( dummyNeighbors ) updateDummyGridCells( index, &dummyGridCells_d ); // if one of the vertices is or has dummy neighbors, keep track of this cell
                if ( dummyVertices ) updateDummyGridCells( index, &dummyGridCells ); // if one of the vertices is a dummy test point, keep track of this cell
#endif
                
#if NO_DIM==3     // end additional z-loop
                z += dx[2];
            }
#endif
            y += dx[1];
        }
        x += dx[0];
    }
    message << "100%\n" << MESSAGE::Flush;
    
    
#ifdef TEST_PADDING
    if (field.density)
        showCellsContainingDummyPoints( &dummyGridCells_d, userOptions, userOptions.outputFilename+"_density", "density" ); // check if there are any cells which may have an error in thedensity estimation due to an incomplete Delaunay tesselation over the region of interest
    if ( field.velocity or field.velocity_gradient or field.scalar or field.scalar_gradient )
        showCellsContainingDummyPoints( &dummyGridCells, userOptions, userOptions.outputFilename+"_fields", "fields" ); // check if there are any cells which may have an error in the fields estimation (all fields except density) due to an incomplete Delaunay tesselation over the region of interest
#endif
}





/* This function interpolates the fields to a redshift cone grid (usefull for galaxy surveys) - it computes the volume averaged value of the fields within the grid cell. To compute the volume average in each grid cell it uses a Monte Carlo algorithm to sample 'noRandPoints' in each grid cell and take the average.
The arguments are:
    'dt' - the periodic Delaunay triangulation; with the density at the vertex points already computed
    'userOptions' - structure which store the number of grid points along each axis and also gives information about the size of the box along each axis. It also contains information about what quantity/qunatities to compute. 
    'quantities' - vectors storing the interpolated grid fields
*/
void interpolateRedshiftCone_averaged_2(DT &dt,
                                        User_options &userOptions,
                                        Quantities *quantities)
{
    // define some variables for easy access of program wide variables
    size_t *nGrid = &(userOptions.gridSize[0]); // size of the density grid along each axis
    int const NN = userOptions.noPoints;        // number of random points in each grid cell used to compute the volume average
    Box boxCoordinates = userOptions.redshiftCone;    // the box coordinates in which the density is computed
    checkAngles( &(boxCoordinates[2]), (NO_DIM-1)*2 );// transform the angles from degrees to radians
    Real *origin = &(userOptions.originPosition[0]);
    vector<Real> boxLength;                     // the size of the particle box of interest along each direction
    for (size_t i=0; i<NO_DIM; ++i)
        boxLength.push_back( boxCoordinates[2*i+1]-boxCoordinates[2*i] );
    
    // define pointers to the vectors storing the output fields
    Field field = userOptions.aField;
    std::vector<Real>                         *density = &(quantities->density);
    std::vector< Pvector<Real,noVelComp> >    *velocity = &(quantities->velocity);
    std::vector< Pvector<Real,noGradComp> >   *velocity_gradient = &(quantities->velocity_gradient);
    std::vector< Pvector<Real,noScalarComp> > *scalar = &(quantities->scalar);
    std::vector< Pvector<Real,noScalarGradComp> > *scalar_gradient = &(quantities->scalar_gradient);
    
    
    // show a message to the user
    MESSAGE::Message message( userOptions.verboseLevel );   // structure used to output messages to the user ( the higher the verboseLevel, the more messages will be showns (check 'message.h' for details) )
    message << "Computing interpolation of the fields volume averaged over the sampling cell on a redshift cone grid with " << (NO_DIM==2 ? "r*psi": "r*theta*psi") << " = " << 
            MESSAGE::printElements( nGrid, NO_DIM, "*" ) << " bins in the region " << (NO_DIM==2 ? "{[r_min,r_max],[psi_min,psi_max]}" : "{[r_min,r_max],[theta_min,theta_max],[psi_min,psi_max]}") << " = " << userOptions.redshiftCone.print() << ".\n"
            << "\t Done:   " << MESSAGE::Flush;
    
    
    //reserve memory for the output fields
    size_t reserveSize = 1;                     // the number of cell in the grid
    for (size_t i=0; i<NO_DIM; ++i)
        reserveSize *= nGrid[i];
    reserveMemory( density, reserveSize, field.density, dt, "density" );    // reserves memory for the density computations (if any) and initializes the density values to 0 if the Delaunay triangulation is not well defined
    reserveMemory( velocity, reserveSize, field.velocity, dt, "velocity" );
    reserveMemory( velocity_gradient, reserveSize, field.velocity_gradient, dt, "velocity gradient" );
    reserveMemory( scalar, reserveSize, field.scalar, dt, "scalar" );
    reserveMemory( scalar_gradient, reserveSize, field.scalar_gradient, dt, "scalar gradient" );
    if ( dt.number_of_vertices()<NO_DIM+1 ) return;     // stop the computation if the Delaunay triangulation is not well defined
    
    
    
    // temporary variables needed for the interpolation
    Point samplePoint;      // the coordinates of the sampling point
    Real dx[NO_DIM];        // grid spacing along each axis
    for (size_t i=0; i<NO_DIM; ++i)
        dx[i] = boxLength[i] / nGrid[i];
    
#ifdef TEST_PADDING
    vector<size_t> incompleteCells_d;// keep track of grid cells where there is an error in the density field computation (because one of the vertices is a dummy point)
    vector<size_t> incompleteCells;  // keep track of grid cells where there is an error in the field (all fields except density) computation (because one of the vertices is a dummy point)
#endif
    int prev = 0, amount100 = 0;    // variable to show the user about the progress of the computation
    
    
    // parameters for searching in the triangulation
    Locate_type lt;
    int li, lj;
    Cell_handle current, previous;
    
    
    // the random number generator
    size_t seed = userOptions.randomSeed;
    if ( userOptions.partNo>=0 )
        seed += 1000*userOptions.partNo; //change seed between different partition computations
    base_generator_type generator( seed );
    boost::uniform_real<> uni_dist(-0.5,0.5);
    boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
    
    
    
    Real r = boxCoordinates[0] + dx[0]/2.;
    for (size_t i=0; i<nGrid[0]; ++i)
    {
        // print information to the user about the ongoing computation
        amount100 = (100 * i)/ nGrid[0];
        if (prev < amount100)
            message.updateProgress( ++prev );
        
        Real theta = boxCoordinates[2] + dx[1]/2.;
        for (size_t j=0; j<nGrid[1]; ++j)
        {
#if NO_DIM==2
            Real tempX = origin[0] + r*cos(theta);
            Real tempY = origin[1] + r*sin(theta);
            samplePoint = Point(tempX,tempY);
            previous = dt.locate( samplePoint, lt, li, previous );
#elif NO_DIM==3
            Real psi = boxCoordinates[4] + dx[2]/2.;
            for (size_t k=0; k<nGrid[2]; ++k)
            {
                Real tempX = origin[0] + r*sin(theta)*cos(psi);
                Real tempY = origin[1] + r*sin(theta)*sin(psi);
                Real tempZ = origin[2] + r*cos(theta);
                samplePoint = Point(tempX,tempY,tempZ);
                previous = dt.locate( samplePoint, lt, li, lj, previous );
#endif
                
                // define some temporary variables used inside the MC loop
                bool dummyNeighbors = false, dummyVertices = false;
                Real tempD = Real(0.);
#ifdef VELOCITY
                Pvector<Real,noVelComp> tempV = Pvector<Real,noVelComp>::zero();
                Pvector<Real,noGradComp> tempVG = Pvector<Real,noGradComp>::zero();
#endif
#ifdef SCALAR
                Pvector<Real,noScalarComp> tempS = Pvector<Real,noScalarComp>::zero();
                Pvector<Real,noScalarGradComp> tempSG = Pvector<Real,noScalarGradComp>::zero();
#endif
                for (int w=0; w<NN; ++w)        // the Monte Carlo loop over the random samples inside each grid cell
                {
                    // get coordinates of the random sample point inside the grid cell and locate the Delaunay traingle/tetrahedron where the random point lies
#if NO_DIM==2
                    Real tempR = r + uni()*dx[0];
                    Real tempTheta = theta + uni()*dx[1];
                    tempX = origin[0] + tempR*cos(tempTheta);
                    tempY = origin[0] + tempR*sin(tempTheta);
                    samplePoint = Point(tempX,tempY);
                    current = dt.locate( samplePoint,lt,li,previous );
#elif NO_DIM==3
                    Real tempR = r + uni()*dx[0];
                    Real tempTheta = theta + uni()*dx[1];
                    Real tempPsi = psi + uni()*dx[2];
                    tempX = origin[0] + tempR*sin(tempTheta)*cos(tempPsi);
                    tempY = origin[1] + tempR*sin(tempTheta)*sin(tempPsi);
                    tempZ = origin[2] + tempR*cos(tempTheta);
                    samplePoint = Point(tempX,tempY,tempZ);
                    current = dt.locate( samplePoint,lt,li,lj, previous );
#endif
                    
#ifdef TEST_PADDING
                    if ( hasDummyNeighbor( current ) ) dummyNeighbors = true; // if one of the vertices is or has dummy neighbors, keep track of this cell
                    if ( hasDummyVertex( current ) ) dummyVertices = true; // if one of the vertices is a dummy test point, keep track of this cell
#endif
                    
                    // check if cell is infinite (one of the vertices is virtual) - don't add anything
                    if ( dt.is_infinite(current) )
                        continue;
                    
                    Real posMatrixInverse[NO_DIM][NO_DIM];  //stores the inverse of the vertex position difference matrix
                    positionMatrix( current, posMatrixInverse );
                    Vertex_handle base = current->vertex(0);  // stores the "base" vertex of the Delaunay cell
                    samplePoint = relativeSamplePoint( base, samplePoint ); // now sample point stores the relative position of the sample point with respect to the base vertex
                    
                    if (field.density)
                    {
                        Real densGrad[NO_DIM];
                        densityGrad(current, posMatrixInverse, densGrad);
                        tempD += densityValue(densGrad,base,samplePoint);
                    }
#ifdef VELOCITY
                    if (field.velocity or field.velocity_gradient)
                    {
                        Real velGrad[NO_DIM][noVelComp];
                        velocityGrad( current, posMatrixInverse, velGrad );
                        if (field.velocity) tempV += velocityValue(velGrad,base,samplePoint);
                        if (field.velocity_gradient) tempVG += velocityGradient(velGrad);
                    }
#endif
#ifdef SCALAR
                    if (field.scalar or field.scalar_gradient)
                    {
#ifdef MY_SCALAR
                        tempS += customScalar(current, posMatrixInverse, base, samplePoint);
#else
                        Real sGrad[NO_DIM][noScalarComp];
                        scalarGrad( current, posMatrixInverse, sGrad );
                        if (field.scalar) tempS += scalarValue(sGrad,base,samplePoint);
                        if (field.scalar_gradient) tempSG += scalarGradient(sGrad);
#endif
                    }
#endif
                }
                
                
                // update the field values
                if (field.density) density->push_back( tempD/NN );
#ifdef VELOCITY
                if (field.velocity) velocity->push_back( tempV/NN );
                if (field.velocity_gradient) velocity_gradient->push_back( tempVG/NN );
#endif
#ifdef SCALAR
                if (field.scalar) scalar->push_back( tempS/NN );
                if (field.scalar_gradient) scalar_gradient->push_back( tempSG/NN );
#endif
                
                // keep track of the cell index if the computation is incomplete
#ifdef TEST_PADDING
#if NO_DIM==2
                size_t k = 0;
#endif
                Real index = gridCellIndex( i,j,k, nGrid );
                if ( dummyNeighbors ) updateDummyGridCells( index, &dummyGridCells_d ); // if one of the vertices is or has dummy neighbors, keep track of this cell
                if ( dummyVertices ) updateDummyGridCells( index, &dummyGridCells ); // if one of the vertices is a dummy test point, keep track of this cell
#endif
                
#if NO_DIM==3     // end additional z-loop
                psi += dx[2];
            }
#endif
            theta += dx[1];
        }
        r += dx[0];
    }
    message << "100%\n" << MESSAGE::Flush;
    
    
#ifdef TEST_PADDING
    if (field.density)
        showCellsContainingDummyPoints( &dummyGridCells_d, userOptions, userOptions.outputFilename+"_density", "density" ); // check if there are any cells which may have an error in thedensity estimation due to an incomplete Delaunay tesselation over the region of interest
    if ( field.velocity or field.velocity_gradient or field.scalar or field.scalar_gradient )
        showCellsContainingDummyPoints( &dummyGridCells, userOptions, userOptions.outputFilename+"_fields", "fields" ); // check if there are any cells which may have an error in the fields estimation (all fields except density) due to an incomplete Delaunay tesselation over the region of interest
#endif
}




/* This function interpolates the fields at user defined sampling points - it computes the volume averaged value of the fields within the grid cell. To compute the volume average in each grid cell it uses a Monte Carlo algorithm to sample 'noRandPoints' in each grid cell and take the average.
The arguments are:
    'dt' - the periodic Delaunay triangulation; with the density at the vertex points already computed
    'userOptions' - structure which store the number of grid points along each axis and also gives information about the size of the box along each axis. It also contains information about what quantity/qunatities to compute. 
    'quantities' - vectors storing the interpolated fields
*/
void interpolateUserSampling_averaged_2(DT &dt,
                                        vector<Sample_point> &samples,
                                        User_options &userOptions,
                                        Quantities *quantities)
{
    int const NN = userOptions.noPoints;        // number of random points in each grid cell used to compute the volume average
    
    // define pointers to the vectors storing the output fields
    Field field = userOptions.aField;
    std::vector<Real>                         *density = &(quantities->density);
    std::vector< Pvector<Real,noVelComp> >    *velocity = &(quantities->velocity);
    std::vector< Pvector<Real,noGradComp> >   *velocity_gradient = &(quantities->velocity_gradient);
    std::vector< Pvector<Real,noScalarComp> > *scalar = &(quantities->scalar);
    std::vector< Pvector<Real,noScalarGradComp> > *scalar_gradient = &(quantities->scalar_gradient);
    
    
    // show a message to the user
    MESSAGE::Message message( userOptions.verboseLevel );   // structure used to output messages to the user ( the higher the verboseLevel, the more messages will be showns (check 'message.h' for details) )
    message << "Computing interpolation of the fields volume averaged over the sampling cell at user given sampling points in the region " << userOptions.region.print() << ".\n"
            << "\t Done:   " << MESSAGE::Flush;
    
    
    //reserve size for the output fields
    size_t reserveSize = samples.size();                   // the number of cell in the grid
    reserveMemory( density, reserveSize, field.density, dt, "density" );    // reserves memory for the density computations (if any) and initializes the density values to 0 if the Delaunay triangulation is not well defined
    reserveMemory( velocity, reserveSize, field.velocity, dt, "velocity" );
    reserveMemory( velocity_gradient, reserveSize, field.velocity_gradient, dt, "velocity gradient" );
    reserveMemory( scalar, reserveSize, field.scalar, dt, "scalar" );
    reserveMemory( scalar_gradient, reserveSize, field.scalar_gradient, dt, "scalar gradient" );
    if ( dt.number_of_vertices()<NO_DIM+1 ) return;     // stop the computation if the Delaunay triangulation is not well defined
    
    
    
    // temporary variables needed for the interpolation
    Point samplePoint;      // the coordinates of the sampling point
    
#ifdef TEST_PADDING
    vector<size_t> incompleteCells_d;// keep track of grid cells where there is an error in the density field computation (because one of the vertices is a dummy point)
    vector<size_t> incompleteCells;  // keep track of grid cells where there is an error in the field (all fields except density) computation (because one of the vertices is a dummy point)
#endif
    int prev = 0, amount100 = 0;    // variable to show the user about the progress of the computation
    
    
    // parameters for searching in the triangulation
    Locate_type lt;
    int li, lj;
    Cell_handle current, previous;
    
    
    // the random number generator
    size_t seed = userOptions.randomSeed;
    if ( userOptions.partNo>=0 )
        seed += 1000*userOptions.partNo; //change seed between different partition computations
    base_generator_type generator( seed );
    boost::uniform_real<> uni_dist(-0.5,0.5);
    boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
    
    
    
    size_t const noElements = samples.size();
    for (size_t i=0; i<noElements; ++i)
    {
        // print information to the user about the ongoing computation
        amount100 = (100 * i)/ noElements;
        if (prev < amount100)
            message.updateProgress( ++prev );
        
        // check if the sample point is inside the region of interest, otherwise continue to next sample point
        if ( not userOptions.region.isPointInBox(samples[i].position()) )
            continue;
        
        // locate the Delaunay cell that contains the center of the grid point
#if NO_DIM==2
        samplePoint = Point( samples[i].position(0), samples[i].position(1) );
        previous = dt.locate( samplePoint, lt, li, previous );
#elif NO_DIM==3
        samplePoint = Point( samples[i].position(0), samples[i].position(1), samples[i].position(2) );
        previous = dt.locate( samplePoint, lt, li, lj, previous );
#endif
        
        
        
        // define some temporary variables used inside the MC loop
        bool dummyNeighbors = false, dummyVertices = false;
        Real tempD = Real(0.);
#ifdef VELOCITY
        Pvector<Real,noVelComp> tempV = Pvector<Real,noVelComp>::zero();
        Pvector<Real,noGradComp> tempVG = Pvector<Real,noGradComp>::zero();
#endif
#ifdef SCALAR
        Pvector<Real,noScalarComp> tempS = Pvector<Real,noScalarComp>::zero();
        Pvector<Real,noScalarGradComp> tempSG = Pvector<Real,noScalarGradComp>::zero();
#endif
        for (int w=0; w<NN; ++w)        // the Monte Carlo loop over the random samples inside each grid cell
        {
            Real tempX = samples[i].position(0) + uni()*samples[i].delta(0);
            Real tempY = samples[i].position(1) + uni()*samples[i].delta(1);
#if NO_DIM==2
            samplePoint = Point(tempX,tempY);
            current = dt.locate( samplePoint,lt,li,previous );
#elif NO_DIM==3
            Real tempZ = samples[i].position(2) + uni()*samples[i].delta(2);
            samplePoint = Point(tempX,tempY,tempZ);
            current = dt.locate( samplePoint,lt,li,lj, previous );
#endif
            
#ifdef TEST_PADDING
            if ( hasDummyNeighbor( current ) ) dummyNeighbors = true; // if one of the vertices is or has dummy neighbors, keep track of this cell
            if ( hasDummyVertex( current ) ) dummyVertices = true; // if one of the vertices is a dummy test point, keep track of this cell
#endif
            
            // check if cell is infinite (one of the vertices is virtual) - don't add anything
            if ( dt.is_infinite(current) )
                continue;
            
            Real posMatrixInverse[NO_DIM][NO_DIM];  //stores the inverse of the vertex position difference matrix
            positionMatrix( current, posMatrixInverse );
            Vertex_handle base = current->vertex(0);  // stores the "base" vertex of the Delaunay cell
            samplePoint = relativeSamplePoint( base, samplePoint ); // now sample point stores the relative position of the sample point with respect to the base vertex
            
            if (field.density)
            {
                Real densGrad[NO_DIM];
                densityGrad(current, posMatrixInverse, densGrad);
                tempD += densityValue(densGrad,base,samplePoint);
            }
#ifdef VELOCITY
            if (field.velocity or field.velocity_gradient)
            {
                Real velGrad[NO_DIM][noVelComp];
                velocityGrad( current, posMatrixInverse, velGrad );
                if (field.velocity) tempV += velocityValue(velGrad,base,samplePoint);
                if (field.velocity_gradient) tempVG += velocityGradient(velGrad);
            }
#endif
#ifdef SCALAR
            if (field.scalar or field.scalar_gradient)
            {
#ifdef MY_SCALAR
                tempS += customScalar(current, posMatrixInverse, base, samplePoint);
#else
                Real sGrad[NO_DIM][noScalarComp];
                scalarGrad( current, posMatrixInverse, sGrad );
                if (field.scalar) tempS += scalarValue(sGrad,base,samplePoint);
                if (field.scalar_gradient) tempSG += scalarGradient(sGrad);
#endif
            }
#endif
        }
        
        
        // update the field values
        if (field.density) density->push_back( tempD/NN );
#ifdef VELOCITY
        if (field.velocity) velocity->push_back( tempV/NN );
        if (field.velocity_gradient) velocity_gradient->push_back( tempVG/NN );
#endif
#ifdef SCALAR
        if (field.scalar) scalar->push_back( tempS/NN );
        if (field.scalar_gradient) scalar_gradient->push_back( tempSG/NN );
#endif
        
        // keep track of the cell index if the computation is incomplete
#ifdef TEST_PADDING
        if ( dummyNeighbors ) updateDummyGridCells( i, &dummyGridCells_d ); // if one of the vertices is or has dummy neighbors, keep track of this cell
        if ( dummyVertices ) updateDummyGridCells( i, &dummyGridCells ); // if one of the vertices is a dummy test point, keep track of this cell
#endif
    }
    message << "100%\n" << MESSAGE::Flush;
    
    
#ifdef TEST_PADDING
    if (field.density)
        showCellsContainingDummyPoints( &dummyGridCells_d, userOptions, userOptions.outputFilename+"_density", "density" ); // check if there are any cells which may have an error in thedensity estimation due to an incomplete Delaunay tesselation over the region of interest
    if ( field.velocity or field.velocity_gradient or field.scalar or field.scalar_gradient )
        showCellsContainingDummyPoints( &dummyGridCells, userOptions, userOptions.outputFilename+"_fields", "fields" ); // check if there are any cells which may have an error in the fields estimation (all fields except density) due to an incomplete Delaunay tesselation over the region of interest
#endif
}










/* This function interpolates the fields to a uniform grid - it computes the volume averaged value of the fields within the grid cell. To compute the volume average in each grid cell it uses equidistant sampling points to sample 'noRandPoints' in each grid cell and take the average.
The arguments are:
    'dt' - the periodic Delaunay triangulation; with the density at the vertex points already computed
    'userOptions' - structure which store the number of grid points along each axis and the number of random sample points in each grid cell. This structure also gives information about the size of the box along each axis.
    'quantities' - structure storing the vectors for the output fields on the grid
*/
void interpolateGrid_averaged_3(DT &dt,
                                User_options &userOptions,
                                Quantities *quantities)
{
    size_t *nGrid = &(userOptions.gridSize[0]); // size of the density grid along each axis
    int const NN = userOptions.noPoints;        // number of equidistant points in each grid cell used to compute the volume average
    int const n = rootN( NN, NO_DIM );          // number of equidistant sampling points along each axis
    Box boxCoordinates = userOptions.region;    // the box coordinates in which the density is computed
    vector<Real> boxLength;                     // the size of the particle box of interest along each direction
    for (size_t i=0; i<NO_DIM; ++i)
        boxLength.push_back( boxCoordinates[2*i+1]-boxCoordinates[2*i] );
    
    
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
            << ". The volume average is done using " << NN << " equidistant sampling points in each grid cell.\n"
            << "\t Done: " << MESSAGE::Flush;
    
    
    //reserve memory for the output fields
    size_t reserveSize = (NO_DIM==2) ? nGrid[0]*nGrid[1] : nGrid[0]*nGrid[1]*nGrid[2];      // the number of cell in the grid
    reserveMemory( density, reserveSize, field.density, dt, "density" );    // reserves memory for the density computations (if any) and initializes the density values to 0 if the Delaunay triangulation is not well defined
    reserveMemory( velocity, reserveSize, field.velocity, dt, "velocity" );
    reserveMemory( velocity_gradient, reserveSize, field.velocity_gradient, dt, "velocity gradient" );
    reserveMemory( scalar, reserveSize, field.scalar, dt, "scalar" );
    reserveMemory( scalar_gradient, reserveSize, field.scalar_gradient, dt, "scalar gradient" );
    if ( dt.number_of_vertices()<NO_DIM+1 ) return;     // stop the computation if the Delaunay triangulation is not well defined
    
    
    // temporary variables needed for the interpolation
    Point samplePoint;      // the coordinates of the MC sampling point
    Real dx[NO_DIM];        // grid spacing along each axis
    for (size_t i=0; i<NO_DIM; ++i)
        dx[i] = boxLength[i] / nGrid[i];
    
#ifdef TEST_PADDING
    vector<size_t> incompleteCells_d;// keep track of grid cells where there is an error in the density field computation (because one of the vertices is a dummy point)
    vector<size_t> incompleteCells;  // keep track of grid cells where there is an error in the field (all fields except density) computation (because one of the vertices is a dummy point)
#endif
    int prev = 0, amount100 = 0;    // variable to show the user about the progress of the computation
    
    
    // parameters for searching in the triangulation
    Locate_type lt;
    int li, lj;
    Cell_handle current, previous;
    
    
    // sequence to store the equidistant sample positions with respect to center of grid cell
    Real samplePosition[NN][NO_DIM];
    for (int i1=0; i1<n; ++i1)
        for (int i2=0; i2<n; ++i2)
#if NO_DIM==2
        {
            int index = i1*n + i2;
#endif
#if NO_DIM==3
            for (int i3=0; i3<n; ++i3)
            {
                int index = i1*n*n + i2*n + i3;
                samplePosition[index][2] = (i3+0.5)*dx[2]/n - dx[2]/2.;
#endif
                samplePosition[index][0] = (i1+0.5)*dx[0]/n - dx[0]/2.;
                samplePosition[index][1] = (i2+0.5)*dx[1]/n - dx[1]/2.;
            }
    
    
    
    
    Real x = boxCoordinates[0] + dx[0]/2.;
    for (size_t i=0; i<nGrid[0]; ++i)
    {
        // print information to the user about the ongoing computation
        amount100 = (100 * i)/ nGrid[0];
        if (prev < amount100)
            message.updateProgress( ++prev );
        
        Real y = boxCoordinates[2] + dx[1]/2.;
        for (size_t j=0; j<nGrid[1]; ++j)
        {
#if NO_DIM==2
            previous = dt.locate( Point(x,y) , lt, li, previous );
#elif NO_DIM==3     // add additional z-loop
            Real z = boxCoordinates[4] + dx[2]/2.;
            for (size_t k=0; k<nGrid[2]; ++k)
            {
                previous = dt.locate( Point(x,y,z) , lt, li, lj, previous );
#endif
                
                
                // define some temporary variables used inside the MC loop
                bool dummyNeighbors = false, dummyVertices = false;
                Real tempD = Real(0.);
#ifdef VELOCITY
                Pvector<Real,noVelComp> tempV = Pvector<Real,noVelComp>::zero();
                Pvector<Real,noGradComp> tempVG = Pvector<Real,noGradComp>::zero();
#endif
#ifdef SCALAR
                Pvector<Real,noScalarComp> tempS = Pvector<Real,noScalarComp>::zero();
                Pvector<Real,noScalarGradComp> tempSG = Pvector<Real,noScalarGradComp>::zero();
#endif
                for (int w=0; w<NN; ++w)        // the loop over the equidistant sample points inside each grid cell
                {
#if NO_DIM==2
                    samplePoint = Point( x + samplePosition[w][0], y + samplePosition[w][1] );
                    current = dt.locate( samplePoint,lt,li,previous );
#elif NO_DIM==3
                    samplePoint = Point( x + samplePosition[w][0], y + samplePosition[w][1], z + samplePosition[w][2] );
                    current = dt.locate( samplePoint,lt,li,lj, previous );
#endif
                    
#ifdef TEST_PADDING
                    if ( hasDummyNeighbor( current ) ) dummyNeighbors = true; // if one of the vertices is or has dummy neighbors, keep track of this cell
                    if ( hasDummyVertex( current ) ) dummyVertices = true; // if one of the vertices is a dummy test point, keep track of this cell
#endif
                    
                    // check if cell is infinite (one of the vertices is virtual) - don't add anything
                    if ( dt.is_infinite(current) )
                        continue;
                    
                    Real posMatrixInverse[NO_DIM][NO_DIM];  //stores the inverse of the vertex position difference matrix
                    positionMatrix( current, posMatrixInverse );
                    Vertex_handle base = current->vertex(0);  // stores the "base" vertex of the Delaunay cell
                    samplePoint = relativeSamplePoint( base, samplePoint ); // now sample point stores the relative position of the sample point with respect to the base vertex
                    
                    if (field.density)
                    {
                        Real densGrad[NO_DIM];
                        densityGrad(current, posMatrixInverse, densGrad);
                        tempD += densityValue(densGrad,base,samplePoint);
                    }
#ifdef VELOCITY
                    if (field.velocity or field.velocity_gradient)
                    {
                        Real velGrad[NO_DIM][noVelComp];
                        velocityGrad( current, posMatrixInverse, velGrad );
                        if (field.velocity) tempV += velocityValue(velGrad,base,samplePoint);
                        if (field.velocity_gradient) tempVG += velocityGradient(velGrad);
                    }
#endif
#ifdef SCALAR
                    if (field.scalar or field.scalar_gradient)
                    {
                        Real sGrad[NO_DIM][noScalarComp];
                        scalarGrad( current, posMatrixInverse, sGrad );
                        if (field.scalar) tempS += scalarValue(sGrad,base,samplePoint);
                        if (field.scalar_gradient) tempSG += scalarGradient(sGrad);
                    }
#endif
                }
                
                
                // update the field values
                if (field.density) density->push_back( tempD/NN );
#ifdef VELOCITY
                if (field.velocity) velocity->push_back( tempV/NN );
                if (field.velocity_gradient) velocity_gradient->push_back( tempVG/NN );
#endif
#ifdef SCALAR
                if (field.scalar) scalar->push_back( tempS/NN );
                if (field.scalar_gradient) scalar_gradient->push_back( tempSG/NN );
#endif
                
                // keep track of the cell index if the computation is incomplete
#ifdef TEST_PADDING
#if NO_DIM==2
                size_t k = 0;
#endif
                Real index = gridCellIndex( i,j,k, nGrid );
                if ( dummyNeighbors ) updateDummyGridCells( index, &dummyGridCells_d ); // if one of the vertices is or has dummy neighbors, keep track of this cell
                if ( dummyVertices ) updateDummyGridCells( index, &dummyGridCells ); // if one of the vertices is a dummy test point, keep track of this cell
#endif
                
#if NO_DIM==3     // end additional z-loop
                z += dx[2];
            }
#endif
            y += dx[1];
        }
        x += dx[0];
    }
    message << "100%\n" << MESSAGE::Flush;
    
    
#ifdef TEST_PADDING
    if (field.density)
        showCellsContainingDummyPoints( &dummyGridCells_d, userOptions, userOptions.outputFilename+"_density", "density" ); // check if there are any cells which may have an error in thedensity estimation due to an incomplete Delaunay tesselation over the region of interest
    if ( field.velocity or field.velocity_gradient or field.scalar or field.scalar_gradient )
        showCellsContainingDummyPoints( &dummyGridCells, userOptions, userOptions.outputFilename+"_fields", "fields" ); // check if there are any cells which may have an error in the fields estimation (all fields except density) due to an incomplete Delaunay tesselation over the region of interest
#endif
}


