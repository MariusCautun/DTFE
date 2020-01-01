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


#ifdef TEST_PADDING
/* This file contains the functions used to test the padding efficiency.
The following functions check if the interpolation values are affected by an incomplete Delaunay tesselation of the region of interest. In the case of density, this means that a Delaunay cell of interest has dummy points as vertices or as neighbors of the vertices.
*/
#include <fstream>



/* Inserts dummy points to test the efficiency of the padding. */
void insertDummyTestParticles(DT &dt,
                            User_options &userOptions)
{
    Box paddedBox = userOptions.paddedBox;
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "Inserting dummy test points on the boundaries of the padded box to test for completness of the Delaunay tesselation inside the box of interest ... " << MESSAGE::Flush;
    Real nprt = (Real) dt.number_of_vertices();
    size_t n = int( pow( nprt, 1./NO_DIM ) ) + 1.;    //number of dummy particles to be inserted along each direction
    
    Real dx[NO_DIM];  //keep track of grid spacing of dummy points along each direction
    for (size_t i=0; i<NO_DIM; ++i)
        dx[i] = (paddedBox[2*i+1] - paddedBox[2*i]) / (n-1);
    Real xMin[NO_DIM], xMax[NO_DIM];
    for (size_t i=0; i<NO_DIM; ++i)
    {
        xMin[i] = paddedBox[2*i];
        xMax[i] = paddedBox[2*i+1];
    }
    
    vector< Point > dummyPoints;

#if NO_DIM==2
    dummyPoints.reserve( 2*NO_DIM*n );
    for (size_t i1=0; i1<n; ++i1)
    {
        dummyPoints.push_back( Point( xMin[0], xMin[1]+dx[1]*i1) );
        dummyPoints.push_back( Point( xMax[0], xMin[1]+dx[1]*i1) );
        
        dummyPoints.push_back( Point( xMin[0]+dx[0]*i1, xMin[1]) );
        dummyPoints.push_back( Point( xMin[0]+dx[0]*i1, xMax[1]) );
    }
#elif NO_DIM==3
    dummyPoints.reserve( 2*NO_DIM*n*n );
    for (size_t i1=0; i1<n; ++i1)
        for (size_t i2=0; i2<n; ++i2)
    {
        dummyPoints.push_back( Point( xMin[0], xMin[1]+dx[1]*i1, xMin[2]+dx[2]*i2) );
        dummyPoints.push_back( Point( xMax[0], xMin[1]+dx[1]*i1, xMin[2]+dx[2]*i2) );
            
        dummyPoints.push_back( Point( xMin[0]+dx[0]*i1, xMin[1], xMin[2]+dx[2]*i2) );
        dummyPoints.push_back( Point( xMin[0]+dx[0]*i1, xMax[1], xMin[2]+dx[2]*i2) );
            
        dummyPoints.push_back( Point( xMin[0]+dx[0]*i1, xMin[1]+dx[1]*i2, xMin[2]) );
        dummyPoints.push_back( Point( xMin[0]+dx[0]*i1, xMin[1]+dx[1]*i2, xMax[2]) );
    }
#endif  // end of NO_DIM condition
    
    
    message << "\n\tNumber of dummy test points to be added: " << dummyPoints.size() << "\n";
    message << "\tPercent of total points: " << 100.*((Real) dummyPoints.size())/nprt << "\n";

    message << "\tUpdating the triangulation.\n\tDone: " << MESSAGE::Flush;
    
    // sort the points to be spatially close - for faster insertion into the triangulation
    size_t prev = 0, amount100 = 0, count = 0;    // variable to show the user about the progress of the computation
    size_t const noPoints = dummyPoints.size();   // total number of dummy points
    CGAL::spatial_sort( dummyPoints.begin(), dummyPoints.end() );
    vector<Point>::iterator it = dummyPoints.begin();
    Vertex_handle vh = dt.nearest_vertex( *it );
    for ( ; it!=dummyPoints.end(); ++it)
    {
        vh = dt.insert( *it, vh );
        vh->info().setDummy();
        
        // show the progress of the computation
        amount100 = (100 * count++)/ noPoints;
        if (prev < amount100)
            message.updateProgress( ++prev );
    }
    dummyPoints.clear();
    message << "100\%\n" << MESSAGE::Flush;
}




/* Check if any of the vertices of a Delaunay cell has dummy neighbors. */
template <typename templateCell>
inline bool hasDummyNeighbor(templateCell &cell)
{
    for (int i=0; i<NO_DIM+1; ++i)
        if ( cell->vertex(i)->info().hasDummyNeighbor() )
            return true;
    return false;
}
/* Check if any of the vertices of a Delaunay cell is a dummy test point. */
template <typename templateCell>
inline bool hasDummyVertex(templateCell &cell)
{
    for (int i=0; i<NO_DIM+1; ++i)
        if ( cell->vertex(i)->info().isDummy() )
            return true;
    return false;
}
/* Keep track of the cells that have wrong values of the density. */
inline void updateDummyGridCells(size_t const gridIndex,
                                 vector<size_t> *dummyGridCells)
{
    if ( dummyGridCells->empty() or gridIndex!=dummyGridCells->back() )
        dummyGridCells->push_back( gridIndex );
}
/* Checks if there are any cells which may have an error in the density estimation due to an incomplete Delaunay tesselation over the region of interest. */
void showCellsContainingDummyPoints(vector<size_t> *dummyGridCells,
                                    User_options &userOptions,
                                    string outputName,
                                    string quantityName,
                                    bool userDefinedSampling = false,
                                    size_t *nGrid = NULL)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    if ( dummyGridCells->size()==0 )
    {
        if ( userOptions.noProcessors!=1 )
            message << "All the '" << quantityName << "' grid cells on this thread have a complete estimation of the value associated to them. There may still be other grid cells on the other processors that have an incomplete estimation of the value associated to them. For that check the output directory for files of the type '" << outputName + ".badGridPoints" << "'\n";
        return;
    }
    if ( nGrid==NULL) nGrid = &(userOptions.gridSize[0]);
    size_t const totalSize = (NO_DIM==2)? (nGrid[0]*nGrid[1]) : (nGrid[0]*nGrid[1]*nGrid[2]);  // the number of cell in the density grid
    char threadId[100] = "";
    if ( userOptions.noProcessors!=1 ) sprintf( threadId, "%d", userOptions.threadId );
    string output = outputName + ".badGridPoints" + threadId;   // name of the output file
    
    sort( dummyGridCells->begin(), dummyGridCells->end() );
    vector<size_t>::iterator it = unique( dummyGridCells->begin(), dummyGridCells->end() );
    dummyGridCells->resize( it - dummyGridCells->begin() );
    message << "The program found that some '" << quantityName << "' grid cells have an incomplete estimation of the value associated to them. The cause of this is the Delaunay tessellation not covering fully all the region where the '" << quantityName << "' is computed. There are " << dummyGridCells->size() << " cells with an incomplete estimation, these represent " << Real(dummyGridCells->size())/totalSize*100 << "\% of the total cells on this thread. The grid cells which have incomplete values will be written to the file '" << output << "'.\n" << MESSAGE::Flush;
    
    std::fstream outputFile;
    outputFile.open( output.c_str(), std::ios::out );
    outputFile << "#Number of cells with an incomplete '" << quantityName << "' value:  " << dummyGridCells->size() << "\n";
    if ( not userDefinedSampling )
    {
        outputFile << "#grid index\tn_x\tn_y\tn_z\n";
        for (int i=0; i<dummyGridCells->size(); ++i)
        {
#if NO_DIM==2
            size_t const tempX = dummyGridCells->at(i)/nGrid[1];
            size_t const tempY = dummyGridCells->at(i) - nGrid[1]*tempX;
            outputFile << "\t" << dummyGridCells->at(i) << "\t" << tempX << "\t" << tempY << "\n";
            
#elif NO_DIM==3
            size_t const tempX = dummyGridCells->at(i)/(nGrid[1]*nGrid[2]);
            size_t const tempY = (dummyGridCells->at(i) - nGrid[1]*nGrid[2]*tempX) / nGrid[2];
            size_t const tempZ = dummyGridCells->at(i) - nGrid[1]*nGrid[2]*tempX - nGrid[2]*tempY;
            outputFile << "\t" << dummyGridCells->at(i) << "\t" << tempX << "\t" << tempY << "\t" << tempZ << "\n";
#endif
        }
    }
    else
    {
        outputFile << "#grid index\n";
        for (int i=0; i<dummyGridCells->size(); ++i)
            outputFile << "\t" << dummyGridCells->at(i) << "\n";
    }
    outputFile.close();
}

#endif  // end " #ifdef TEST_PADDING "
