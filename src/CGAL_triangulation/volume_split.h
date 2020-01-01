#ifndef VOLUME_SPLIT_HEADER
#define VOLUME_SPLIT_HEADER

#include "../Pvector.h"
#include <cmath>


#ifndef VOLUME_TOL
    #define VOLUME_TOL Real(1.e-3)
#endif

#if NO_DIM==2
    #define NO_PAIRS 3
#else
    #define NO_PAIRS 6
#endif


namespace VOLUME_SPLIT
{

struct Vertex
{
    Pvector<Real,NO_DIM> coords;
    Pvector<int,NO_DIM>  nMin;
    Pvector<int,NO_DIM>  nMax;
    
    inline Real& operator[](int i)
    {
        return this->coords[i];
    }
    
    static Vertex zero()
    {
        return Vertex();
    }
    
    inline void update_nMinMax(int const i, Real const p, Real const tol)
    {
        this->nMin[i]  = (int) (p-tol);
        this->nMax[i]  = (int) (p+tol);
    }
};




struct Simplex
{
    std::vector< Vertex > pos;
    Real vol;   // keeps the volume of the simplex
    
    Simplex()
    {
        pos.assign( NO_DIM+1, Vertex::zero() );
    }
    
    /* Overload the [] operator. */
    Vertex& operator[](int i)
    {
        return this->pos[i];
    }
    
    /* Assign a set of vertices. */
    template <typename T>
    void assign(T const cell)
    {
        for (int i=0; i<NO_DIM+1; ++i)
            for (int j=0; j<NO_DIM; ++j)
                pos[i][j] = cell->vertex(i)->point()[j];
    }
    
    // returns the volume of a simplex
    Real volume()
    {
        Real posDiff[NO_DIM][NO_DIM];
        for (int i=0; i<NO_DIM; ++i)
            for (int j=0; j<NO_DIM; ++j)
                posDiff[i][j] = this->pos[i+1][j] - this->pos[0][j];
        
#if NO_DIM==2
        return Real(.5) * std::fabs( posDiff[0][0]*posDiff[1][1] - posDiff[0][1]*posDiff[1][0] );
#elif NO_DIM==3
        return Real(1./6.) * std::fabs( posDiff[0][0]*posDiff[1][1]*posDiff[2][2] + posDiff[0][1]*posDiff[1][2]*posDiff[2][0] + posDiff[0][2]*posDiff[1][0]*posDiff[2][1] - posDiff[0][2]*posDiff[1][1]*posDiff[2][0] - posDiff[0][0]*posDiff[1][2]*posDiff[2][1] - posDiff[0][1]*posDiff[1][0]*posDiff[2][2] );
#endif
    }
    
    // returns the center of mass of the simplex
    void massCenter(Vertex &cm)
    {
        for (int i=0; i<NO_DIM; ++i)
        {
            cm[i] = Real(0.);
            for (int j=0; j<NO_DIM+1; ++j)
                cm[i] += this->pos[j][i];
            cm[i] /= Real(NO_DIM+1);
        }
    }
};








struct VolumeSplit
{
    // variable that remain the same for all simplices
    Real start[NO_DIM];     // keeps the origin of the full grid
    size_t  grid[NO_DIM];   // keeps the dimensions of the grid
    Real dx[NO_DIM];        // the grid spoacing along each dimension
    int  pairs[NO_PAIRS][2];// the posible permutations of vertice pair that give all the edges of the triangle/tetrahedron
    Real tol;               // tolerance factor that deals with comparing real numbers
    Real tolVolume;
    Real cellVolume;        // the volume of a grid cell
    
    // variables that change with simplices
    Real currentStart[NO_DIM];  // the origin of a new, smaller grid, that fully encompasses the current simplex
    int  newGrid[NO_DIM];       // the dimensions of this new grid (it has the same dx as the main grid)
    std::vector< Pvector<Real,NO_DIM+1> > *results;    //will store the contributions for each cell of the smaller grid
    
    std::vector< Simplex > simplices;
    
    
    VolumeSplit(Real *startPos, size_t *gridDims, Real *dxValues)
    {
        tol = Real(1.e-5);
        tolVolume = 1.e-3;
        cellVolume = Real(1.);
        for (int i=0; i<NO_DIM; ++i)
        {
            start[i] = startPos[i];
            grid[i] = gridDims[i];
            dx[i] = dxValues[i];
            cellVolume *= dx[i];
        }
         simplices.assign( 50000, Simplex() );
        
#if NO_DIM==2
        pairs[0][0] = 0; pairs[0][1] = 1;
        pairs[1][0] = 1; pairs[1][1] = 2;
        pairs[2][0] = 2; pairs[2][1] = 0;
#elif  NO_DIM==3
        pairs[0][0] = 0; pairs[0][1] = 1;
        pairs[1][0] = 0; pairs[1][1] = 2;
        pairs[2][0] = 0; pairs[2][1] = 3;
        pairs[3][0] = 1; pairs[3][1] = 2;
        pairs[4][0] = 1; pairs[4][1] = 3;
        pairs[5][0] = 2; pairs[5][1] = 3;
#endif
    }
    
    
    // sorts the edges (or pairs of points) in descending order of their lengths
    void sortPairs(Simplex &pos,
                   int pairs[NO_PAIRS][2],
                   int newPairs[NO_PAIRS][2])
    {
        Real len[NO_PAIRS];     // variable that stores the length of each edge of the simplex
        int  count[NO_PAIRS];   // variable to keep track of the position of the edge when edges are ordered in descending order of their lengths
        // find the length of each edge
        for (int i=0; i<NO_PAIRS; ++i)
        {
            len[i] = Real(0.);
            int i1 = pairs[i][0], i2 = pairs[i][1];
            for (int j=0; j<NO_DIM; ++j)
            {
                Real temp = pos[i1][j] - pos[i2][j];
                len[i] += temp*temp;
            }
            count[i] = 0;
        }
        // sort the edges according to their length
        for (int i=0; i<NO_PAIRS; ++i)
            for (int j=i+1; j<NO_PAIRS; ++j)
                if ( len[i]>=len[j] )
                    count[j] += 1;
                else
                    count[i] += 1;
        // copy the edges according to the order stored in 'count'
        for (int i=0; i<NO_PAIRS; ++i)
            for (int j=0; j<2; ++j)
               newPairs[count[i]][j] = pairs[i][j];
    }
    
    
    // returns the grid index for a cell of the small grid
    int index_smallGrid(int *n)
    {
        if (n[0]<0 or n[0]>=this->newGrid[0]) return -1;
        if (n[1]<0 or n[1]>=this->newGrid[1]) return -1;
#if NO_DIM==2
        return n[0]*newGrid[1] + n[1];
#elif NO_DIM==3
        if (n[2]<0 or n[2]>=this->newGrid[2]) return -1;
        return n[0]*newGrid[1]*newGrid[2] + n[1]*newGrid[2] + n[2];
#endif
    }
    
    // int returns the grid index of a point associated to the smaller grid
    int index_smallGrid(Vertex pos)
    {
        int n[NO_DIM];
        for (int i=0; i<NO_DIM; ++i)
            n[i] = (int) pos[i];
        return index_smallGrid( n );
    }
    
    
    // computes the contribution of a simplex that is fully contained within a grid cell
    void simplexContribution(Simplex &pos)
    {
        this->simplexContribution( pos, pos.volume() );
    }
    void simplexContribution(Simplex &pos, Real vol)
    {
        // find the center of mass to know in which cell this simplex lies
        Vertex cm;
        pos.massCenter( cm );
        int index = this->index_smallGrid( cm );
        if (index<0)
            return;
        Real volume = vol * cellVolume;
        
        // store the values in the output matrix
        (*this->results)[index][NO_DIM] += volume;
        for (int i=0; i<NO_DIM; ++i)
            (*this->results)[index][i] += ( cm[i]*dx[i] + currentStart[i] ) * volume;
    }
    
    
    // tests if the simplex intersects with the grid wall. If true, splits the simplex in 2, if false, calls 'simplexContribution' function.
    void simplexIteration(Simplex &pos, Real vol)
    {
        int sortedPairs[NO_PAIRS][2];
        this->sortPairs( pos, this->pairs, sortedPairs );
        
        
        // check if the simplex splits into smaller units
        for (int i=0; i<NO_PAIRS; ++i)
        {
            int i1 = sortedPairs[i][0]; //label of 1st vertex
            int i2 = sortedPairs[i][1]; //label of 2nd vertex
            
            // check intersection with x-axis, y-axis and z-axis for 3D case (only x- and y- for 2D)
            for (int j=0; j<NO_DIM; ++j)
            {
                if ( pos[i1].nMin[j]>pos[i2].nMax[j] or pos[i1].nMax[j]<pos[i2].nMin[j] ) // this expression is true only where there is an intersection of the j-axis (0=x,1=y,2=z) with this edge
                {
                    int const j2 = (j+1) % NO_DIM;// j2=1 for x-axis, j2=2 for y-axis, j2=0 for z-axis
                    int const j3 = (j+2) % NO_DIM;// j3=2 for x-axis, j3=0 for y-axis, j3=1 for z-axis
                    
                    // check wich vertex has the smaller coordinate, this determines the value of 'nMin'
                    int nMin = pos[i1].nMax[j];
                    if ( pos[i1].nMin[j]>pos[i2].nMax[j] )
                        nMin = pos[i2].nMax[j];
                    
                    // get the slope of the edge along the j2 and j3 axes
                    Real temp = pos[i1][j]-pos[i2][j];    // the difference in coordinates along the axis in question
                    Real slope2 = (pos[i1][j2]-pos[i2][j2]) / temp;
                    Real constant2 = pos[i1][j2] - slope2*pos[i1][j];
#if NO_DIM==3
                    Real slope3 = (pos[i1][j3]-pos[i2][j3]) / temp;
                    Real constant3 = pos[i1][j3] - slope3*pos[i1][j];
#endif
                    
                    // get the new point along which to divide the simplex
                    Vertex newPoint;
                    newPoint[j]  = Real(nMin+1);
                    newPoint[j2] = slope2*newPoint[j] + constant2;
                    newPoint.nMin[j]  = nMin;
                    newPoint.nMax[j]  = nMin+1;
                    newPoint.update_nMinMax( j2, newPoint[j2], this->tol );
#if NO_DIM==3
                    newPoint[j3] = slope3*newPoint[j] + constant3;
                    newPoint.update_nMinMax( j3, newPoint[j3], this->tol );
#endif
                    
                    // reiterate over two new simplices
                    Simplex pos2 = pos;
                    pos2[i2] = newPoint;
                    Real vol1 = pos2.volume();
                    // do next step for first simplex
                    if ( vol1<this->tolVolume )
                        this->simplexContribution( pos, vol1 );
                    else
                        this->simplexIteration( pos2, vol1 );
                    // do next step for second simplex
                    pos2[i1] = pos[i2];
                    vol1 = vol - vol1;
                    if ( vol1<this->tolVolume )
                        this->simplexContribution( pos, vol1 );
                    else
                        this->simplexIteration( pos2, vol1 );
                    return;
                }
            }
        }
        
        // the simplex does not split any more, so compute it's contribution
        this->simplexContribution( pos, vol );
    }
    
    
    void simplexIteration2(Simplex &initSimplex, Real volPrev)
    {
        simplices[0] = initSimplex;
        int current = 0, total = 1;
        simplices[0].vol = volPrev;
        
        // loop until all simplices are fully contained in a single grid cell
        while ( true )
        {
            int sortedPairs[NO_PAIRS][2];
            Simplex *pos = &(simplices[current]);
            this->sortPairs( *pos, this->pairs, sortedPairs );
            bool simplexFinished = false;
            bool simplexSplit = false;
            
            
            // check if the simplex splits into smaller simplices
            for (int i=0; i<NO_PAIRS; ++i)
            {
                int i1 = sortedPairs[i][0]; //label of 1st vertex
                int i2 = sortedPairs[i][1]; //label of 2nd vertex
                
                // check intersection with x-axis, y-axis and z-axis for 3D case (only x- and y- for 2D)
                for (int j=0; j<NO_DIM; ++j)
                {
                    if ( (*pos)[i1].nMin[j]>(*pos)[i2].nMax[j] or (*pos)[i1].nMax[j]<(*pos)[i2].nMin[j] ) // this expression is true only where there is an intersection of the j-axis (0=x,1=y,2=z) with this edge
                    {
                        int const j2 = (j+1) % NO_DIM;// j2=1 for x-axis, j2=2 for y-axis, j2=0 for z-axis
                        int const j3 = (j+2) % NO_DIM;// j3=2 for x-axis, j3=0 for y-axis, j3=1 for z-axis
                        
                        // check wich vertex has the smaller coordinate, this determines the value of 'nMin'
                        int nMin = (*pos)[i1].nMax[j];
                        if ( (*pos)[i1].nMin[j]>(*pos)[i2].nMax[j] )
                            nMin = (*pos)[i2].nMax[j];
                        
                        // get the slope of the edge along the j2 and j3 axes
                        Real temp = (*pos)[i1][j]-(*pos)[i2][j];    // the difference in coordinates along the axis in question
                        Real slope2 = ((*pos)[i1][j2]-(*pos)[i2][j2]) / temp;
                        Real constant2 = (*pos)[i1][j2] - slope2*(*pos)[i1][j];
#if NO_DIM==3
                        Real slope3 = ((*pos)[i1][j3]-(*pos)[i2][j3]) / temp;
                        Real constant3 = (*pos)[i1][j3] - slope3*(*pos)[i1][j];
#endif
                        
                        // get the new point along which to divide the simplex
                        Vertex newPoint;
                        newPoint[j]  = Real(nMin+1);
                        newPoint[j2] = slope2*newPoint[j] + constant2;
                        newPoint.nMin[j]  = nMin;
                        newPoint.nMax[j]  = nMin+1;
                        newPoint.update_nMinMax( j2, newPoint[j2], this->tol );
#if NO_DIM==3
                        newPoint[j3] = slope3*newPoint[j] + constant3;
                        newPoint.update_nMinMax( j3, newPoint[j3], this->tol );
#endif
                        
                        // reiterate over two new simplices
                        simplices[total] = *pos;
                        simplices[total][i1] = newPoint;
                        simplices[total].vol = simplices[total].volume();
                        // do next step for first simplex
                        (*pos)[i2] = newPoint;
                        pos->vol -= simplices[total].vol;
                        if ( pos->vol<this->tolVolume )
                        {
                            this->simplexContribution( *pos, pos->vol ); // the simplex does not split any more, so compute it's contribution since below volume tolerance
                            simplexFinished = true;
                        }
                        // do next step for second simplex
                        if ( simplices[total].vol<this->tolVolume )
                            this->simplexContribution( simplices[total], simplices[total].vol ); // the simplex does not split any more, so compute it's contribution since below volume tolerance
                        else
                            ++total;    // added an additional simplex that needs to be split
                        
                        
                        // continue to a different simplex since this one was split already
                        simplexSplit = true;
                        break;
                    }
                }
                if ( simplexSplit ) break;
            }
            
            if ( not simplexSplit ) // the simplex does not split any more, so compute it's contribution
            {
                this->simplexContribution( *pos, pos->vol );
                simplexFinished = true;
            }
            
            if ( simplexFinished ) // go to the next simplex => increase current by 1
                if ( ++current==total )   // there are no more simplices left
                    return;
        }
    }
    
    
    // returns the grid index along a coordinate for the main grid
    inline int mainGridIndex(Real const x, int const axis)
    {
        Real temp = (x - this->start[axis]) / this->dx[axis];
        if (temp>=Real(0.))
            return int( temp );
        else
            return int( temp ) - 1;
    }
    
    // finds the grid that fully encompasses the simplex
    void findBoundingGrid(Simplex &pos,
                          vector<size_t> &indices)
    {
        int nMin[NO_DIM], nMax[NO_DIM];
        for (int j=0; j<NO_DIM; ++j)
        {
            nMin[j] = this->mainGridIndex( pos[0][j], j );
            nMax[j] = nMin[j];
        }
        
        // find the min and max extensions of the grid along each dimension
        for (int i=1; i<NO_DIM+1; ++i)
            for (int j=0; j<NO_DIM; ++j)
            {
                int temp = this->mainGridIndex( pos[i][j], j );
                if (temp<nMin[j]) nMin[j] = temp;
                if (temp>nMax[j]) nMax[j] = temp;
            }
        
        // add an extra cell for the upper bound
        for (int j=0; j<NO_DIM; ++j)
            nMax[j] += 1;
        
        // initialize the new grid
        int newSize = 1;
        for (int j=0; j<NO_DIM; ++j)
        {
            this->currentStart[j] = this->start[j] + nMin[j] * this->dx[j];
            this->newGrid[j] =  nMax[j] -  nMin[j];
            newSize *= this->newGrid[j];
        }
        
        // compute the indices of the large grid corresponding to the smaller grid
        indices.assign( newSize, size_t(-1) );
        for (int i0=nMin[0], n0=0; i0<nMax[0]; ++i0, ++n0)
        {
            if ( i0<0 or i0>=grid[0] ) continue;
            for (int i1=nMin[1], n1=0; i1<nMax[1]; ++i1, ++n1)
            {
                if ( i1<0 or i1>=grid[1] ) continue;
#if NO_DIM==2
                indices[ n0*this->newGrid[1]+n1 ] = i0*this->grid[1]+i1;
#elif NO_DIM==3
                for (int i2=nMin[2], n2=0; i2<nMax[2]; ++i2, ++n2)
                {
                    if ( i2<0 or i2>=grid[2] ) continue;
                    indices[ n0*this->newGrid[1]*this->newGrid[2]+n1*this->newGrid[2]+n2 ] = i0*this->grid[1]*this->grid[2]+i1*this->grid[2]+i2;
                }
#endif
            }
        }
        
        //change the vertex positions to the new start of the grid
        for (int i=0; i<NO_DIM+1; ++i)
            for (int j=0; j<NO_DIM; ++j)
            {
                pos[i][j] = (pos[i][j] - this->currentStart[j]) / this->dx[j];
                pos[i].update_nMinMax( j, pos[i][j], this->tol );
            }
    }
    
    
    // this function takes care to initialize the computations for the start of each simplex
    void findIntersection(Simplex &pos,
                          vector<size_t> &indices,
                          vector< Pvector<Real,NO_DIM+1> > &contributions)
    {
        this->results = &contributions;
        
        // find the size of the smallest grid that fully encompasses the input simplex
        this->findBoundingGrid( pos, indices );
        contributions.assign( indices.size(), Pvector<Real,NO_DIM+1>::zero() );
        
        // call the function that computes the volume intersection
        Real vol = pos.volume();
        this->tolVolume = Real(1.)<vol ? VOLUME_TOL : VOLUME_TOL*vol;
        if ( this->tolVolume<Real(1.e-6) )  this->tolVolume = Real(1.e-6);
        this->simplexIteration2( pos, vol );
    }
};


}   //end of namespace VOLUME_SPLIT
#endif
