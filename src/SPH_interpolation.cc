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


#include <vector>
#include <cmath>

// the next line uses boost multi_array library
#include "kdtree/kdtree2.hpp"
// #include "kdtree/kdtree2.cpp"

#include "define.h"
#include "particle_data.h"
#include "quantities.h"
#include "box.h"
#include "user_options.h"
#include "message.h"




void SPH_interpolation(vector<Particle_data> &p,
                        User_options &userOptions,
                        kdtree2_array &gridPoints,
                        size_t const totalGrid,
                        Quantities *q);




/* This function interpolates the density and velocity to grid using the SPH (smoothed particle hydrodynamics) method. */
void SPH_interpolation(vector<Particle_data> *particles,
                       vector<Sample_point> &samples,
                       User_options &userOptions,
                       Quantities *q)
{
    // compute the positions of the grid sampling points
    size_t totalGrid = 0;
    if ( not samples.empty() )
        totalGrid = samples.size();
    else
    {
        size_t const *grid = &(userOptions.gridSize[0]);
        totalGrid = (NO_DIM==2) ? grid[0]*grid[1] : grid[0]*grid[1]*grid[2];
    }
    
    kdtree2_array gridPoints(extents[totalGrid][NO_DIM]);
    if ( not samples.empty() )      // user defined sample points
    {
        for (size_t i=0; i<totalGrid; ++i)
            for (size_t j=0; j<NO_DIM; ++j)
                gridPoints[i][j] = samples[i].position(j);
    }
    else if ( not userOptions.redshiftConeOn )     // sample points on a rectangular grid
    {
        size_t const *grid = &(userOptions.gridSize[0]);
        Box box = userOptions.region;    // the box coordinates in which the density is computed
        Real y[NO_DIM];
        Real dy[NO_DIM];
        for (size_t i=0; i<NO_DIM; ++i)
            dy[i] = (box[2*i+1] - box[2*i]) / grid[i];
        
        size_t index = 0;
        y[0] = box[0] + 0.5*dy[0];
        for (size_t i1=0; i1<grid[0]; ++i1)
        {
            y[1] = box[2] + 0.5*dy[1];
            for (size_t i2=0; i2<grid[1]; ++i2)
            {
#if NO_DIM==3
                y[2] = box[4] + 0.5*dy[2];
                for (size_t i3=0; i3<grid[2]; ++i3)
                {
#endif
                    for (size_t j=0; j<NO_DIM; ++j)
                        gridPoints[index][j] = y[j];
                    ++index;
#if NO_DIM==3
                    y[2] += dy[2];
                }
#endif
                y[1] += dy[1];
            }
            y[0] += dy[0];
        }
    }
    else if ( userOptions.redshiftConeOn )     //sample points on a light cone grid
    {
        size_t const *grid = &(userOptions.gridSize[0]);
        Box box = userOptions.redshiftCone;    // the box coordinates in which the density is computed
        Real const *origin = &(userOptions.originPosition[0]);
        Real dx[NO_DIM];
        for (size_t i=0; i<NO_DIM; ++i)
            dx[i] = (box[2*i+1] - box[2*i]) / grid[i];
        
        size_t index = 0;
        Real r = box[0] + dx[0]/2.;
        for (size_t i=0; i<grid[0]; ++i)
        {
            Real theta = box[2] + dx[1]/2.;
            for (size_t j=0; j<grid[1]; ++j)
            {
#if NO_DIM==2
                gridPoints[index][0] = origin[0] + r*std::cos(theta);
                gridPoints[index++][1] = origin[1] + r*std::sin(theta);
#elif NO_DIM==3
                Real phi = box[4] + dx[2]/2.;
                for (size_t k=0; k<grid[2]; ++k)
                {
                    gridPoints[index][0] = origin[0] + r*sin(theta)*cos(phi);
                    gridPoints[index][1] = origin[1] + r*sin(theta)*sin(phi);
                    gridPoints[index++][2] = origin[2] + r*cos(theta);
                    phi += dx[2];
                }
#endif
                theta += dx[1];
            }
            r += dx[0];
        }
    }
    
    
    
    // now interpolate the results to a grid
    SPH_interpolation( *particles, userOptions, gridPoints, totalGrid, q );
    particles->clear();
}



/* Computes the values of the SPH smoothing kernel (except the h^-2 (in 2D) and h^-3 (in 3D) normalizeation factor). */
Real SPH_smoothingKernel(Real x)
{
#if NO_DIM==2
    Real factor = Real(10./(7.*PI));
#elif NO_DIM==3
    Real factor = Real(1./PI);
#endif
    
    if ( std::fabs(x)<=1. )
        return factor * Real(1. - 1.5*x*x + 0.75*x*x*x);
    else if ( std::fabs(x)<=2. )
    {
        Real temp = Real( 2.-x );
        return factor * Real(0.25 * temp*temp*temp);
    }
    else
        return Real(0.);
}

/* Computes the derivative of the SPH smoothing kernel (except the h^-2 (in 2D) and h^-3 (in 3D) normalizeation factor). */
Real SPH_smoothingKernelDerivative(Real x)
{
#if NO_DIM==2
    Real factor = Real(10./(7.*PI));
#elif NO_DIM==3
    Real factor = Real(1./PI);
#endif
    
    if ( std::fabs(x)<=1. )
        return factor * Real(-3.*x + 2.25*x*x);
    else if ( std::fabs(x)<=2. )
    {
        Real temp = Real( 2.-x );
        return factor * Real(-0.75 * temp*temp);
    }
    else
        return Real(0.);
}


/* Function that returns the constant in front of W(r,h), constant that depends only on h. */
inline Real hFactor(Real h)
{
#if NO_DIM==2
    return Real(1.)/(h*h);
#elif NO_DIM==3
    return Real(1.)/(h*h*h);
#endif
}


/* This function computes the derivative of a vector in the SPH method. 
It returns:
    1/2 m W_deriv(r,h) * vec * \vec{r}/r
where:
    factor = 1/2 m W_deriv(r,h) .
*/
template <typename T, size_t Nvector, size_t Nderiv>
Pvector<T,Nderiv> getSPH_derivative(T const factor,
                                    Pvector<T,Nvector> &vec,
                                    T *pos1,
                                    T *pos2)
{
    Pvector<T,Nderiv> result;
    T vecR[NO_DIM], res1 = T(0.), res2 = T(0.);
    for (int i=0; i<NO_DIM; ++i)
    {
        vecR[i] = pos1[i] - pos2[i];
        res1 += pos1[i]*pos1[i];
        res2 += vecR[i]*vecR[i];
    }
    
    T res = std::sqrt( res2 );
    if ( res2/res1<T(1.e-6) )   // if res2==0
        for (int i=0; i<NO_DIM; ++i)
            vecR[i] = T(0.);
    else
        for (int i=0; i<NO_DIM; ++i)
            vecR[i] /= res;
        
    for (int i=0; i<Nvector; ++i)
        for (int j=0; j<NO_DIM; ++j)
            result[i*NO_DIM+j] = factor * vec[i] * vecR[j];
    return result;
}




/* This function uses the SPH method to interpolate quantities to a grid.
NOTE: The algorithm implemented here is described at: http://ciera.northwestern.edu/StarCrash/manual/html/node7.html
It has the following steps:
    
    1) Loop over all particles - for each particle (here denoted by i) define h_i=R/2 where R is the distance to the N-th closest neighbor. Within the same loop assign the W(r_ij,h_i) part of the contribution to the particle density of neighbor j. (The full expression for density is: \rho = \sum_j m_j 1/2 [W(r_ij,h_i) + W(r_ij,h_j)]. )
    
    2) Loop over all the density grid points and assign a smoothing length to each cell (h_cell=R/2 where R is the sphere in which we find the N-th closest particle neighbors). Use this loop to also add the first part of the SPH smoothing kernel (similar to above).
    
    3) Loop again over all particles to assign the 2nd part of the smoothing kernel to the grid quantity. This can be done via two methods:
        - for a regular grid, loop over all the grid points that are with r<h_i for particle i.
        - construct a tree from the grid point distribution and use tree searching algorithms to find all the grid points with h_i distance from particle i.
*/
void SPH_interpolation(vector<Particle_data> &p,
                        User_options &userOptions,
                        kdtree2_array &gridPoints,
                        size_t const totalGrid,
                        Quantities *q)
{
    MESSAGE::Message message(userOptions.verboseLevel);
    message << "\nInterpolating the fields to a grid using the SPH method with " << userOptions.SPH_neighbors << " neighbors. ";
    if ( userOptions.gridSize.empty() )
        message << "The interpolation takes place at " << totalGrid << " user defined sampling points:\n" << MESSAGE::Flush;
    else if ( not userOptions.redshiftConeOn )
        message << "The interpolation takes place inside the box of coordinates " << userOptions.region.print() << " on a " << MESSAGE::printElements( userOptions.gridSize, "*" ) << " grid:\n" << MESSAGE::Flush;
    else if ( userOptions.redshiftConeOn )
        message << "The interpolation takes place inside the spherical coordinates " << userOptions.region.print() << " on a " << MESSAGE::printElements( userOptions.gridSize, "*" ) << " grid:\n" << MESSAGE::Flush;
    
    
    // first construct the kdtree
    message << "Computing the kdtree for the SPH interpolation ... " << MESSAGE::Flush;
    boost::timer t;
    t.restart();
    size_t const noParticles = p.size();
    kdtree2_array dataPoints(extents[noParticles][NO_DIM]);
    for (size_t i=0; i<noParticles; ++i)
        for (size_t j=0; j<NO_DIM; ++j)
            dataPoints[i][j] = p[i].position(j);
    kdtree2* tree;
    tree = new kdtree2( dataPoints, false );
    tree->sort_results = true;
    
    message << "Done.\n" << MESSAGE::Flush;
    printElapsedTime( &t, &userOptions, "kdtree construction" );
    
    
    // construct a table with the values of the smoothing kernel
    int const NTable = 1000;
    Real dx = Real( 2. / NTable);
    Real W[NTable+1], W_deriv[NTable+1];
    for (int i=0; i<=NTable; ++i)
    {
        W[i] = SPH_smoothingKernel( dx*(i+0.5) );
        W_deriv[i] = SPH_smoothingKernelDerivative( dx*(i+0.5) );
    }
    
    
    // find the smoothing length and SPH density associated to each particle
    message << "Computing the smoothing scale and density at each particle position.\n\tDone:  " << MESSAGE::Flush;
    t.restart();
    kdtree2_result_vector result;   //stores the result of the nearest neighbors
    int N = int( userOptions.SPH_neighbors );
    std::vector<Real> smoothingLength( noParticles, Real(0.) ); // stores the smoothing length for each particle
    std::vector<Real> density( noParticles, Real(0.) );         // stores the density at each particle position
    Real *h = &(smoothingLength[0]);
    Real *d = &(density[0]);
    size_t prev = 0, amount100 = 0;
    // first find the smoothing length and density associated to each particle
    for (size_t i=0; i<noParticles; ++i)
    {
        amount100 = (100 * i)/ noParticles;
        if (prev < amount100)
            message.updateProgress( ++prev );
        
        tree->n_nearest_around_point(i,0,N,result); // returns the neighbors in ascending order according to the distance they are at
        h[i] = Real(std::sqrt( result[N-1].dis ) / 2.); // choose h such that within 2h there are N neighbors
        
        // get the contribution of this particle to the other particles' density - use only the W(r.h_i) part, not the full W_ij
        Real const c1 = hFactor(h[i]);
        for (size_t j=0; j<result.size(); ++j)
        {
            Real const temp = std::sqrt( result[j].dis );
            int const id = result[j].idx;
            int bin1 = int( temp / (h[i]*dx));
            
            d[i] += Real(0.5) * p[id].weight() * c1*W[bin1];
            d[id] += Real(0.5) * p[i].weight() * c1*W[bin1];
        }
    }
    message << "100%\n" << MESSAGE::Flush;
    printElapsedTime( &t, &userOptions, "SPH particle density" );
    
    
    
    // computes the quantities of interest at grid points
    message << "Computing the interpolated fields on the grid.\n\tDone:  " << MESSAGE::Flush;
    t.restart();
    // reserve memory for the output
    q->density.reserve( totalGrid );    //always need memory for the density
    if ( userOptions.aField.velocity ) q->velocity.reserve( totalGrid );
    if ( userOptions.aField.velocity_gradient ) q->velocity_gradient.reserve( totalGrid );
    if ( userOptions.aField.scalar ) q->scalar.reserve( totalGrid );
    if ( userOptions.aField.scalar_gradient ) q->scalar_gradient.reserve( totalGrid );
    
    
    // Compute the grid smoothing scale and the 1st contribution to the grid interpolated quantities (the W(r_ij,h_i) part)
    std::vector<float> y( NO_DIM, float(0.) );
    prev = 0; amount100 = 0;
    for (size_t i=0; i<totalGrid; ++i)
    {
        amount100 = (50 * i) / totalGrid;
        if (prev < amount100)
            message.updateProgress( ++prev );
        
        for (size_t j=0; j<NO_DIM; ++j)
            y[j] = gridPoints[i][j];
        tree->n_nearest( y, N, result );
        Real const tempH = Real(std::sqrt( result[N-1].dis ) / 2.);
        Real const c1 = hFactor(tempH);
        
        // temporary variable for the different quantities that are computed
        Real resDens = Real(0.);
        Pvector<Real,noVelComp> resVel = Pvector<Real,noVelComp>::zero();
        Pvector<Real,noScalarComp> resIntensive = Pvector<Real,noScalarComp>::zero();
        Pvector<Real,noScalarComp> resExtensive = Pvector<Real,noScalarComp>::zero();
//         Pvector<Real,noGradComp> resVelGrad = Pvector<Real,noGradComp>::zero();
//         Pvector<Real,noScalarGradComp> resGradIntensive = Pvector<Real,noScalarGradComp>::zero();
//         Pvector<Real,noScalarGradComp> resGradExtensive = Pvector<Real,noScalarGradComp>::zero();
        
        for (size_t j=0; j<result.size(); ++j)
        {
            int const id = result[j].idx;
            Real tempR = std::sqrt( result[j].dis );
            int bin1 = int( tempR / (tempH*dx));
            
            Real temp1 = Real(0.5) * p[id].weight() * c1*W[bin1];
            resDens += temp1;
            resVel += p[id].velocity() * temp1;
            resIntensive += p[id].scalar() * temp1;
            resExtensive += p[id].scalar() * temp1 / d[id];
            
//             Real temp2 = Real(0.5) * p[id].weight() * c1/tempH*W_deriv[bin1];
//             resVelGrad += getSPH_derivative<Real,noVelComp,NO_DIM*noVelComp>( temp2, p[id].velocity(), &(y[0]), &(p[id].position()[0]) );
//             resGradIntensive += getSPH_derivative<Real,noScalarComp,NO_DIM*noScalarComp>( temp2, p[id].scalar(), &(y[0]), &(p[id].position()[0]) );
//             resGradExtensive += resGradIntensive / d[id];
        }
        
        q->density.push_back( resDens );
        if ( userOptions.aField.velocity ) q->velocity.push_back( resVel );
//         if ( userOptions.field.velocity_gradient ) q->velocity_gradient.push_back( resVelGrad );
        if ( userOptions.extensive )
        {
            if ( userOptions.aField.scalar ) q->scalar.push_back( resExtensive );
//             if ( userOptions.field.scalar_gradient ) q->scalar_gradient.push_back( resGradExtensive );
        }
        else
        {
            if ( userOptions.aField.scalar ) q->scalar.push_back( resIntensive );
//             if ( userOptions.field.scalar_gradient ) q->scalar_gradient.push_back( resGradIntensive );
        }
    }
    
    // compute the grid points tree
    delete tree;
    tree = new kdtree2( gridPoints, false );
    tree->sort_results = true;
    for (size_t i=0; i<noParticles; ++i)
    {
        amount100 = 50 + (50 * i) / noParticles;
        if (prev < amount100)
            message.updateProgress( ++prev );
        
        for (size_t j=0; j<NO_DIM; ++j)
            y[j] = p[i].position(j);
        Real const distance = Real(4.) * h[i]*h[i]; //(2h)^2
        tree->r_nearest( y, distance, result );    //now locate all the grid points within the smoothing radius of the given particle
        
        Real const c1 = Real(0.5) * p[i].weight() * hFactor(h[i]);
        for (size_t j=0; j<result.size(); ++j)
        {
            int const id = result[j].idx;
            Real y2[NO_DIM];
            for (int i1=0; i1<NO_DIM; ++i1)
                y2[i1] = gridPoints[id][i1];
            Real tempR = std::sqrt( result[j].dis );
            int bin1 = int( tempR / (h[i]*dx));
            
            Real temp1 = c1*W[bin1];
//            Real temp2 = c1/h[i]*W_deriv[bin1];
            q->density[id] += temp1;
            if ( userOptions.aField.velocity ) q->velocity[id] += p[i].velocity() * temp1;
//             if ( userOptions.field.velocity_gradient ) q->velocity_gradient[id] += getSPH_derivative<Real,noVelComp,NO_DIM*noVelComp>( temp2, p[i].velocity(), y2, &(y[0]) );
            if ( userOptions.extensive )
            {
                if ( userOptions.aField.scalar ) q->scalar[id] += p[i].scalar() * temp1;
//                 if ( userOptions.field.scalar_gradient ) q->scalar_gradient[id] += getSPH_derivative<Real,noScalarComp,NO_DIM*noScalarComp>( temp2, p[i].scalar(), y2, &(y[0]) );
            }
            else
            {
                if ( userOptions.aField.scalar ) q->scalar[id] += p[i].scalar() * temp1 / d[i];
//                 if ( userOptions.field.scalar_gradient ) q->scalar_gradient[id] += getSPH_derivative<Real,noScalarComp,NO_DIM*noScalarComp>( temp2/d[i], p[i].scalar(), y2, &(y[0]) );
            }
        }
    }
    delete tree;
    
    
    // if compute velocity, divide by density at grid cell
    if ( userOptions.aField.velocity )
        for (size_t i=0; i<totalGrid; ++i)
            q->velocity[i] /= q->density[i];
//     if ( userOptions.field.velocity_gradient )
//         for (size_t i=0; i<totalGrid; ++i)
//             q->velocity_gradient[i] /= q->density[i];
    
    // if compute scalar field and sclara field is intensive, divide by density at grid cell
    if ( userOptions.aField.scalar and not userOptions.extensive )
        for (size_t i=0; i<totalGrid; ++i)
            q->scalar[i] /= q->density[i];
//     if ( userOptions.field.scalar_gradient and not userOptions.extensive )
//         for (size_t i=0; i<totalGrid; ++i)
//             q->scalar_gradient[i] /= q->density[i];
    
    // if density computation - multiply the density by the normalization factor
    if ( userOptions.aField.density )
    {
        Real factor = Real( 1. / userOptions.averageDensity );    //normalization factor for the density
        for (std::vector<Real>::iterator it=q->density.begin(); it!=q->density.end(); ++it)
            (*it) *= factor;
    }
    else q->density.clear();
    
    
    message << "100%\n" << MESSAGE::Flush;
    printElapsedTime( &t, &userOptions, "SPH grid interpolation" );
}














// /* This function uses the SPH method to interpolate quantities to a grid.
// NOTE: this is an incomplete implementation of the SPH method.
// */
// void SPH_interpolation_regular_grid(vector<Particle_data> &p,
//                                     User_options &userOptions,
//                                     Quantities *q)
// {
//     MESSAGE::Message message(userOptions.verboseLevel);
//     message << "\nInterpolating the fields to a grid using the SPH method with " << userOptions.SPH_neighbors << " neighbors. The interpolation takes place inside the box of coordinates " << userOptions.region.print() << " on a " << MESSAGE::printElements( userOptions.gridSize, "*" ) << " grid ... " << MESSAGE::Flush;
//     
//     
//     // first construct the kdtree
//     size_t const noParticles = p.size();
//     kdtree2_array dataPoints(extents[noParticles][NO_DIM]);
//     for (size_t i=0; i<noParticles; ++i)
//         for (size_t j=0; j<NO_DIM; ++j)
//             dataPoints[i][j] = p[i].position(j);
//     kdtree2* tree;
//     tree = new kdtree2( dataPoints, true );
//     tree->sort_results = true;
//     
//     
//     // construct a table with the values of the smoothing kernel
//     int const NTable = 1000;
//     Real dx = Real( 2. / NTable);
//     Real W[NTable+1];
//     for (int i=0; i<=NTable; ++i)
//         W[i] = SPH_smoothingKernel( dx*(i+0.5) );
//     
//     
//     message << "\nComputing the smoothing scale and density at each particle position.\n\tDone:  " << MESSAGE::Flush;
//     boost::timer t;
//     t.restart();
//     // find the smoothing length and SPH density associated to each particle
//     kdtree2_result_vector result;   //stores the result of the nearest neighbors
//     int N = int( userOptions.SPH_neighbors );
//     std::vector<Real> smoothingLength( noParticles, Real(0.) ); // stores the smoothing length for each particle
//     std::vector<Real> density( noParticles, Real(0.) );         // stores the density at each particle position
//     Real *h = &(smoothingLength[0]);
//     Real *d = &(density[0]);
//     size_t prev = 0, amount100 = 0;
//     // first find the smoothing length associated to each particle
//     for (size_t i=0; i<noParticles; ++i)
//     {
//         amount100 = (50 * i)/ noParticles;
//         if (prev < amount100)
//             message.updateProgress( ++prev );
//         
//         tree->n_nearest_around_point(i,0,N,result); // returns the neighbors in ascending order according to the distance they are at
//         h[i] = Real(std::sqrt( result[N-1].dis ) / 2.); // choose h such that within 2h there are N neighbors
//     }
//     // find the density at each particle position
//     for (size_t i=0; i<noParticles; ++i)
//     {
//         amount100 = 50 + (50 * i)/ noParticles;
//         if (prev < amount100)
//             message.updateProgress( ++prev );
//         
//         tree->n_nearest_around_point(i,0,N,result);
//         Real const c1 = hFactor(h[i]);
//         for (size_t j=0; j<result.size(); ++j)
//         {
//             Real const c2 = hFactor(h[j]);
//             
//             Real const temp = std::sqrt( result[j].dis );
//             int const id = result[j].idx;
//             int bin1 = int( temp / (h[i]*dx));
//             if ( bin1>NTable ) bin1 = NTable;   // if this is true, than 'W[bin1]=0', but also 'W[NTable]=0', hence we can have 'bin1=NTable'
//             int bin2 = int( temp / (h[id]*dx));
//             if ( bin2>NTable ) bin2 = NTable;
//             
//             d[i] += Real(0.5) * p[id].weight() * (c1*W[bin1] + c2*W[bin2]);
//         }
//     }
//     message << "100%\n" << MESSAGE::Flush;
//     printElapsedTime( &t, &userOptions, "SPH particle density" );
//     
//     
//     message << "Computing the interpolated fields on the grid.\n\tDone:  " << MESSAGE::Flush;
//     t.restart();
//     // computes the quantities of interest at grid points
//     size_t const *grid = &(userOptions.gridSize[0]);
//     Box box = userOptions.region;    // the box coordinates in which the density is computed
//     std::vector<float> y( NO_DIM, float(0.) );
//     Real dy[NO_DIM];
//     for (size_t i=0; i<NO_DIM; ++i)
//         dy[i] = (box[2*i+1] - box[2*i]) / grid[i];
//     // reserve memory for the output
//     size_t totalGrid = (NO_DIM==2) ? grid[0]*grid[1] : grid[0]*grid[1]*grid[2];
//     if ( userOptions.field.density ) q->density.reserve( totalGrid );
//     if ( userOptions.field.velocity ) q->velocity.reserve( totalGrid );
//     if ( userOptions.field.scalar ) q->scalar.reserve( totalGrid );
//     Real factor = Real( 1. / userOptions.averageDensity );    //normalization factor for the density
//     
//     prev = 0; amount100 = 0;
//     y[0] = 0.5*dy[0];
//     for (size_t i1=0; i1<grid[0]; ++i1)
//     {
//         amount100 = (100 * i1)/ grid[0];
//         if (prev < amount100)
//             message.updateProgress( ++prev );
//         
//         y[1] = 0.5*dy[1];
//         for (size_t i2=0; i2<grid[1]; ++i2)
//         {
// #if NO_DIM==3
//             y[2] = 0.5*dy[2];
//             for (size_t i3=0; i3<grid[2]; ++i3)
//             {
// #endif
//                 tree->n_nearest( y, N, result );
//                 Real const tempH = Real(std::sqrt( result[N-1].dis ) / 2.);
//                 Real const c1 = hFactor(tempH);
//                 // temporary variable for the different quantities that are computed
//                 Real resDens = Real(0.);
//                 Pvector<Real,noVelComp> resVel = Pvector<Real,noVelComp>::zero();
//                 Pvector<Real,noScalarComp> resIntensive = Pvector<Real,noScalarComp>::zero();
//                 Pvector<Real,noScalarComp> resExtensive = Pvector<Real,noScalarComp>::zero();
//                 
//                 for (size_t j=0; j<result.size(); ++j)
//                 {
//                     int const id = result[j].idx;
//                     Real const c2 = hFactor(h[id]);
//             
//                     Real temp = std::sqrt( result[j].dis );
//                     int bin1 = int( temp / (tempH*dx));
//                     if ( bin1>NTable ) bin1 = NTable;
//                     int bin2 = int( temp / (h[id]*dx));
//                     if ( bin2>NTable ) bin2 = NTable;
//                     
//                     temp = Real(0.5) * p[id].weight() * (c1*W[bin1] + c2*W[bin2]);
//                     resDens += temp;
//                     resVel += p[id].velocity() * temp;
//                     resIntensive += p[id].scalar() * temp;
//                     resExtensive += p[id].scalar() * temp / d[id];
//                 }
//                 
//                 if ( userOptions.field.density ) q->density.push_back( resDens * factor );
//                 if ( userOptions.field.velocity ) q->velocity.push_back( resVel / resDens );
//                 if ( userOptions.field.scalar and userOptions.extensive ) q->scalar.push_back( resExtensive );
//                 else if ( userOptions.field.scalar ) q->scalar.push_back( resIntensive / resDens );
//                 
// #if NO_DIM==3
//                 y[2] += dy[2];
//             }
// #endif
//             y[1] += dy[1];
//         }
//         y[0] += dy[0];
//     }
//     message << "100%\n" << MESSAGE::Flush;
//     printElapsedTime( &t, &userOptions, "SPH grid interpolation" );
// }


