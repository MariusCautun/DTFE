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

#include "define.h"
#include "particle_data.h"
#include "user_options.h"
#include "message.h"

#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
typedef boost::mt19937 base_generator_type;

typedef std::vector<Particle_data>::iterator   vectorIterator;


/* This function selects from the input data a random subsample, with the subsample size given by 'User_options::randomSample'. */
void randomSample(std::vector<Particle_data> particles,
                  std::vector<Particle_data> *subSample,
                  User_options userOptions )
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "\nGenerating a random subsample of the data set of size " << userOptions.randomSample << " using the radom generator seed " << userOptions.randomSeed << " ... " << MESSAGE::Flush;
    
    
    size_t noParticles = size_t( particles.size() * (userOptions.randomSample+.001) );
    subSample->reserve(noParticles);
    
    size_t const seed = userOptions.randomSeed;
    base_generator_type generator( seed );
    boost::uniform_real<> uni_dist(0.,1.);
    boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
    double threshold = double(userOptions.randomSample);
    
    for (vectorIterator it=particles.begin(); it!=particles.end(); ++it)
    {
        if ( uni()<=threshold )
            subSample->push_back( *it );
    }
    
    message << "Done.\n\tKept " << subSample->size() << " particles which represent " << std::setprecision(4) << 100.*Real(subSample->size())/ Real(particles.size()) << "\% of the initial data size.\n\n" << MESSAGE::Flush;
}



/* This function generates random particles in a box of unit length. */
void randomParticles(std::vector<Particle_data> *particles,
                     User_options *userOptions )
{
    size_t noParticles = 1;
    for (int i=0; i<NO_DIM; ++i)
        noParticles *= userOptions->poisson;
    particles->reserve( noParticles );
    MESSAGE::Message message( userOptions->verboseLevel );
    message << "\nGenerating a Poisson distribution of " << noParticles << " particles of equal weight in a box of size unity ... " << MESSAGE::Flush;
    
    size_t const seed = userOptions->randomSeed + size_t(1.e8);
    base_generator_type generator( seed );
    boost::uniform_real<> uni_dist(0.,1.);
    boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
    
    for (size_t i=0; i<noParticles; ++i)
    {
        Particle_data temp;
        temp.weight() = 1.;
        for (int j=0; j<NO_DIM; ++j)
            temp.pos[j] = uni();
        
        particles->push_back( temp );
    }
    
    for (int i=0; i<NO_DIM; ++i)
    {
        userOptions->boxCoordinates[2*i] = 0.;
        userOptions->boxCoordinates[2*i+1] = 1.;
    }
    
    message << "Done.\n\n" << MESSAGE::Flush;
}


