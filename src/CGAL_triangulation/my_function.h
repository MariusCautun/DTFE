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



/* You can use this file to implement your own function (quantity) that you want to be volume averaged using the DTFE method. */


#ifdef SCALAR

//uncommenting the following uses the function below to assign values to the variable 'scalar', values that will be volume averaged afterwards
//#define MY_SCALAR


/* This function takes as input the density value 'density', the density gradient 'densityGradient' (has NO_DIM components), the velocity 'velocity' (has NO_DIM components) and the velocity gradient 'velocityGradient' (has NO_DIM*NO_DIM components). Compute the new quantity you want and save it in the 'scalar' variable (remeber to choose the correct value for the '-DNO_SCALARS' macro in the Makefile since  'noScalarComp=NO_SCALARS').

The density, density gradient, velocity and velocity gradient are the values of those quantities at the point 'pointPosition'.

The program will then volume average the each component of the 'scalar' variable. You must run the program with the '--field scalar / scalar_a' options for this computation to take part. 
*/
void inline personalizedFunction(Point &pointPosition,
                                 Real density,
                                 Real *densityGradient,
                                 Pvector<Real,noVelComp> &velocity,
                                 Real velocityGradient[][noVelComp],
                                 Pvector<Real,noScalarComp> &scalar)
{
    /* example of possible text */
    scalar[0] = density * densityGradient[1];   // density times the y-component of the density gradient
    scalar[1] = densityGradient[2] * velocity[2]; // product of z-components of density gradient and velocity
    // velocityGradient[i][j] - the entry (i,j) of the velocity gradient matrix = the derivative along direction i (0=x, 1=y, 2=z) of velocity component j (0=x, 1=y, 2=z)
}

#endif
