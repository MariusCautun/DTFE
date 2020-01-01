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


/** The following functions are defined in "density_interpolation.cc" file. */


 /* Computes the minor determinant coresponding to entry (row,column) of a 3x3 matrix. */
inline Real minorDeterminant(Real matrix[][NO_DIM],
                             size_t const row,
                             size_t const column);

/* Computes the inverse of a 2x2 or 3x3 matrix. */
void matrixInverse(Real matrix[][NO_DIM],
                   Real result[][NO_DIM]);

/* Computes the matrix multiplication of two matrices. */
template <size_t N>
inline void matrixMultiplication(Real M1[][NO_DIM],
                                 Real M2[][N],
                                 Real result[][N])
{
    for (int i=0; i<NO_DIM; ++i)
        for (int j=0; j<N; ++j)
        {
            result[i][j] = 0.;
            for (int k=0; k<NO_DIM; ++k)
                result[i][j] += M1[i][k] * M2[k][j];
        }
}
inline void matrixMultiplication(Real M1[][NO_DIM],
                                 Real *M2,
                                 Real *result)
{
    for (int i=0; i<NO_DIM; ++i)
    {
        result[i] = 0.;
        for (int j=0; j<NO_DIM; ++j)
            result[i] += M1[i][j] * M2[j];
    }
}

/* Computes the vertex position difference matrix for a given Delaunay cell. */
void vertexPositionMatrix(Cell_handle &cell,
                          Real vertexMatrix[][NO_DIM]);

