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


/*
  This class implements a physical vector. The components of the vector can be access as the elements on an array (i.e. using []) with indices from 0 to N-1, where N is the number of components of the vector.
The following mathematical operations are supported:
        Pvector + Pvector       = vector addition 
        Pvector += Pvector      = vector addition with the left side of += (same as 'int i+= 2;')
        Pvector - Pvector       = vector substraction
        Pvector -= Pvector      = vector substraction
        Pvector * scalar        = vector multiplication with scalar value
        Pvector *= scalar       = vector multiplication with scalar value
        Pvector / scalar        = vector division via scalar value
        Pvector /= scalar       = vector division via scalar value
        
        Pvector == Pvector      = check if the two vectors are the same (NOTE: doesnt not make sense for float and double since one should add a precision level - this is not included)
        Pvector != Pvector      = not ( Pvector == Pvector )
        zero()                  = returns a null vector (all entries are 0.)
*/


#ifndef PVECTOR_HEADER
#define PVECTOR_HEADER

template <typename T, size_t N>
struct Pvector
{
    T x[N];
    
    Pvector() {
        for (size_t i=0; i<N; ++i)
            x[i] = T(0); }

#if NO_DIM==2
    Pvector(T x1, T x2) {
        x[0] = x1;
        x[1] = x2;}
#elif NO_DIM==3
    Pvector(T x1, T x2, T x3) {
        x[0] = x1;
        x[1] = x2;
        x[2] = x3;}
#endif
    Pvector(T *xVector) {
        for (size_t i=0; i<N; ++i)
            x[i] = xVector[i]; }
    
    /* Overload the [] operator. */
    T& operator[](size_t i)
    {
        return x[i];
    }
    
   
    Pvector operator +(Pvector<T,N> other)    // Overload the + operator.
    {
        Pvector<T,N> result;
        for (size_t i=0; i<N; ++i )
            result[i] = this->x[i] + other[i];
        return result;
    }
    Pvector& operator +=(Pvector<T,N> other)    // Overload the += operator.
    {
        for (size_t i=0; i<N; ++i )
            x[i] += other[i];
        return *this;
    }
    Pvector operator -(Pvector<T,N> other)    // Overload the - operator.
    {
        Pvector<T,N> result;
        for (size_t i=0; i<N; ++i )
            result[i] = this->x[i] - other[i];
        return result;
    }
    Pvector& operator -=(Pvector<T,N> other)    // Overload the += operator.
    {
        for (size_t i=0; i<N; ++i )
            x[i] -= other[i];
        return *this;
    }
    Pvector operator *(T other)    // Overload the * operator.
    {
        Pvector<T,N> result;
        for (size_t i=0; i<N; ++i )
            result[i] = this->x[i] * other;
        return result;
    }
    Pvector& operator *=(T other)    // Overload the *= operator.
    {
        for (size_t i=0; i<N; ++i )
            x[i] *= other;
        return *this;
    }
    Pvector operator /(T other)    // Overload the / operator.
    {
        Pvector<T,N> result;
        for (size_t i=0; i<N; ++i )
            result[i] = this->x[i] / other;
        return result;
    }
    Pvector& operator /=(T other)    // Overload the /= operator.
    {
        for (size_t i=0; i<N; ++i )
            x[i] /= other;
        return *this;
    }
    
    bool operator ==(Pvector<T,N> other)    // Overload the == operator.
    {
        bool equal = true;
        for (size_t i=0; i<N; ++i )
            if( not x[i] == other[i] )
            {
                equal = false; 
                break;
            }
        return equal;
    }
    bool operator !=(Pvector<T,N> other)    // Overload the != operator.
    {
        return !( (*this)==other );
    }
    
    
    static Pvector zero()      // Creates a null Pvector
    {
        Pvector<T,N> result;
        for (size_t i=0; i<N; ++i )
            result[i] = T(0.);
        return result;
    }
};



#endif
