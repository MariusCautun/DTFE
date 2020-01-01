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



#ifndef MISCELLANEOUS_HEADER
#define MISCELLANEOUS_HEADER

#include <string>
#include <cmath>

#include "message.h"



/* Checks if a number is within a given interval. */
template <typename T1>
void intervalCheck(T1 const target,
                   T1 const minValue, T1 const maxValue,
                   std::string errorMessage)
{
    if ( target<minValue or target>maxValue )
    {
        MESSAGE::Error error;
        error << "Some program variable failed a consitency check. The variable " << errorMessage << " has the value " << target << ", but it should be between " << minValue << " to " << maxValue << " ." << MESSAGE::EndError;
    }
}
/* Checks if a number is >= than a lower bound. */
template <typename T1>
void lowerBoundCheck(T1 const target,
                    T1 const minValue,
                    std::string errorMessage)
{
    if ( target<minValue )
    {
        MESSAGE::Error error;
        error << "Some program variable failed a consitency check. The variable " << errorMessage << " has the value " << target << ", but it should be larger or equal than " << minValue << " ." << MESSAGE::EndError;
    }
}
/* Checks if a number is <= than an upper bound. */
template <typename T1>
void upperBoundCheck(T1 const target,
                    T1 const maxValue,
                    std::string errorMessage)
{
    if ( target>maxValue )
    {
        MESSAGE::Error error;
        error << "Some program variables failed a consitency check. The variable " << errorMessage << " has the value " << target << ", but it should be smaller or equal than " << maxValue << " ." << MESSAGE::EndError;
    }
}


// The following 4 functions can be use only in conjunction with the 'boost/program_options.hpp' library.
/* Function used to check that 'opt1' and 'opt2' are not specified at the same time. */
template <typename T>
void conflicting_options(const T &vm,
                        const char* opt1, const char* opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted() && vm.count(opt2) && !vm[opt2].defaulted())
    {
        MESSAGE::Error error;
        error << "Conflicting options '" << opt1 << "' and '" << opt2 << "'!\n" << MESSAGE::EndError;
    }
}

/* Function used to check that of 'for_what' is specified, then 'required_option' is specified too. */
template <typename T>
void option_dependency(const T &vm,
                        const char* for_what, const char* required_option)
{
    if (vm.count(for_what) && !vm[for_what].defaulted())
        if (vm.count(required_option) == 0 || vm[required_option].defaulted())
        {
            MESSAGE::Error error;
            error << "Option '" << for_what << "' requires option '" << required_option << "'!\n" << MESSAGE::EndError;
        }
}


/* Function used to check that 'opt1' is not supplied in the absence of 'opt2'. */
template <typename T>
void superfluous_options(const T &vm,
                        const char* opt1, const char* opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted())
        if (vm.count(opt2) == 0 || vm[opt2].defaulted())
        {
            MESSAGE::Error error;
            error << "The option '" << opt1 << "' can be used only in the presence of '" << opt2 << "'. It does not make sense to use '" << opt1 << "' otherwise!\n" << MESSAGE::EndError;
        }
}
/* Function used to check that 'opt1' is not supplied if 'optionsOn' is false.
NOTE: 'OptionOn' should be true if one or several options were suplied. Give the name of the options in 'optionsName'. */
template <typename T>
void superfluous_options(const T &vm,
                        const char* opt1,
                        const bool optionsOn,
                        std::string const optionsName)
{
    if (vm.count(opt1) && !vm[opt1].defaulted())
        if ( not optionsOn )
        {
            MESSAGE::Error error;
            error << "The option '" << opt1 << "' can be used only in the presence of option/s: " << optionsName << ". It does not make sense to use '" << opt1 << "' otherwise!\n" << MESSAGE::EndError;
        }
}




/* Returns the minimum and maximum values of an array. */
template <typename T, typename T_INT> inline T minimum(T *x, T_INT const size )
{
    T temp = x[0];
    for (T_INT i=1; i<size; ++i)
        if ( temp>x[i] )
            temp = x[i];
    return temp;
}
template <typename T, typename T_INT> inline T maximum(T *x, T_INT const size )
{
    T temp = x[0];
    for (T_INT i=1; i<size; ++i)
        if ( temp<x[i] )
            temp = x[i];
    return temp;
}




/* Quicksort algorithm for sorting an array in increasing order according to the elements. It sorts the array only between elements iMin to iMax.
For further details see: http://en.wikipedia.org/wiki/Quicksort
*/
template <typename T1, typename T2>
void quicksort( T1 values[], T2 iMin, T2 iMax )
{
    if ( iMin>=iMax )   //condition to stop the iterative computation
        return;
    
    T1 temp = values[iMax]; //the value against which we compare = the pivot
    T2 i = iMin, i1=iMin, i2=iMax;
    do
    {
        if ( values[i]>temp )   //if value larger than pivot
        {
            if ( i==i1 )    //move value to the right of the pivot
                values[i2] = values[i];
            //now the value is already to the right of the pivot
            --i2;
            i = i2;
        }
        else
        {
            if ( i==i2 )    //move value to the left of the pivot
                values[i1] = values[i];
            ++i1;
            i = i1;
        }
    } while ( i1!=i2 );
    
    //now copy the pivot to i=i1=i2 since that array element is free
    values[i] = temp;
    
    //we still have to order again the elements:
    quicksort( values, iMin, i-1);    //at the left of the pivot
    quicksort( values, i+1, iMax);    //at the right of the pivot
}


/* This function takes root "rootPower" from an integer type and returns also an integer type. If the input integer is input = result^rootPower, than it returns result, otherwise it terminates the program.
*/
template <typename T>
T rootN(T const input,
        T const rootPower)
{
    lowerBoundCheck( input, 1, "'input' in function 'rootN'" );
    lowerBoundCheck( rootPower, 0, "'rootPower' in function 'rootN'" );
    
    T result;
    double temp = pow( input-1., 1./rootPower ); //since pow returns a double, we take the root from 'input-1' to be sure that we get a double smaller than the integer we hope to find
    result = T( temp ) + 1;	//returns integer part of temp + 1
    
    
    // check that indeed input = result^rootPower
    T temp2 = 1;
    for (T i=0; i<rootPower; ++i)
        temp2 *= result;
    if (temp2==input)
        return result;
    else
    {
        MESSAGE::Error error;
        error << "The first argument of the function 'rootN' is " << input << " which cannot be written as " << result << "^" << rootPower << "=" << temp <<"." << MESSAGE::EndError;
    }
    return T(0.);
}

/* This function checks if the root "rootPower" from an integer is an integer (returns true), otherwise returns false. */
template <typename T>
T isRootN(T const input,
        T const rootPower)
{
    lowerBoundCheck( input, 1, "'input' in function 'rootN'" );
    lowerBoundCheck( rootPower, 0, "'rootPower' in function 'rootN'" );
    
    T result;
    double temp = pow( input-1., 1./rootPower ); //since pow returns a double, we take the root from 'input-1' to be sure that we get a double smaller than the integer we hope to find
    result = T( temp ) + 1;	//returns integer part of temp + 1
    
    
    // check that indeed input = result^rootPower
    T temp2 = 1;
    for (T i=0; i<rootPower; ++i)
        temp2 *= result;
    if (temp2==input)
        return true;
    return false;
}


#endif

