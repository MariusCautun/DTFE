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
    These classes implement messages for the user during the runtime of the program.
    Depending on the 'verboseLevel', it does the following:
        verboseLevel = 3 - prints all messages; including 'progressMessages'
        verboseLevel = 2 - prints only user messages, warnings and errors
        verboseLevel = 1 - prints only Warnings and Error messages
        verboseLevel = 0 - prints only Error messages
        
NOTE: This class has been design such that the messages are sent as 'boots::tuple' of elements such that it can output complex messages. Example of usage:
    Message m( verboseLevel );
    m.error( make_tuple("Test number ", 10, ". Testing the message interface ", "etc...", ... ) );
*/


#ifndef MESSAGE_HEADER
#define MESSAGE_HEADER

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cstdlib>


#define OUT std::cout


namespace MESSAGE
{


/* Structure to let the program know the end of an error message.
NOTE: Must be always called to stop the program. */
struct EndErrorMessage
{
    inline friend std::ostream& operator << (std::ostream& out, EndErrorMessage &end)
    {
        out << " The program ended unsucessfully!\n\n";
        exit( EXIT_FAILURE );
        return out;
    }
};


/* Structure to let the program know the end of a warning message.
NOTE: Can easily be modified to end the program - see class 'EndErrorMessage'. */
struct EndWarningMessage
{
    inline friend std::ostream& operator << (std::ostream& out, EndWarningMessage &end)
    {
        out << "\n\n" << std::flush;
        return out;
    }
};


/* Flushes a message from the buffer. */
struct FlushMessage
{
    inline friend std::ostream& operator << (std::ostream& out, FlushMessage &end)
    {
        out << std::flush;
        return out;
    }
};



/* Sends an error message to the user. This class will always send messages, indifferently of the value of the verbose level.
NOTE: Must end the error message with 'EndError' or 'MESSAGE::EndErrorMessage' for the program to terminate. */
struct Error
{
    int _verboseLevel;  //this has no effect on the class behavior, offered just for consistency
    int iterationNo;
    
    // class constructors
    Error(void)
    {
        _verboseLevel = 0;
        iterationNo = 0;
    }
    Error(int const verboseLevel)
    {
        if ( verboseLevel<0 )
            OUT << "\n\n~~~ ERROR ~~~ " << "The constructor of the 'Error' class must receive an integer value >= " << 0 << ". This represents the verbose level of the program - check 'message.h' for additional details." << " The program ended unsucessfully!\n";
        _verboseLevel = verboseLevel;
        iterationNo = 0;
    }
    
    // Show error message and stop program execution
    //NOTE: Must end the error message with 'EndError' or 'MESSAGE::EndErrorMessage' for the program to terminate
    template <typename T>
    inline Error& operator << (T right)
    {
        if ( (iterationNo++)==0 )
            OUT << "\n\n~~~ ERROR ~~~ ";
        OUT << right;
        return *this;
    }
};


/* Sends a warning message to the user. This class will send messages only if verbose level >= 1.
NOTE: Must end the warning message with 'EndWarning' or 'MESSAGE::EndWarningMessage' for the warning message to end accordingly. */
struct Warning
{
    int _verboseLevel;
    int iterationNo;
    
    // class constructor
    Warning(int const verboseLevel)
    {
        if ( verboseLevel<0 )
        {
            Error error;
            error << "The constructor of the 'Warning' class must receive an integer value >= " << 0 << ". This represents the verbose level of the program - check 'message.h' for additional details." << EndErrorMessage();
        }
        _verboseLevel = verboseLevel;
        iterationNo = 0;
    }
    
    // Show warning message
    // NOTE: Must end the warning message with 'EndWarning' or 'MESSAGE::EndWarningMessage' for the warning message to end accordingly.
    template <typename T>
    inline Warning& operator << (T right)
    {
        if ( _verboseLevel<1 )
            return *this;
        if ( (iterationNo++)==0 )
            OUT << "\n\n~~~ WARNING ~~~ ";
        OUT << right;
        return *this;
    }
};


/* Sends a message to the user. */
struct Message
{
    int _verboseLevel;
    int _percentangeDone;
    
    
    // Class constructor
    Message(int const verboseLevel)
    {
        if ( verboseLevel<0 )
        {
            Error error;
            error << "The constructor of the 'Message' class must receive an integer value >= " << 0 << ". This represents the verbose level of the program - check 'message.h' for additional details." << EndErrorMessage();
        }
        _verboseLevel = verboseLevel;
        _percentangeDone = 0;
    }
    
    // Overaloads the << which is used to send messages to the user
    template <typename T>
    inline Message& operator << (T right)
    {
        if ( _verboseLevel>=2 )
            OUT << right;
        return *this;
    }
    
    // Show progress to the user
    inline void updateProgress(int const percentageDone)
    {
        if ( _verboseLevel<3 or _percentangeDone==percentageDone )
            return;
        _percentangeDone = percentageDone;
        char const del[4] = "\b\b\b";
        OUT << std::setw(2) << _percentangeDone << "\%" << std::flush << del;
    }
};



/* Prints to a string the components of an array/vector using the second string to separate the values. */
template <typename T>
std::string printElements(T *ptr, size_t const size, std::string separation=" ")
{
    if ( size<1 ) return "";
    
    std::ostringstream buffer;
    buffer << ptr[0];
    for (size_t i=1; i<size; ++i)
        buffer << separation << ptr[i];
    return buffer.str();
}
template <typename T>
std::string printElements(T &vec, std::string separation=" ")
{
    return printElements( &(vec[0]), vec.size(), separation );
}


static const EndErrorMessage   EndError   = EndErrorMessage();
static const EndWarningMessage EndWarning = EndWarningMessage();
static const FlushMessage      Flush      = FlushMessage();
} // end namespace MESSAGE




//! Template functions for outputing errors
/* Throw error and stop the program.
NOTE: The function 'throwError' can take from 1 to 5 parameters, otherwise one will need to call directly the MESSAGE::Error class. */
template <typename T1, typename T2, typename T3, typename T4, typename T5>
void throwError(T1 message_1,
                T2 message_2,
                T3 message_3,
                T4 message_4,
                T5 message_5)
{
    MESSAGE::Error error;
    error << message_1 << message_2 << message_3 << message_4 <<  message_5 << MESSAGE::EndError;
}
template <typename T1, typename T2, typename T3, typename T4>
void throwError(T1 m1, T2 m2, T3 m3, T4 m4)
{
    throwError( m1, m2, m3, m4, "" );
}
template <typename T1, typename T2, typename T3>
void throwError(T1 m1, T2 m2, T3 m3)
{
    throwError( m1, m2, m3, "", "" );
}
template <typename T1, typename T2>
void throwError(T1 m1, T2 m2)
{
    throwError( m1, m2, "", "", "" );
}
template <typename T1>
void throwError(T1 m1)
{
    throwError( m1, "", "", "", "" );
}


#endif
