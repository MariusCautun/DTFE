
#ifndef DEFINE_HEADER
#define DEFINE_HEADER

#include <cstddef>



// choose the number of dimensions
#ifndef NO_DIM
     #define NO_DIM 3
#endif
#if NO_DIM!=2
     #if NO_DIM!=3
          #error "The preprocessor macro 'NO_DIM' can take only the values 2 or 3!"
     #endif
#endif
#define NO_DIM2 (NO_DIM*NO_DIM)



// define some shorthand notations
#ifdef DOUBLE
    typedef double                         Real;
#else
    typedef float                          Real;
#endif



/* define program constants giving the number of components for the velocity and scalar point/particle properties
NOTE: Do not modify the following definitions. */
static const size_t noVelComp = NO_DIM;                     // number of velocity components
static const size_t noGradComp = NO_DIM * NO_DIM;           // number of velocity gradient components
static const size_t noShearComp = (NO_DIM*(NO_DIM+1))/2-1;  /* number of independent velocity shear components (2 for 2D and 5 for 3D) - '-1' comes because the shear is traceless (no need to store component yy/zz for 2D/3D) */
static const size_t noVortComp = ((NO_DIM-1)*NO_DIM)/2;     // number of independent velocity vorticity components (1 for 2D and 3 for 3D)

#ifndef NO_SCALARS
    #define NO_SCALARS 1
#endif
static const size_t noScalarComp = NO_SCALARS;               // the number of scalar components
static const size_t noScalarGradComp = noScalarComp * NO_DIM;// the number of derivatives of the scalar field (NO_DIM derivatives for each component of the scalar field)



// pi
#define PI   3.14159265


// define conversion factor
#ifndef MPC_UNIT
#define MPC_UNIT 1.
#endif



#endif
