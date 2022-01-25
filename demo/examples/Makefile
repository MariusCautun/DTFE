# Makefile for compiling the DTFE code on Linux systems


# Path to the GSL, Boost C++ and CGAL libraries - must be set by user (only if they aren't installed in the default system path); C++ compiler - preferably a version that supports OpenMP
GSL_PATH = 
BOOST_PATH = 
CGAL_PATH = 
CC = g++44


# paths to where to put the object files and the executables files. If you build the DTFE library than you also need to specify the directory where to put the library and the directory where to copy the header files needed by the library (choose an empty directory for the header files).
OBJ_DIR = ./o
BIN_DIR = ./
LIB_DIR = ./
INC_DIR = ./DTFE_include


# Choose the compiler directives
OPTIONS = 
#------------------------ set the number of spatial dimensions (2 or 3 dimensions)
OPTIONS += -DNO_DIM=3 
#------------------------ set type of variables - float (comment the next line) or double (uncomment the next line)
# OPTIONS += -DDOUBLE 
#------------------------ set which quantities can be computed (can save memory by leaving some out)
# do compute velocity and velocity related components 
OPTIONS += -DVELOCITY 
# do compute scalar fields
OPTIONS += -DSCALAR 
# number of scalar fields
OPTIONS += -DNO_SCALARS=1 
#------------------------ additional compiler options
# enable to check if the padding gives a complete Delaunay Tesselation of the region of interest
# OPTIONS += -DTEST_PADDING 
# enable this option if to use OpenMP (share the workload between CPU cores sharing the same RAM)
OPTIONS += -DOPEN_MP 
# enable this if the input data is not in Mpc units (value=what is 1 MPc in the units of the data - in this example the data is in kpc)
OPTIONS += -DMPC_UNIT=1000. 
#------------------------ options usefull when using DTFE as a library
# uncomment the line to get access to a function that returns the Delaunay triangulation of the point set
# OPTIONS += -DTRIANGULATION 


#------------------------ compiler directive that affect only the help messages when using the '-h / --help' option (it does not affect the program in any other way)- if the option is uncommented, than it will show that set of options
OPTIONS += -DFIELD_OPTIONS 
OPTIONS += -DREGION_OPTIONS 
# OPTIONS += -DPARTITION_OPTIONS 
# OPTIONS += -DPADDING_OPTIONS 
OPTIONS += -DDENSITY_OPTIONS 
# OPTIONS += -DREDSHIFT_CONE_OPTIONS 
OPTIONS += -DADDITIONAL_OPTIONS 







###############  DO NOT MODIFY BELOW THIS LINE  ###########################
# do not modify below this line
SRC = ../../src

INCLUDES = 
LIBRARIES = 

ifneq ($(strip $(GSL_PATH)),)
	INCLUDES += -I/$(GSL_PATH)/include 
	LIBRARIES += -L/$(GSL_PATH)/lib 
endif
ifneq ($(strip $(BOOST_PATH)),)
	INCLUDES += -I/$(BOOST_PATH)/include 
	LIBRARIES += -L/$(BOOST_PATH)/lib 
endif
ifneq ($(strip $(CGAL_PATH)),)
	INCLUDES += -I/$(CGAL_PATH)/include 
	LIBRARIES += -L/$(CGAL_PATH)/lib 
endif

COMPILE_FLAGS = -frounding-math -O3 -fopenmp -DNDEBUG $(OPTIONS)
DTFE_INC = $(INCLUDES)
DTFE_LIB = -rdynamic $(LIBRARIES) -lCGAL -lboost_thread -lboost_filesystem -lboost_program_options -lgsl -lgslcblas -lm




MAIN_SOURCES = DTFE.h input_output.cc input_output.h message.h math_functions.h
DTFE_SOURCES = DTFE.cpp define.h particle_data.h user_options.h box.h quantities.h user_options.cc quantities.cc subpartition.h random.cc TSC_interpolation.cc SPH_interpolation.cc kdtree/kdtree2.hpp Pvector.h message.h
TRIANG_SOURCES = triangulation_interpolation.cpp define.h CGAL_include_2D.h CGAL_include_3D.h vertexData.h particle_data.h particle_data_traits.h user_options.h box.h quantities.h density_interpolation.cc velocity_interpolation.cc scalar_interpolation.cc Pvector.h message.h

FILES = DTFE.h DTFE.cpp triangulation_interpolation.cpp density_interpolation.cc velocity_interpolation.cc CGAL_include_2D.h  CGAL_include_3D.h scalar_interpolation.cc define.h CGAL_include_3D.h vertexData.h particle_data.h particle_data_traits.h user_options.h user_options.cc box.h quantities.h subpartition.h quantities.cc Pvector.h message.h TSC_interpolation.cc SPH_interpolation.cc random.cc kdtree/kdtree2.hpp kdtree/kdtree2.cpp


DTFE_particle_density: set_directories $(OBJ_DIR)/DTFE.o $(OBJ_DIR)/triangulation_interpolation.o $(OBJ_DIR)/main_particle_density.o $(OBJ_DIR)/kdtree2.o Makefile
	$(CC) $(COMPILE_FLAGS) $(OBJ_DIR)/DTFE.o $(OBJ_DIR)/triangulation_interpolation.o $(OBJ_DIR)/main_particle_density.o $(OBJ_DIR)/kdtree2.o $(DTFE_LIB) -o $(BIN_DIR)/DTFE_particle_density


$(OBJ_DIR)/main_particle_density.o: $(addprefix $(SRC)/, $(MAIN_SOURCES)) main_particle_density.cpp Makefile
	$(CC) $(COMPILE_FLAGS) $(DTFE_INC) -o $(OBJ_DIR)/main_particle_density.o -c main_particle_density.cpp

$(OBJ_DIR)/DTFE.o: $(addprefix $(SRC)/, $(DTFE_SOURCES)) Makefile
	$(CC) $(COMPILE_FLAGS) $(DTFE_INC) -o $(OBJ_DIR)/DTFE.o -c $(SRC)/DTFE.cpp

$(OBJ_DIR)/kdtree2.o: $(SRC)/kdtree/kdtree2.hpp $(SRC)/kdtree/kdtree2.cpp Makefile
	$(CC) -O3 -ffast-math -fomit-frame-pointer $(BOOST_INC) -o $(OBJ_DIR)/kdtree2.o -c $(SRC)/kdtree/kdtree2.cpp

$(OBJ_DIR)/triangulation_interpolation.o: $(addprefix $(SRC)/, $(TRIANG_SOURCES)) Makefile
	$(CC) $(COMPILE_FLAGS) $(DTFE_INC) -o $(OBJ_DIR)/triangulation_interpolation.o -c $(SRC)/triangulation_interpolation.cpp


clean:
	rm -f $(BIN_DIR)/DTFE_particle_density $(OBJ_DIR)/*.o


set_directories:
	@ if !( test -d $(OBJ_DIR) ); \
	then mkdir $(OBJ_DIR); \
	fi
	@ if !( test -d $(BIN_DIR) ); \
	then mkdir $(BIN_DIR); \
	fi
