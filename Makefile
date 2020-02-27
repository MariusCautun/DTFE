# Makefile for compiling the DTFE code on Linux systems


# Path to the GSL, Boost C++ and CGAL libraries - must be set by user (only if they aren't installed in the default system path) -- (NOTE: You must add only the directory where the libraries are installed, the program will add the '/lib' and '/include' parts automatically); C++ compiler - preferably a version that supports OpenMP
GSL_PATH = /cosma/local/gsl/2.4
BOOST_PATH = /cosma/local/boost/gnu_7.3.0/1_67_0
CGAL_PATH = /cosma/home/dphlss/cautun/Programs/stow
CC = g++
# set the following if you have installed the HDF5 library and would like to read in HDF5 gadget files
HDF5_PATH = 


# paths to where to put the object files and the executables files. If you build the DTFE library than you also need to specify the directory where to put the library and the directory where to copy the header files needed by the library (choose an empty directory for the header files).
OBJ_DIR = ./o
BIN_DIR = ./
LIB_DIR = ./
INC_DIR = ./DTFE_include



############################# Choose the compiler directives ##################################

############################# Overall options ##################################
OPTIONS = 
#------------------------ set the number of spatial dimensions (2 or 3 dimensions)
OPTIONS += -DNO_DIM=3 
#------------------------ set type of variables - float (comment the next line) or double (uncomment the next line)
# OPTIONS += -DDOUBLE 

############################# Quantities to be computed ##################################
#------------------------ set which quantities can be computed (can save memory by leaving some out)
# Comment this line if you don't need to compute velocity and velocity related components 
OPTIONS += -DVELOCITY 
# Comment this line if you don't need to interpolate additional fields stored in the scalar variable
OPTIONS += -DSCALAR 
# number of components of the scalar variable
OPTIONS += -DNO_SCALARS=1 

############################# Input and output operations default settings ##################################
#------------------------ set which are the default input and output functions for doing data io
# default function to read the input data (101-multiple gadget file, 102-single gadget file, 105-HDF5 gadget file, 111-text file, ... see documentation for more options). The input file type can be set during runtime using the option '--input'. This makefile option only sets a default input file in the case none is given via the program options.
OPTIONS += -DINPUT_FILE_DEFAULT=101 
# default value for the units of the input data (value=what is 1 Mpc in the units of the data - in this example the data is in kpc). You can change this also during runtime using the program option '--MpcUnit'.
OPTIONS += -DMPC_UNIT=1000. 
# default function to write the output data (101-binary file, 111-text file, ... see documentation for more options). The output file type can be set during runtime using the option '--output'. This makefile option only sets a default output file in the case none is given via the program options.
OPTIONS += -DOUTPUT_FILE_DEFAULT=101  
#101 for binary file, 100 my density file

############################# additional compiler options ##################################
# enable this option if to use OpenMP (share the workload between CPU cores sharing the same RAM)
OPTIONS += -DOPEN_MP 
# enable to check if the padding gives a complete Delaunay Tesselation of the region of interest
# OPTIONS += -DTEST_PADDING 
# enable this option to shift from position space to redshift space; You also need to activate this option during run-time using '--redshift-space arguments'
OPTIONS += -DREDSHIFT_SPACE

#------------------------ options usefull when using DTFE as a library
# uncomment the line to get access to a function that returns the Delaunay triangulation of the point set
# OPTIONS += -DTRIANGULATION 


############################# Help menu messages options ##################################
#------------------------ compiler directive that affect only the help messages when using the '-h / --help' option (it does not affect the program in any other way)- if the option is uncommented, than it will show that set of options in the help menu
OPTIONS += -DFIELD_OPTIONS 
OPTIONS += -DREGION_OPTIONS 
# OPTIONS += -DPARTITION_OPTIONS 
# OPTIONS += -DPADDING_OPTIONS 
OPTIONS += -DAVERAGING_OPTIONS 
# OPTIONS += -DREDSHIFT_CONE_OPTIONS 
OPTIONS += -DADDITIONAL_OPTIONS 








###############  DO NOT MODIFY BELOW THIS LINE  ###########################
# do not modify below this line
SRC = ./src
INCLUDES = 
LIBRARIES = 

ifneq ($(strip $(GSL_PATH)),)
	INCLUDES += -I/$(strip $(GSL_PATH))/include 
	LIBRARIES += -L/$(strip $(GSL_PATH))/lib 
endif
ifneq ($(strip $(BOOST_PATH)),)
	INCLUDES += -I/$(strip $(BOOST_PATH))/include 
	LIBRARIES += -L/$(strip $(BOOST_PATH))/lib 
endif
ifneq ($(strip $(CGAL_PATH)),)
	INCLUDES += -I/$(strip $(CGAL_PATH))/include 
	LIBRARIES += -L/$(strip $(CGAL_PATH))/lib 
endif
ifneq ($(strip $(HDF5_PATH)),)
	INCLUDES += -I/$(strip $(HDF5_PATH))/include 
	LIBRARIES += -L/$(strip $(HDF5_PATH))/lib -lhdf5 -lhdf5_cpp
	OPTIONS += -DHDF5
endif



COMPILE_FLAGS = -frounding-math -O3 -fopenmp -DNDEBUG $(OPTIONS)
DTFE_INC = $(INCLUDES)
# the following libraries should work in most cases
DTFE_LIB = -rdynamic $(LIBRARIES) -lCGAL -lboost_thread -lboost_filesystem -lboost_program_options -lgmp -lgsl -lgslcblas -lm -lboost_system



IO_SOURCES = $(addprefix io/, input_output.h gadget_reader.cc text_io.cc binary_io.cc my_io.cc)
MAIN_SOURCES = main.cpp DTFE.h message.h user_options.h input_output.cc $(IO_SOURCES)
DTFE_SOURCES = DTFE.cpp define.h particle_data.h user_options.h box.h quantities.h user_options.cc quantities.cc subpartition.h random.cc CIC_interpolation.cc TSC_interpolation.cc SPH_interpolation.cc kdtree/kdtree2.hpp Pvector.h message.h miscellaneous.h
TRIANG_SOURCES = $(addprefix CGAL_triangulation/, triangulation.cpp triangulation_miscellaneous.cc unaveraged_interpolation.cc averaged_interpolation_1.cc averaged_interpolation_2.cc padding_test.cc CGAL_include_2D.h CGAL_include_3D.h vertexData.h particle_data_traits.h) define.h particle_data.h user_options.h box.h quantities.h Pvector.h message.h math_functions.h

ALL_FILES = $(DTFE_SOURCES) $(TRIANG_SOURCES) $(MAIN_SOURCES) kdtree/kdtree2.hpp kdtree/kdtree2.cpp
LIB_FILES = $(DTFE_SOURCES) $(TRIANG_SOURCES)

HEADERS_1 = DTFE.h define.h user_options.h particle_data.h quantities.h Pvector.h math_functions.h  message.h box.h miscellaneous.h
HEADERS_2 = $(addprefix CGAL_triangulation/, CGAL_include_2D.h CGAL_include_3D.h vertexData.h particle_data_traits.h)



DTFE: set_directories $(OBJ_DIR)/DTFE.o $(OBJ_DIR)/triangulation.o $(OBJ_DIR)/main.o $(OBJ_DIR)/kdtree2.o Makefile
	$(CC) $(COMPILE_FLAGS) $(OBJ_DIR)/DTFE.o $(OBJ_DIR)/triangulation.o $(OBJ_DIR)/main.o $(OBJ_DIR)/kdtree2.o $(DTFE_LIB) -o $(BIN_DIR)/DTFE


$(OBJ_DIR)/main.o: $(addprefix $(SRC)/, $(MAIN_SOURCES)) Makefile
	$(CC) $(COMPILE_FLAGS) $(DTFE_INC) -o $(OBJ_DIR)/main.o -c $(SRC)/main.cpp

$(OBJ_DIR)/DTFE.o: $(addprefix $(SRC)/, $(DTFE_SOURCES)) Makefile
	$(CC) $(COMPILE_FLAGS) $(DTFE_INC) -o $(OBJ_DIR)/DTFE.o -c $(SRC)/DTFE.cpp

$(OBJ_DIR)/kdtree2.o: $(SRC)/kdtree/kdtree2.hpp $(SRC)/kdtree/kdtree2.cpp Makefile
	$(CC) -O3 -ffast-math -fomit-frame-pointer $(DTFE_INC) -o $(OBJ_DIR)/kdtree2.o -c $(SRC)/kdtree/kdtree2.cpp

$(OBJ_DIR)/triangulation.o: $(addprefix $(SRC)/, $(TRIANG_SOURCES)) Makefile
	$(CC) $(COMPILE_FLAGS) $(DTFE_INC) -o $(OBJ_DIR)/triangulation.o -c $(SRC)/CGAL_triangulation/triangulation.cpp


library: set_directories set_directories_2 $(addprefix $(SRC)/, $(LIB_FILES) ) copy_headers Makefile
	$(CC) $(COMPILE_FLAGS) -fPIC $(DTFE_INC) -o $(OBJ_DIR)/DTFE_l.o -c $(SRC)/DTFE.cpp
	$(CC) -O3 -ffast-math -fomit-frame-pointer -fPIC $(DTFE_INC) -o $(OBJ_DIR)/kdtree2_l.o -c $(SRC)/kdtree/kdtree2.cpp
	$(CC) $(COMPILE_FLAGS) -fPIC $(DTFE_INC) -o $(OBJ_DIR)/triangulation_l.o -c $(SRC)/CGAL_triangulation/triangulation.cpp
	$(CC) $(COMPILE_FLAGS) -shared $(OBJ_DIR)/DTFE_l.o $(OBJ_DIR)/triangulation_l.o $(OBJ_DIR)/kdtree2_l.o $(DTFE_LIB) -o $(LIB_DIR)/libDTFE.so


clean:
	rm -f $(BIN_DIR)/DTFE $(OBJ_DIR)/*.o

copy_headers:
	cp $(addprefix $(SRC)/, $(HEADERS_1)) $(INC_DIR)
	cp $(addprefix $(SRC)/, $(HEADERS_2)) $(INC_DIR)/CGAL_triangulation

set_directories:
	@ if !( test -d $(OBJ_DIR) ); \
	then mkdir $(OBJ_DIR); \
	fi
	@ if !( test -d $(BIN_DIR) ); \
	then mkdir $(BIN_DIR); \
	fi

set_directories_2:
	@ if !( test -d $(LIB_DIR) ); \
	then mkdir $(LIB_DIR); \
	fi
	@ if !( test -d $(INC_DIR) ); \
	then mkdir $(INC_DIR); \
	fi
	@ if !( test -d $(INC_DIR)/CGAL_triangulation ); \
	then mkdir $(INC_DIR)/CGAL_triangulation; \
	fi
