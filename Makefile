# ----------------------------------------------------------------
# environment
CC = mpicxx

# ----------------------------------------------------------------
# options
#CFLAGS = #-std=c++11
CFLAGS = -I /opt/homebrew/Cellar/hypre/2.31.0/include -I /opt/homebrew/Cellar/open-mpi/5.0.3_1/include#-std=c++11
CLIBS = -L /opt/homebrew/Cellar/hypre/2.31.0/lib -lHYPRE
OPTIMIZE = -Ofast#-fopenmp
EXE_NAME = a.out

# ----------------------------------------------------------------
# executables
Normal: common.hpp functions.hpp main.cpp common.cpp support.cpp initial.cpp input.cpp EMfield.cpp fluid.cpp solver.cpp output.cpp
	$(CC) $(CFLAGS) $(CLIBS) $(OPTIMIZE) main.cpp common.cpp support.cpp initial.cpp input.cpp EMfield.cpp fluid.cpp solver.cpp output.cpp -o $(EXE_NAME)
