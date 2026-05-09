# ----------------------------------------------------------------
# environment
CC = mpicxx

# ----------------------------------------------------------------
# options
#CFLAGS = #-std=c++11
CFLAGS = -I ./include -I /opt/homebrew/Cellar/hypre/2.31.0/include -I /opt/homebrew/Cellar/open-mpi/5.0.3_1/include#-std=c++11
CLIBS = -L /opt/homebrew/Cellar/hypre/2.31.0/lib -lHYPRE
OPTIMIZE = -Ofast#-fopenmp
EXE_NAME = a.out

# ----------------------------------------------------------------
# executables
Normal: ./src/main.cpp ./src/common.cpp ./src/support.cpp ./src/initial.cpp ./src/input.cpp ./src/EMfield.cpp ./src/fluid.cpp ./src/solver.cpp ./src/output.cpp
	$(CC) $(CFLAGS) $(CLIBS) $(OPTIMIZE) ./src/main.cpp ./src/common.cpp ./src/support.cpp ./src/initial.cpp ./src/input.cpp ./src/EMfield.cpp ./src/fluid.cpp ./src/solver.cpp ./src/output.cpp -o $(EXE_NAME)
