SYSTEM = x86-64_linux
LIBFORMAT = static_pic

CPLEXDIR = /home/mdelorme/data/volume_2/ibm/ILOG/CPLEX_Studio2211/cplex
CPDIR = /home/mdelorme/data/volume_2/ibm/ILOG/CPLEX_Studio2211/cpoptimizer
CONCERTDIR = /home/mdelorme/data/volume_2/ibm/ILOG/CPLEX_Studio2211/concert
CPLEXBINDIR = $(CPLEXDIR)/bin/$(SYSTEM)
CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLIBDIR = $(CPDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEXFLAGS = -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include -I$(CPDIR)/include  -L$(CPLIBDIR) -lcp -L$(CPLEXLIBDIR) -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread -ldl 

GUROBI_DIR = /data/volume_2/gurobi1001/linux64
CLIB     = -L$(GUROBI_DIR)/lib/ -lgurobi100
INC      = $(GUROBI_DIR)/include/
CPPLIB   = -L$(GUROBI_DIR)/lib/ -lgurobi_c++ $(CLIB)
GUROBI_OPTS = -I$(INC) $(CPPLIB) -lpthread -lm -m64 

CC = g++ 
CFLAGS = -O2 -Wall -ansi -pedantic -DIL_STD -std=c++11
DEBUG = -pg -g -Wall -ansi -pedantic -DIL_STD 
OBJECTS = main.o Allocation.o time.o 
 
exec : $(OBJECTS) 
	$(CC) $(CFLAGS) -o DEC3 $(OBJECTS) $(CPLEXFLAGS) $(GUROBI_OPTS)
 
.cpp.o : 
	$(CC) $(CFLAGS) $(CPLEXFLAGS) $(GUROBI_OPTS) -c $< -o $@ 
 
clean : 
	rm -f $(OBJECTS) DEC3 