FORT = gfortran
CC = g++
FFLAGS = -O3 -fopenmp
CFLAGS = -O3 -fopenmp
TSG_INC = -I/panfs/panasas1/users/jburkardt/include
TSG_LIB = -L/panfs/panasas1/users/jburkardt/libcpp/Linux -ltsg

all: tsg_prb

tsg_prb: tsg_prb.o tsgFCWrapper.o tsgFBridge.o
	$(FORT) $(FFLAGS) tsg_prb.o tsgFModule.mod tsgFBridge.o tsgFCWrapper.o -o tsg_prb $(TSG_LIB) -lstdc++
	
tsg_prb.o: tsg_prb.f90 tsgFModule.mod
	$(FORT) $(FFLAGS) -c tsg_prb.f90 -o tsg_prb.o
	
tsgFModule.mod: tsgFModule.f90
	$(FORT) $(FFLAGS) -c tsgFModule.f90 -o tsgFModule.mod
	
tsgFCWrapper.o: tsgFCWrapper.cpp
	$(CC) $(CFLAGS) $(TSG_INC) -c tsgFCWrapper.cpp -o tsgFCWrapper.o
	
tsgFBridge.o: tsgFBridge.f90 tsgFModule.mod 
	$(FORT) $(FFLAGS) -c tsgFBridge.f90 -o tsgFBridge.o
	
clean:
	rm -fr tsg_prb
	rm -fr *.o
	rm -fr *.mod
