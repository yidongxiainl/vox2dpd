CC = g++

####flags
CFLAGS = -O3 -Wall -Wno-unused-result
#CFLAGS = -Og -g -pg -fPIC
CFLAGS += -std=c++0x

####Source files
EXE = vox2dpd
OBJ = vox2dpd.o

####Executable
VOXCEL2ATOM:$(EXE)
$(EXE):$(OBJ)
	$(CC) $(CFLAGS) -o $(EXE) $(OBJ)
.cpp.o:
	$(CC) $(CFLAGS) -c $*.cpp

#### clean up old builds
clean:
	rm -f *.o
	rm -f $(EXE)

