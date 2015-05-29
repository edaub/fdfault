CC=mpic++
CFLAGS=-c -O3
IFLAGS=
LFLAGS=
EFLAGS=-O3
EXEC=../fdfault

fdfault : block.o boundary.o cartesian.o coord.o domain.o fd.o fields.o friction.o front.o frontlist.o interface.o load.o main.o material.o outputlist.o outputunit.o pert.o problem.o rk.o slipweak.o swparam.o surface.o utilities.o
	$(CC) $(EFLAGS) $(LFLAGS) $(IFLAGS) -o $(EXEC) block.o boundary.o cartesian.o coord.o domain.o \
		fd.o fields.o friction.o front.o frontlist.o interface.o load.o main.o material.o \
		outputlist.o outputunit.o pert.o problem.o rk.o slipweak.o swparam.o surface.o utilities.o

block.o : block.hpp boundary.hpp cartesian.hpp coord.hpp fd.hpp material.hpp surface.hpp block.cpp
	$(CC) $(CFLAGS) block.cpp

boundary.o : boundary.hpp coord.hpp fd.hpp fields.hpp material.hpp surface.hpp boundary.cpp
	$(CC) $(CFLAGS) boundary.cpp

cartesian.o : cartesian.hpp coord.hpp cartesian.cpp
	$(CC) $(CFLAGS) cartesian.cpp

coord.o : coord.hpp coord.cpp
	$(CC) $(CFLAGS) coord.cpp

domain.o : block.hpp cartesian.hpp domain.hpp fd.hpp fields.hpp interface.hpp rk.hpp slipweak.hpp domain.cpp
	$(CC) $(CFLAGS) domain.cpp

fd.o : fd.hpp coord.hpp fd.cpp
	$(CC) $(CFLAGS) fd.cpp

fields.o : fields.hpp coord.hpp cartesian.hpp fields.cpp
	$(CC) $(CFLAGS) fields.cpp

friction.o : block.hpp cartesian.hpp fd.hpp fields.hpp friction.hpp interface.hpp load.hpp friction.cpp
	$(CC) $(CFLAGS) friction.cpp

front.o : cartesian.hpp domain.hpp fields.hpp front.hpp interface.hpp utilities.h front.cpp
	$(CC) $(CFLAGS) front.cpp

frontlist.o : domain.hpp front.hpp frontlist.hpp frontlist.cpp
	$(CC) $(CFLAGS) frontlist.cpp

interface.o : block.hpp boundary.hpp cartesian.hpp coord.hpp fields.hpp interface.hpp surface.hpp interface.cpp
	$(CC) $(CFLAGS) interface.cpp

load.o : load.hpp pert.hpp load.cpp
	$(CC) $(CFLAGS) load.cpp

main.o : problem.hpp main.cpp
	$(CC) $(CFLAGS) main.cpp

material.o : material.hpp material.cpp
	$(CC) $(CFLAGS) material.cpp

outputlist.o : outputlist.hpp outputunit.hpp outputlist.cpp
	$(CC) $(CFLAGS) outputlist.cpp

outputunit.o : outputunit.hpp utilities.h outputunit.cpp
	$(CC) $(CFLAGS) outputunit.cpp

pert.o : pert.hpp pert.cpp
	$(CC) $(CFLAGS) pert.cpp

problem.o : domain.hpp outputlist.hpp frontlist.hpp problem.hpp rk.hpp problem.cpp
	$(CC) $(CFLAGS) problem.cpp

rk.o : rk.hpp rk.cpp
	$(CC) $(CFLAGS) rk.cpp

slipweak.o : block.hpp cartesian.hpp fd.hpp fields.hpp friction.hpp interface.hpp swparam.hpp slipweak.hpp slipweak.cpp
	$(CC) $(CFLAGS) slipweak.cpp

swparam.o : pert.hpp swparam.hpp swparam.cpp
	$(CC) $(CFLAGS) swparam.cpp

surface.o : coord.hpp surface.hpp surface.cpp
	$(CC) $(CFLAGS) surface.cpp

utilities.o : utilities.h utilities.cpp
	$(CC) $(CFLAGS) utilities.cpp

clean:
	rm *.o $(EXEC)