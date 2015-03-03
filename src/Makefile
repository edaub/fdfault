CC=mpic++
CFLAGS=-c -O2
EXEC = ../fdfault

block.o : block.hpp boundary.hpp cartesian.hpp coord.hpp fd.hpp material.hpp surface.hpp block.cpp
	$(CC) $(CFLAGS) block.cpp

boundary.o : boundary.hpp coord.hpp fd.hpp fields.hpp material.hpp surface.hpp boundary.cpp
	$(CC) $(CFLAGS) boundary.cpp

cartesian.o : cartesian.hpp coord.hpp cartesian.cpp
	$(CC) $(CFLAGS) cartesian.cpp

coord.o : coord.hpp coord.cpp
	$(CC) $(CFLAGS) coord.cpp

domain.o : block.hpp cartesian.hpp domain.hpp fd.hpp fields.hpp rk.hpp domain.cpp
	$(CC) $(CFLAGS) domain.cpp

fd.o : fd.hpp coord.hpp fd.cpp
	$(CC) $(CFLAGS) fd.cpp

fields.o : fields.hpp coord.hpp cartesian.hpp fields.cpp
	$(CC) $(CFLAGS) fields.cpp

friction.o: friction.hpp interface.hpp friction.cpp
	$(CC) $(CFLAGS) friction.cpp

main.o : problem.hpp main.cpp
	$(CC) $(CFLAGS) main.cpp

material.o : material.hpp material.cpp
	$(CC) $(CFLAGS) material.cpp

outputlist.o : outputlist.hpp outputunit.hpp outputlist.cpp
	$(CC) $(CFLAGS) outputlist.cpp

outputunit.o : outputunit.hpp outputunit.cpp
	$(CC) $(CFLAGS) outputunit.cpp

problem.o : domain.hpp outputlist.hpp problem.hpp rk.hpp problem.cpp
	$(CC) $(CFLAGS) problem.cpp

rk.o : rk.hpp rk.cpp
	$(CC) $(CFLAGS) rk.cpp

surface.o : coord.hpp surface.hpp surface.cpp
	$(CC) $(CFLAGS) surface.cpp

fdfault : block.o boundary.o cartesian.o coord.o domain.o fd.o fields.o main.o material.o outputlist.o outputunit.o problem.o rk.o surface.o
	$(CC) -o $(EXEC) block.o boundary.o cartesian.o coord.o domain.o \
		fd.o fields.o main.o material.o outputlist.o outputunit.o \
		problem.o rk.o surface.o