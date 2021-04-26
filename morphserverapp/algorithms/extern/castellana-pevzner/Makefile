# Use "make all" to make "oe" and "oeac" files
# Use "make oe/oeac" to make only oe/oeac files

CC = g++

all:			oe oea oeac
					echo "done."

oe:			oe.o protein.o lattice.o common.o
					$(CC) oe.o protein.o lattice.o common.o -o oe

oe.o:			algo.c 
					$(CC) -DDISTANCE_CONSISTENT -o oe.o -c algo.c

oea:			oea.o protein.o lattice.o common.o 
					$(CC) oea.o protein.o lattice.o common.o -o oea

oea.o:			algo.c 
					$(CC) -DDISTANCE_CONSISTENT -DHEURISTIC -o oea.o -c algo.c

oeac:			oeac.o protein.o lattice.o common.o 
					$(CC) oeac.o protein.o lattice.o common.o -o oeac

oeac.o:			algo.c 
					$(CC) -DDISTANCE_CONSISTENT -DANGLE_CONSISTENT -o oeac.o -c algo.c

protein.o:		protein.c
					$(CC) -c protein.c

common.o:		common.c
					$(CC) -c common.c

lattice.o:		lattice.c
					$(CC) -c lattice.c

clean:
					rm *.o oe oeac oea
