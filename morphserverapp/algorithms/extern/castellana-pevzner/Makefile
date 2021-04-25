# Use "make all" to make "oe" and "oeac" files
# Use "make oe/oeac" to make only oe/oeac files

OPT = -ffast-math -O3 -funroll-loops -finline-functions -Winline
WARN = -Wall -pedantic -std=c99
CFLAGS = $(WARN)
CC = gcc

all:			oe oea oeac
					echo "done."

oe:			oe.o protein.o lattice.o common.o
					$(CC) $(OPT) $(CFLAGS) oe.o protein.o lattice.o common.o $(LIBS) -o oe

oe.o:			algo.c 
					$(CC) $(OPT) $(CFLAGS) $(INCDIR) -DDISTANCE_CONSISTENT -o oe.o -c algo.c

oea:			oea.o protein.o lattice.o common.o 
					$(CC) $(OPT) $(CFLAGS) oea.o protein.o lattice.o common.o $(LIBS) -o oea

oea.o:			algo.c 
					$(CC) $(OPT) $(CFLAGS) $(INCDIR) -DDISTANCE_CONSISTENT -DHEURISTIC -o oea.o -c algo.c

oeac:			oeac.o protein.o lattice.o common.o 
					$(CC) $(OPT) $(CFLAGS) oeac.o protein.o lattice.o common.o $(LIBS) -o oeac

oeac.o:		algo.c 
					$(CC) $(OPT) $(CFLAGS) $(INCDIR) -DDISTANCE_CONSISTENT -DANGLE_CONSISTENT -o oeac.o -c algo.c

protein.o:		protein.c
					$(CC) $(OPT) $(CFLAGS) $(INCDIR) -c protein.c

common.o:			common.c
					$(CC) $(OPT) $(CFLAGS) $(INCDIR) -c common.c

lattice.o:		lattice.c
					$(CC) $(OPT) $(CFLAGS) $(INCDIR) -c lattice.c

clean:
					rm *.o oe oeac
