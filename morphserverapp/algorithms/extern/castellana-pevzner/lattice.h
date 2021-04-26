#ifndef LATTICE_H_
#define LATTICE_H_

#include <stdio.h>
#include <stdlib.h>

#include "common.h"

#define LATTICE_DENSITY 6

#define MIN_LATTICE_EDGE_LEN 1.
#define MAX_LATTICE_EDGE_LEN 2.
#define LATTICE_EDGE_INCREM  0.5


typedef struct {

   double edgeLen;

   int Q;
   
   double **p;      // stores coordinates of lattice points w.r.t. the lattice center    

   double *vscores; // stores VScores of lattice points

} Lattice;


extern void _set_Q(Lattice *pLat);
extern void _set_p(Lattice *pLat);
extern void _set_vscores(Lattice *pLat);

extern Lattice init_lattice(double edgeLen);
extern void free_lattice(Lattice lat);

extern int validate_lattice(Lattice lat);

extern Lattice init_lattice_around(Lattice lat);
extern void set_lattice_around(Lattice latAround, Lattice lat, double *v);

extern void get_vertex(double* vi, Lattice lat, double* center, int i); // stores in 'vi' the i-th vertex of the lattice 'lat' built around 'center' 

extern void log_lattice(Lattice lat);

#endif
