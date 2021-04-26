#ifndef RUN_H_ 
#define RUN_H_

#include <math.h>
#include <stdlib.h>

#include "common.h"
#include "protein.h"
#include "lattice.h"

const int RES_OK   =  0;
const int RES_FAIL = -1;

// initialization
extern void init();

// scoring
extern inline double dist_2(double* v1, double* v2);
extern inline double escore(Protein proInter, int j, double* v1, double* v2);
extern inline double ascore(Protein proInter, int j, double* v1, double* v2, double* v3);
extern void adjust_vscores(Lattice *pLatAround, Protein proInter, int j); 

// morphing
extern double conv(double a, double b, double alpha);
extern Protein intermediate(Protein proPrev, Protein proFinal, double alpha);

#ifdef ANGLE_CONSISTENT
extern void solveOESP(Lattice lat, Protein proInter, short*** back, int *indFinalI, int *indFinalH);
extern void adjustIntermediate(Lattice lat, int indFinalI, int indFinalH, Protein proInter, short*** back);
extern int proteinize(Protein proInter, double latticeEdgeLen); 
#else
extern void solveOESP(Lattice lat, Protein proInter, short**back, int *indFinal);
extern void adjustIntermediate(Lattice lat, int indFinal, Protein proInter, short** back);
extern int proteinize(Protein proInter, double latticeEdgeLen);
#endif

extern int morph(Protein *p, int k);

extern void free_protein_sequence(Protein* p, int k);

//extern int main();

#endif  // RUN_H_
