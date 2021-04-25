#ifndef PROTEIN_H_
#define PROTEIN_H_

#include "common.h"

typedef struct {
    int n;
    double ** rs;
    int * idx;
} Protein;


// allocate memory for a vertex
extern double *new_vertex();

// allocate memory for a protein comprising n residues
extern Protein new_protein(int n);

// release memory used by a protein 
extern void free_protein(Protein p);

// scan a protein comprising n residues
extern void scan_protein(Protein p, int n);

// copy a protein 
extern Protein copy_protein(Protein p);

// output a protein
extern void print_protein(Protein p);
// output a protein + return symbol 
extern void println_protein(Protein p);

// validate protein 
// output: 1 if is valid, 0 otherwise
extern int validate_protein(Protein p);

extern void log_protein(Protein p); 

#endif // PROTEIN_H_


