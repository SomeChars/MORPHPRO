#include "protein.h"

#ifndef STDIO
    #define STDIO
    #include <stdio.h>
#endif
#ifndef STDLIB
    #define STDLIB
    #include <stdlib.h>
#endif
#ifndef MATH
    #define MATH
    #include <math.h>
#endif


double *new_vertex() {

    double *v = (double *)malloc(sizeof(double) * 3);

    assert(v != NULL);
  
    return v;
}

Protein new_protein(int n) {

    assert(n > 0);

    Protein ans;
    ans.n = n;       

    double **d = (double **) malloc(sizeof(double *) * n);
    assert(d != NULL);

    for(int i = 0; i < n; i ++) {
        d[i] = new_vertex();
        assert(d[i] != NULL);
    }

    ans.rs = d;
    ans.idx = (int *) malloc(sizeof(int) * n);

    return ans;
}

void free_protein(Protein p) {

    for(int i = 0; i < p.n; i ++) {
        if (p.rs[i] != NULL) free(p.rs[i]);
    }

    if (p.idx != NULL) free(p.idx);
}

void scan_protein(Protein p, int n) {
    assert(n > 0);
    assert(p.n == n);

    for(int i = 0; i < n; i ++) {
        for(int j = 0; j < 3; j ++) {
            scanf("%lf", &p.rs[i][j]);
        }
    }

    for(int i = 0; i < n; i ++) {
        scanf("%d", &p.idx[i]);
    }
}

Protein copy_protein(Protein p){ 

    assert(p.n > 0);
    assert(p.rs != NULL);
    for (int i=0; i < p.n; i++) 
        assert(p.rs[i] != NULL);
    assert(p.idx != NULL);

    int n = p.n;    
    
    Protein pNew = new_protein(n);

    for(int i = 0; i < n; i ++) {
        
	for (int j = 0; j < 3; j++) 
            pNew.rs[i][j] = p.rs[i][j];

        pNew.idx[i] = p.idx[i];
    }

    return pNew;
}

void print_protein(Protein p) {
    for(int i = 0; i < p.n; i++) {
        for(int j = 0; j < 3; j++) {
            printf("%lf ", p.rs[i][j]);
        }
    }
}

void println_protein(Protein p) {
    print_protein(p);
    printf("\n");
}

int validate_protein(Protein p) {

   if (p.n < 1) return 0;

   if (p.rs == NULL) return 0;

   for (int i = 0; i < p.n; i++) 
       if (p.rs[i] == NULL) return 0;
  
   if (p.idx == NULL) return 0;

   return 1;
}

void log_protein(Protein p) {
 
    char info[80];

    for (int i = 0; i < p.n; i++) {
        sprintf(info, "%lf %lf %lf, index = %d\n", p.rs[i][0], p.rs[i][1], p.rs[i][2], p.idx[i]);
        write_to_log(info);
    }
}

