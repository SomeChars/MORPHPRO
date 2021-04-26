#include "algo.h"

const double MIN_NEIGHBOR_DIST = 3.7;
const double MAX_NEIGHBOR_DIST = 3.9;
const double MIN_NEIGHBOR_DIST_2 = 13.69;
const double MAX_NEIGHBOR_DIST_2 = 15.21;

double COS_70;
double COS_120; 
double COS_70_2; 
double COS_120_2; 


char status[80];


void init()
{
    COS_70    = cos(70. * M_PI / 180.);
    COS_120   = -0.5; 
    COS_70_2  = COS_70 * COS_70; // ~0.11
    COS_120_2 = 0.25; 
}



//---------------------------------------------------------------------------------------
//
//  Scores
//
//---------------------------------------------------------------------------------------

inline double dist_2(double* v1, double* v2) {

   assert(v1 != NULL);
   assert(v2 != NULL);

   return (v2[0]-v1[0])*(v2[0]-v1[0]) + (v2[1]-v1[1])*(v2[1]-v1[1]) + (v2[2]-v1[2])*(v2[2]-v1[2]);
}

inline double escore(Protein proInter, int j, double* v1, double* v2) {

   assert((1 <= j) && (j < proInter.n));
   assert(proInter.idx != NULL);
   assert(v1 != NULL);
   assert(v2 != NULL);

   if (proInter.idx[j-1] != proInter.idx[j]-1) 
       //the current and the previous residues were non-consecutive in the initial protein
       return 0.;

   double d_2 = dist_2(v1, v2);

   if ((MIN_NEIGHBOR_DIST_2 < d_2) && (d_2 < MAX_NEIGHBOR_DIST_2))
       return 0.;

   return INF;
}

inline double ascore(Protein proInter, int j, double* v1, double* v2, double* v3) {

   assert((2 <= j) && (j < proInter.n));
   assert(proInter.idx != NULL);
   assert(v1 != NULL);
   assert(v2 != NULL);
   assert(v3 != NULL);

   if ((proInter.idx[j-1] != proInter.idx[j]-1) || (proInter.idx[j-2] != proInter.idx[j-1]-1)) 
       //the current and the previous residues, or the previous and its preceeding residue were non-consecutive in the initial protein
       return 0.;

   double aa = dist_2(v1, v2);
   double bb = dist_2(v3, v2);
   double cc = dist_2(v1, v3);

   double dd = 2*sqrt(aa*bb);

   double cos_v1v2v3 = (aa + bb - cc)/dd;

   if ((COS_70 > cos_v1v2v3) && (cos_v1v2v3 > COS_120))
      return 0.;

   return INF;
}

void adjust_vscores(Lattice *pLatAround, Protein proInter, int j) {
// *pLatAround is a lattice around the j-th residue of proInter
// vscores for its points are set equal to squared distances from the center

#ifdef DISTANCE_CONSISTENT

    assert(validate_lattice(*pLatAround));
    assert(proInter.rs != NULL);
    assert((0 <= j) && (j < proInter.n));

    for (int i = 0; i < pLatAround->Q; i++) {

       assert(pLatAround->vscores[i] > 0.);

       if (pLatAround->vscores[i] < 4.0)  
          for (int m = 0; m < j-2; m++) 
              pLatAround->vscores[i] += 100./dist_2(pLatAround->p[i], proInter.rs[m]);
    }

#endif
}

//---------------------------------------------------------------------------------------
//
//  Morphing proteins
//
//---------------------------------------------------------------------------------------

double conv(double a, double b, double alpha) {

   assert((0 < alpha) && (alpha < 1)); 

   return ((1-alpha)*a + alpha*b);
}

Protein intermediate(Protein proPrev, Protein proFinal, double alpha) {

   assert(proPrev.n == proFinal.n);

   assert(validate_protein(proPrev));
   assert(validate_protein(proFinal));
   
   assert((0 < alpha) && (alpha < 1)); 

   int n = proPrev.n;

   Protein proInter = new_protein(n);

   for (int i = 0; i < n; i++) {

       for (int j = 0; j < 3; j++)
            proInter.rs[i][j] = conv(proPrev.rs[i][j], proFinal.rs[i][j], alpha);

       proInter.idx[i] = proPrev.idx[i];
   }

  return proInter;
}


#ifdef ANGLE_CONSISTENT

void solveOESP(Lattice lat, Protein proInter, short*** back, int *indFinalI, int *indFinalH) {

      assert(lat.Q > 0);
      assert(validate_protein(proInter));
      assert(back != NULL);
      for (int j = 0; j < proInter.n; j++) {
         assert(back[j] != NULL);
         for (int i = 0; i < lat.Q; i++)
             assert(back[j][i] != NULL);
      }
      assert(indFinalI != NULL);
      assert(indFinalH != NULL);

      
      Lattice* pLatPrevPrev = (Lattice*)malloc(sizeof(Lattice));
      assert(pLatPrevPrev != NULL);
      (*pLatPrevPrev) = init_lattice_around(lat);
      set_lattice_around(*pLatPrevPrev, lat,  proInter.rs[0]);
      adjust_vscores(pLatPrevPrev, proInter, 0); 

      Lattice* pLatPrev = (Lattice*)malloc(sizeof(Lattice));
      assert(pLatPrev != NULL);
      (*pLatPrev) = init_lattice_around(lat);
      set_lattice_around(*pLatPrev, lat, proInter.rs[1]);
      adjust_vscores(pLatPrev, proInter, 1); 

      Lattice* pLatCurr = (Lattice*)malloc(sizeof(Lattice));
      assert(pLatCurr != NULL);
      (*pLatCurr) = init_lattice_around(lat);


      int Q = lat.Q; //lattice size

      double** pathPrev = (double**)malloc(sizeof(double*) * Q);
      assert(pathPrev != NULL);
       for (int h = 0; h < Q; h++) {
           pathPrev[h] = (double*)malloc(sizeof(double) * Q);
           assert(pathPrev[h] != NULL);
       }

       // initialization: PATH(v_i2h1)
       for (int i = 0; i < Q; i++) 
           for (int h = 0; h < Q; h++) {
              pathPrev[i][h] = INF;
              double eScore = escore(proInter, 1, pLatPrevPrev->p[h], pLatPrev->p[i]);
              if (eScore < INF - 1)
                  pathPrev[i][h] = pLatPrev->vscores[i] + eScore + pLatPrevPrev->vscores[h];
       }

       double** pathCurr = (double**)malloc(sizeof(double*) * Q);
       assert(pathCurr != NULL);
       for (int h = 0; h < Q; h++) {
           pathCurr[h] = (double*)malloc(sizeof(double) * Q);
           assert(pathCurr[h] != NULL);
       }

       int n = proInter.n; // protein length

       int hasSolution = RES_FAIL;

       for (int j = 2; j < n; j++) {

           hasSolution = RES_FAIL;

           set_lattice_around(*pLatCurr, lat, proInter.rs[j]);
           adjust_vscores(pLatCurr, proInter, j); 
       
           for (int i = 0; i < Q; i++) {

               for (int h = 0; h < Q; h++) {

                   pathCurr[i][h] = INF;

                   double eScore = escore(proInter, j, pLatPrev->p[h], pLatCurr->p[i]); 
              
                   if (eScore > INF - 1) 
                       continue;

                   double min = INF;

                   for (int g = 0; g < Q; g++) {

                       if (pathPrev[h][g] > INF - 1)
                           continue;

                       double score = pathPrev[h][g] + ascore(proInter, j, pLatCurr->p[i], pLatPrev->p[h], pLatPrevPrev->p[g]);

                       if (score < min) {
                            min = score; 
                            back[j][i][h] = (short)g;           
                            hasSolution = RES_OK;    
                       }
                   }

                   if (min < INF - 1) 
                       pathCurr[i][h] = pLatCurr->vscores[i] + eScore + min;
                }

            }

            if (j < n-1) {
               Lattice* pLatTmp = pLatPrevPrev;
               pLatPrevPrev = pLatPrev;
               pLatPrev = pLatCurr;
               pLatCurr = pLatTmp;

               double** pathTmp = pathPrev; 
               pathPrev = pathCurr;    
               pathCurr = pathTmp;
            } 

         if (hasSolution == RES_FAIL) {
sprintf(status, "could not handle the %d-th residue\n", j);
// write_to_log(status);
             break;
         }

//sprintf(status, "handled the %d-th residue\n", j);
//write_to_log(status);
      }

// if (hasSolution == RES_OK)
//     write_to_log("computed PATH values\n");
// else 
//     write_to_log("no solution for the current lattice\n");

 
     int indFinI = -1;
     int indFinH = -1;

     if (hasSolution == RES_OK) {
  
         // retrieve OESP solution 
         double minTotal = INF - 1;
       
           for (int i = 0; i < Q; i++) {
               for (int h = 0; h < Q; h++) {
                   if (pathCurr[i][h] < minTotal) {
                        minTotal = pathCurr[i][h];
                        indFinI = i;
                        indFinH = h;
                   }  
               }
           }
     }

     *indFinalI = indFinI;
     *indFinalH = indFinH;

     if (pathCurr != NULL) {
        for (int h = 0; h < Q; h++) 
           if (pathCurr[h] != NULL) free(pathCurr[h]);
        free(pathCurr);
     }
     if (pathPrev != NULL) {
       for (int h = 0; h < Q; h++) 
          if (pathPrev[h] != NULL) free(pathPrev[h]);
       free(pathPrev);
     }
         

     free_lattice(*pLatCurr);
     if (pLatCurr != NULL) free(pLatCurr);

     free_lattice(*pLatPrev);
     if (pLatPrev != NULL) free(pLatPrev);

     free_lattice(*pLatPrevPrev);
     if (pLatPrevPrev != NULL) free(pLatPrevPrev);
}

void adjustIntermediate(Lattice lat, int indFinalI, int indFinalH, Protein proInter, short*** back) {

// write_to_log("starting adjustment...\n");

     assert(indFinalI >= 0);
     assert(validate_protein(proInter));
     assert(back != NULL);
     for (int j = 0; j < proInter.n; j++) {
        assert(back[j] != NULL);
        for (int i = 0; i < lat.Q; i++)
            assert(back[j][i] != NULL);
     }

     int n = proInter.n;

     int indCurrI = indFinalI;
     int indCurrH = indFinalH;

     for (int j = n-1; j >= 0; j--) {
           	
        double vCorr[3];

        get_vertex(vCorr, lat, proInter.rs[j], indCurrI);

        for (int k=0; k < 3; k++) 
             proInter.rs[j][k] = vCorr[k];

       int indCurrH_tmp = indCurrH;
       if (j > 0)
           indCurrH = back[j][indCurrI][indCurrH];
       indCurrI = indCurrH_tmp;

       if (j > 1)
           assert((0 <= indCurrH) && (indCurrH < lat.Q));
       else
           assert(indCurrH == -1);
    }

// write_to_log("completed adjustment\n");
}

int proteinize(Protein proInter, double latticeEdgeLen) {

       assert(validate_protein(proInter));

       assert((latticeEdgeLen > MIN_LATTICE_EDGE_LEN - EPS) && (latticeEdgeLen < MAX_LATTICE_EDGE_LEN + EPS));

       int res = RES_FAIL;

       Lattice lat = init_lattice(latticeEdgeLen);

//log_lattice(lat);

       int Q = lat.Q; //lattice size

       int n = proInter.n; // protein length

       // allocate memory for an index array for backtracking, and initialize it with '-1's
       short*** back = (short***)malloc(sizeof(short**) * n);
       assert(back != NULL);

       for (int j = 0; j < n; j++) {

           back[j] = (short**)malloc(sizeof(short*) * Q);
           assert(back[j] != NULL);

           for (int i = 0; i < Q; i++) { 
                back[j][i] = (short*)malloc(sizeof(short) * Q);
                assert(back[j][i] != NULL);

                for (int h = 0; h < Q; h++) 
                    back[j][i][h] = -1.;                          
           } 
       }

       int indFinalI = -1;
       int indFinalH = -1;

       solveOESP(lat, proInter, back, &indFinalI, &indFinalH);

sprintf(status, "indFinalI = %d, indFinalH = %d\n", indFinalI, indFinalH);
// write_to_log(status);

       if (indFinalI >= 0) {
            assert(indFinalH >= 0);
            adjustIntermediate(lat, indFinalI, indFinalH, proInter, back);
            res = RES_OK;
// write_to_log("retrieved OESP solution\n");
       }
// 	else 
// write_to_log("no OESP solution\n");

       if (back != NULL) {
           for (int j = 0; j < n; j++) {
               if (back[j] != NULL) 
                   for (int i = 0; i < Q; i++) 
                        if (back[j][i] != NULL) free(back[j][i]); 
               free(back[j]);
           }
       }
       free(back);

       free_lattice(lat);

       return res;
}

#else


void solveOESP(Lattice lat, Protein proInter, short**back, int *indFinal) {

    assert(lat.Q > 0);
    assert(validate_protein(proInter));
    assert(back != NULL);
    for (int j = 0; j < proInter.n; j++) 
       assert(back[j] != NULL);
   assert(indFinal != NULL);

     Lattice* pLatPrev = (Lattice*)malloc(sizeof(Lattice));
     assert(pLatPrev != NULL);
     (*pLatPrev) = init_lattice_around(lat);
     set_lattice_around(*pLatPrev, lat, proInter.rs[0]);
     adjust_vscores(pLatPrev, proInter, 1); 

     Lattice* pLatCurr = (Lattice*)malloc(sizeof(Lattice));
     assert(pLatCurr != NULL);
     (*pLatCurr) = init_lattice_around(lat);


    int Q = lat.Q; //lattice size

    double* pathPrev = (double*)malloc(sizeof(double)*Q);
    assert(pathPrev != NULL);

    // initialization: PATH(v_i1)
    for (int i = 0; i < Q; i++) 
        pathPrev[i] = pLatPrev->vscores[i];

     double* pathCurr = (double*)malloc(sizeof(double)*Q);
     assert(pathCurr != NULL);


     int n = proInter.n; // protein length

     int hasSolution = RES_FAIL;

     for (int j = 1; j < n; j++) {

         hasSolution = RES_FAIL;
           
         set_lattice_around(*pLatCurr, lat, proInter.rs[j]);
         adjust_vscores(pLatCurr, proInter, j); 

         for (int i = 0; i < Q; i++) {

             pathCurr[i] = INF;

             double min = INF;

             for (int h = 0; h < Q; h++) {

                 if (pathPrev[h] > INF - 1)
                    continue;

                 double score = pathPrev[h] + escore(proInter, j, pLatPrev->p[h], pLatCurr->p[i]);

#ifdef HEURISTIC
                 if (j >= 2) {

                    int indPrevPrev = back[j-1][h];
                    assert((0 <= indPrevPrev) && (indPrevPrev < Q));

                    double vPrevPrev[3];
                    get_vertex(vPrevPrev, lat, proInter.rs[j-2], indPrevPrev);

                    score += ascore(proInter, j, vPrevPrev, pLatPrev->p[h], pLatCurr->p[i]);                   
                 }

#endif

                 if (score < min) {
                     min = score; 
                     back[j][i] = (short)h;   
                     hasSolution = RES_OK;            
                 }
             }

             if (min < INF - 1)
                 pathCurr[i] = pLatCurr->vscores[i] + min;
         }

         if (j < n-1) {
            Lattice* pLatTmp = pLatPrev;
            pLatPrev = pLatCurr;
            pLatCurr = pLatTmp;

            double* pathTmp = pathPrev; 
            pathPrev = pathCurr;    
            pathCurr = pathTmp; 
         }  

         if (hasSolution == RES_FAIL) {
sprintf(status, "could not handle the %d-th residue\n", j);
// write_to_log(status);
             break;
         }

//sprintf(status, "handled the %d-th residue\n", j);
//write_to_log(status);
      }

// if (hasSolution == RES_OK)
//     write_to_log("computed PATH values\n");
// else 
//     write_to_log("no solution for the current lattice\n");


     int indFin = -1; 

     if (hasSolution == RES_OK) {

        // retrieve OESP solution 
        double minTotal = INF - 1;

        for (int i = 0; i < Q; i++) {
           if (pathCurr[i] < minTotal) {
                minTotal = pathCurr[i];
                indFin = i;
           }  
        }
     }

     *indFinal = indFin;
   
      
     if (pathCurr != NULL) free(pathCurr);
     if (pathPrev != NULL) free(pathPrev);


     free_lattice(*pLatCurr);
     if (pLatCurr != NULL) free(pLatCurr);

     free_lattice(*pLatPrev);
     if (pLatPrev != NULL) free(pLatPrev);
}

void adjustIntermediate(Lattice lat, int indFinal, Protein proInter, short** back) {

     assert(indFinal >= 0);
     assert(validate_protein(proInter));
     assert(back != NULL);
     for (int j = 0; j < proInter.n; j++) 
         assert(back[j] != NULL);

     int n = proInter.n;

     int indCurr = indFinal;

     for (int j = n-1; j >= 0; j--) {
           	
        double vCorr[3];

        get_vertex(vCorr, lat, proInter.rs[j], indCurr);

        for (int k=0; k < 3; k++) 
             proInter.rs[j][k] = vCorr[k];

       indCurr = back[j][indCurr];

       if (j > 0)
           assert((0 <= indCurr) && (indCurr < lat.Q));
       else
           assert(indCurr == -1);
    }
}

int proteinize(Protein proInter, double latticeEdgeLen) {

       assert(validate_protein(proInter));

       assert((latticeEdgeLen > MIN_LATTICE_EDGE_LEN - EPS) && (latticeEdgeLen < MAX_LATTICE_EDGE_LEN + EPS));

       int res = RES_FAIL;

       Lattice lat = init_lattice(latticeEdgeLen);

//log_lattice(lat);

       int Q = lat.Q; //lattice size

       int n = proInter.n; // protein length

       // allocate memory for an index array for backtracking, and initialize it with '-1's
       short** back = (short**)malloc(sizeof(short*) * n);
       assert(back != NULL);

       for (int j = 0; j < n; j++) {

           back[j] = (short*)malloc(sizeof(short) * Q);
           assert(back[j] != NULL);

           for (int i = 0; i < Q; i++) 
                back[j][i] = -1;
       }

       int indFinal = -1;

       solveOESP(lat, proInter, back, &indFinal);

       if (indFinal >= 0) {
            adjustIntermediate(lat, indFinal, proInter, back);
            res = RES_OK;
// write_to_log("retrieved OESP solution\n");
       }
// 	else 
// write_to_log("no OESP solution\n");

       for (int j = 0; j < n; j++) 
           if (back[j] != NULL) free(back[j]);
       free(back);

       free_lattice(lat);

       return res;
}

#endif

int morph(Protein *p, int k) {

    assert(p != NULL);
    assert(k >= 0);

    int res = RES_FAIL;

    for(int i = 1; i <= k; i++) {

        double alpha = 1.0 / (k + 2 - i);

sprintf(status, "computing %d-th intermediate, alpha=%f... \n", i, alpha);
// write_to_log(status); 

        Protein proInter = intermediate(p[i - 1], p[k+1], alpha);

sprintf(status, "done\n");
// write_to_log(status); 

sprintf(status, "proteinizing %d-th intermediate... \n", i);
// write_to_log(status); 

         double latticeEdgeLen = MIN_LATTICE_EDGE_LEN;

         while (latticeEdgeLen < MAX_LATTICE_EDGE_LEN + EPS) {

  		res = proteinize(proInter, latticeEdgeLen);

              if (res == RES_OK) 
                  break;

              latticeEdgeLen += LATTICE_EDGE_INCREM;  
	  }

         if (res == RES_OK) {
              p[i] = copy_protein(proInter);
              free_protein(proInter); 
sprintf(status, "done\n");
// write_to_log(status); 
	  }
         else { 
sprintf(status, "failed\n");
// write_to_log(status); 
              free_protein(proInter); 
              break;
	  }

    }

    return res;     
}

//---------------------------------------------------------------------------------------
//
//  Cleaning up
//
//---------------------------------------------------------------------------------------

void free_protein_sequence(Protein* p, int k) {

    assert(k >= 0);    // k - number of intermediate states

    assert(p != NULL); 

    for (int i=0; i < k+2; i++) 
        free_protein(p[i]);
    
    free(p);  
}

//---------------------------------------------------------------------------------------
//
//  Main function
//
//---------------------------------------------------------------------------------------

int main(int _n, char** args) {

init_log();
// write_to_log("starting...\n");

    init();

    int n = -1; // number of residues in the chains being morphed
    int k = -1; // number of intermediate states

    scanf("%d", &n);

    Protein proInitial = new_protein(n);   // initial conformation
    scan_protein(proInitial, n);

    Protein proFinal   = new_protein(n);   // final conformation
    scan_protein(proFinal, n);

    scanf("%d", &k);

    assert( n >  0); 
    assert( k >= 0);

sprintf(status, "Morphing will be performed between two chains with %d residues each.\n", n);   
// write_to_log(status);
sprintf(status, "The morphing path will pass through %d intermediate states.\n", k);
// write_to_log(status);

/*
double distMax = -1.;
for (int j = 0; j < proInitial.n; j++) {
   double dist = sqrt(dist_2(proInitial.rs[j], proFinal.rs[j]));
      if (dist > distMax)
            distMax = dist;
   sprintf(status, "distance between the %d-th residues: %8.4f\n", j, dist);
   write_to_log(status);
}
sprintf(status, "max distance between residues: %8.4f\n", distMax);
write_to_log(status);
*/


    Protein *p = (Protein*)malloc(sizeof(Protein)*(k+2));
    assert(p != NULL); 

    p[0]   = copy_protein(proInitial);
    p[k+1] = copy_protein(proFinal);

    free_protein(proInitial);
    free_protein(proFinal);
    
    int res = morph(p, k);

    if (res == RES_FAIL) {
// write_to_log("Morphing failed.\n");       
       printf("FAIL");
       free_protein_sequence(p, k);
       return 1;
    }

    assert(res == RES_OK);

// write_to_log("Morphing completed.\n");

    printf("SUCC ");
    printf("%d\n", k + 2);
    for(int i = 0; i < k + 2; ++i ) {
        println_protein(p[i]);
    }

/*
for (int i=0; i < k+2; i++) {
  sprintf(status, "logging the %d-th conormation:\n", i);
  write_to_log(status);
  log_protein(p[i]); 
}
*/

    free_protein_sequence(p, k);

    return 0;
}

