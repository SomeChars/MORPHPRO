/*******************************************************************************
 *  -/_|:|_|_\- 
 *
 *  File:           main.c
 *
 *  Function:       Rapid calculation of the least-squares rotation using a 
 *                  quaternion-based characteristic polynomial and 
 *                  a cofactor matrix
 *
 *  Author(s):      Douglas L. Theobald
 *                  Department of Biochemistry
 *                  MS 009
 *                  Brandeis University
 *                  415 South St
 *                  Waltham, MA  02453
 *                  USA
 *
 *                  dtheobald@brandeis.edu
 *                  
 *                  Pu Liu
 *                  Johnson & Johnson Pharmaceutical Research and Development, L.L.C.
 *                  665 Stockton Drive
 *                  Exton, PA  19341
 *                  USA
 *
 *                  pliu24@its.jnj.com
 * 
 *
 *    If you use this QCP rotation calculation method in a publication, please
 *    reference:
 *
 *      Douglas L. Theobald (2005)
 *      "Rapid calculation of RMSD using a quaternion-based characteristic
 *      polynomial."
 *      Acta Crystallographica A 61(4):478-480.
 *
 *      Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2009)
 *      "Fast determination of the optimal rotational matrix for weighted 
 *      superpositions."
 *      in press, Journal of Computational Chemistry
 *
 *
 *  Copyright (c) 2009, Pu Liu and Douglas L. Theobald
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without modification, are permitted 
 *  provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this list of 
 *    conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice, this list 
 *    of conditions and the following disclaimer in the documentation and/or other materials 
 *    provided with the distribution.
 *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to 
 *    endorse or promote products derived from this software without specific prior written 
 *    permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *
 *  Source:         started anew.
 *
 *  Change History:
 *    2009/04/13      Started source
 *  
 ******************************************************************************/

 /* Sample code to use the routine for fast RMSD & rotational matrix calculation
 For the example provided below, the minimum least-squares RMSD for the two
    7-atom fragments should be 0.719106 A.

    The rotation quaterion should be:

    -8.620063e-01   3.435505e-01   1.242953e-01  -3.513814e-01 

    And the corresponding 3x3 rotation matrix is:

    [     0.72216358     0.69118937    -0.02714790 ]
    [    -0.52038257     0.51700833    -0.67963547 ]
    [    -0.45572112     0.50493528     0.73304748 ]

	Compile instruction:

	 make

	 How to incorporate the code into your own package.

	 1. copy the qcprot.h and qcprot.c into the source directory
	 2. change your own code to call the fast rotational routine and include the qcprot.h
	 3. change your make file to include qcprot.c
*/

#include "qcprot.h"

double **MatInit(const int rows, const int cols);

void MatDestroy(double ***matrix);

static void Mat3Print(double *matrix);

void PrintCoords(const double **coords, const int len);


int main()
{  
    double          rmsd;
    double        **frag_a, **frag_b;
    int             len;
    scanf("%d", &len);
    double        rotmat[9];

    frag_a = (double **)malloc(3 * sizeof(double *));
    frag_b = (double **)malloc(3 * sizeof(double *));
    for(int i = 0; i < 3; i ++) {
        frag_a[i] = (double *) malloc(len * sizeof(double));
        frag_b[i] = (double *) malloc(len * sizeof(double));
    }
    for(int i = 0; i < len; i ++) {
        for(int j = 0; j < 3; j ++)
        scanf("%lf", &frag_a[j][i]);
    }
    for(int i = 0; i < len; i ++) {
        for(int j = 0; j < 3; j ++)
        scanf("%lf", &frag_b[j][i]);
    }
    rmsd = CalcRMSDRotationalMatrix((double **) frag_a, (double **) frag_b, len, rotmat, NULL);
	/*
    PrintCoords((const double **) frag_a, len);
    PrintCoords((const double **) frag_b, len);
	*/
    //printf("\nrmsd = %f\n", rmsd);

    printf("%lf ", rmsd);
    for(int i = 0; i < 9; i ++) printf("%lf ", rotmat[i]);
//	Mat3Print(rotmat);

    free(frag_a);
    free(frag_b);

    return 0;
}


double **MatInit(const int rows, const int cols)
{
    int             i;
    double        **matrix = NULL;
    double         *matspace = NULL;

    matspace = (double *) calloc((rows * cols), sizeof(double));
    if (matspace == NULL)
    {
        perror("\n ERROR");
        printf("\n ERROR: Failure to allocate matrix space in MatInit(): (%d x %d)\n", rows, cols);
        exit(EXIT_FAILURE);
    }

    /* allocate room for the pointers to the rows */
    matrix = (double **) malloc(rows * sizeof(double *));
    if (matrix == NULL)
    {
        perror("\n ERROR");
        printf("\n ERROR: Failure to allocate room for row pointers in MatInit(): (%d)\n", rows);
        exit(EXIT_FAILURE);
    }

    /*  now 'point' the pointers */
    for (i = 0; i < rows; i++)
        matrix[i] = matspace + (i * cols);

    return(matrix);
}


void MatDestroy(double ***matrix_ptr)
{
    double **matrix = *matrix_ptr;

    if (matrix != NULL)
    {
        if (matrix[0] != NULL)
        {
            free(matrix[0]);
            matrix[0] = NULL;
        }

        free(matrix);
        *matrix_ptr = NULL;
    }
}


static void Mat3Print(double *matrix)
{
    int             i;

    printf("\n");
    for (i = 0; i < 3; ++i)
    {
        printf(" [ % 14.8f % 14.8f % 14.8f ]\n",
               matrix[3 * i],
               matrix[3 * i + 1],
               matrix[3 * i + 2]);
    }

    fflush(NULL);
}


void PrintCoords(const double **coords, const int len)
{
    int             i;

    for (i = 0; i < len; ++i)
        printf("\n % 8.3f % 8.3f % 8.3f", coords[0][i], coords[1][i], coords[2][i]);
    putchar('\n');
}


