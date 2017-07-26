/**************************************************************************
% Fuzzy Random-valued Impulse Noise Reduction Method 
%
%  The paper of the FIDRM is proposed in: 
%
%  Stefan Schulte, Mike Nachtegael, Valérie De Witte, 
%  Dietrich Van Der Weken and  Etienne E. Kerre:
%  A Fuzzy Impulse Noise Detection and Reduction Method.
%  IEEE Transactions on Image Processing 15(5), 2006, 1153-1162
%  
% Stefan Schulte (stefan.schulte@Ugent.be):
% Last modified: 03/01/06
%
%**************************************************************************/

#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double LARGE (double x, double p1, double p2) {
   double res = 0.0;
   if ((x > p1) && (x < p2))   res = (x-p1)/(p2-p1);
   else if (x > p2)            res = 1.0;
   else                        res = 0.0;
   return res;
}

double absol(double a) {
   double b;
   if(a<0)   b=-a;
   else   b=a;
   return b;
}

double minimum(double a, double b) {
   double x;
   if(a<=b)   x = a;
   else   x=b;
   return x;
}

double maximum(double a, double b) {
   double x;
   if(a<=b)   x = b;
   else   x=a;
   return x;
}


/**************************************************************************************
*  The main function for the calculation of the FIDRM output
*   for the random-valued impulse noise case 
***************************************************************************************/
void callDenoise2(double **A1,double **mem,double *filt,int M, int N,int W) { 
   int i,j,k,l;   
   double som, wei;

   int rand1a = 0;
   int rand1b = 0;
   int rand2a = 0;
   int rand2b = 0;
        
   for(i=0; i<M; i++){
      for(j=0; j<N; j++){
         /* step 1. Determine the local window*/      
         if(i < W) {
            rand1a = i;
            rand1b = W;
         }
         else {
              if (i>M-W-1){
               rand1a = W;
               rand1b = M-i-1;
            }
            else{
               rand1a = W;
               rand1b = W;
            }
         }

         
         if(j < W) {
            rand2a = j;
            rand2b = W;
         }
         else {
            if (j > N-W-1){
               rand2a = W;
               rand2b = N-j-1;
            }
            else{
               rand2a = W;
               rand2b = W;
            }
         }
          
         if(mem[i][j] > 0.2){
             som = 0; wei = 0;
             for (k=-rand1a; k<=rand1b; k++){
                 for (l=-rand2a; l<=rand2b; l++){
                     som += (1.0-mem[i+k][j+l])*A1[i+k][j+l];
                     wei += (1.0-mem[i+k][j+l]);
                 }
             }
             if (wei != 0)
                 filt[i+j*M] = (som/wei);
             else filt[i+j*M] = A1[i][j];
         }
         else filt[i+j*M] =  A1[i][j];
      }
   }
}  /* End of callFuzzyShrink */


#define Im1      prhs[0]
#define MEM      prhs[1]
#define WSIZE    prhs[2]

#define OUT plhs[0]

/**
*  The interaction with Matlab (mex):
*        nlhs = amount of output arguments (= 1)
*        nrhs = amount of input arguments (= 3)
*     *plhs[] = link to the output 
*     *prhs[] = link to the input 
*
**/
void mexFunction( int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[] ) {
    int row, col, i, M, N, M2, N2,W;
    double **mem, **A1, **filt;
    
    if (nlhs!=1)
        mexErrMsgTxt("It requires one output arguments only [OUT].");
    if (nrhs!=3)
       mexErrMsgTxt("this method requires three input argument [Im1, FIXED, WSIE]");

    /* Get input values */  
    M = mxGetM(Im1);
    N = mxGetN(Im1);

    M2 = mxGetM(MEM);
    N2 = mxGetN(MEM);
    
    W = mxGetScalar(WSIZE);

    /**
    * Allocate memory for return matrices 
    **/
    OUT = mxCreateDoubleMatrix(M, N, mxREAL);  
    filt = mxGetPr(OUT);

    /**
    * Dynamic allocation of memory for the input array
    **/
    A1 = malloc(M*sizeof(int));
    for(i=0;i<M;i++)
      A1[i] = malloc(N*sizeof(double));

    mem = malloc(M2*sizeof(int));
    for(i=0;i<M2;i++)
      mem[i] = malloc(N2*sizeof(double));

     /**
     * Convert ARRAY_IN and INPUT_MASK to 2x2 C arrays (MATLAB stores a two-dimensional matrix 
     * in memory as a one-dimensional array) 
     ***/
     for (col=0; col < N; col++)
         for (row=0; row < M; row++) {
             A1[row][col] = (mxGetPr(Im1))[row+col*M];
	      }
	      
     for (col=0; col < N2; col++)
         for (row=0; row < M2; row++) {
             mem[row][col] = (mxGetPr(MEM))[row+col*M2];
	      }
	      
	callDenoise2(A1,mem,filt,M,N,W);

    for(i=0;i<M;i++)  free(A1[i]);
    free(A1); 
    for(i=0;i<M2;i++)  free(mem[i]);
    free(mem); 
}
/* end mexFcn*/