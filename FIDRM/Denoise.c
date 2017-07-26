/**************************************************************************
% Fuzzy Fixed-valued Impulse Noise Reduction Method 
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

double NOISE(double x, double ** noise, int amount) {
   double ret, tmp;
   int i;
   ret = 0.0; tmp = 0.0;
   
   for (i = 0; i<amount; i++){
      if((x == noise[i][0]) || ((x < noise[i][0]) && (x >= noise[i][1])) || ((x > noise[i][0]) && (x<=noise[i][3])))
          tmp = 1.0;
      else if((x<noise[i][1]) && (x>noise[i][2]))  tmp = ((noise[i][1]-x)/(noise[i][1]-noise[i][2]));
      else if ((x>noise[i][3])&&(x<noise[i][4])) tmp = ((noise[i][4]-x)/(noise[i][4]-noise[i][3]));
      else tmp = 0.0;
      ret = maximum(ret,tmp);
   }
   return ret;
}


/**************************************************************************************
*  The main function for the calculation of the FIDRM output
*   for the fixed-valued impulse noise case 
***************************************************************************************/
void callDenoise(double **A1,double **fixed,double *filt,int M, int N,int W) { 
   int i,j,k,l, amount,it, loop;   
   double som, wei, cen, noise, tmp, res, tel;
   int *over_i, *over_j, OVER, OVER2, ok;

   int rand1a = 0;
   int rand1b = 0;
   int rand2a = 0;
   int rand2b = 0;
        
   over_i = malloc(M*N*sizeof(double));
   over_j = malloc(M*N*sizeof(double));

   amount = fixed[10][0];
   it = 1;
   OVER = 0;
   OVER2 = 0;
   while((it == 1)||(OVER != 0)) {
      if(it==1){
          for(i=2; i<M-2; i++){
              for(j=2; j<N-2; j++){
                   ok = 0;
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
                 /* end step 1. */      
         
                 cen = A1[i][j];
                 noise = NOISE(cen, fixed, amount);
                 if(noise > 0){
                    som = 0.0; wei = 0.0;
                    for (k=i-rand1a; k<=i+rand1b; k++){
                       for (l=j-rand2a; l<=j+rand2b; l++){
                          tmp = 1.0 - NOISE(A1[k][l], fixed, amount);
                          wei += tmp;
                          som += tmp*A1[k][l];
                       }
                    }
                    if (wei == 0) res = A1[k][l];
                    else   res = som/wei;
                    noise = NOISE(res, fixed, amount);
                    if (noise > 0){
                       over_i[OVER2] = i;
                       over_j[OVER2] = j;
                       OVER2++;
                       ok = 1;
                    }
                 }
                 else res = A1[i][j];
                 filt[i+j*M] = res;
              }
          }
      }
      else {
         OVER2 = 0;
         ok = 0;
         for (loop = 0; loop<OVER;loop++){
            i = over_i[loop];
            j = over_j[loop];
 
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
            /* end step 1. */      

            if(NOISE(A1[i][j], fixed, amount) >0){
                som = 0.0; wei = 0.0; tel = 0;
                for (k=i-rand1a; k<=i+rand1b; k++){
                   for (l=j-rand2a; l<=j+rand2b; l++){
                      tmp = 1.0 - NOISE(filt[k+l*M], fixed, amount);
                      wei += tmp;
                      som += tmp*filt[k+l*M];
                      tel++;
                   }
                }
                if(tel!=0){
                   if (wei ==0) res = filt[i+j*M];
                   else res = som/wei;
                }
                else res = filt[i+j*M];
                
                if (NOISE(res, fixed, amount) > 0){
                   over_i[OVER2] = i;
                   over_j[OVER2] = j;
                   OVER2++;
                   ok = 1;
                }
            }
            else res = filt[i+j*M];
            filt[i+j*M] = res;
         }
      }
      it++;
      W++;
      if ((OVER == OVER2)||(it==15)) break;
      OVER = OVER2;
   }
   free(over_i); 
   free(over_j); 
} 


#define Im1      prhs[0]
#define FIXED    prhs[1]
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
    double **fixed, **A1, **filt;
    
    if (nlhs!=1)
        mexErrMsgTxt("It requires one output arguments only [OUT].");
    if (nrhs!=3)
       mexErrMsgTxt("this method requires three input argument [Im1, FIXED, WSIE]");

    /* Get input values */  
    M = mxGetM(Im1);
    N = mxGetN(Im1);

    M2 = mxGetM(FIXED);
    N2 = mxGetN(FIXED);
    
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

    fixed = malloc(M2*sizeof(int));
    for(i=0;i<M2;i++)
      fixed[i] = malloc(N2*sizeof(double));

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
             fixed[row][col] = (mxGetPr(FIXED))[row+col*M2];
	      }
	      

    callDenoise(A1,fixed,filt,M,N,W);

    for(i=0;i<M;i++)  free(A1[i]);
    free(A1); 
    for(i=0;i<M2;i++)  free(fixed[i]);
    free(fixed); 
}
/* end mexFcn*/