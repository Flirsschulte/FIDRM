/**************************************************************************
% Fuzzy Fixed-valued Impulse Noise Detection Method 
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

/********************************************************************************************
*  The main function for the calculation of the membership degrees in the fuzzy set noise
*********************************************************************************************/
void callMem(double **A1,double *spv,int M, int N) { 
   int i,j,k, loop;   
   double *med, *histo,*sorted,*maxi;
   double som, tel, hlp,res, tmp, slope1, diff1, diff2, diff3;
   double l1, l2, l3, minslop1, minslop2, slope2;
   double a = 80.0, b = 160.0;
   int **xy, **adj;

   int rand1a = 0;
   int rand1b = 0;
   int rand2a = 0;
   int rand2b = 0;
   
   xy = malloc(9*sizeof(int));
   adj = malloc(9*sizeof(int));
   for(i=0;i<9;i++){
      xy[i] = malloc(2*sizeof(int));
      adj[i] = malloc(2*sizeof(int));
   }
   
   med = malloc(8*sizeof(double));
   histo = malloc(256*sizeof(double));
   sorted = malloc(256*sizeof(double));
   maxi = malloc(10*sizeof(double));
      
   xy[0][0] = -1; xy[0][1] = -1; xy[1][0] = -1; xy[1][1] = 0; 
   xy[2][0] = -1; xy[2][1] =  1; xy[3][0] =  0; xy[3][1] = -1; 
   xy[4][0] =  0; xy[4][1] =  0; xy[5][0] =  0; xy[5][1] = 1; 
   xy[6][0] =  1; xy[6][1] = -1; xy[7][0] =  1; xy[7][1] = 0; 
   xy[8][0] =  1; xy[8][1] =  1; 

   adj[0][0] = 1; adj[0][1] = -1; adj[1][0] =  0; adj[1][1] = -1; 
   adj[2][0] = 1; adj[2][1] =  1; adj[3][0] =  1; adj[3][1] =  0; 
   adj[4][0] = 0; adj[4][1] =  0; adj[5][0] =  1; adj[5][1] =  0; 
   adj[6][0] = 1; adj[6][1] =  1; adj[7][0] =  0; adj[7][1] = -1; 
   adj[8][0] = 1; adj[8][1] = -1; 
               
      
   for(i=2; i<M-2; i++){
      for(j=2; j<N-2; j++){
         /**************************************
         *       Gradient values method
         **************************************/
         tel = 0;
         for (k = 0; k<9; k++){
            if (k==4) continue;
            diff1 = A1[i+xy[k][1]][j+xy[k][0]] - A1[i][j];
            diff2 = A1[i+xy[k][1]+adj[k][1]][j+xy[k][0]+adj[k][0]] - A1[i+adj[k][1]][j+adj[k][0]];
            diff3 = A1[i+xy[k][1]-adj[k][1]][j+xy[k][0]-adj[k][0]] - A1[i-adj[k][1]][j-adj[k][0]];
            
            l1 = LARGE(absol(diff1),a,b);
            l2 = LARGE(absol(diff2),a,b);
            l3 = LARGE(absol(diff3),a,b);
            
            res = 1.0-  maximum( maximum(minimum(1-l1,1-l2), minimum(1-l1,1-l3)), maximum(minimum(l1,l2),minimum(l1,l3)));
            res = l1*res;
            
            for (loop = 0; loop<tel; loop++){
                if (res < med[loop]){
                   hlp = med[loop];
                   med[loop] = res;
                   res = hlp;
                }
            }
            med[(int)tel] = res;
            tel++;
         }
         if (minimum( minimum(med[6],med[7]),med[5]) >= 0.5) histo[(int)A1[i][j]]++;
      }
   }

   
   for (i=0; i<256; i++) sorted[i] =  histo[i];

   for (i=0; i<256; i++) 
      for (j=i+1; j<256; j++)
         if(sorted[j]>sorted[i]){
            tmp = sorted[j];
            sorted[j] = sorted[i];
            sorted[i] = tmp;
         }
   
   som = 0.0;
   for (j = 0; j<10; j++)
      for (i = 0; i<256; i++){
         if (j==0) som +=  histo[i];
         if( histo[i]==sorted[j]) maxi[j] = i;
      }
   
   // Test if fixed-valued impulse noise is detected or not
   if ((sorted[0]/som) < 0.1){ 
     printf("\t NO FIXED VALUE IMPULSE NOISE DETECTED \n");
     spv[0] = -1;
   }
   else{
      minslop1 = 20;
      minslop2 = 20;
      for (i = 0; i<10; i++){
         slope1 = 0;
         for (j = maxi[i]+1; j< minimum(256,maxi[i]+31);j++){
            if( histo[j] <=  histo[(int)maxi[i]]*0.015) break;
            slope1++;
         }
         slope1 = minimum(slope1, minslop1);
         minslop1 = minimum(slope1,minslop1);

         slope2 = 0;
         for (j = maxi[i]-1; j >= maximum(0,maxi[i]-30);j--){
            if( histo[j] <=  histo[(int)maxi[i]]*0.015) break;
            slope2++;
         }
         slope2 = minimum(slope2,minslop2);
         minslop2 = slope2;
      
         spv[i+0*11] = maxi[i];
         spv[i+1*11] = maxi[i] - minimum(slope2,30)-(minimum(slope2,30)*(1.0/3.0));
      
         spv[i+2*11] = maxi[i] - minimum(slope2,30);
         spv[i+3*11] = maxi[i] + minimum(slope1,30);
         spv[i+4*11] = maxi[i] + minimum(slope1,30)+(minimum(slope1,30)*(1.0/3.0));
      }

   
      spv[10+0*5] = 0;
      for (i =0; i<10; i++){
         if ((sorted[i]/som) < 0.005) break;
         spv[10+0*5] = spv[10+0*5] + 1;
      }
   }
   
   for(i=0;i<9;i++)  free(xy[i]);
   free(xy); 
   for(i=0;i<9;i++)  free(adj[i]);
   free(adj); 
   free(med); 
   free(histo); 
   free(sorted); 
   free(maxi); 
   
}  


#define Im1      prhs[0]

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
    int row, col, i, M, N;
    double **spv, **A1;
    
    if (nlhs!=1)
        mexErrMsgTxt("It requires one output arguments only [M1].");
    if (nrhs!=1)
       mexErrMsgTxt("this method requires one input argument [Im1]");

    /* Get input values */  
    M = mxGetM(Im1);
    N = mxGetN(Im1);

    /**
    * Allocate memory for return matrices 
    **/
    OUT = mxCreateDoubleMatrix(11, 5, mxREAL);  
    spv = mxGetPr(OUT);

    /**
    * Dynamic allocation of memory for the input array
    **/
    A1 = malloc(M*sizeof(int));
    for(i=0;i<M;i++)
      A1[i] = malloc(N*sizeof(double));

     /**
     * Convert ARRAY_IN and INPUT_MASK to 2x2 C arrays (MATLAB stores a two-dimensional matrix 
     * in memory as a one-dimensional array) 
     ***/
     for (col=0; col < N; col++)
         for (row=0; row < M; row++) {
             A1[row][col] = (mxGetPr(Im1))[row+col*M];
	      }
	      
    callMem(A1,spv,M,N);

    for(i=0;i<M;i++)  free(A1[i]);
    free(A1); 
}
/* end mexFcn*/