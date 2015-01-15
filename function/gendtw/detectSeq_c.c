/**
 * (C) Copyright 2014 Víctor Ponce-López <vponcel@uoc.edu>
 */

/** 
 * This is the C/MEX code to detect start-end from the optimal dynamic time
 * warping path from two multi-dimensional signals
 *
 * compile: 
 *     mex detectSeq_c.c
 *
 * usage:
 *     se=detectSeq_c(W)  
 *     where W is the Cost Matrix obtained from the dtw_c function which 
 *     computes the DTW matrix between two multi-dimensional signals
 */

#include "matrix.h"
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

int * detectSeq_c(float *w,int ns, int m)
{
    int i,j,k;
    int * se;
    float last,minw;
    bool foundPath;
    
    se = (int *)malloc(2*sizeof(int));
    last = w[m*(ns-1) +(m-1)];
       
    /* Obtain start-end from the optimal path */
    se[1] = m;
    j=m;
    k=j-1;    
    foundPath = false;
    /* printf("%f ",w[j*(ns-1) +(k)]); */
    while (!foundPath && j > 0)
    {
        se[0] = j;
        --k;        
        if (w[(j-1)*(ns-1) +(k-1)] <= w[(j-1)*(ns-1) +(k)])
        {
            if (w[(j-1)*(ns-1) +(k-1)] <= w[j*(ns-1) +(k)])
            {
                --j;
                --k;
            }
            else if (w[(j-1)*(ns-1) +(k)] <= w[j*(ns-1) +(k)])
            {
                --j;
            }                
        }
        else
        {
            if (w[(j-1)*(ns-1) +(k)] <= w[j*(ns-1) +(k)])
            {
                --j;
            }
        }
        if (w[j*(ns-1) +(k)] == 0 || j==1)
        {
            foundPath = true;
            if (j==1) se[0]=1;
        }        
        
        /*
        printf("%f ",w[j*(ns-1) +(k)]);
        printf("%f ",w[(j-1)*(ns-1) +(k)]);
        printf("%f ",w[(j-1)*(ns-1) +(k-1)]);
         */        
    }  
    
    
    /* Print se 
    for (i=0;i<2;i++)
    {
        printf("%d ",se[i]);
    }
    */
    
    return se;
}

/* the gateway function */

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    int ns,m;
    int dims[2];
	int i,j;
    float *w;
    int *dp,*se;
	
    /*  check for proper number of arguments */

    if(nrhs!=1)
    {
        mexErrMsgIdAndTxt( "MATLAB:aligngesture_c:invalidNumInputs",
                "One unique argument required.");
    }
    if(nlhs>1)
    {
        mexErrMsgIdAndTxt( "MATLAB:dtw_c:invalidNumOutputs",
                "aligngesture_c: One output required.");
    }
    
    /*  create a pointer to the input matrix w */

    w = mxGetPr(prhs[0]);
    
    /*  get the dimensions of the matrix input w */

    ns = mxGetM(prhs[0]);
    m = mxGetN(prhs[0]);
      
    /*
    for(i=0;i<ns;i++)
    {
        for(j=0;j<m;j++)
        {
            printf("%f  ",w[j*(ns) +i]);
        }
        printf("\n");
    }
    */
    
    /*  set the output pointer to the output vector */
   
	se = detectSeq_c(w,ns,m);
    /*
    for (i=0;i<2;i++)
    {
        printf("%d ",se[i]);
    }
    */

    dims[0] = 1;
    dims[1] = 2;
    plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    
	/*  create a C pointer to a copy of the output matrix */

    dp = mxGetPr(plhs[0]);
    
    for (i=0;i<2;i++)
    {
        memcpy(&dp[i], &se[i], sizeof(int));
    }
    
    /* Release . Also tested with mxFree() when allocating using mex calls */
    /*free(se); */

    
    return;
    
}
