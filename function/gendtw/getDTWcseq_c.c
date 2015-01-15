/**
 * (C) Copyright 2014 Víctor Ponce-López <vponcel@uoc.edu>
 */

/** 
 * This is the C/MEX code that computes the updated minimum costs given a 
 * DTW matrix
    % input:
    %   W: Warping cost sequence
    % output    
    %   s: array of updated minimum costs
 *
 * compile: 
 *     mex getDTWcseq_c.c
 *
 * usage:
 *     s=getDTWcseq_c(W)  
 *     where W is the Cost Matrix obtained from the dtw_c function which 
 *     computes the DTW matrix between two multi-dimensional signals
 */

#include "matrix.h"
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

float * getDTWcseq_c(float *w, int ns, int m)
{
    int i,j,k,se[2],a;
    float * s;
    bool foundPath;
    
    if (m<=1) return NULL;
    /*printf("%d\n",m);*/
    s=(float *)malloc((m-1)*sizeof(float));
    for (i=0;i<m-1;i++)
    {
        s[i]=FLT_MAX;
    }
    
    /* Obtain start-end from the optimal path */
    /*printf("%f ",w[j*(ns-1) +k]);*/
    for (i=1;i<=m;i++)
    {
        se[1]=i;
        j=i;
        se[0]=j;
        k=j-1;
        foundPath = false;
        while (!foundPath && j>1)
        {   
            --k;        
            if (w[(j-1)*(ns-1) +(k-1)] <= w[(j-1)*(ns-1) +k])
            {
                if (w[(j-1)*(ns-1) +(k-1)] <= w[j*(ns-1) +k])
                {
                    --j;
                    --k;
                }
                else if (w[(j-1)*(ns-1) +k] <= w[j*(ns-1) +k])
                {
                    --j;
                }                
            }
            else
            {
                if (w[(j-1)*(ns-1) +k] <= w[j*(ns-1) +k])
                {
                    --j;
                }
            }
            if (w[j*(ns-1) +k] == 0 && j>=1)
            {
                foundPath = true;
                se[0] = j;
            }
            
            /*            
            printf("%f ",w[j*(ns-1) +k]);
            printf("%f ",w[(j-1)*(ns-1) +k]);
            printf("%f ",w[(j-1)*(ns-1) +(k-1)]);
            */
        }
        /*
        printf("%d ",se[0]);
        printf("%d ",se[1]);
        printf("\n");
        */
        k=se[1]-1;
        for (a=se[0];a<se[1];a++)
        {
            if (w[se[1]*(ns-1) +k] < s[a-1])
            {
                s[a-1]=w[se[1]*(ns-1) +k];
                /*printf("%.2f ",s[a-1]);*/
            }                
        }
    }
    
    return s;
}

/* the gateway function */

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    int ns,m;
    int i,j,dims[2];
    float *w,*s;
    float *dp;
	
    /*  check for proper number of arguments */
    if(nrhs!=1)
    {
        mexErrMsgIdAndTxt( "MATLAB:aligngesture_c:invalidNumInputs",
                "One unique argument required.");
    }
    if(nlhs!=1)
    {
        mexErrMsgIdAndTxt( "MATLAB:dtw_c:invalidNumOutputs",
                "getDTWcseq_c: One output required.");
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
   
	s = (float *)getDTWcseq_c(w,ns,m);
    /*
    for (i=0;i<ns;i++)
    {
        printf("%d ",s[i]);
    }
    */

    dims[0] = 1;
    dims[1] = m-1;
    plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    
	/*  create a C pointer to a copy of the output matrix */

    dp = mxGetPr(plhs[0]);
    
    for (i=0;i<m-1;i++)
    {
        memcpy(&dp[i], &s[i], sizeof(float));
    }

    free(s);
    
    return;    
}