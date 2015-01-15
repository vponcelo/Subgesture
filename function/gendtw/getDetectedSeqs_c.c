/**
 * (C) Copyright 2014 Víctor Ponce-López <vponcel@uoc.edu>
 */

/** 
 * This is a very fast C/MEX code that returns a logical vector with the 
 * detected frames obtained from the optimal dynamic time warping paths. 
 * Thus, given 1a DTW cost matrix between two multi-dimensional signals, 
 * the start-end positions of the optimal path are obtained using the 
 * detectSeq_c function. 
 *
 * compile: 
 *     mex getDetectedSeqs_c.c
 *
 * usage:
 *     detSeqLog=getDetectedSeqs_c(W,int32(idx),detSeqLog,maxWLen);
 *     where W is the Cost submatrix of maxWLen length, idx is a vector 
 *     with several threshold values where matrix W must be cut at each 
 *     time for evaluation, and detSeqLog is the input logical vector to be
 *     modified.  
 */

#include "matrix.h"
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

int * detectSeq_c(float *w, int ns, int m, int maxWLen)
{
    int i,j,k;
    int * se;
    float last;
    bool foundPath;
    
    se = (int *)malloc(2*sizeof(int));
    last = w[m*(ns-1) +(m-1)];
       
    /* Obtain start-end from the optimal path */
    se[1] = m;
    j=m;
    k=j-1;    
    foundPath = false;
    /*printf(" %d %f ",j,w[j*(ns-1) +(k)]);*/
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
        printf("%d %f\n",j,w[j*(ns-1) +(k)]);
        printf("%f ",w[(j-1)*(ns-1) +(k)]);
        printf("%f ",w[(j-1)*(ns-1) +(k-1)]);
         */        
    }
    /*printf("%d %d\n",foundPath,j);*/
    
    /* Print se 
    for (i=0;i<2;i++)
    {
        printf("%d ",se[i]);
    }
    */
    
    return se;
}

bool * getDetectedSeqs_c(bool * detSeqLog, float *w, int nw, int m, int *idx, int k, int maxWLen)
{
    int i,j,l,n,a;
    int * se;
    float * w_cut;
    bool maxAllocated;
    
    w_cut=NULL;
    se=NULL;
    maxAllocated=false;
    
    for (l=0;l<k;l++)
    {
        n = idx[l]-1;
        
        /* This if statement may be uncommented to perform different 
         * thresholding evaluation. Then, lines 49,57,79 of the g.m file 
         * must be also uncommented. Otherwise very large (or untractable)
         * computation may happen when W becomes huge */
        /*if (detSeqLog[n]==0)
        {*/
            /* print matrix values */
            /*printf("Matrix cutting point: %d\nk=%d\n",n,k);*/
            if (n>=maxWLen) 
            {
                a=n-maxWLen;
                n=maxWLen;
            }
            else
            {
                maxAllocated=false;
                a=0;
            }
            if (!maxAllocated)
            {   
                w_cut = (float *)realloc(w_cut,nw*n*sizeof(float));
                if (n==maxWLen) maxAllocated=true;
            }
            /*printf("a: %d\nmaxAllocated: %d\n",a,maxAllocated);*/
            for (i=0;i<nw;i++)
            {            
                for(j=0;j<n;j++)
                {                    
                    w_cut[j*(nw) +i] = w[(j+a+1)*(nw) +i];                    
                }
            }
            /*
            for(i=0;i<nw;i++)
            {
                for(j=0;j<n;j++)
                {
                    printf("%f ",w_cut[j*(nw) +i]);
                }
                printf("\n");
            }
            */
            
            se = (int *)detectSeq_c(w_cut,nw,n,maxWLen);
            
            /*
            for (i=0;i<2;i++)
            {
                printf("%d ",se[i]);
            }*/

            for (j=a+se[0]-1;j<a+se[1];j++)
            {
                detSeqLog[j] = true;
                /*printf("%d,%d,%d,%d\n",se[0],se[1],j+1,detSeqLog[j]);*/
            }
        /*}*/
    }
    
    free(w_cut);
    free(se);
    
    return detSeqLog;
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    int nw,m,dims[2];
    int i,k,l,maxWLen;
    /*int j;*/
    float *w;
    int *idx;
    bool *dp,*detSeqLog;
	
    /*  check for proper number of arguments */
    if(nrhs!=4)
    {
        mexErrMsgIdAndTxt( "MATLAB:aligngesture_c:invalidNumInputs",
                "Three arguments required.");
    }
    if(nlhs!=1)
    {
        mexErrMsgIdAndTxt( "MATLAB:dtw_c:invalidNumOutputs",
                "aligngesture_c: One output required.");
    }
    
    /*  create a pointer to the input matrix w */
    w = mxGetPr(prhs[0]);
    
    /*  get the dimensions of the matrix input w */
    nw = mxGetM(prhs[0]);
    m = mxGetN(prhs[0]);
    
    /*  create a pointer to the input vector idx */
    idx = mxGetPr(prhs[1]);
    k = mxGetN(prhs[1]);
    
    /*  create a pointer to the input vector detSeqLog */
    detSeqLog = mxGetPr(prhs[2]);
    l = mxGetN(prhs[2]);
      
    /*  input scalar maxWLen */
    maxWLen = mxGetScalar(prhs[3]);
    /*printf("%d\n",maxWLen);*/
           
    /*
    for(i=0;i<nw;i++)
    {
        for(j=1;j<m;j++)
        {
            printf("%f  ",w[j*(nw) +i]);
        }
        printf("\n");
    }
    */
    
    /*  set the output pointer to the output vector 
    for (i=0;i<l;i++)
    {
        printf("%d ",detSeqLog[i]);
    }
    printf("\n");
    */
     
    /*
    for (i=0;i<k;i++)
    {
        printf("%d ",idx[i]);
    }
     */
    
    detSeqLog = getDetectedSeqs_c(detSeqLog,w,nw,m,idx,k,maxWLen);
    
    dims[0] = 1;
    dims[1] = l;
    plhs[0] = mxCreateLogicalArray(2,dims);
    
    /*  create a C pointer to a copy of the output matrix */
    dp = mxGetPr(plhs[0]);
    
    for (i=0;i<2;i++)
    {
        memcpy(&dp[i],&detSeqLog[i],l*sizeof(bool));
    }
    
    /*
    for (i=0;i<l;i++)
    {
        printf("%d ",dp[i]);
    }
    printf("\n");
    */
    
    /*free(detSeqLog);*/
    
    return;
    
}
