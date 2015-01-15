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

int * detectSeq_c(float *w, int ns, int m, int maxWLen)
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
    /*printf("%f ",w[j*(ns-1) +(k)]);*/
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
