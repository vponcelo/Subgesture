/**
 * 2014 modified version by Víctor Ponce-López <vponcel@uoc.edu>
 * Originally from Copyright (C) 2013 Quan Wang <wangq10@rpi.edu> 
 */

/** 
 * This is the C/MEX code of dynamic time warping of two multi-dimensional
 * signals
 *
 * compile: 
 *     mex dtw_c.c
 *
 * usage:
 *     d=dtw_c(s,t,a)  or  d=dtw_c(s,t,a,w)
 *     where s is signal 1, t is signal 2, a is a flag to align (1) or 
 *     detect (0) and w is window parameter 
 */

#include "matrix.h"
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

double **dtw_c(double *s, double *t, int a, int w, int ns, int nt, int m, double *DS, double *KM, double *KT, int nDS, int mDS, int nKM, int mKM, int nKT, int mKT)
{
    double d=0;
    int sizediff=ns-nt>0 ? ns-nt : nt-ns;    
    double ** D;
    double minKM,minKT;
    int i,j;
    int j1,j2,jin;
    int c,locKM,locKT;
    double cost,temp;
    
    /* printf("ns=%d, nt=%d, m=%d, w=%d, s[0]=%f, t[0]=%f\n",ns,nt,m,w,s[0],t[0]); */
    
    
    if(w!=-1 && w<sizediff) w=sizediff; /* adapt window size */
    
    /* create D */
    /* D=(double **)calloc((ns),sizeof(double *));	*/
    D=(double **)malloc((ns)*sizeof(double *));     /* mxMalloc and mxCalloc have been also tested */
    for(i=0;i<ns;i++)
    {
        D[i]=(double *)malloc((nt)*sizeof(double));
        /* D[i]=(double *)calloc((nt),sizeof(double)); */
    }
    
    
    /* initialization */
    for(i=0;i<ns;i++)
    {
        for(j=0;j<nt;j++)
        {
            if (j==0)
            {
                if (a==1)
                {
                    D[i][j]=-1;
                } else {
                    D[i][j]=0;
                }
            } else {
                D[i][j]=-1;
            }
        }
    }
    D[0][0]=0;
    
    /* dynamic programming */
    for(i=1;i<ns;i++)
    {
        if(w==-1)
        {
            j1=1;
            j2=nt;
        }
        else
        {
            j1= i-w>1 ? i-w : 1;
            j2= i+w<nt ? i+w : nt;
        }
        for(j=j1;j<j2;j++)
        {
            if (nDS==-1)
            {
                double partCosts = 0;
                for (jin=0;jin<m;jin++)
                {   /* ns i nt refers to the size of D. s and t have size (ns-1,nt-1) */
                    partCosts += (s[jin*(ns-1) +(i-1)]-t[jin*(nt-1) +(j-1)])*(s[jin*(ns-1) +(i-1)]-t[jin*(nt-1) +(j-1)]);  /* Euclidean distance */
                    /* partCosts += fabs(s[jin*(ns-1) +(i-1)]-t[jin*(nt-1) +(j-1)]);   */ /* Absolute distance */
                    /*    printf("(%f - %f)^2+",s[jin*(ns-1) +(i-1)],t[jin*(nt-1) +(j-1)]);          */
                }
                /*printf("\n"); */

                cost = sqrt(partCosts);     /* Euclidean distance */
            }
            else
            {
                minKM = KM[(j-1)*nKM];                
                locKM = 0;
                minKT = KT[(i-1)*nKT];
                locKT = 0;
                for ( c = 1 ; c < nKM ; c++ ) 
                {
                    if ( KM[(j-1)*nKM+c] < minKM ) 
                    {
                       minKM = KM[(j-1)*nKM+c];
                       locKM = c;
                    }
                    if ( KT[(i-1)*nKT+c] < minKT ) 
                    {
                       minKT = KT[(i-1)*nKT+c];
                       locKT = c;
                    }
                    /* printf("%f , %f\n",minKM,minKT); */
                    /* printf("%d , %d\n",locKM,locKT); */
                }     
                
                cost = DS[locKM*nKM+locKT];
            }
            /* printf("%f , %d, %d\n",cost,i,j); */            
                
            temp=D[i-1][j];
            if(D[i][j-1]!=-1)
            {
                /*printf("(%f,%f):\t",temp,D[i][j-1]); */
                if(temp==-1 || D[i][j-1]<temp) temp=D[i][j-1];                
                /*printf("%f\t",temp); */
            }
            if(D[i-1][j-1]!=-1) 
            {
                /*printf("\n(%f,%f):\t",temp,D[i-1][j-1]); */
                if(temp==-1 || D[i-1][j-1]<temp) temp=D[i-1][j-1];
                /*printf("%f\t",temp); */
            }
            
            D[i][j]=cost+temp;
            /*printf("accumulated cost: %f\n",D[i][j]); */
        }
    }
    
    
    /* view matrix D */
    /*
    for(i=0;i<ns;i++)
    {
        for(j=0;j<nt;j++)
        {
            printf("%f  ",D[i][j]);
        }
        printf("\n");
    }    
    */
    
    /* free D */
    /*for(i=0;i<ns+1;i++)
    {
        free(D[i]);
    }
    // free(D);
              // Now, matrix D must be returned and cannot be freed here */
    
    /*d=D[ns][nt]; */
    
    return D;
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *s,*t,*DS,*KM,*KT;
    int a,w;    
    int ns,nt,m,nDS,mDS,nKM,mKM,nKT,mKT;
	int i,j;
    double *dp;
	double **D;
	
    /*  check for proper number of arguments */
    if(nrhs!=3 && nrhs!=4 && nrhs!=7)
    {
        mexErrMsgIdAndTxt( "MATLAB:dtw_c:invalidNumInputs",
                "Two or three inputs required.");
    }
    if(nlhs>1)
    {
        mexErrMsgIdAndTxt( "MATLAB:dtw_c:invalidNumOutputs",
                "dtw_c: One output required.");
    }
    
    /* check to make sure w is a scalar */
    if(nrhs==3)
    {
        w=-1;
        DS = NULL;
        nDS = -1;
        mDS = -1;
        KM = NULL;
        nKM = -1;
        mKM = -1;
        KT = NULL;        
        nKT = -1;
        mKT = -1;
    }
    else if(nrhs==4)
    {
        if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
                mxGetN(prhs[2])*mxGetM(prhs[2])!=1 )
        {
            mexErrMsgIdAndTxt( "MATLAB:dtw_c:wNotScalar",
                    "dtw_c: Input w must be a scalar.");
        }
        
        /*  get the scalar input w */
        w = (int) mxGetScalar(prhs[3]);
    }
    else if(nrhs==7)
    {
        DS = mxGetPr(prhs[4]);
        nDS = mxGetM(prhs[4]);
        mDS = mxGetN(prhs[4]);
        if (nDS!=mDS)
        {
            mexErrMsgTxt("dtw_c: Input DS must be of dimension nxn.");
        }
        /*printf("%d , %d \n ",nDS,mDS); */
        KM = mxGetPr(prhs[5]);
        nKM = mxGetM(prhs[5]);
        if (nDS!=nKM)
        {
            mexErrMsgTxt("dtw_c: Input DS and KM must have the same number of rows.");
        }
        /*printf("%d , %d \n ",nDS,nKM); */
        mKM = mxGetN(prhs[5]);
        KT = mxGetPr(prhs[6]);
        nKT = mxGetM(prhs[6]);        
        if (nDS!=nKT)
        {
            mexErrMsgTxt("dtw_c: Input DS and KT must have the same number of rows.");
        }
        /*printf("%d , %d \n ",nDS,nKT); */
        mKT = mxGetN(prhs[6]);
        
    }
    /*
    for(i=0;i<nKM;i++)
    {
        for(j=0;j<mKM;j++)
        {
            printf("%f  ",KM[j*(nKM) +i]);
        }
        printf("\n");
    }
    */
    
    /*  create a pointer to the input matrix s */
    s = mxGetPr(prhs[0]);
    
    /*  create a pointer to the input matrix t */
    t = mxGetPr(prhs[1]);
    
    /*  create a pointer to the input flag a */
    a = (int)mxGetScalar(prhs[2]);
    
    /*  get the dimensions of the matrix input s */
    ns = mxGetM(prhs[0]);
    if (mKT!=-1 && ns!=mKT)
    {
        mexErrMsgTxt("dtw_c: Number of rows of s must have the same number of columns of KT.");
    }
    /*printf("%d , %d \n ",ns,mKT); */
    ns++;
    m = mxGetN(prhs[0]);
      
    /*
    for(i=0;i<ns;i++)
    {
        for(j=0;j<m;j++)
        {
            printf("%f  ",s[j*(ns) +i]);
        }
        printf("\n");
    }
    */
    
    /*  get the length of the matrix input t */
    nt = mxGetM(prhs[1]);
    if (mKM!=-1 && nt!=mKM)
    {
        mexErrMsgTxt("dtw_c: Number of rows of t must have the same number of columns of KM.");
    }
    /*printf("%d , %d \n ",nt,mKM); */
    nt++;
    /*
    for(i=0;i<nt;i++)
    {
        for(j=0;j<m;j++)
        {
            printf("%f  ",t[j*(nt) +i]);
        }
        printf("\n");
    }
    */
    
    /*  set the output pointer to the output matrix */
    
	D = dtw_c(s,t,a,w,ns,nt,m,DS,KM,KT,nDS,mDS,nKM,mKM,nKT,mKT);
    
    plhs[0] = mxCreateDoubleMatrix( nt, ns, mxREAL);
	
	/*  create a C pointer to a copy of the output matrix */
    dp = mxGetPr(plhs[0]);	
    
    for(i=0; i<ns; i++)
    {
        memcpy(&dp[(nt)*i], D[i], sizeof(double)*(nt));
    }    
    
    /* These are the possible releases. Also tested with mxFree() when allocating using mex calls    */
    for(i=0;i<ns;i++)
    {
        free(D[i]);     
    }
    free(D);        
    
    return;
    
}
