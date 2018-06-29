/* QR factorization with column pivoting, written by Yuanzhe Xi*/
/* Modified from Jiachen Xue's rrqr.c code for real matrices   */
/* [Q,R,rk] = rrqr(A,tol), A will be overwritten by Q after factorization */
/* In order to use this code, type mex rrqr.c in your matlab */
/* 2014 Feb version */

#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
double rrqr(double *A, double *R, int rows, int cols, int *p, double tol,double *nflops) {
    double rk;
    double iii;
    double *vn = (double *)calloc(cols, sizeof(double));
    int *piv = (int *)calloc(rows, sizeof(int));
    double vmax = 0.0;
    int imax = -1;
    int i, j, k, small, i1;
    
    *nflops = 0.0;
    for (i=0; i<cols; i++) {
        vn[i] = 0.0;
        for (j=0; j<rows; j++) {
            vn[i] += A[i*rows+j]*A[i*rows+j];
        }
        p[i] = i+1;
    }
    *nflops += (2*rows-1)*cols;
    
    small = (rows < cols ? rows:cols);
/*    small = cols;
    if (small > rows)
    {
        small = rows;
    }*/
    rk = small-1;
    iii=small-1;
    for (i=0; i<small; i++) {
        for (j=i; j<cols; j++) {
            if (fabs(vn[j]) > vmax) {
                vmax = vn[j];
                imax = j;
            }
        }
        vmax = 0.0;
        
        if (imax > i) {
            double tmp;
            /*exchange cols in A*/
            for (j=0; j<rows; j++) {
                tmp = A[i*rows + j];
                A[i*rows+j] = A[imax*rows + j];
                A[imax*rows + j] = tmp;
            }
            
            /*exchange cols in R*/
            for (j=0; j<i; j++) {
                tmp = R[j*cols + i];
                R[j*cols + i] = R[j*cols + imax];
                R[j*cols + imax] = tmp;
            }
            
            piv[i] = imax;
            vn[imax] = vn[i];
            
            i1 = p[i]; p[i] = p[imax]; p[imax] = i1;
        }
        
        
        /*R(i,i) = norm(A(1:m,i))*/
        for (j=0; j<rows; j++) {
            R[i*cols + i] += A[i*rows + j] * A[i*rows + j];
        }
        R[i*cols+i] = sqrt(R[i*cols + i]);
        

        if (R[i*cols + i] < 1e-20) {
            rk = (double)(i-1);
	        iii=i;
            break;
        }
        
        /*A(:,i) = A(:,i)/R(i,i)*/
        for (j=0; j<rows; j++) {
            A[i*rows + j] = A[i*rows + j]/R[i*cols+i];
        }
        
        /*vn(i) = R(i,i)*/
        vn[i] = R[i*cols + i];
        
        /*R(i,i+1:n) = A(1:m,i)'*A(1:m,i+1:n)*/
        for (j=i+1; j<cols; j++) {
            for (k=0; k<rows; k++) {
                R[i*cols+j] += A[i*rows + k] * A[j*rows + k];
            }
        }
        
        /*A(1:m,i+1:n) = A(1:m,i+1:n)-A(1:m,i)*R(i,i+1:n)*/
        for (j=i+1; j<cols; j++) {
            for (k=0; k<rows; k++) {
                A[j*rows+k] -= A[i*rows+k] * R[i*cols+j];
            }
        }
        
       if (R[i*cols + i]/R[0] < tol) {
            rk = (double)(i-1);
	        iii=i;
            break;
        }
        
        /*vn(i+1:n) = vn(i+1:n)-abs(R(i,i+1:n)).^2*/
        for (j=i+1; j<cols; j++) {
            vn[j] -= R[i*cols + j]*R[i*cols + j];
        }
        *nflops += (3+4*(cols-i-1))*rows+2*(cols-i-1);
    }
    
   //*nflops += 2*cols*rows-cols + 2*rk*rows + (2*cols-rk+1)*rk/2 + rk*rows + (4*rows-1)/2*(2*cols-1-rk)*rk;
    /* important */
    /*	if ((int)rk==rows) {
		iii=(int)rk-1;
	}
	else {
		iii=(int)rk;
	}
    */
    for (i=(int)iii; i>=0; i--) {
      /*        int mn = (int)rk * piv[i]-1;
        if (mn > rk) {
            mn = rk;
	    }  */
      /* */
      int mn=(int)rk;
      /**/
        if (piv[i] > 0) {
            
            double tmp;
            for (j=0; j<=(int)mn; j++) {
                tmp = R[j*cols + i];
                R[j*cols + i] = R[j*cols + piv[i]];
                R[j*cols + piv[i]] = tmp;
            }
        }
    }
    
    free(vn);
    free(piv);
    return rk;
}


double rrqrc(double *A, double *Ai, double *R, double *Ri,int rows, int cols, int *p, double tol,double *nflops) {
    double rk = (double)rows;
    double iii;
    double *vn = (double *)calloc(cols, sizeof(double));
    int *piv = (int *)calloc(rows, sizeof(int));
    double vmax = 0.0;
    int imax = -1;
    int i, j, k, small, i1;

    small = cols;
    if (small > rows)
    {
        small = rows;
    }
    rk = small-1;
    iii=small-1;

    *nflops = 0.0;
    for (i=0; i<cols; i++) {
        vn[i] = 0.0;
        for (j=0; j<rows; j++) {
            vn[i] += A[i*rows+j]*A[i*rows+j]+Ai[i*rows+j]*Ai[i*rows+j];
        }
        p[i] = i+1;
    }
    *nflops += (2*rows-1)*cols;
    
    for (i=0; i<small; i++) {
        for (j=i; j<cols; j++) {
            if (fabs(vn[j]) > vmax) {
                vmax = vn[j];
                imax = j;
            }
        }
        vmax = 0.0;
        
        if (imax > i) {
            double tmp;
            /*exchange cols in A*/
            for (j=0; j<rows; j++) {
                tmp = A[i*rows + j];
                A[i*rows+j] = A[imax*rows + j];
                A[imax*rows + j] = tmp;
                tmp = Ai[i*rows + j];
                Ai[i*rows+j] = Ai[imax*rows + j];
                Ai[imax*rows + j] = tmp;
            }
            
            /*exchange cols in R*/
            for (j=0; j<i; j++) {
                tmp = R[j*cols + i];
                R[j*cols + i] = R[j*cols + imax];
                R[j*cols + imax] = tmp;
                tmp = Ri[j*cols + i];
                Ri[j*cols + i] = Ri[j*cols + imax];
                Ri[j*cols + imax] = tmp;
            }
            
            piv[i] = imax;
            vn[imax] = vn[i];
            
            i1 = p[i]; p[i] = p[imax]; p[imax] = i1;
        }
        
        
        /*R(i,i) = norm(A(1:m,i))*/
        for (j=0; j<rows; j++) {
            R[i*cols + i] += A[i*rows + j] * A[i*rows + j]+Ai[i*rows + j] * Ai[i*rows + j];
        }
        R[i*cols+i] = sqrt(R[i*cols + i]);
        Ri[i*cols+i] = 0;
        
         if (R[i*cols + i] < 1e-20) {
            rk = (double)(i-1);
	    iii=i;
            break;
        }
        
        /*A(:,i) = A(:,i)/R(i,i)*/
        for (j=0; j<rows; j++) {
            A[i*rows + j] = A[i*rows + j]/R[i*cols+i];
            Ai[i*rows + j] = Ai[i*rows + j]/R[i*cols+i];
        }
        
        /*vn(i) = R(i,i)*/
        vn[i] = R[i*cols + i];
        
        /*R(i,i+1:n) = A(1:m,i)'*A(1:m,i+1:n)*/
        for (j=i+1; j<cols; j++) {
            for (k=0; k<rows; k++) {
                R[i*cols+j] += A[i*rows + k] * A[j*rows + k]+Ai[i*rows + k] * Ai[j*rows + k];
                Ri[i*cols+j] += A[i*rows + k] * Ai[j*rows + k]-Ai[i*rows + k] * A[j*rows + k];
            }
        }   
        
        /*A(1:m,i+1:n) = A(1:m,i+1:n)-A(1:m,i)*R(i,i+1:n)*/
        for (j=i+1; j<cols; j++) {
            for (k=0; k<rows; k++) {
	      /*  A[j*rows+k] -= A[i*rows+k] * R[i*cols+j]-Ai[i*rows+k] * Ri[i*cols+j];
		  Ai[j*rows+k] -= A[i*rows+k] * Ri[i*cols+j]+Ai[i*rows+k] * R[i*cols+j]; */
	      double temp,tempi;
	      temp=A[i*rows+k]*R[i*cols+j]-Ai[i*rows+k]*Ri[i*cols+j];
	      tempi=A[i*rows+k]*Ri[i*cols+j]+Ai[i*rows+k]*R[i*cols+j];
	      A[j*rows+k]-=temp;
	      Ai[j*rows+k]-=tempi;
            }
        }
        
        if (R[i*cols + i]/R[0] < tol) {
            rk = (double)(i-1);
	    iii=i;
            break;
        }
        
        /*vn(i+1:n) = vn(i+1:n)-abs(R(i,i+1:n)).^2*/
        for (j=i+1; j<cols; j++) {
            vn[j] -= R[i*cols + j]*R[i*cols + j]+Ri[i*cols + j]*Ri[i*cols + j];
        }
        
        *nflops += (3+4*(cols-i-1))*rows+2*(cols-i-1);
    }
    /* important */
    /*	if ((int)rk==rows) {
		iii=(int)rk-1;
	}
	else {
		iii=(int)rk;
		} */
    //*nflops += 2*cols*rows-cols + 2*rk*rows + (2*cols-rk+1)*rk/2 + rk*rows + (4*rows-1)/2*(2*cols-1-rk)*rk;
    for (i=(int)iii; i>=0; i--) {
      /*  int mn = (int)rk * piv[i];
        if (mn > rk) {
            mn = rk;
	    }*/
      int mn=(int)rk;
      
        if (piv[i] > 0) {
            double tmp;
            for (j=0; j<=mn; j++) {
                tmp = R[j*cols + i];
                R[j*cols + i] = R[j*cols + piv[i]];
                R[j*cols + piv[i]] = tmp;
                tmp = Ri[j*cols + i];
                Ri[j*cols + i] = Ri[j*cols + piv[i]];
                Ri[j*cols + piv[i]] = tmp;
            }
        }
    }
    
    
    free(vn);
    free(piv);
    return rk;
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if (mxIsComplex(prhs[0]))
    {
        double *A0, *A0i;
        A0 = mxGetPr(prhs[0]);
        A0i = mxGetPi(prhs[0]);
        int m = mxGetM(prhs[0]);
        int n = mxGetN(prhs[0]);
        
        double *A = (double *) calloc (m*n, sizeof(double));
        double *Ai = (double *) calloc (m*n, sizeof(double));
        int *p1 = (int *) calloc (n, sizeof(int));

        int i,j;
        for (i=0; i<n; i++) {
            memcpy(A+i*m, A0+i*m, sizeof(double)*m);
            memcpy(Ai+i*m, A0i+i*m, sizeof(double)*m);
        }
        
        double tol = mxGetScalar(prhs[1]);
        
        double *R_local = (double *) calloc (m*n, sizeof(double));
        double *R_locali = (double *) calloc (m*n, sizeof(double));
        double nflops = 0;
        double rk = rrqrc(A, Ai, R_local, R_locali, m, n, p1, tol,&nflops);
        plhs[0] = mxCreateDoubleMatrix(m,(int)rk+1,mxCOMPLEX);
        plhs[1] = mxCreateDoubleMatrix((int)rk+1,n, mxCOMPLEX);
        plhs[2] = mxCreateDoubleScalar(nflops);
        plhs[3] = mxCreateDoubleScalar(rk+1);
        plhs[4] = mxCreateDoubleMatrix(1, n, mxREAL);
       
        double *Q = mxGetPr(plhs[0]);
        double *Qi = mxGetPi(plhs[0]);
        double *R = mxGetPr(plhs[1]);
        double *Ri = mxGetPi(plhs[1]);
        double *p = mxGetPr(plhs[4]);

        for (i=0; i<(int)rk+1; i++) {
            memcpy(Q+i*m, A+i*m, sizeof(double)*m);
            memcpy(Qi+i*m, Ai+i*m, sizeof(double)*m);
        }
        
        for (i=0; i<n; i++) {
            for (j=0; j<(int)rk+1; j++) {
	              R[i*((int)rk+1) + j] = R_local[j*n+i];
	              Ri[i*((int)rk+1) + j] = R_locali[j*n+i];
            }
            p[i] = (int)p1[i];
        }
        
        free(p1);
        free(R_local);
        free(R_locali);
        free(A);
        free(Ai);
    }
    else
    {
        double *A0;
        A0 = mxGetPr(prhs[0]);
        int m = mxGetM(prhs[0]);
        int n = mxGetN(prhs[0]);
        double *A = (double *) calloc (m*n, sizeof(double));
        int *p1 = (int *) calloc (n, sizeof(int));
       
        int i,j;
        for (i=0; i<n; i++) {
            memcpy(A+i*m, A0+i*m, sizeof(double)*m);
        }
        
        double tol = mxGetScalar(prhs[1]);
        
        double *R_local = (double *) calloc (m*n, sizeof(double));
        double nflops = 0;
        double rk = rrqr(A, R_local, m, n, p1, tol,&nflops);
        plhs[0] = mxCreateDoubleMatrix(m,(int)rk+1,mxREAL);
        plhs[1] = mxCreateDoubleMatrix((int)rk+1,n, mxREAL);
        plhs[2] = mxCreateDoubleScalar(nflops);
        plhs[3] = mxCreateDoubleScalar(rk+1);
        plhs[4] = mxCreateDoubleMatrix(1, n, mxREAL);
        
        double *Q = mxGetPr(plhs[0]);
        double *R = mxGetPr(plhs[1]);
        double *p = mxGetPr(plhs[4]);

        for (i=0; i<(int)rk+1; i++) {
            memcpy(Q+i*m, A+i*m, sizeof(double)*m);
        }
        for (i=0; i<n; i++) {
            for (j=0; j<(int)rk+1; j++) {
	              R[i*((int)rk+1) + j] = R_local[j*n+i];
            }
            p[i] = (int)p1[i];
        }

        free(p1);
        free(R_local);
        free(A);   
    }
}