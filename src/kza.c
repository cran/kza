/*
    kza.c
    Copyright (C) 2005  Brian D. Close

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <errno.h>
#include <R.h>
#include <Rdefines.h>
#include "kz.h"

typedef struct {
    int x,y;
} point;

/* upper left is 0,0 */
typedef struct{
    point upperLeft, lowerRight;
} boundingBox ;

/* adaptive curvature function */
static double adaptive(double d, double m)
{
	return( 1 - (d/m) );
}

static double mavg1d(SEXP x, int a, int b)
{
	double s=0.00;
	long i, z;
  
	for(i=a, z=0; i<b; i++)
	{
		if (R_FINITE(REAL(x)[i])) {
			z++;
			s += REAL(x)[i];
		}
	}
	if (z==0) return R_NaN;
	return s/z;
}

SEXP kza1d(SEXP v, SEXP y, SEXP window, SEXP iterations, SEXP minimum_window_length, SEXP tolerance)
{
	int i;
	int n,q;
	SEXP d, dprime;
	double m;
	long t, qh, qt;
	int min_window_length;
	SEXP tmp, ans;
	double eps;

    eps = REAL(tolerance)[0];
	q = INTEGER_VALUE(window);
	min_window_length = INTEGER_VALUE(minimum_window_length);
    n = LENGTH(y);
    
	PROTECT(d = allocVector(REALSXP, n));
	PROTECT(dprime = allocVector(REALSXP, n));
    R_differenced(y, d, dprime, q);
 
    m = R_maximum(d);

    PROTECT(tmp = allocVector(REALSXP, n));
    copyVector(tmp, v);
    
    PROTECT(ans = allocVector(REALSXP, n));
    for(i=0; i<INTEGER_VALUE(iterations); i++) {
    	for (t=0; t<n; t++) {
		    if (fabs(REAL(dprime)[t]) < eps) { /* dprime[t] = 0 */
			    qh = (int) floor(q*adaptive(REAL(d)[t], m));
			    qt = (int) floor(q*adaptive(REAL(d)[t], m));
		    } else if (REAL(dprime)[t] < 0) {
		    	qh = q;
    		    qt = (int) floor(q*adaptive(REAL(d)[t], m));
	    	} else {
		    	qh = (int) floor(q*adaptive(REAL(d)[t], m));
			    qt = q;
			}
			qt = ((qt) < min_window_length) ? min_window_length : qt;
			qh = ((qh) < min_window_length) ? min_window_length : qh;
	
	        /* check bounds */
        	qh = (qh > n-t-1) ? n-t-1 : qh; /* head past end of series */
            qt = (qt > t) ? t : qt;  	        		
   		    REAL(ans)[t] = mavg1d(tmp, t-qt, t+qh+1); 
	    }
	    copyVector(tmp, ans);
	}
	UNPROTECT(4);
	return ans;
}

static double mavg2d(SEXP v, boundingBox b)
{
	double s;
	int i, j, z;
	int nr;

    nr = nrows(v);
    s=0.0;
	for(i=b.upperLeft.x, z=0; i<=b.lowerRight.x; i++) {
	    for(j=b.upperLeft.y; j<=b.lowerRight.y; j++) {
		    if (R_FINITE(REAL(v)[i+nr*j])) {
			    z++;
			    s += REAL(v)[i+nr*j];
		    }
		}
	}
	if (z==0) return R_NaN;
	return s/z;
}

SEXP kza2d(SEXP v, SEXP kz, SEXP window, SEXP iterations, SEXP minimum_window_length, SEXP tolerance)
{
    int q, min_window_length;
	int i, j, k;
	SEXP dx, dy, dprimex, dprimey;
	double mx, my;
	int qx_head, qx_tail, qy_head, qy_tail;
	int col_head, col_tail, row_head, row_tail;
	SEXP dim, tmp, ans;
	int nr, nc;
	boundingBox b;
	double epsilon;

    epsilon = REAL(tolerance)[0];
	q = INTEGER_VALUE(window);
	min_window_length = INTEGER_VALUE(minimum_window_length);
	
    dim = GET_DIM(v);
	nr = INTEGER(dim)[0]; nc = INTEGER(dim)[1];
	
	PROTECT(tmp = allocMatrix(REALSXP, nr, nc));
    copyMatrix(tmp, v, 1);

	/* calculate d1 = |Z(x+q,y) - Z(x-q,y)| , d2 = |Z(x,y+q) - Z(x,y-q)| */
	PROTECT(dx = allocMatrix(REALSXP, nr, nc));
	PROTECT(dy = allocMatrix(REALSXP, nr, nc));

    /*  origin is upper left (0,0)
        positive directions are down (y) and right (x) 
     */
	for(i=0; i<nr; i++) { /* ith row (y) */
	    for(j=0;j<nc;j++) { /* jth column (x) */
	        row_tail = (i-q<0 ? 0 : i-q);
	        col_tail = (j-q<0 ? 0 : j-q);
	        row_head = (i+q>=nr ? nr-1 : i+q);
	        col_head = (j+q>=nc ? nc-1 : j+q);
	        REAL(dx)[i+nr*j] = fabs(REAL(kz)[i+nr*col_head] - REAL(kz)[i+nr*col_tail]);
	        REAL(dy)[i+nr*j] = fabs(REAL(kz)[row_head+nr*j] - REAL(kz)[row_tail+nr*j]);
	    }
	}

	/* d'(t) = d(t+1)-d(t) */
	PROTECT(dprimex = allocMatrix(REALSXP, nr, nc));
	PROTECT(dprimey = allocMatrix(REALSXP, nr, nc));

	for(i=0; i<nr; i++) { /* ith row  (y)*/
	    for(j=0;j<nc-1;j++) { /* jth column (x) */
	        REAL(dprimex)[i+j*nr] = REAL(dx)[i+nr*(j+1)]-REAL(dx)[i+nr*j];
	    }
	    REAL(dprimex)[i+nr*j] = REAL(dprimex)[i+nr*(j-1)];
	}
    for(j=0;j<nc;j++) {
    	for(i=0; i<nr-1; i++) {
	        REAL(dprimey)[i+j*nr] = REAL(dy)[i+1+nr*j]-REAL(dy)[i+nr*j];
	    }
	    REAL(dprimey)[nr-1+nr*j] = REAL(dprimex)[nr-2+nr*j];
    }
    
    mx = R_maximum(dx);
    my = R_maximum(dy);

	PROTECT(ans = allocMatrix(REALSXP, nr, nc));

    for(k=0; k<INTEGER_VALUE(iterations); k++) {

        /* determine bounding box for i,j coordinate */
    	for (i=0; i<nr; i++) {
    	    for(j=0; j<nc; j++) { /* coordinate (x (col),y (row) is [i+nr*j] */
		        if (fabs(REAL(dprimex)[i+nr*j]) < epsilon) { /* dprimex[i,j] = 0 */
			        qx_head = (int) floor(q*adaptive(REAL(dx)[i+j*nr], mx));
			        qx_tail = (int) floor(q*adaptive(REAL(dx)[i+j*nr], mx));
		        } else if (REAL(dprimex)[i+nr*j] < 0) {
		        	qx_head = q;
    		        qx_tail = (int) floor(q*adaptive(REAL(dx)[i+j*nr], mx));
	    	    } else {
		        	qx_head = (int) floor(q*adaptive(REAL(dx)[i+j*nr], mx));
			        qx_tail = q;
			    }

		        if (fabs(REAL(dprimey)[i+nr*j]) < epsilon) { /* dprimey[i,j] = 0 */
			        qy_head = (int) floor(q*adaptive(REAL(dy)[i+j*nr], my));
			        qy_tail = (int) floor(q*adaptive(REAL(dy)[i+j*nr], my));
		        } else if (REAL(dprimey)[i+nr*j] < 0) {
		        	qy_head = q;
    		        qy_tail = (int) floor(q*adaptive(REAL(dy)[i+j*nr], my));
	    	    } else {
		        	qy_head = (int) floor(q*adaptive(REAL(dy)[i+j*nr], my));
			        qy_tail = q;
			    }

			    qx_tail = ((qx_tail) < min_window_length) ? min_window_length : qx_tail;
			    qx_head = ((qx_head) < min_window_length) ? min_window_length : qx_head;
			    qy_tail = ((qy_tail) < min_window_length) ? min_window_length : qy_tail;
			    qy_head = ((qy_head) < min_window_length) ? min_window_length : qy_head;
	            
	            /* check bounds */
        	    qy_head = (qy_head > nc) ? nc : qy_head; /* head past end of series */
                qy_tail = (qy_tail > i+nr*j) ? i+nr*j : qy_tail;

                b.upperLeft.x = (j-qx_tail<0 ? 0 : j-qx_tail);
                b.upperLeft.y = (i-qy_tail<0 ? 0 : i-qy_tail);
                b.lowerRight.x = (j+qx_head>=nc ? nc-1 : j+qx_head);
                b.lowerRight.y = (i+qy_head>=nr ? nr-1 : i+qy_head);
                if ((b.lowerRight.y - b.upperLeft.y) + (b.lowerRight.x-b.upperLeft.x))
                    REAL(ans)[i+j*nr]=mavg2d(tmp, b);
                else /* zero sized box - just use avg of one item */
   		            REAL(ans)[i+j*nr]=REAL(tmp)[i+j*nr];
   		    }
	    }
	}

    UNPROTECT(6);
	return ans;
}

/* x : matrix including class(its) */
SEXP kza(SEXP x, SEXP kz, SEXP window, SEXP iterations, SEXP min_window, SEXP tol)
{
    R_len_t nx, ny;
    SEXP ans=R_NilValue, dim, class;
       
    dim = getAttrib(x, R_DimSymbol);
    class = getAttrib(x, R_ClassSymbol);
       
    nx = INTEGER(dim)[0]; ny = INTEGER(dim)[1];

    if (isMatrix(x)) {
        PROTECT(ans = kza2d(x, kz, window, iterations, min_window, tol)); 
        UNPROTECT(1); 
    }
    else if (isVector(x)) {
        PROTECT(ans = kza1d(x, kz, window, iterations, min_window, tol)); 
        UNPROTECT(1); 
    }
    else error("Input is not a Matrix.");
    
    /* its package 
    if (strcmp("its", CHAR(STRING_ELT(class, 0)) == 0)  {
    
    }
    */
    
    return(ans);
}

/* x : matrix including class(its) */
SEXP kz(SEXP x, SEXP window, SEXP iterations)
{
    SEXP ans=R_NilValue, dim, class;
       
    dim = GET_DIM(x);
       
    if (isMatrix(x)) {
        PROTECT(ans = kz2d(x, window, iterations)); 
        UNPROTECT(1);
    }
    else if (isVector(x)) {
        PROTECT(ans = kz1d(x, window, iterations)); 
        UNPROTECT(1);
    }
    else error("Too many dimensions.");
    
    
    /* its package 
    class = getAttrib(x, R_ClassSymbol);
    if (strcmp("its", CHAR(STRING_ELT(class, 0)) == 0)  {
    
    }
    */
    return(ans);
}

SEXP kzsv(SEXP kza_data, SEXP kz_data, SEXP window, SEXP minimum_window_length, SEXP iterations, SEXP tolerance)
{
    SEXP ans=R_NilValue;
    
    PROTECT(ans = R_kzsv(kza_data, kz_data, window, minimum_window_length, tolerance));
    UNPROTECT(1);
    return(ans);
}
