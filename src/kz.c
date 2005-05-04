/*
    kz.c
    Copyright (C) 2005 Brian D. Close

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

#include <R.h> 
#include <Rdefines.h>
#include "kz.h"


static double mavg2d(SEXP v, int row, int col, int w)
{
	double s=0.00;
	int i, j, z=0;
	int startrow, startcol, endrow, endcol;
	int nr=nrows(v);

    if (isMatrix(v)) {
        
        startrow = (row-w>0 ? row-w : 0);
        startcol = (col-w>0 ? col-w : 0);

        endrow = (row+w<nr ? row+w+1 : nrows(v));
        endcol = (col+w<ncols(v) ? col+w+1 : ncols(v));
        
	    for(i=startrow, z=0, s=0.00; i<endrow; i++) {
	        for(j=startcol; j<endcol; j++) {
	    	    if (R_FINITE(REAL(v)[i+nr*j])) { 
	    		    z++;
	    		    s += REAL(v)[i+nr*j];
	    	    }
	    	}
	    }
	} else error("Input is not a vector or Matrix.");
	if (z == 0) return R_NaN;
	return s/z;
}

static double mavg1d(SEXP v, int col, int w)
{
	double s=0.00;
	int i, z=0;
	int startcol, endcol;

	if (isVector(v)) {
        startcol = (col-w>0 ? col-w : 0);
        endcol = (col+w<LENGTH(v) ? col+w+1 : LENGTH(v));
        
	    for(i=startcol, z=0, s=0.00; i<endcol; i++) {
    	    if (R_FINITE(REAL(v)[i])) {
    		    z++;
    		    s += REAL(v)[i];
	    	}
	    }
	} else error("Input is not a vector or Matrix.");
	if (z == 0) return R_NaN;
	return s/z;
}

SEXP kz1d(SEXP x, SEXP window, SEXP iterations)
{
	int p;
	int i, k;
	int m;
	SEXP ans, tmp, dim;
	int n;

	m = (2 * INTEGER_VALUE(window)) + 1; 
	p = (m-1)/2;
	
	dim = GET_DIM(x);
	n = LENGTH(x);
	
	PROTECT(tmp = allocVector(REALSXP, LENGTH(x)));
	PROTECT(ans = allocVector(REALSXP, LENGTH(x)));
    copyVector(tmp, x);

	for(k=0; k<INTEGER_VALUE(iterations); k++) {
        for(i=0;i<LENGTH(x);i++) {
            REAL(ans)[i] = mavg1d(tmp, i, p);
        }
        /* copyMatrix (destination, source, byrow) */
   	    copyVector(tmp, ans); 
	}
	UNPROTECT(2);
	return ans;
}

SEXP kz2d(SEXP x, SEXP window, SEXP iterations)
{
	int p;
	int i, j, k;
	int m;
	SEXP ans, tmp, dim;
	int nr, nc;

	m = (2 * INTEGER_VALUE(window)) + 1; 
	p = (m-1)/2;
	
	dim = GET_DIM(x);
	nr = INTEGER(dim)[0]; nc = INTEGER(dim)[1];
	
	PROTECT(ans = allocMatrix(REALSXP, nr, nc));
	PROTECT(tmp = allocMatrix(REALSXP, nr, nc));
    copyMatrix(tmp, x, 1);

	for(k=0; k<INTEGER_VALUE(iterations); k++) {
	    for(i=0;i<nr; i++) {
	        for(j=0;j<nc;j++) {
	            REAL(ans)[i+nr*j] = mavg2d(tmp, i, j, p);
	        }
        }
        /* copyMatrix (destination, source, byrow) */
   	    copyMatrix(tmp, ans, 1); 
	}
	UNPROTECT(2);
	return ans;
}
