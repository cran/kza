/*
    kza.c
    Copyright (C) 2003  Brian D. Close

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
#include "kz.h"

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif


static double maximum(double *v, int length) 
{
	double m;
	int i;
	m = v[0];

	for(i=0; i<length; i++) m = MAX(v[i], m);

	return m;
}

/* adaptive curvature function */
static double adaptive(double d, double m)
{
	return( 1 - (d/m) );
}

/* v = data vector 
 * n = size of data vector 
 * q = half window size 
 * k = iterations 
 * build filter, apply to raw data 
 */
void kza(double *v, long *n_arg, long *q_arg, long *k_arg)
{
	int i,j;
	long n,q,k;
	double *d, *dprime;
	double m;
	int t, qh, qt;
	double *y;

	n = *n_arg;	q = *q_arg; k = *k_arg;
	y = malloc( n*sizeof(double) );
	if (y == NULL) error("malloc failed: errno\n", errno);
	memcpy(y, v, n*sizeof(double));
	kz(y, n_arg, q_arg, k_arg);

	/* calculate d = |Z(t+q) - Z(t-q)| */
	d = malloc( n*sizeof(double) );
	if (d == NULL) error("malloc failed: errno\n", errno);
	for (i=0; i<q; i++) d[i] = fabs(y[i+q] - y[0]);
	for (i=q; i<n-q; i++) {d[i] = fabs(y[i+q] - y[i-q]);}
	for (i=n-q; i<n; i++) {d[i] = fabs(y[n-1] - y[i]);}

	/* d'(t) = d(t+1)-d(t) */
	dprime = malloc( n*sizeof(double) );
	if (dprime == NULL) error("malloc failed: errno\n", errno);
	for(i=0; i<n-1; i++) dprime[i] = d[i+1]-d[i];
	dprime[n-1] = 0;

    m = maximum(d, n);
	memset(y, '\0', n*sizeof(double));

    for(i=0; i<k; i++) {
        for(t=0; t<q; t++) { y[t] = NaN; }
    	for (t=q; t<n-q; t++) {
	    	if (dprime[t] < 0) {
		    	qh = q;
    		    qt = (int) floor(q*adaptive(d[t], m)+0.5);
	    	} else if (dprime[t] > 0) {
		    	qh = (int) floor(q*adaptive(d[t], m)+0.5);
			    qt = q;
		    } else { /* dprime[t] = 0 */
			    qh = (int) floor(q*adaptive(d[t], m)+0.5);
			    qt = (int) floor(q*adaptive(d[t], m)+0.5);
		    }
		    for (j=t-qt; j<=t+qh; j++) { y[t] += v[j]; }
		    y[t] /= (qh+qt+1);
	    }
	    for (t=n-q; t<n; t++) {
	        y[t] = NaN;
        }
        
	    memcpy(v,y, n*sizeof(double));
	}

	/* clean up */
	free(d);
	free(dprime);
	free(y);
}
