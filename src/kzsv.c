/*
    kzsv.c
    kz sample variance
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
#include <errno.h>
#include <math.h>
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

/*  sum((y[i] - mean(y))^2)
 * y is filtered data from kza
 * n_arg is size of y
 * q is half window size
 * d is difference values
 * sigma^2 = sum((y[i] - mean(y))^2)/(qh + qt + 1)
 */
void kzsv(double *y, long *n_arg, long *q_arg, double *d)
{
	double *v;
	long n, q;
	double m;
	long t;
	int i, cnt;
	long qh,qt;
	double *dprime;
	double dmax;

	n = *n_arg;
	q = *q_arg;
	qt = qh = q;

	v = calloc(sizeof(double),n);
	if (v == NULL) error("malloc failed: errno\n", errno);

	/* filter */
	dprime = malloc( n*sizeof(double) );
	if (dprime == NULL) error("malloc failed: errno\n", errno);
	for(i=0; i<n-1; i++) dprime[i] = d[i+1]-d[i];
	dprime[n-1] = 0;
	
   	dmax = maximum(d, n);

    for (t=0; t<q; t++) {
        v[t] = NaN;
    }
	for (t=q; t<n-q; t++) {
	    /* set head and tail size of filter */
    	if (dprime[t] < 0) {
	    	qh = q;
   		    qt = (int) MAX(floor(q*adaptive(d[t], dmax)+0.5), 1);
    	} else if (dprime[t] > 0) {
	    	qh = (int) MAX(floor(q*adaptive(d[t], dmax)+0.5), 1);
		    qt = q;
	    } else { /* dprime[t] = 0 */
		    qh = (int) MAX(floor(q*adaptive(d[t], dmax)+0.5), 1);
		    qt = (int) MAX(floor(q*adaptive(d[t], dmax)+0.5), 1);
	    }
	    m = avg(y+t-qt,qt+qh+1);
		for(cnt=0, i=t-qt;i<=t+qh;i++) {
    		if (!isnan(y[i])) {
                v[t] += (y[i] - m)*(y[i] - m);
                cnt++;
            }
		}
		if (cnt>0) v[t] /= cnt;
	}

	for (t=n-q; t<n; t++) {
	    v[t] = NaN;
	}
	 memcpy(y, v, n*sizeof(double));

	free(v);
	free(dprime);
}
