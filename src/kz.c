/*
    kz.c
    Copyright (C) 2003 Brian D. Close

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
#include <errno.h>
#include <memory.h>
#include <float.h>
#include <math.h>
#include "kz.h"
#include <R.h>

double avg(double x[], int m)
{
	double s=0.00;
	int i, z;
  
	if (m==0) return NaN;

	for(i=0, z=0; i<m; i++)
	{
		if (!isnan(x[i])) {
			z++;
			s += x[i];
		}
	}
	if (z==0) return NaN;
	return s/z;
}


/* n = size of input vector */
/* m = size of filter */
/* k = iterations */
void kz(double *data_vector, long *n_arg, long *m_arg, long *k_arg)
{
	int p;
	int i, j, k;
	double *x, *y;
	long m, n, iterations;

	m = (2 * *m_arg) + 1; n = *n_arg; iterations = *k_arg;

	x = malloc(m*sizeof(double));
	if (x == NULL) error("malloc failed: errno\n", errno);
	y = malloc( n*sizeof(double) );
	if (y == NULL) error("malloc failed: errno\n", errno);
	for(i=0; i<n; i++) y[i] = data_vector[i];

	p = (m-1)/2;

	for(k=0; k<iterations; k++) 
	{
		memset(x, '\0', m*sizeof(double));
		for (i=0; i<p; i++)
			x[i] = NaN;
		for(i=p; i<m; i++)
			x[i] = y[i-p];

		for(i=0; i<n; i++)
		{
			y[i] = avg(x, m);

			/* setup x for next iteration */
			for(j=0; j<m-1; j++)
				x[j] = x[j+1];
			if (i+p+1<n)
				x[m-1] = y[i+p+1];
			else
				x[m-1] = NaN;
		}
	}
	for(i=0; i<n; i++) data_vector[i] = y[i];
		
	/* clean up */
	free(x);
	free(y);
}


long min_distance(double *x, long n)
{
    long d=0;
    long i;
    
    for(i=1; i<n; i++) {
        d = (x[i]-x[i-1]>d ? x[i]-x[i-1] : d);
    }
    return d;
}
