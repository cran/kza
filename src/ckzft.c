/*
    kz.c
    Copyright (C) 2013 Brian D. Close

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#include <complex.h>
#include <math.h>

#define notNaN(x)   ((x)==(x))
#define isNaN(x)  (!((x)==(x)))
#define MIN(y,x) ((x)<(y) && (x)==(x) ? (x) : (y))
#define MAX(y,x) ((x)>(y) && (x)==(x) ? (x) : (y))
#define SQR(x) ((x)*(x))

double setPrecision(double x, int prec) {
    return 
    	ceil( x * pow(10,(double)prec) - .4999999999999)
    	/ pow(10,(double)prec);
}

double complex cexp(double complex z)
{
	double complex w;
	double r, x, y;
	x = creal(z);
	y = cimag(z);
	r = exp(x);
	w = r * cos(y) + r * sin(y) * I;
	return w;
}

// inY is indexed starting from 1 to end of time sequence

void ckzft(double *outR, double *outImg, double *inX, const int *vectorSize, double *inY, double *inWindow, double *inScale, double *inLambda)
{
	double *r;
	double *img;
	double *x;
	double *t;
	double m;
	double lambda;
	int xLength;
	int start;
	int i,j, k;
	int n;
	double complex c;
	double pi;
	double complex theta;
	double NaN = (0.0/0.0);
	int here=0;
	
	pi = M_PI;
	x  = inX;
	xLength = *vectorSize;
	t = inY;
	m = *inWindow;
	lambda = *inLambda;
	r = outR; img = outImg;

	// i is center k is index in time
	for (i=0, start=0; i<=t[xLength-1]; i++, r++, img++) {
		for(j=-(m-1)/2, n=0, k=start; j<=(m-1)/2; j++) {
			if (i+j < 0) {continue;}
			if (i+j > t[xLength-1]) {break;}
			if (t[k] < i-(m-1)/2) {start++; k++; --j; continue;}
			if (t[k]==i+j)
			{
				if (notNaN(x[k])) {
					theta = 0-(2*pi*lambda*j)*I;
					c = x[k] * cexp(theta);
					if (isNaN(*r)) { *r=0; *img=0; }
					*r += creal(c);
					*img += cimag(c);
					n++;
				}
				k++;
			}
		}
		if (n) {
			here++;
			*r /= n;
			*r = setPrecision(*r,14);
			*img /= n;
		} else {
			*r = NaN;
			*img = NaN;
		}
	}
}

