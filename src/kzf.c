#include <R.h>
#include "kz.h"

/* y is raw data
 * smooth it with kz
 * find adaptive filter
 */
void kzf(double *y, long *n_arg, long *q_arg, double *d, long *k_arg)
{
    int n, q, k;
    int i;

	n = *n_arg;	q = *q_arg; k = *k_arg;
	kz(y, n_arg, q_arg, k_arg);
	
	/* calculate d = |Z(t+q) - Z(t-q)| */
	for (i=0; i<q; i++) {d[i] = fabs(y[i+q] - y[0]);}
	for (i=q; i<n-q; i++) {d[i] = fabs(y[i+q] - y[i-q]);}
	for (i=n-q; i<n; i++) {d[i] = fabs(y[n-1] - y[i]);}
}    
