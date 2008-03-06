/*
 *	pol.c
 *
 * 	Copyright 2006 Johan de Jong
 *
 *	This file is part of Frobenius
 *
 *	Frobenius is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; either version 2 of the License, or
 *	(at your option) any later version.
 *
 *	Frobenius is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with Frobenius; if not, write to the Free Software Foundation, 
 *	Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *									*/

#include <stdio.h>
#include <stdlib.h>

#include "data.h"
#include "scalar.h"
#include "pol.h"

/* Assumes g is reduced modulo f. */
void prime_power(polynomial h, polynomial g, polynomial f)
{
	int i;
	polynomial tmp;

#ifdef KIJKEN
	test_pol(h);
	test_pol(g);
	test_pol(f);
#endif

	make_pol(&tmp);
	copy_pol(tmp, g);

	h->degree = 0;
	h->coeffs[0] = 1;

	i = 1;
	do {
		if (prime & i) {
			pol_mult(h, h, tmp);
			r_reduce(h, h, f);
		}
		pol_mult(tmp, tmp, tmp);
		r_reduce(tmp, tmp, f);
		i = 2*i;
	} while (i <= prime);

	free_pol(&tmp);
}

/* Assumes f is square free, has no divisors of degree <= d and
 * h is the (p^d)th power of x modulo f.
 * Returns the first integer e >= d such that f has factors of degree e,
 * and sets g equal to the product of those factors, f equal to
 * f/g and h equal to x^(p^e) modulo f/g.
 * Returns 0 if f irreducible, in which case f is unchanged but h and
 * g do not have meaning. */
int next_degree(polynomial g, polynomial f, unsigned int d, polynomial h)
{
	polynomial tmp1, tmp2;

#ifdef KIJKEN
	if ((f == g) || (f == h) || (g == h)) {
		printf("Should not have same pols in next_degree.");
		exit(1);
	}
#endif

	d++;
	while (2*d <= f->degree) {
		prime_power(h, h, f);
		if (h->degree == 1 && h->coeffs[1] == 1 && h->coeffs[0] == 0) {
			copy_pol(g, f);
			f->degree = 0;
			f->coeffs[0] = 1;
			resize_pol(f, f->degree);
			h->degree = 0;
			h->coeffs[0] = 0;
			resize_pol(h, h->degree);
			return(d);
		}
		h->coeffs[1] = (h->coeffs[1] + prime - 1) % prime;
		gcd(g, h, f);
		h->coeffs[1] = (h->coeffs[1] + 1) % prime;
		if (g->degree > 0) {
			make_pol(&tmp1);
			make_pol(&tmp2);
			qr_reduce(tmp1, f, tmp2, g);
			copy_pol(f, tmp2);
			r_reduce(h, h, f);
			free_pol(&tmp1);
			free_pol(&tmp2);
			return(d);
		}
		d++;
	}
	return(0);
}

void print_degrees(polynomial f)
{
	int e;
	polynomial h;
	polynomial g;

	make_pol(&h);
	make_pol(&g);

	h->degree = 1;
	h->coeffs[1] = 1;

	e = 0;
	do {
		e = next_degree(g, f, e, h);
		if (e) {
			printf("We have factors of degree %d"
			" with multiplicity %d.\n", e, g->degree/e);
			print_pol(g);
		} else {
			printf("Irreducible factor of degree %d.\n",
				f->degree);
			print_pol(f);
		}
	} while (e);
	free_pol(&h);
	free_pol(&g);
}
