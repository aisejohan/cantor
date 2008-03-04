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

/* Assumes f is square free, and h is the (p^d)th power of x modulo f.
 * Returns the first integer e >= d such that f has factors of degree e,
 * and sets g equal to the product of those factors, f equal to
 * f/g and h equal to x^(p^e) modulo f/g.
int next_degree(polynomial g, polynomial f, int d, polynomial h)
{
}
*/
