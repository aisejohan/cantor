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
#include "utils.h"

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

/* Can improve by introducing squaring function. */
scalar *frobs_mod_f(polynomial f)
{
	int d, e, i;
	polynomial h;
	polynomial g;
	scalar *frobs;

	d = f->degree;
	frobs = (scalar *) calloc(d*d, sizeof(scalar));

	make_pol(&h);
	h->degree = 0;
	h->coeffs[0] = 1;

	frobs[0] = 1;

	make_pol(&g);
	g->degree = prime;
	resize_pol(g, g->degree);
	i = 0;
	while (i < prime) {
		g->coeffs[i] = 0;
		i++;
	}
	g->coeffs[prime] = 1;
	r_reduce(g, g, f);

	e = 1;
	while (e < d) {
		pol_mult(h, h, g);
		r_reduce(h, h, f);
		i = 0;
		while (i <= h->degree) {
			frobs[e*d + i] = h->coeffs[i];
			i++;
		}
		e++;
	}

	free_pol(&h);
	free_pol(&g);
	return(frobs);
}

void fast_prime_power_mod(polynomial h, int d, scalar *frobs)
{
	int i, j;
	scalar c;
	unsigned int *tmp;

	tmp = calloc(d, sizeof(unsigned int));
	i = 0;
	while (i <= h->degree) {
		j = 0;
		while (j < d) {
			tmp[j] += sc_mul(h->coeffs[i], frobs[i*d + j]);
			j++;
		}
		i++;
	}
	i = d - 1;
	while (((c = tmp[i] % prime) == 0) && (i > 0))  i--;
	h->degree = i;
	resize_pol(h, h->degree);
	j = 0;
	while (j < i) {
		h->coeffs[j] = tmp[j] % prime;
		j++;
	}
	h->coeffs[h->degree] = c;
	free(tmp);
}

int *easy_case(polynomial f)
{
	int *list;

	if (f->degree == 0) {
		list = (int *)malloc(1*sizeof(int));
		list[0] = 0;
		return(list);
	}	

	if (f->degree == 1) {
		list = (int *)malloc(2*sizeof(int));
		list[0] = 1;
		list[1] = 1;
		return(list);
	}

	list = NULL;
	return(list);
}


/* f has to be square free and prime to x. */
int *list_degrees_sq_x_free(polynomial f)
{
	int a, d, e, i, nr;
	int *list;
	scalar *frobs;
	polynomial h;
	polynomial g;
	polynomial tmp1;
	polynomial tmp2;

	list = easy_case(f);
	if (list) return(list);

	make_pol(&h);
	h->degree = 1;
	h->coeffs[1] = 1;

	make_pol(&g);
	make_pol(&tmp1);
	make_pol(&tmp2);

	list = (int *)malloc(f->degree*sizeof(int));

	a = 0;
	e = 1;
	d = f->degree;
	frobs = frobs_mod_f(f);
	while (f->degree >= 2*e) {
		fast_prime_power_mod(h, d, frobs);
		h->coeffs[1] = (h->coeffs[1] + prime - 1) % prime;
		gcd(g, h, f);
		h->coeffs[1] = (h->coeffs[1] + 1) % prime;
		if (g->degree > 0) {
			nr = g->degree / e;
			i = 0;
			while (i < nr) {
				a++;
				list[a] = e;
				i++;
			}
			qr_reduce(tmp1, f, tmp2, g);
			copy_pol(f, tmp2);
		}
		e++;
	}
	if (f->degree) {
		a++;
		list[a] = f->degree;
	}
	list[0] = a;

	free_pol(&h);
	free_pol(&g);
	free_pol(&tmp1);
	free_pol(&tmp2);
	free(frobs);

	return(list);
}

int *list_degrees(polynomial f)
{
	int a, i, j;
	int *t1, *t2;
	int *list;
	polynomial h;
	polynomial g;
	polynomial g1;
	polynomial f1;
	polynomial tmp1;
	polynomial tmp2;

	list = easy_case(f);
	if (list) return(list);

	make_pol(&h);

	i = 0;
	while (f->coeffs[i] == 0) i++;
	if (i) {
		t1 = malloc((i+1)*sizeof(int));
		a = 0;
		while (a < i) {
			a++;
			t1[a] = 1;
		}
		t1[0] = a;
		h->degree = f->degree - i;
		resize_pol(h, h->degree);
		j = 0;
		while (j <= h->degree) {
			h->coeffs[j] = f->coeffs[j + i];
			j++;
		}
		t2 = list_degrees(h);
		list = merge(t1, t2);
		free(t1);
		free(t2);
		free_pol(&h);
		return(list);
	}	

	deriv(h, f);

	/* The case where f is a pth power. */
	if ((h->degree == 0) && ( h->coeffs[0] == 0)) {
		h->degree = f->degree / prime;
		resize_pol(h, h->degree);
		i = 0;
		while (i <= h->degree) {
			h->coeffs[i] = f->coeffs[prime*i];
			i++;
		}
		t1 = list_degrees(h);
		list = (int *)malloc((1 + t1[0]*prime)*sizeof(int));
		a = 0;
		i = 1;
		while (i <= t1[0]) {
			j = 0;
			while (j < prime) {
				a++;
				list[a] = t1[i];
				j++;
			}
			i++;
		}
		list[0] = a;
		free_pol(&h);
		free(t1);
		return(list);
	}

	make_pol(&g);
	gcd(g, h, f);

	if (g->degree > 0) {
		make_pol(&tmp1);
		make_pol(&tmp2);
		make_pol(&g1);
		make_pol(&f1);
		qr_reduce(tmp1, f, tmp2, g);
		gcd(g1, tmp2, g);
		qr_reduce(tmp1, tmp2, f1, g1);
		list = list_degrees_sq_x_free(f1);
		t2 = list_degrees_sq_x_free(g1);
		t1 = merge(t2, list);
		free(list);
		free(t2);
		t2 = list_degrees(g);
		list = merge(t1, t2);
		free(t1);
		free(t2);
		free_pol(&h);
		free_pol(&g);
		free_pol(&tmp1);
		free_pol(&tmp2);
		free_pol(&g1);
		free_pol(&f1);
		return(list);
	}

	free_pol(&h);
	free_pol(&g);
	return(list_degrees_sq_x_free(f));
}


		
		

		
