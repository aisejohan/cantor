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

void make_pol(polynomial *f)
{
	*f = (polynomial) malloc(sizeof(struct pol));
	(*f)->degree = 0;
	(*f)->length = min_length;
	(*f)->coeffs = (scalar *) malloc(((*f)->length+1)*sizeof(scalar));
	if (!(*f)->coeffs) {
		perror("Malloc failed in make_pol!");
		exit(1);
	}
}

void free_pol(polynomial *f)
{
	free((*f)->coeffs);
	free(*f);
}

/* This resizes the storage space of the polynomial, but
 * does not allow you to change the polynomial. */
void resize_pol(polynomial f, unsigned int new_length)
{
	void *ptr;
	
	if (new_length > f->length) {
		if (new_length < 2*f->length) new_length = 2*f->length;
		ptr = realloc(f->coeffs, (new_length+1)*sizeof(scalar));
		if (!ptr) {
			perror("Realloc failed in resize_pol!");
			exit(1);
		}
		f->coeffs = (scalar *) ptr;
		f->length = new_length;
	} else {
		if (new_length < min_length) new_length = min_length;
		if (2*new_length > f->length) return;
		if (new_length < f->degree) return;
		ptr = realloc(f->coeffs, (new_length+1)*sizeof(scalar));
		if (!ptr) {
			perror("Realloc failed in resize_pol!");
			exit(1);
		}
		f->coeffs = (scalar *) ptr;
		f->length = new_length;
	}
}

/* Sets g equal to a copy of f. Will truncate for you.		*/
void copy_pol(polynomial g, polynomial f)
{
	int i;
	
	g->degree = f->degree;
	resize_pol(g, f->degree);
	for(i = 0; i <= f->degree; i++) g->coeffs[i] = f->coeffs[i];
}

/* Prints a polynomial. 				*/
void print_pol(polynomial f)
{
	int i,y=0;

	if (f->coeffs[0] % prime) {
			y=1;
			print_scalar(f->coeffs[0]);
	} else if(f->degree == 0) printf("0");

	for(i=1; i<=f->degree; i++) {
		if (f->coeffs[i] % prime) {
			if (y) printf(" + ");
			y=1;
			print_scalar(f->coeffs[i]);
			printf("*x^%d", i);
		}
	}
	printf("\n");
}

/* Makes a random monic polynomial. */
void random_pol(polynomial f, unsigned int d)
{
	int i;

	f->degree = d;
	resize_pol(f, d);

	for(i = 0; i < d; i++) f->coeffs[i] = rand() % prime;

	f->coeffs[d] = 1;
}

/* This should also work if some of these are the same. */
void pol_add(polynomial h, polynomial g, polynomial f)
{
	int i, d;
	scalar c;
	
	if (f->degree > g->degree) {
		d = f->degree;
		i = f->degree;
		resize_pol(h, i);
		while (i > g->degree) {
			h->coeffs[i] = f->coeffs[i];
			i--;
		}
	} else if (f->degree < g->degree) {
		d = g->degree;
		i = g->degree;
		resize_pol(h, i);
		while (i > f->degree) {
			h->coeffs[i] = g->coeffs[i];
			i--;
		}
	} else {
		i = f->degree;
		while ((c = (f->coeffs[i] + g->coeffs[i]) % prime) == 0 &&
			i > 0) i--;
		/* Note that with this convention the zero polynomial
		 * has degree 0.*/
		d = i;
		resize_pol(h, i);
		h->coeffs[i] = c;
		i--;
	}
	while (i >= 0) {
		h->coeffs[i] = (g->coeffs[i] + f->coeffs[i]) % prime;
		i--;
	}
	h->degree = d;
}

void times_int(polynomial f, int a, polynomial g)
{
	int i;
	scalar c;

	i = g->degree;
	while ((c = (a * g->coeffs[i]) % prime) == 0 && i > 0) i--;
	f->degree = i;
	resize_pol(f, f->degree);
	f->coeffs[i] = c;
	i--;
	while (i >= 0) {
		f->coeffs[i] = (a * g->coeffs[i]) % prime;
		i--;
	}
}

void times_scalar(polynomial f, scalar a, int power, polynomial g)
{
	int i;
	scalar c;

	i = g->degree;
	while ((c = (a * g->coeffs[i]) % prime) == 0 && i > 0) i--;
	f->degree = i + power;
	resize_pol(f, f->degree);
	f->coeffs[i + power] = c;
	i--;
	while (i >= 0) {
		f->coeffs[i + power] = (a * g->coeffs[i]) % prime;
		i--;
	}
	i = power - 1;
	while (i >= 0) {
		f->coeffs[i] = 0;
		i--;
	}
}

void pol_mult(polynomial h, polynomial g, polynomial f)
{
	int i,j;
	polynomial tmp;

	make_pol(&tmp);
	tmp->degree = f->degree + g->degree;
	resize_pol(tmp, f->degree + g->degree);

	i = f->degree + g->degree;
	while (i >= 0) {
		tmp->coeffs[i] = 0;
		i--;
	}

	i = 0;
	while (i <= f->degree) {
		j = 0;
		while (j <= g->degree) {
			tmp->coeffs[i+j] = (tmp->coeffs[i+j] +
				f->coeffs[i] * g->coeffs[j]) % prime;
			j++;
		}
		i++;
	}

	i = f->degree + g->degree;
	while (tmp->coeffs[i] == 0 && i > 0) i--;
	h->degree = i;
	resize_pol(h, h->degree);
	while (i >= 0) {
		h->coeffs[i] = tmp->coeffs[i];
		i--;
	}
	free_pol(&tmp);
}

/* r = g + qf */
void qr_reduce(polynomial r, polynomial g, polynomial q, polynomial f)
{
	int i,j;
	scalar c;
	polynomial tmp;

	if ((q == f) || (r == f) || (g == f) || (q == r)) {
		printf("reduce: equal pols not allowed.\n");
		exit(1);
	}

	if (r != g) copy_pol(r, g);

	i = r->degree - f->degree;
	if (i >= 0) {
		q->degree = i;
		resize_pol(q, q->degree);
	} else {
		q->degree = 0;
		resize_pol(q, 0);
		q->coeffs[0] = 0;
	}
	make_pol(&tmp);
	c = sc_inv(f->coeffs[f->degree]);
	c = (prime - c) % prime;
	while (i >= 0) {
		q->coeffs[i] = (c * r->coeffs[r->degree]) % prime;
		times_scalar(tmp, q->coeffs[i], i, f);
		pol_add(r, r, tmp);
		j = r->degree - f->degree;
		i--;
		while ((i > j) && (i >= 0)) {
			q->coeffs[i] = 0;
			i--;
		}
	}
	free_pol(&tmp);
}

void r_reduce(polynomial r, polynomial g, polynomial f)
{
	int i;
	scalar c, q;
	polynomial tmp;

	if ((r == f) || (g == f)) {
		printf("reduce: equal pols not allowed.\n");
		exit(1);
	}

	if (r != g) copy_pol(r, g);

	i = r->degree - f->degree;
	make_pol(&tmp);
	c = sc_inv(f->coeffs[f->degree]);
	c = (prime - c) % prime;
	while (i >= 0) {
		q = (c * r->coeffs[r->degree]) % prime;
		times_scalar(tmp, q, i, f);
		pol_add(r, r, tmp);
		i = r->degree - f->degree;
	}
	free_pol(&tmp);
}
