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

void make_pol(polynomial f)
{
	int i;

	f->degree = 0;
	f->length = min_length;
	f->coeffs = (mscalar *) malloc((f->length+1)*sizeof(mscalar));
	if (!f->coeffs) {
		perror("Malloc failed in make_pol!");
		exit(1);
	}
	for(i = 0; i <= f->length; i++) make_scalar(f->coeffs[i]);
}

void free_pol(polynomial f)
{
	int i;
	
	for(i = 0; i <= f->length; i++) free_scalar(f->coeffs[i]);
	free(f->coeffs);
}

/* This resizes the storage space of the polynomial, but
 * does not allow you to change the polynomial. */
void resize_pol(polynomial f, unsigned int new_length)
{
	int i;
	void *ptr;
	
	if (new_length > f->length) {
		if (new_length < 2*f->length) new_length = 2*f->length;
		ptr = realloc(f->coeffs, (new_length+1)*sizeof(mscalar));
		if (!ptr) {
			perror("Realloc failed in resize_pol!");
			exit(1);
		}
		f->coeffs = (mscalar *) ptr;
		for(i = f->length + 1; i <= new_length; i++)
			make_scalar(f->coeffs[i]);
		f->length = new_length;
	} else {
		if (2*new_length > f->length) return;
		if (new_length < min_length) new_length = min_length;
		if (new_length < f->degree) {
			printf("Truncating too much in resize_pol! STOP.");
			exit(1);
		}
		for(i = f->length; i > new_length; i--)
			free_scalar(f->coeffs[i]);
		ptr = realloc(f->coeffs, (new_length+1)*sizeof(mscalar));
		if (!ptr) {
			perror("Realloc failed in resize_pol!");
			exit(1);
		}
		f->coeffs = (mscalar *) ptr;
		f->length = new_length;
	}
}

/* Copies a polynomial. Will truncate for you.		*/
void copy_pol(polynomial f, polynomial g)
{
	int i;
	
	g->degree = f->degree;
	resize_pol(g, f->degree);
	for(i = 0; i <= f->degree; i++) sc_copy(f->coeffs[i], g->coeffs[i]);
}

/* Same, but copies negative. 				*/
void negate_pol(polynomial f, polynomial g)
{
	int i;

	g->degree = f->degree;
	resize_pol(g, f->degree);
	for(i = 0; i <= f->degree; i++) sc_negate(f->coeffs[i], g->coeffs[i]);
}

/* Prints a polynomial. 				*/
void print_pol(polynomial f)
{
	int i,y=0;

	if(!sc_is_zero(f->coeffs[0])) {
			y=1;
			printmscalar(f->coeffs[0]);
	} else if(f->degree == 0) printf("0");

	for(i=1;i<=f->degree;i++) {
		if(!sc_is_zero(f->coeffs[i])) {
			if(y) printf(" + ");
			y=1;
			printmscalar(f->coeffs[i]);
			printf("*x^%d", i);
		}
	}
	printf("\n");
}

/* Makes a random monic polynomial. */
void random_pol(polynomial f, unsigned int d)
{
	int i,c;

	f->degree = d;
	resize_pol(f, d);

	for(i = 0; i < d; i++) {
		c = rand() % p;
		ito_sc(c, f->coeffs[i]);
	}

	sc_one(f->coeffs[d]);
}

/* This should also work if some of these are the same. */
void pol_add(polynomial f, polynomial g, polynomial h)
{
	mscalar c;
	int i,d;

	make_scalar(c);
	
	if (f->degree > g->degree) {
		d = f->degree;
		resize_pol(h, d);
		i=f->degree;
		while (i > g->degree) {
			sc_copy(f->coeffs[i], h->coeffs[i]);
			i--;
		}
	} else if (f->degree < g->degree) {
		d = g->degree;
		resize_pol(h, d);
		i=g->degree;
		while (i > f->degree) {
			sc_copy(g->coeffs[i], h->coeffs[i]);
			i--;
		}
	} else {
		/* degrees are equal, maybe cancellation */
		i=f->degree;
		sc_add(f->coeffs[i], g->coeffs[i], c);
		while ((sc_is_zero(c)) && (i > 0)) {
			i--;
			sc_add(f->coeffs[i], g->coeffs[i], c);
		}
		/* Note that with this convention the zero polynomial
		 * has degree 0.*/
		d = i;
		h -> degree = d;
		resize_pol(h, d);
	}
	while (i >= 0) {
		sc_add(f->coeffs[i], g->coeffs[i], h->coeffs[i]);
		i--;
	}

	h->degree = d;

	free_scalar(c);
}

/* Replaces g by the product of a and f (possibly zero).	*/
void times_int(int a, polynomial f, polynomial g)
{
	int i;
	mscalar c;

	make_scalar(c);

	i = f->degree;
	sc_imult(a, f->coeffs[i], c);
	while (sc_is_zero(c) && i > 0) {
		i--;
		sc_imult(a, f->coeffs[i], c);
	}
	g->degree = i;
	resize_pol(g, g->degree);
	while (i >= 0) {
		sc_imult(a, f->coeffs[i], g->coeffs[i]);
		i--;
	}

	free_scalar(c);
}

/* Replaces g by the product of a and f (possibly zero).	*/
void times_scalar(mscalar a, polynomial f, polynomial g)
{
	int i;
	mscalar c;

	make_scalar(c);
	
	i = f->degree;
	sc_mult(a, f->coeffs[i], c);
	while (sc_is_zero(c) && i > 0) {
		i--;
		sc_mult(a, f->coeffs[i], c);
	}
	g->degree = i;
	resize_pol(g, g->degree);
	while (i >= 0) {
		sc_mult(a, f->coeffs[i], g->coeffs[i]);
		i--;
	}

	free_scalar(c);
}


/* Multiply f and g and put the result in h. 
 * This is optimized a little bit so uses
 * how scalars work. This is bad... */
void pol_mult(polynomial f, polynomial g, polynomial h)
{
	int i,j;
	polynomial t;

	make_pol(t);
	resize_pol(t, f->degree + g->degree);

	for(i = 0; i <= f->degree; i++) {
		for(j = 0; j <= g->degree; j++) {
			mpz_addmul(t->coeffs[i+j], f->coeffs[i],
				g->coeffs[j]);
		}
	}

	i = f->degree + g->degree;
	while (sc_is_zero(t->coeffs[i]) && i > 0) i--;
	t->degree = i;

	h->degree = t->degree;
	resize_pol(h, h->degree);
	for(i = 0; i<= t->degree; i++)
		mpz_mod(h->coeffs[i], t->coeffs[i], modulus);

	free_pol(t);
}
