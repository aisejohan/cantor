/*
 *	emgcd.c
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
 #include <string.h>

#include "data.h"
#include "scalar.h"
#include "pol.h"


void make_emgcd(polynomial **emgcd)
{
	int i;

	*emgcd = (polynomial *)malloc(6*sizeof(polynomial));
	i = 0;
	while (i <= 5) {
		make_pol(&((*emgcd)[i]));
		i++;
	}
}


polynomial *easy_emgcd(polynomial *u)
{
	polynomial *emgcd;

	if (u[0]->degree <= u[1]->degree) {
		printf("Wrong degrees in emgcd routine.\n");
		exit(1);
	}

	if (2*u[1]->degree < u[0]->degree) {
		make_emgcd(&emgcd);
		copy_pol(emgcd[0], u[0]);
		copy_pol(emgcd[1], u[1]);
		emgcd[2]->degree = 0;
		emgcd[3]->degree = 0;
		emgcd[4]->degree = 0;
		emgcd[5]->degree = 0;
		emgcd[2]->coeffs[0] = 1;
		emgcd[3]->coeffs[0] = 0;
		emgcd[4]->coeffs[0] = 0;
		emgcd[5]->coeffs[0] = 1;
		return(emgcd);
	}

	emgcd = NULL;
	return(emgcd);
}

/* Assumes that deg(f) >= k > 0, and writes
 * b x^k + c = f with deg(c) < k. */
void split_polynomial(polynomial b, polynomial c, polynomial f, int k)
{
	int i;

	b->degree = f->degree - k;
	resize_pol(b, b->degree);
	i = b->degree;
	while (i >= 0) {
		b->coeffs[i] = f->coeffs[i + k];
		i--;
	}
	i = k - 1;
	while ((f->coeffs[i] == 0) && (i > 0)) i--;
	c->degree = i;
	resize_pol(c, c->degree);
	while (i >= 0) {
		c->coeffs[i] = f->coeffs[i];
		i--;
	}
}

/* Assumes that k >= 0.
 * f = b * x^k + c	 */
void put_back(polynomial f, polynomial b, polynomial c, int k)
{
	int i;

	if (b->degree == 0 && b->coeffs[0] == 0) {
		copy_pol(f, c);
		return;
	}
	f->degree = b->degree + k;
	resize_pol(f, f->degree);
	memset(f->coeffs, 0, (f->degree + 1)*sizeof(scalar));
	i = 0;
	while (i <= b->degree) {
		f->coeffs[i + k] = b->coeffs[i];
		i++;
	}
	pol_add(f, f, c);
}

/* Notation for the matrix ENGCD is:
 * 	emgcd[0]	emgcd[2]	emgcd[4]
 * 	emgcd[1]	emgcd[3]	emgcd[5]
 * for the vector U is
 * 	u[0]
 * 	u[1]
 * Assumes that the degree of u[0] > degree of u[1]. */
polynomial *do_emgcd(polynomial *u)
{
	int k, m;
	polynomial t;
	polynomial q;
	polynomial b[2];
	polynomial c[2];
	polynomial d[2];
	polynomial *emgcd1, *emgcd2;


	emgcd1 = easy_emgcd(u);
	if (emgcd1) {
		return(emgcd1);
	}

	/* Note that m > 0. */
	m = (u[0]->degree + 1)/2;

	make_pol(&b[0]);
	make_pol(&b[1]);
	make_pol(&c[0]);
	make_pol(&c[1]);
	split_polynomial(b[0], c[0], u[0], m);
	/* Note that deg(u[1]) >= m. */
	split_polynomial(b[1], c[1], u[1], m);

	/* Note that deg(b[0]) > deg(b[1]) */
	emgcd1 = do_emgcd(b);
	
	make_pol(&d[0]);
	make_pol(&d[1]);

	/*
	d[0] = emgcd1[0]*x^m + emgcd1[2]*c[0] + emgcd1[4]*c[1]
	d[1] = emgcd1[1]*x^m + emgcd1[3]*c[0] + emgcd1[5]*c[1]
	*/
	
	make_pol(&t);

	pol_mult(t, emgcd1[2], c[0]);
	put_back(d[0], emgcd1[0], t, m);
	pol_mult(t, emgcd1[4], c[1]);
	pol_add(d[0], d[0], t);

	pol_mult(t, emgcd1[3], c[0]);
	put_back(d[1], emgcd1[1], t, m);
	pol_mult(t, emgcd1[5], c[1]);
	pol_add(d[1], d[1], t);

	if (2*d[1]->degree < u[0]->degree) {
		copy_pol(emgcd1[0], d[0]);
		copy_pol(emgcd1[1], d[1]);
		free_pol(&b[0]);
		free_pol(&b[1]);
		free_pol(&c[0]);
		free_pol(&c[1]);
		free_pol(&d[0]);
		free_pol(&d[1]);
		free_pol(&t);
		return(emgcd1);
	}

	make_pol(&q);

	/* Note this q is the negative of the q from the paper. */
	qr_reduce(t, d[0], q, d[1]);
	k = 2*m - d[1]->degree;

	if (t->degree < k) {
		make_emgcd(&emgcd2);
		emgcd2[2]->coeffs[0] = 1;
		emgcd2[5]->coeffs[0] = 1;
		if (d[1]->degree < k) {
			copy_pol(c[0], d[1]);
		} else {
			split_polynomial(b[0], c[0], d[1], k);
			copy_pol(emgcd2[0], b[0]);
		}
		copy_pol(c[1], t);
	} else {
		split_polynomial(b[0], c[0], d[1], k);
		split_polynomial(b[1], c[1], t, k);

		emgcd2 = do_emgcd(b);
	}

	/*
	emgcd1[0] = emgcd2[0]*x^m + emgcd2[2]*c[0] + emgcd2[4]*c[1]
	emgcd1[1] = emgcd2[1]*x^m + emgcd2[3]*c[0] + emgcd2[5]*c[1]
	*/

	pol_mult(t, emgcd2[2], c[0]);
	put_back(emgcd1[0], emgcd2[0], t, k);
	pol_mult(t, emgcd2[4], c[1]);
	pol_add(emgcd1[0], emgcd1[0], t);

	pol_mult(t, emgcd2[3], c[0]);
	put_back(emgcd1[1], emgcd2[1], t, k);
	pol_mult(t, emgcd2[5], c[1]);
	pol_add(emgcd1[1], emgcd1[1], t);

	/*
	 * Matrix multiplication done by hand
	 * 0  1 			emgcd1[2] emgcd1[4]
	 * 		times
	 * 1  q				emgcd1[3] emgcd1[5]
	 * put into
	 * 	emgcd1[3] emgcd1[5]
	 * 	d[1]      c[1]
	 */

	pol_mult(t, q, emgcd1[3]);
	pol_add(d[1], emgcd1[2], t);
	pol_mult(t, q, emgcd1[5]);
	pol_add(c[1], emgcd1[4], t);

	/*
	 * Matrix multiplication done by hand
	 * emgcd2[2] emgcd2[4]		emgcd1[3] emgcd1[5]
	 * 			times
	 * emgcd2[3] emgcd2[5]		d[1]      c[1]
	 * put into
	 * 	emgcd1[2] emgcd1[4]
	 * 	emgcd1[3] emgcd1[5]
	 */

	pol_mult(t, emgcd2[2], emgcd1[3]);
	pol_mult(q, emgcd2[4], d[1]);
	pol_add(emgcd1[2], t, q);

	pol_mult(t, emgcd2[2], emgcd1[5]);
	pol_mult(q, emgcd2[4], c[1]);
	pol_add(emgcd1[4], t, q);

	pol_mult(t, emgcd2[3], emgcd1[3]);
	pol_mult(q, emgcd2[5], d[1]);
	pol_add(emgcd1[3], t, q);

	pol_mult(t, emgcd2[3], emgcd1[5]);
	pol_mult(q, emgcd2[5], c[1]);
	pol_add(emgcd1[5], t, q);

	free_pol(&b[0]);
	free_pol(&b[1]);
	free_pol(&c[0]);
	free_pol(&c[1]);
	free_pol(&d[0]);
	free_pol(&d[1]);
	free_pol(&t);
	free_pol(&q);
	free_pol(&emgcd2[0]);
	free_pol(&emgcd2[1]);
	free_pol(&emgcd2[2]);
	free_pol(&emgcd2[3]);
	free_pol(&emgcd2[4]);
	free_pol(&emgcd2[5]);
	return(emgcd1);
}

void variant_gcd(polynomial g, polynomial h, polynomial f)
{
	polynomial a[2];
	polynomial *emgcd;

	make_pol(&a[0]);
	make_pol(&a[1]);
	copy_pol(a[0], h);
	copy_pol(a[1], f);

	do {
		r_reduce(a[1], a[1], a[0]);
		emgcd = do_emgcd(a);
		free_pol(&a[1]);
		a[1] = emgcd[0];
		a[0] = emgcd[1];
		free_pol(&emgcd[2]);
		free_pol(&emgcd[3]);
		free_pol(&emgcd[4]);
		free_pol(&emgcd[5]);
	} while (a[0]->degree > 0 || a[0]->coeffs[0] != 0);

	copy_pol(g, a[1]);
}

