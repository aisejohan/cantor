/*
 *	modp_pol.c
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
#include "modp_scalar.h"
#include "pol.h"
#include "utils.h"

unsigned int do_mod(polynomial f, unsigned int **a)
{
	int i,d;
	struct exponents e;

	e = take_exponents(f);

	if (e.val > 0) 
		d = 0;
	else 
		d = e.piv;

	*a = (unsigned int *) malloc((d+1)*sizeof(unsigned int));
	if (!*a) {
		perror("Malloc failed in do_mod.");
		exit(1);
	}

	for(i = 0; i <= d; i++) (*a)[i] = mpz_fdiv_ui(f->coeffs[i], p);

	return(d);
}

void lift(unsigned int first, unsigned int last, unsigned int *a, polynomial f)
{
	int i;

	f->degree = last - first ;
	resize_pol(f, f->degree);

	for(i = first; i <= last; i++) ito_sc(a[i], f->coeffs[i - first]);

	i = last;
	while (!a[i] && i > first) i--;
	f->degree = i - first;
}

int euclidean_mod_p(polynomial f, polynomial g, polynomial a, polynomial b)
{
	int i;
	unsigned int c;
	unsigned int d1,d2,dA,dB;
	unsigned int *q1,*q2,*A,*B,*T;

	d1 = do_mod(f, &q1);
	d2 = do_mod(g, &q2);

	/* Deal with di = 0 first, so later we have di > 0. */
	if (d1 == 0 || d2 == 0) {
		if ((d1 == 0 && d2 > 0 && q1[0] == 0) ||
			(d1 > 0 && d2 == 0 && q2[0] == 0) ||
			(d1 == 0 && d2 ==0 && q1[0] == 0 && q2[0] == 0)) {
				printf("Impossible in euclidean_mod_p.\n");
				exit(1);
		}
		if (d1 == 0 && q1[0] != 0) {
			q1[0] = INV[q1[0]];
			lift(0, 0, q1, a);
			q2[0] = 0;
			lift(0, 0, q2, b);
			free(q2);
			free(q1);
			return(1);
		}
		q2[0] = INV[q2[0]];
		lift(0, 0, q2, b);
		q1[0] = 0;
		lift(0, 0, q1, a);
		free(q2);
		free(q1);
		return(1);
	}

	dA = d1 + d2 + d2;
	A = (unsigned int *) malloc((dA+1)*sizeof(unsigned int));
	if (!A) {
		perror("Not again! euclidean_mod_p.\n");
		exit(1);
	}
	for(i = 0; i <= dA; i++) A[i] = 0;

	dB = d2 + d1 + d2;
	B = (unsigned int *) malloc((dB+1)*sizeof(unsigned int));
	if (!B) {
		perror("Not again! euclidean_mod_p.\n");
		exit(1);
	}
	for(i = 0; i <= dB; i++) B[i] = 0;
	
	for(i = 0; i <= d1; i++) A[i+d1+d2] = q1[i];
	A[d1] = 1;

	for(i = 0; i <= d2; i++) B[i+d1+d2] = q2[i];
	B[0] = 1;

	while (dB >= d1 + d2) {
		while (dA >= dB) {
			c = MUL[A[dA]][NEG[INV[B[dB]]]];
			for(i = 0; i <= dB - 1; i++)
				A[i + dA - dB] = 
					ADD[A[i + dA - dB]][MUL[c][B[i]]];
			dA = dA - 1;
			while (A[dA] == 0 && dA > 0) dA--;
		}
		/* Swap. */
		T = A;
		A = B;
		B = T;
		c = dA;
		dA = dB;
		dB = c;
	}

	c = INV[A[dA]];
	for(i = 0; i <= dA - 1; i++) A[i] = MUL[c][A[i]];

	lift(d1, d1 + d2 - 1, A, a);
	lift(0, d1 - 1, A, b);

	if (dA != d1 + d2) return(0);
	return(1);

	free(A);
	free(B);
	free(q2);
	free(q1);
}
