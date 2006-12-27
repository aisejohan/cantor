/*
 *	reduce.c
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
#include "modp_pol.h"

/* Reduces f by g. This assumes that
 * 	g is prepared, and
 * 	f is not equal g.	 */
void reduce_by_g(polynomial f, polynomial g)
{
	int i;
	int d;
	mscalar u,c,t;

	make_scalar(u);
	make_scalar(c);
	make_scalar(t);

	d = g->degree;
	sc_copy(g->coeffs[d],u);
	if (valuation(u) != 0) {
		printf("In reduce_by_g, g not prepared. STOP.");
		exit(1);
	}
	sc_inv(u, u);
	sc_negate(u, u);

	while (f->degree >= d) {
		sc_mult(u, f->coeffs[f->degree], c);
		for(i = d - 1; i >= 0; i--) {
			sc_mult(c, g->coeffs[i], t);
			sc_add(t, f->coeffs[f->degree - d + i],
				f->coeffs[f->degree - d + i]);
		}
		f->degree = f->degree - 1;
	}
	resize_pol(f,f->degree);

	free_scalar(u);
	free_scalar(c);
	free_scalar(t);
}

/* Reduces f by g. This assumes that
 * 	g is prepared,
 * 	all three polynomials are different.
 * Puts the quotient in q so that fold = fnew + q*g. */
void reduce_by_g_q(polynomial f, polynomial g, polynomial q)
{
	int i;
	int d;
	mscalar u,c,t;

	make_scalar(u);
	make_scalar(c);
	make_scalar(t);

	d = g->degree;
	sc_copy(g->coeffs[d],u);
	if (valuation(u) != 0) {
		printf("In reduce_by_g, g not prepared. STOP.");
		exit(1);
	}
	sc_inv(u, u);

	if (f->degree < d) {
		q->degree = 0;
		sc_zero(q->coeffs[0]);
		resize_pol(q, q->degree);
		goto einde;
	}

	q->degree = f->degree - g->degree;
	resize_pol(q, q->degree);
	while (f->degree >= d) {
		sc_mult(u, f->coeffs[f->degree], q->coeffs[f->degree - d]);
		for(i = d - 1; i >= 0; i--) {
			sc_mult(q->coeffs[f->degree - d], g->coeffs[i], t);
			sc_sub(f->coeffs[f->degree - d + i], t,
				f->coeffs[f->degree - d + i]);
		}
		f->degree = f->degree - 1;
	}
	resize_pol(f,f->degree);

einde:
	free_scalar(u);
	free_scalar(c);
	free_scalar(t);
}

/* Assumes g is prepared. Reduces f by g and sets fi to the inverse. 
 * Also assumes f is invertible mod g. */
void invert_unit_mod_g(polynomial f, polynomial g, polynomial fi)
{
	int i;
	polynomial t,s;

	reduce_by_g(f, g);

	make_pol(t);
	make_pol(s);

	/* First do euclid mod p. */
	i = euclidean_mod_p(f, g, fi, s);
	if (!i || g->degree == 0) {
		printf("Not invertible! invert_unit_mod_g\n");
		exit(1);
	}

	pol_mult(f, fi, t);
	reduce_by_g(t, g);
	while (t->degree > 0) {
		negate_pol(t, s);
		sc_ui_add(2, s->coeffs[0], s->coeffs[0]);
		pol_mult(t, s, t);
		reduce_by_g(t, g);
		pol_mult(s, fi, fi);
		reduce_by_g(fi, g);
	}

	sc_one(s->coeffs[0]);
	sc_sub(t->coeffs[0], s->coeffs[0], s->coeffs[0]);
	if (!sc_is_zero(s->coeffs[0])) {
		sc_inv(t->coeffs[0], t->coeffs[0]);
		times_scalar(t->coeffs[0], fi, fi);
	}

	free_pol(s);
	free_pol(t);
}
