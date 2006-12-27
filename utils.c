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

/* Computes the exponents of f. */
struct exponents take_exponents(polynomial f)
{
	int i;
	unsigned int val;
	struct exponents e;
	
	e.val = r;
	e.piv = 0;

	for(i = f->degree; i >= 0; i--) {
		if (!sc_is_zero(f->coeffs[i])) {
			val = valuation(f->coeffs[i]);
			if (val < e.val) {
				e.val = val;
				e.piv = i;
				if (val == 0) return(e);
			}
		}
	}

	return(e);
}

/* Puts the derivative of f in g. */
void deriv(polynomial f, polynomial g)
{
	int i;

	resize_pol(g, f->degree);
	for(i = 1; i <= f->degree; i++)
		sc_imult(i, f->coeffs[i], g->coeffs[i-1]);
	
	i = f->degree - 1;
	while (sc_is_zero(g->coeffs[i]) && i > 0) i--;
	g->degree = i;

}

/* Puts the inverse of u into ui. Assumptions:
 * 	u is a unit, and
 * 	u and ui are not the same. 	*/
void invert_unit(polynomial u, polynomial ui)
{
	polynomial t,s;
	struct exponents e;

	e = take_exponents(u);
	if (e.val != 0 || e.piv != 0) {
		printf("Polynomial not invertible in invert_unit.\n");
		exit(1);
	}

	make_pol(t);
	make_pol(s);

	ui->degree = 0;
	sc_inv(u->coeffs[0], ui->coeffs[0]);
	pol_mult(u, ui, t);
	while (t->degree > 0) {
		negate_pol(t, s);
		sc_ui_add(2, s->coeffs[0], s->coeffs[0]);
		pol_mult(t, s, t);
		pol_mult(s, ui, ui);
	}

	sc_one(s->coeffs[0]);
	sc_sub(t->coeffs[0], s->coeffs[0], s->coeffs[0]);
	if (!sc_is_zero(s->coeffs[0])) {
		sc_inv(t->coeffs[0], t->coeffs[0]);
		times_scalar(t->coeffs[0], ui, ui);
	}

	free_pol(s);
	free_pol(t);
}

/* Substitutes the variable x in f by g and puts the result in h. 
 * Assumes that all three pols are different. */
void substitute(polynomial f, polynomial g, polynomial h)
{
	int i;

	h->degree = 0;
	sc_copy(f->coeffs[f->degree], h->coeffs[0]);

	if (f->degree == 0) return;

	for(i = f->degree - 1; i >= 0; i--) {
		pol_mult(g, h, h);
		sc_add(f->coeffs[i], h->coeffs[0], h->coeffs[0]);
	}
}
