/*
 *	frob.c
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
#include "data.h"
#include "scalar.h"
#include "pol.h"
#include "utils.h"
#include "reduce.h"

void frob_lift(polynomial f, polynomial g)
{
	int i, k;
	polynomial fx, fp, rest, fxg, b, c;
	
	make_pol(fx);
	make_pol(fp);
	make_pol(rest);
	make_pol(fxg);
	make_pol(b);
	make_pol(c);

	deriv(f, fx);

	copy_pol(f, fp);
	for(i = 2; i <= p; i++)	pol_mult(f, fp, fp);

	g->degree = p;
	resize_pol(g, p);
	for(i = 0; i <= p-1; i++) sc_zero(g->coeffs[i]);
	sc_one(g->coeffs[p]);

	k = 1;

	while (k < r) {
		printf("This is stage %d.\n", k);
		substitute(f, g, rest);
printf("subs\n");
		reduce_by_g(rest, fp);
printf("red\n");
		substitute(fx, g, fxg);
printf("subs\n");
		invert_unit_mod_g(fxg, fp, c);
printf("inv\n");
		pol_mult(rest, c, c);
printf("mul\n");
		reduce_by_g(c, fp);
printf("red\n");
		negate_pol(c, c);
printf("add\n");
		pol_add(c, g, g);
		k = 2*k;
	}
}
