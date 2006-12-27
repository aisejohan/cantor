/*
 *	inout.c
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

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>

#include "data.h"
#include "scalar.h"
#include "pol.h"

void read_in(polynomial f)
{
	int i;
	int uit;
	FILE *fd;
	mscalar c;

	make_scalar(c);

	fd = fopen("input", "r");
	if (!fd) {
		perror("Not opening file <input>. Does it exist?");
		exit(1);
	}

	i=0;
	uit = mpz_inp_str(c, fd, 10);
	while (uit) {
		f->degree = i;
		resize_pol(f, i);
		sc_copy(c, f->coeffs[i]);
		i++;
		uit = mpz_inp_str(c, fd, 10);
	}

	fclose(fd);
	free_scalar(c);
}

/* Error message on writing. */
void error_msg(int uit)
{
	if(uit < 0) {
		perror("Not writing correctly to file! uit = %d.");
		exit(1);
	}
}

/* For use in pari... */
void put_out(polynomial f)
{
	int i;
	int uit;
	FILE *fd;

	fd = fopen("output", "w");
	if (!fd) {
		perror("Not opening file <output>!?");
		exit(1);
	}

	uit = fprintf(fd, "g = \\\n");
	error_msg(uit);

	i = 0;
	while (i <= f->degree) {
		uit = fprintf(fd, "Mod(");
		error_msg(uit);
		uit = mpz_out_str(fd, 10, f->coeffs[i]);
		error_msg(uit);
		uit = fprintf(fd, ", %d^%d)*x^%d + \\\n", p, r, i);
		error_msg(uit);
		i++;
	}
	uit = fprintf(fd,"0;\n");
	error_msg(uit);

	fclose(fd);
}
