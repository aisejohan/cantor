/*
 *	scalar.c
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

	fprintf(fd, "g = \\\n");

	i = 0;
	uit = mpz_out_str(fd, 10, f->coeffs[i]);
	fprintf(fd, " + \\\n");
	while (uit && i <= f->degree) {
		uit = mpz_out_str(fd, 10, f->coeffs[i]);
		fprintf(fd, "*x^%d + \\\n", i);
		i++;
	}
	fprintf(fd,"0;\n");

	fclose(fd);
}

