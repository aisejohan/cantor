/*
 *	fl.c
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
#include <unistd.h>

#include "data.h"
#include "scalar.h"
#include "modp_scalar.h"
#include "pol.h"
#include "utils.h"
#include "reduce.h"
#include "modp_pol.h"
#include "frob.h"
#include "inout.h"

int main()
{
	polynomial f,g;

	setup_scalars();
	setup_modp_scalars();

	printf("/* The prime is: %d. */\n",p);
	printf("/* The power is: %d. */\n",r);

	make_pol(f);
	make_pol(g);

	read_in(f);
	print_pol(f);
	frob_lift(f, g);
	print_pol(g);
	printf("\n");
	put_out(g);

	free_pol(g);
	free_pol(f);

	exit(0);
}
