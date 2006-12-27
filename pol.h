/*
 *	pol.h
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
void make_pol(polynomial f);
void free_pol(polynomial f);
void copy_pol(polynomial f, polynomial g);
void negate_pol(polynomial f, polynomial g);
void print_pol(polynomial f);
void random_pol(polynomial f, unsigned int d);
void pol_add(polynomial f, polynomial g, polynomial h);
void times_int(int mult, polynomial f, polynomial g);
void times_scalar(mscalar a, polynomial f, polynomial g);
void pol_mult(polynomial f, polynomial g, polynomial h);
void resize_pol(polynomial f, unsigned int new_length);
