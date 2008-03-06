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


struct xu_pol {
	unsigned int dx;
	unsigned int du;
	int **coeffs; /* coeffs[i][j] is coeff of x^i u^j */
};

struct sparse_pol {
	unsigned int degree;
	unsigned int length; /* nonzero terms i=0,...,length */
	unsigned int *exps;
	scalar *coeffs;
};

typedef struct sparse_pol *sparse_polynomial;
typedef struct xu_pol *xu_polynomial;

void make_sparse_pol(sparse_polynomial *f, unsigned int len);
void convert_to_sparse(sparse_polynomial *sA, polynomial A);
void free_sparse_pol(sparse_polynomial *f);
void r_reduce_sparse(polynomial r, polynomial g, sparse_polynomial f);
void random_xu(xu_polynomial *f, unsigned int dx, unsigned int du);
void read_xu(xu_polynomial *f, unsigned int dx, unsigned int du);
void free_xu_pol(xu_polynomial *f);
int xu_to_sparse(sparse_polynomial *sA, xu_polynomial A);
void print_xu_pol(xu_polynomial A);
void sparse_to_pol(polynomial A, sparse_polynomial sA);
void print_degrees_sparse(polynomial f, sparse_polynomial sf);
