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

#include <stdio.h>
#include <stdlib.h>

#include "data.h"
#include "scalar.h"
#include "pol.h"
#include "xu_and_sparse.h"

void make_sparse_pol(sparse_polynomial *f, unsigned int len)
{
	*f = (sparse_polynomial) malloc(sizeof(struct sparse_pol));
	(*f)->degree = 0;
	(*f)->length = len;
	(*f)->exps = (unsigned int *)
		malloc(((*f)->length+1)*sizeof(unsigned int));
	(*f)->coeffs = (scalar *)
		malloc(((*f)->length+1)*sizeof(scalar));
	if ((!(*f)->coeffs) || (!(*f)->exps)) {
		perror("Malloc failed in make_sparse_pol!");
		exit(1);
	}
}

void free_sparse_pol(sparse_polynomial *f)
{
	free((*f)->exps);
	free((*f)->coeffs);
	free(*f);
}

void convert_to_sparse(sparse_polynomial *sA, polynomial A)
{
	int len, i, j;

	len = 0;
	i = 0;
	while (i < A->degree) {
		if (A->coeffs[i]) len++;
		i++;
	}
	make_sparse_pol(sA, len);
	(*sA)->degree = A->degree;
	i = 0;
	j = 0;
	while (i <= A->degree) {
		if (A->coeffs[i]) {
			(*sA)->exps[j] = i;
			(*sA)->coeffs[j] = A->coeffs[i];
			j++;
		}
		i++;
	}
}


/* r = g mod f */
void r_reduce_sparse(polynomial r, polynomial g, sparse_polynomial f)
{
	int i, j;
	scalar c, q;

	if (r != g) copy_pol(r, g);

	i = r->degree - f->degree;
	c = sc_neg_inv(f->coeffs[f->length]);
	while (i >= 0) {
		q = sc_mul(c, r->coeffs[r->degree]);
		j = 0;
		while (j < f->length) {
			r->coeffs[f->exps[j]+i] =
				sc_sum(r->coeffs[f->exps[j]+i],
					sc_mul(q, f->coeffs[j]));
			j++;
		}
		r->degree--;
		while (!r->coeffs[r->degree] && r->degree > 0) r->degree--;
		i = r->degree - f->degree;
	}
}

void random_xu(xu_polynomial *f, unsigned int dx, unsigned int du)
{
	int i,j;

	*f = (xu_polynomial) malloc(sizeof(struct xu_pol));
	(*f)->dx = dx;
	(*f)->du = du;
	(*f)->coeffs = (int **)malloc((dx + 1)*sizeof(int *));
	i = 0;
	while (i <= dx) {
		(*f)->coeffs[i] = (int *)malloc((du + 1)*sizeof(int));
		j = 0;
		while (j <= du) {
			(*f)->coeffs[i][j] = (rand() % 3) - 1;
/*			rand() - (RAND_MAX > 1);  */
			j++;
		}
		i++;
	}
}

void read_xu(xu_polynomial *f, unsigned int dx, unsigned int du)
{
	int i,j,c;

	printf("Give the coefficients of f: \n");

	*f = (xu_polynomial) malloc(sizeof(struct xu_pol));
	(*f)->dx = dx;
	(*f)->du = du;
	(*f)->coeffs = (int **)malloc((dx + 1)*sizeof(int *));
	i = 0;
	while (i <= dx) {
		(*f)->coeffs[i] = (int *)malloc((du + 1)*sizeof(int));
		j = 0;
		while (j <= du) {
			printf("Coeff of x^%d u^%d: ",i, j);
			scanf("%d", &c);
			(*f)->coeffs[i][j] = c;
			j++;
		}
		i++;
	}
}

void free_xu_pol(xu_polynomial *f)
{
	int i;

	i = 0;
	while (i <= (*f)->dx) {
		free((*f)->coeffs[i]);
		i++;
	}
	free((*f)->coeffs);
	free(*f);
}


/* Assumes prime > dx. */
int xu_to_sparse(sparse_polynomial *sA, xu_polynomial A)
{
	int nr, degree, i, j, k;
	int a;

	degree = 0;
	nr = 0;
	j = 0;
	while (j <= A->du) {
		i = 0;
		while (i <= A->dx) {
			if (A->coeffs[i][j] % prime) {
				nr++;
				if (degree < i + prime*j) degree = i + j*prime;
			}
			i++;
		}
		j++;
	}
	if ((nr == 0) || (prime <= A->dx)) return(0);

	make_sparse_pol(sA, nr - 1);
	(*sA)->degree = degree;

	k = 0;
	j = 0;
	while (j <= A->du) {
		i = 0;
		while (i <= A->dx) {
			a = A->coeffs[i][j] % prime;
			if (a) {
				if (a < 0) a = prime + a;
				(*sA)->exps[k] = i + prime*j;
				(*sA)->coeffs[k] = a;
				k++;
			}
			i++;
		}
		j++;
	}
	return(nr);
}

void print_xu_pol(xu_polynomial A)
{
	int i, j, y, a;

	y = 0;
	j = 0;
	while (j <= A->du) {
		i = 0;
		while (i <= A->dx) {
			a = A->coeffs[i][j];
			if (a) {
				if (y) printf("+ ");
				y = 1;
				printf("%d ", a);
				if (i) printf("* x^%d ", i);
				if (j) printf("* u^%d ", j);
			}
			i++;
		}
		j++;
	}
	printf("\n");
}


void sparse_to_pol(polynomial A, sparse_polynomial sA)
{
	int i;

	A->degree = sA->degree;
	resize_pol(A, sA->degree);
	
	i = 0;
	while (i <= sA->degree) {
		A->coeffs[i] = 0;
		i++;
	}

	i = 0;
	while (i <= sA->length) {
		A->coeffs[sA->exps[i]] = sA->coeffs[i];
		i++;
	}
}

void
prime_power_sparse(
	polynomial h,
	polynomial g,
	polynomial f,
	sparse_polynomial sf)
{
	int i;
	polynomial tmp;

#ifdef KIJKEN
	test_pol(h);
	test_pol(g);
	test_pol(f);
#endif

	make_pol(&tmp);
	copy_pol(tmp, g);

	h->degree = 0;
	h->coeffs[0] = 1;

	i = 1;
	do {
		if (prime & i) {
			pol_mult(h, h, tmp);
			r_reduce_sparse(h, h, sf);
			r_reduce(h, h, f);
		}
		pol_mult(tmp, tmp, tmp);
		r_reduce_sparse(tmp, tmp, sf);
		r_reduce(tmp, tmp, f);
		i = 2*i;
	} while (i <= prime);

	free_pol(&tmp);
}

int next_degree_sparse(
	polynomial g,
	polynomial f,
	sparse_polynomial sf,
	unsigned int d,
	polynomial h)
{
	polynomial tmp1, tmp2;

#ifdef KIJKEN
	if ((f == g) || (f == h) || (g == h)) {
		printf("Should not have same pols in next_degree.");
		exit(1);
	}
#endif

	d++;
	while (2*d <= f->degree) {
		prime_power_sparse(h, h, f, sf);
		if (h->degree == 1 && h->coeffs[1] == 1) {
			if (h->coeffs[0] == 0) {
				copy_pol(g, f);
				f->degree = 0;
				f->coeffs[0] = 1;
				resize_pol(f, f->degree);
				h->degree = 0;
				h->coeffs[0] = 0;
				resize_pol(h, h->degree);
				return(d);
			}
			/* This is the case where h = x + a, a not 0. */
			g->degree = 0;
			g->coeffs[0] = 1;
		} else {
			h->coeffs[1] = (h->coeffs[1] + prime - 1) % prime;
			gcd(g, h, f);
			h->coeffs[1] = (h->coeffs[1] + 1) % prime;
		}
		if (g->degree > 0) {
			make_pol(&tmp1);
			make_pol(&tmp2);
			qr_reduce(tmp1, f, tmp2, g);
			copy_pol(f, tmp2);
			if (f->degree) r_reduce(h, h, f); /* FIXME */
			free_pol(&tmp1);
			free_pol(&tmp2);
			return(d);
		}
		d++;
	}
	return(0);
}

void print_degrees_sparse(polynomial f, sparse_polynomial sf)
{
	int e, i, y;
	polynomial h;
	polynomial g;

	make_pol(&h);
	make_pol(&g);

	h->degree = 1;
	h->coeffs[1] = 1;

	y = 1;
	e = 0;
	do {
		e = next_degree_sparse(g, f, sf, e, h);
		if (e) {
			i = 0;
			while (i < g->degree/e) {
				if (y) {
					printf("%d", e);
					y = 0;
				} else {
					printf(", %d", e);
				}
				fflush(stdout);
				i++;
			}
			if (!f->degree) e = 0;
		} else {
			if (f->degree) {
				if (y) {
					printf("%d", f->degree);
				} else {
					printf(", %d", f->degree);
				}
			}
		}
	} while (e);
	free_pol(&h);
	free_pol(&g);
}
