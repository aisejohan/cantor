/*
 *	scalar.h
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

scalar *neg_invs;
scalar **sums;
scalar **muls;
int prime;
void change_prime(int p);
void print_scalar(scalar a);

/*
scalar sc_neg_inv(scalar i);
scalar sc_sum(scalar i, scalar j);
scalar sc_mul(scalar i, scalar j);
*/

#define sc_neg_inv(i)		neg_invs[i]
#define sc_sum(i,j)		sums[i][j]
#define sc_mul(i,j)		muls[i][j]
