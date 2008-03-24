default(timer, 1);

/* Only works for lists of positive integers. */
merge(v, w) =
{
	local(a,b,u,uu,i,j,k);
	a = #v;
	b = #w;
	u = vector(a+b+1);
	u[1] = -1;
	i = 1;
	j = 1;
	k = 2;
	while ((i <= a) && (j<= b),
		if (v[i] < w[j],
			if (u[k-1] < v[i],
				u[k] = v[i];
				k++
			);
			i++
		,
			if (u[k-1] < w[j],
				u[k] = w[j];
				k++
			);
			j++
		)
	);
	while (i <= a,
		if(u[k-1] < v[i], u[k] = v[i]; k++);
		i++
	);
	while (j <= b,
		if(u[k-1] < w[j], u[k] = w[j]; k++);
		j++
	);
	uu = vector(k-2, i, u[i+1]);
	return(uu)
}


do_list(f, n) =
{
local(p,i,j,k,l,fp,lijst,sums,translate,nr,d);
if(n,,n=20);
i = 1;
while(i <= n,
	p = prime(i);
	print1(p,": ");
	fp = subst(f, u, x^p);
	sums = [0];
	lijst = factormod(fp, p);
	j = 1;
	while (j <= matsize(lijst)[1],
		d = poldegree(lijst[j,1]);
		k = 1;
		while (k <= lijst[j,2],
			nr = #sums;
			translate = vector(nr,l,sums[l]+d);
			sums = merge(sums, translate);
			print1(d, ", ");
			k++
		);
		j++
	);
	print("");
	closest = 1;
	distance_squared = (p - sums[1])^2;
	j = 2;
	while (j <= #sums,
		d = (p - sums[j])^2;
		if (d < distance_squared, closest = j; distance_squared = d);
		j++
	);
	if(distance_squared < 4*p,
		print("Closest is ",sums[closest],
			" with dist^2 = ",distance_squared)
	);
	i++
);
return
}
