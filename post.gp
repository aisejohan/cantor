default(timer, 1);

list(f, n) =
{
local(p,i);
if(n,,n=20);
i = 1;
while(i <= n,
	p = prime(i);
	print1(p,": ");
	fp = subst(f, u, x^p);
	lijst = factormod(fp, p);
	for(j = 1, matsize(lijst)[1],
		print1(poldegree(lijst[j, 1]),", ")
	);
	print("");
	i++
);
return
}
