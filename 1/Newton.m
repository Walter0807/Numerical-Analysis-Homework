f = @(x) log(1+x)-cos(x);
fd = @(x) 1/(1+x)+sin(x);
p = 1/2;
q = 0;
epsilon = 1e-6;
t = 0;

while (1)
	t = t + 1;
	fprintf('p_%d = %f, f(p_%d) = %f\n', t, p, t, f(p));
	if (abs(p-q)<=epsilon)
        break
    end
    q = p;
	p = p - f(p)/fd(p);
end



