solve(0,1)
global t;
t = 0;

function x = solve(a,b)
    f = @(x) log(1+x)-cos(x);
    global t;
	epsilon = 1e-6;
	p = (a+b)/2;
    t = t + 1;
    fprintf('p_%d = %f, f(p_%d) = %f\n', t, p, t, f(p));
	if (p-a<=epsilon)
		x = p;
    elseif (f(a)*f(p)<0)
        x = solve(a,p);
    else
        x = solve(p,b);
	end
end        

