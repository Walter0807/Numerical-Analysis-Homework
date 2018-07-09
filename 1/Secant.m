f = @(x) log(1+x)-cos(x);
p = zeros(1,1000);
epsilon = 1e-6;
p(2) = 1;
i = 2;
while (abs(p(i)-p(i-1))>epsilon)
  p(i+1) = p(i) - (f(p(i))*(p(i)-p(i-1)))/(f(p(i))-f(p(i-1)));    
  fprintf('p_%d = %f, f(p_%d) = %f\n', i, p(i+1), i, f(p(i+1)));
  i = i + 1;
end

