f = @(x) log(1+x)-cos(x);
p = [0,1,0,0,0];
for i = 2:4 
  p(i+1) = p(i) - (f(p(i))*(p(i)-p(i-1)))/(f(p(i))-f(p(i-1)));    
  fprintf('p_%d = %f, f(p_%d) = %f\n', i, p(i+1), i, f(p(i+1)));
end
