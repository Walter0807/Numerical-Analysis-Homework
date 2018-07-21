
# Permuted LU

Initialization & find solution:


```matlab
A = [2 5 -9 3; 5 6 -4 2; 3 -4 2 7; 11 7 4 -8];
b = [151 103 16 -32]';
b1 = b;
ans = A\b
n = size(A,1);
P = eye(n);
```


    ans =
    
        3.0000
        5.0000
      -11.0000
        7.0000



LU factorization with partial pivoting:


```matlab
for k = 1:n-1
    [v,pos]=max(A(k:n,k));
    % swap rows
    A([k;pos+k-1],:) = A([pos+k-1;k],:);
    b([k;pos+k-1],:) = b([pos+k-1;k],:);
    P([k;pos+k-1],:) = P([pos+k-1;k],:);
    for i = k+1:n
        r = A(i,k) / A(k,k);
        b(i,1) = b(i,1) - r*b(k,1);
        A(i,k+1:n) = A(i,k+1:n) - r*A(k,k+1:n);
        A(i,k) = r;
    end
end
% Print L
fprintf('L:')
eye(n) + tril(A,-1)
% Print U
fprintf('U:')
triu(A)
% Print P
P
% Print LU
fprintf('Check LU:')
(eye(n) + tril(A,-1))*triu(A)
```

    L:
    ans =
    
        1.0000         0         0         0
        0.1818    1.0000         0         0
        0.4545    0.7561    1.0000         0
        0.2727   -1.5854   -9.4444    1.0000
    
    U:
    ans =
    
       11.0000    7.0000    4.0000   -8.0000
             0    3.7273   -9.7273    4.4545
             0         0    1.5366    2.2683
             0         0         0   37.6667


​    
    P =
    
         0     0     0     1
         1     0     0     0
         0     1     0     0
         0     0     1     0
    
    Check LU:
    ans =
    
       11.0000    7.0000    4.0000   -8.0000
        2.0000    5.0000   -9.0000    3.0000
        5.0000    6.0000   -4.0000    2.0000
        3.0000   -4.0000    2.0000    7.0000



Solve `x`:


```matlab
% Print x
x = zeros(n,1);
A = triu(A);
for i = n: -1: 1
    x(i,1) = (b(i,1) - A(i,:)*x)/A(i,i);
end
x
```


    x =
    
        3.0000
        5.0000
      -11.0000
        7.0000



# Tridiagonal LU

**a**

Using Gaussian Elimination, 

1. `row 2` +=  `row 1` * 1/2
2. `row 3` +=  `row 2` * 2/3
3. `row 4` +=  `row 3` * 3/4


$$
L = \left[
 \begin{matrix}
   1 &  & & \\
   -\frac{1}{2} & 1  \\
    & -\frac{2}{3} & 1 \\
   & & -\frac{3}{4} & 1
  \end{matrix}
  \right]
$$

$$
U = \left[
 \begin{matrix}
   2 & -1 & & \\
    & \frac{3}{2} & -1 &  \\
    & & \frac{4}{3} & -1 \\
   & & & \frac{5}{4}
  \end{matrix}
  \right]
$$

$$
Ux = L^{-1}b = \left[
 \begin{matrix}
   1 \\
   \frac{3}{2} \\
   2 \\
   \frac{5}{2}
  \end{matrix}
  \right]
$$

By backward substitution,
$$
x = \left[
 \begin{matrix}
   2 \\
   3 \\
   3 \\
   2
  \end{matrix}
  \right]
$$


**b**



```matlab
A = [2 -1 0 0; -1 2 -1 0; 0 -1 2 -1; 0 0 -1 2];
b = [1 1 1 1]';
n = size(A,1);

for k = 1:n-1
    bound = min(k+2,n);
    r = A(k+1,k) / A(k,k);
    b(k+1,1) = b(k+1,1) - r*b(k,1);
    A(k+1,k+1:bound) = A(k+1,k+1:bound) - r*A(k,k+1:bound);
    A(k+1,k) = r;
end
L = eye(n) + tril(A,-1)
U = triu(A)
x = zeros(n,1);
A = triu(A);
for i = n: -1: 1
    bound = min(k+2,n);
    x(i,1) = (b(i,1) - U(i,i:bound)*x(i:bound,1))/U(i,i);
end
x
```


    L =
    
        1.0000         0         0         0
       -0.5000    1.0000         0         0
             0   -0.6667    1.0000         0
             0         0   -0.7500    1.0000


​    
    U =
    
        2.0000   -1.0000         0         0
             0    1.5000   -1.0000         0
             0         0    1.3333   -1.0000
             0         0         0    1.2500


​    
    x =
    
        2.0000
        3.0000
        3.0000
        2.0000



Define -, *, \ operations of matrix elements as basic operations:

Computing rounds and basic operations per loop,
$$
Number of computations = 7(n-2) + 5 +3(n-1) + 1 = 10n-11
$$
and
$$
N = 3n-2
$$

Therefore, number of computations roughly equals to $\frac{10}{3}N$

**c**

As demonstrated, the results are same.


# Facts about LTs

**a**


$n\times n$ matrix $L$ is invertible $\Leftrightarrow$ $det(L)\neq 0$
$$
\begin{eqnarray}
det(L) &=& 
\sum\limits_{i=1}^{n}l_{i,n}(-1)^{i+n} det(L_{n-1})
\\ &=& l_{n,n} \cdot det(L_{n-1})
\\ &=& l_{n,n} \cdot l_{n-1,n-1} \cdot det(L_{n-2})
\\ &=& \cdots
\\ &=& \prod\limits_{i=1}^n l_{i,i}
\end{eqnarray}
$$
where
$$
L_n = L(1 : n, 1 : n)
$$
Therefore, if $L$ is any $n \times n$ lower triangular invertible matrix, then $l_{jj} \neq 0$ for every $1 \leq j \leq n$.

**b**

Define 
$$
L(i,j) = l_{i,j}, L^{-1}(i,j) = p_{i,j}
$$

For $1<k\leq n$, 
$$
\begin{eqnarray}
\sum\limits_{j=1}^{n}l_{1,j}p_{j,k} &=& I_{1,k} \\
&=& l_{1,1} \cdot p_{1,k} \\
&=& 0
\end{eqnarray}
$$
From *(a)*, $l_{1,1} \neq 0$. So $p_{1,k} = 0$.

Suppose for $i \leq i_0, i<k \leq n$, We have $p_{i,k} = 0$. Now consider $i = i_0+1, k > i$

$$
\sum\limits_{j=1}^{n}l_{i+1,j}p_{j,k} = I_{i+1,k} = 0
$$
When $j<i-1 = i_0$, $p_{j,k} = 0$; When $j>i$, $l_{i+1,j} = 0$.

Therefore,

$$
\sum\limits_{j=1}^{n}l_{i+1,j}p_{j,k} = l_{i,i} \cdot p_{i,k} = 0 
$$

From *(a)*, $l_{i,i} \neq 0$. So $p_{i,k} = 0$.

Then for $i \leq i_0+1, k>i$, We have $p_{i,k} = 0$.

Using *complete induction*, We have $p_{i,k} = 0$ for $1 \leq i < k \leq n$.

So, if $L$ is an $n×n$ invertible lower triangular matrix, then $L^{−1}$ is also lower triangular.

**c**

Following the conventions in Problem 3 of HW2, we have
$$
e_k^{T}g_{j} = 0+0+\cdots+0 = 0
$$
where
$$
1 \leq k < j \leq n-1
$$
So,
$$
\begin{eqnarray}
L_1L_2\cdots L_{n-1} &=& (I + g_1e_1^{T})(I + g_2e_2^{T}) \cdots (I + g_{n-1}e_{n-1}^{T})\\
&=& [I + (g_1e_1^{T} + g_2e_2^{T}) + g_1e_1^{T}g_2e_2^{T}] (I + g_{3}e_3^{T}) \cdots (I + g_{n-1}e_{n-1}^{T})\\
&=& [I + (g_1e_1^{T} + g_2e_2^{T})](I + g_{3}e_k^{T}) \cdots (I + g_{n-1}e_k^{T})\\
&=& [I + (g_1e_1^{T} + g_2e_2^{T}+ g_3e_3^{T})](I + g_{4}e_k^{T}) \cdots (I + g_{n-1}e_k^{T})\\
&=& \cdots \\
&=& I + \sum\limits_{k=1}^{n-1}g_ke_k^{T}\\
\end{eqnarray}
$$




