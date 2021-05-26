As = [ 1 2 3 ; 4 5 6 ; 7 8 9] 
bs = [1 1 4]' ;

nx = 128
ny = 100 
nz = 50 
n = [nx-1 ny-1 0]' 

S = diag([-1 -1 +1]) 

A = As * S
b = As*n + bs

i = [0 0 0]'
is = S*i + n

xs = As*is + bs
x = A*i + b




i2 = [nx-1 ny-1 nz-1]'
is2 = S*i2 + n

xs2 = As*is2 + bs
x2 = A*i2 + b





