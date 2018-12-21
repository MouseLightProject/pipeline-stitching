function [c,ceq] = edgeconstraint(x,model,pinit,dimcent)

% init_vals = [feval(model,pinit,1) feval(model,pinit,dimcent) feval(model,pinit,2*dimcent)];
% iter_vals = [feval(model,x,1) feval(model,x,dimcent) feval(model,x,2*dimcent)]; 
init_vals = [feval(model,pinit,1) feval(model,pinit,2*dimcent)];
iter_vals = [feval(model,x,1) feval(model,x,2*dimcent)]; 

% below is a nice trick to compansate for any stage movement errors.
% Sometimes, stage move x-um, which results in X-pixels, but if you look at
% images, they are X +/- eps, where eps can be around 3-4 pixels. By taking
% out, p(3) which corresponds to shift, we get better estimate

% difs = ( init_vals-pinit(3) - (iter_vals-x(3))); 
difs = ( init_vals-pinit(3) - (iter_vals-x(3))); 

c = norm(difs);
ceq = [];
