function [c,ceq] = edgeconstraint(parameters, model, initial_parameters, dimcent)
    % init_vals = [feval(model,pinit,1) feval(model,pinit,dimcent) feval(model,pinit,2*dimcent)];
    % iter_vals = [feval(model,x,1) feval(model,x,dimcent) feval(model,x,2*dimcent)];
    initial_extremal_y_values = [feval(model,initial_parameters,1) feval(model,initial_parameters,2*dimcent)];
    extremal_y_values = [feval(model,parameters,1) feval(model,parameters,2*dimcent)];
    
    % below is a nice trick to compansate for any stage movement errors.
    % Sometimes, stage move x-um, which results in X-pixels, but if you look at
    % images, they are X +/- eps, where eps can be around 3-4 pixels. By taking
    % out, p(3) which corresponds to shift, we get better estimate
    
    initial_shift = initial_parameters(3) ;
    shift = parameters(3) ;
    
    initial_extremal_y_values_relative_to_shift = initial_extremal_y_values-initial_shift ;
    extremal_y_values_relative_to_shift = extremal_y_values-shift ;
    
    delta = initial_extremal_y_values_relative_to_shift - extremal_y_values_relative_to_shift ;
    
    c = norm(delta) ;
    ceq = [] ;
%     ceq = norm(delta) ;
%     c = [] ;
end
