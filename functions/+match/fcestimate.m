function [X_, Y_, final_parameters, valid] =  fcestimate(X_, Y_, axis_index, params, pinit_model, dims)
    %FCESTIMATE Summary of this function goes here
    %   Detailed explanation goes here
    viz = params.viz;
    % model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
    model = params.model;

%     if isfield(params,'dims')
%         dims = params.dims;
%     else
%         dims = [1024 1536 251];
%     end
    other_axis_index = 3 - axis_index ;
    dimcent = dims(other_axis_index)/2 ;  % center of image along curvature axis

    % polynomial coeeficients (p3-p2(y-p1)^2):
    % p(1) : imaging center ~ dims/2 +/- %10
    % p(2) : curvature: -/+[1e-5 1e-6] for x & y, use initialization if
    % avaliable, curvature might flip sign based on objective or reduce to 0 as
    % medium changes (due to temperature or adding liquid solution, etc.)
    % p(3): avarage displacement: between [[1-%overlap]*dims dims],
    % initialization might not be useful, as this reduces to mean descriptor
    % displacement
    initial_parameters = pinit_model(axis_index,:);
%     dispvec = X_-Y_;
%     y = dispvec(:,iadj);
%     if iadj==1 % along x
%         pinit = [dimcent dimcent^-2 median(y)]; % median might be off for p(3), as curvature will bend the displacement. doing center weighted median will be more accurate.
%     else % along y
%         pinit = [dimcent -dimcent^-2 median(y)]; % median might be off for p(3), as curvature will bend the displacement. doing center weighted median will be more accurate.
%     end


    %%
    [final_parameters,x,y,yest] = fit2disp(X_,Y_,axis_index,model,initial_parameters,dimcent);
    is_outlier_from_match_index = abs(y-yest)>2; % reject anything more than 2 pix away

    %%
    if viz ,
        util.debug.vizCurvature(x, y, model, final_parameters, initial_parameters, is_outlier_from_match_index, axis_index) ;
    end
    %if viz, figno=303; util.debug.vizCurvature; end %#ok<NASGU>

    %%
    % if percentage of outliers is large, don't do correction!!
    fraction_outlier = mean(is_outlier_from_match_index) ;
    if fraction_outlier < .25 ,
        X_ = X_(~is_outlier_from_match_index,:);
        Y_ = Y_(~is_outlier_from_match_index,:);
    %         % fit to inliers to improve estimation. use previous estimation as
    %         % initialization
    %         out1 = out;
    %         [out,x,y,yest] = fit2disp(X_,Y_,iadj,model,out,dimcent);
    %         x_range = min(x):max(x);
    %         y_range = feval(model,out,x_range);
    %         y_range_init = feval(model,pinit,x_range);
    %         outliers = abs(y-yest)>2; % reject anything more than 2 pix away
    %         yest1 = feval(model,out1,x)
    %         yest2 = feval(model,out,x)
        valid=1;
    else
        final_parameters = nan(1,3);
        valid = 0;
    end
end



function [final_parameters, x, y, y_est] = fit2disp(X_, Y_, axis_index, model, initial_parameters, dimcent)
    % fits polynomial on displacement fields for a given direction
    %% set upper and lower boundaries
    %p(1): imaging center is around image center (ideally). pinit(1)*.05 rougly
    % corresponds to 2.5% of image size and set based on previous imaging that
    % we have done. On '2018-08-18' sample, imaging center is around [465 733].
    % If we use image center as initialization, it corresponds to .1 shift on
    % 'x', if we use beads or previous imaging as initialization, we get .01
    % shift which significantly improves.
    lb1 = initial_parameters(1)-initial_parameters(1)*0.1;
    ub1 = initial_parameters(1)+initial_parameters(1)*0.1;
    % p(2): curvature should rely on initialization as magnitude and sign might change
    % percent ratios do not make sense here as this number get squared, so
    % provide a large range
    % (dims/2).^2*1e-5 ~ [2.5 6] pixels in x & y. we have "-" cancave
    % curvature for x overlap and "+" convex curvature for y overlap. by
    % setting [-2e-5 1.5e-5], we roughly force "maximum" of 5 pixel warp along x,
    % and 9 pixel warp along y overlap regions
    lb2 = -2.0e-5 ;
    ub2 = +1.5e-5 ;
    
    % p(3): mean displacement is initialized based on descriptors, and stage is
    % mostly accurate, use a tight bound on stage displacement
%     lb3 = initial_parameters(3)-3;
%     ub3 = initial_parameters(3)+3;
    % Not anymore! ALT, 2021-06-03
    lb3 = initial_parameters(3)-40 ;
    ub3 = initial_parameters(3)+40 ;
    
    ub = [ub1 ub2 ub3];
    lb = [lb1 lb2 lb3];
    %%
    dispvec = X_-Y_;
    y = dispvec(:,axis_index);
    
    %%% for non focus axis reject outliers based on vector norm. This should be
    %%% (roughly) constant for non curvature directions
    %     vcomp = dispvec(:,setdiff(1:3,iadj));
    %     medvcomp = median(vcomp);
    %     normvcomp = vcomp-ones(size(vcomp,1),1)*medvcomp;
    %     normvcomp = sqrt(sum(normvcomp.*normvcomp,2));
    %     validinds = normvcomp<util.get1DThresh(normvcomp,20,.95);
    
    validinds = 1:length(y);
    X_ = X_(validinds,:);
    % Y_ = Y_(validinds,:);
    y = y(validinds,:);
    
    other_axis_index = 3 - axis_index ;
    x = X_(:,other_axis_index) ;
    
    
    error_function = @(p) sum((y-feval(model,p,x)).^2);
    % sqerr = @(p) sum((y-feval(model,p,x)).^2);
    
    options = optimoptions('fmincon', 'Display', 'none') ;
    %options.ConstraintTolerance = 1e-9 ;
    options.OptimalityTolerance = 1e-9 ;
    %nonlfun = @(parameters)(match.edgeconstraint(parameters, model, initial_parameters, dimcent)) ;
    %initial_error = feval(error_function, initial_parameters)
    %[initial_cneq, initial_ceq] = nonlfun(initial_parameters)
    [final_parameters, final_error, exit_code] = fmincon(error_function,initial_parameters,[],[],[],[],lb,ub,[],options) ;  %#ok<ASGLU>
    %[final_cneq, final_ceq] = nonlfun(final_parameters)
    %final_error
    %exit_code
    y_est = feval(model,final_parameters,x) ;
end



