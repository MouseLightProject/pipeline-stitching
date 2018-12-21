function [X_,Y_,out,valid] =  fcestimate(X_,Y_,iadj,params,pinit_model)
%FCESTIMATE Summary of this function goes here
%   Detailed explanation goes here
viz = params.viz;
% model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
model = params.model;

if isfield(params,'dims')
    dims = params.dims;
else
    dims = [1024 1536 251];
end
dimcent = dims(setdiff([1:2],iadj))/2; %#ok<NBRAK> % center of image along curvature axis

% polynomial coeeficients (p3-p2(y-p1)^2):
% p(1) : imaging center ~ dims/2 +/- %10
% p(2) : curvature: -/+[1e-5 1e-6] for x & y, use initialization if
% avaliable, curvature might flip sign based on objective or reduce to 0 as
% medium changes (due to temperature or adding liquid solution, etc.)
% p(3): avarage displacement: between [[1-%overlap]*dims dims],
% initialization might not be useful, as this reduces to mean descriptor
% displacement
if nargin==5
    %pinit = params.init(iadj,:);
    pinit = pinit_model(iadj,:);
else
    dispvec = X_-Y_;
    y = dispvec(:,iadj);
    if iadj==1 % along x
        pinit = [dimcent dimcent^-2 median(y)]; % median might be off for p(3), as curvature will bend the displacement. doing center weighted median will be more accurate.
    else % along y
        pinit = [dimcent -dimcent^-2 median(y)]; % median might be off for p(3), as curvature will bend the displacement. doing center weighted median will be more accurate.
    end
end


%%
[out,x,y,yest] = fit2disp(X_,Y_,iadj,model,pinit,dimcent);
outliers = abs(y-yest)>2; % reject anything more than 2 pix away

%%
% if percentage of otliers is large, dont do correction!!
if sum(outliers)/length(outliers) < .25
    X_ = X_(~outliers,:);
    Y_ = Y_(~outliers,:);
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
    out = nan(1,3);
    valid = 0;
end

%%
if viz, figno=303; util.debug.vizCurvature; end %#ok<NASGU>

%%
end

function [out,x,y,yest] = fit2disp(X_,Y_,iadj,model,pinit,dimcent)
% fits polynomial on displacement fields for a given direction
%% set upper and lower boundaries
%p(1): imaging center is around image center (ideally). pinit(1)*.05 rougly
% corresponds to 2.5% of image size and set based on previous imaging that
% we have done. On '2018-08-18' sample, imaging center is around [465 733].
% If we use image center as initialization, it corresponds to .1 shift on
% 'x', if we use beads or previous imaging as initialization, we get .01
% shift which significantly improves.
lb1 = pinit(1)-pinit(1)*0.1;
ub1 = pinit(1)+pinit(1)*0.1;
% p(2): curvature should rely on initialization as magnitude and sign might change
% percent ratios do not make sense here as this number get squared, so
% provide a large range
% (dims/2).^2*1e-5 ~ [2.5 6] pixels in x & y. we have "-" cancave
% curvature for x overlap and "+" convex curvature for y overlap. by
% setting [-2e-5 1.5e-5], we roughly force "maximum" of 5 pixel warp along x,
% and 9 pixel warp along y overlap regions
lb2 = -2e-5;
ub2 = 1.5e-5;
% p(3): mean displacement is initialized based on descriptors, and stage is
% mostly accurate, use a tight bound on stage displacement
lb3 = pinit(3)-3;
ub3 = pinit(3)+3;

ub = [ub1 ub2 ub3];
lb = [lb1 lb2 lb3];
%%
dispvec = X_-Y_;
y = dispvec(:,iadj);

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

x = X_(:,setdiff([1 2],iadj));


fun = @(p) sum((y-feval(model,p,x)).^2);
% sqerr = @(p) sum((y-feval(model,p,x)).^2);

options = optimoptions('fmincon','Display', 'off');
options.ConstraintTolerance = 1e-9;
options.OptimalityTolerance = 1e-9;
nonlfun = @(x) match.edgeconstraint(x,model,pinit,dimcent);
out = fmincon(fun,pinit,[],[],[],[],lb,ub,nonlfun,options);
% modd = pinit_model(:,:);
% modd(iadj,:) = out;
% match.vizCurvature(modd)
%%
% outlier rejection based on parametric model
yest = feval(model,out,x);
end



