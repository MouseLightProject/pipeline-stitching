function [paireddescriptor, curvemodel] = xymatch(descriptors, neigs, scopeloc, params, model)
%ESTIMATESCOPEPARAMETERS Summary of this function goes here
%
% [OUTPUTARGS] = ESTIMATESCOPEPARAMETERS(INPUTARGS) Explain usage here
%
% Inputs:
%
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2016/09/12 10:38:28 $	$Revision: 0.1 $
% Copyright: HHMI 2016
%%
%addpath(genpath('./thirdparty'))
neigs_row_count = size(neigs,1);
if nargin<5
    model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
end

debug = 0;
%res = 0;
viz=0;
fignum = 101;

projectionThr = 20; % distance between target and projected point has to be less than this number
dims = params.imagesize;
imsize_um = params.imsize_um;

% slid = [[75 960];[0 dims(2)];[0 dims(3)]];
% expensionshift = [0 0 20]; % HEURISTICS:: tissue expends, so overlap is bigger between tiles

%%
optimopts = statset('nlinfit');
optimopts.RobustWgtFun = 'bisquare';
% opt.method='nonrigid_lowrank';
opt.method='nonrigid';
opt.beta=6;            % the width of Gaussian kernel (smoothness), higher numbers make transformation more stiff
opt.lambda=16;          % regularization weight
opt.viz=0;              % show every iteration
opt.outliers=0.9;       % use 0.7 noise weight
opt.fgt=0;              % do not use FGT (default)
opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=1;          % compute correspondence vector at the end of registration (not being estimated by default)
%     opt.max_it=100;         % max number of iterations
%     opt.tol=1e-10;          % tolerance
matchparams.model = model;
matchparams.optimopts = optimopts;
matchparams.opt = opt;
matchparams.projectionThr = projectionThr;
matchparams.debug = debug;
matchparams.viz = viz;
matchparams.fignum = fignum;
matchparams.opt.beta=2;
% matchparams.opt.method = 'nonrigid';

%%
pixshiftpertile = nan(neigs_row_count,3,2);
for ii=1:neigs_row_count
    for iadj = 1:2
        if isfinite(neigs(ii,iadj+1))
            stgshift = 1000*(scopeloc.loc(neigs(ii,iadj+1),:)-scopeloc.loc(neigs(ii,1),:));
            pixshift = round(stgshift.*(dims-1)./(imsize_um));
            pixshiftpertile(ii,:,iadj) = pixshift;
        end
    end
end
if viz
    these = isfinite(pixshiftpertile(:,1,1));onx = pixshiftpertile(these,1,1);medonx = median(onx);  %#ok<UNRCH>
    these = isfinite(pixshiftpertile(:,2,2));ony = pixshiftpertile(these,2,2);medony = median(ony);
    pixshiftpertile(isnan(pixshiftpertile(:,1,1)),:,1) = ones(sum(isnan(pixshiftpertile(:,1,1))),1)*[medonx,0,0];
    pixshiftpertile(isnan(pixshiftpertile(:,2,2)),:,2) = ones(sum(isnan(pixshiftpertile(:,2,2))),1)*[0,medony,0];
    figure(123), cla,subplot(121),hist(onx,100),hold on,subplot(122),hist(ony,100) %#ok<HIST>
end
% replace any nans (boundary tiles without an adjacent tile) with median values
meds = squeeze(median(pixshiftpertile,1,'omitnan'))';
for iadj = 1:2
    % replace any nan rows with median
    these = isnan(pixshiftpertile(:,1,iadj));
    pixshiftpertile(:,:,iadj) = util.rowreplace(pixshiftpertile(:,:,iadj),these,meds(iadj,:));
end
%%
% p(1): imaging center, ideally dims/2
% p(2): polynomial multiplier, "+" for diverging corner (hourglass),
%       "-" for converging (blob). For scope one, p(2) multipliers are "+"
%       for x direction, "-" for y direction
% p(3): displacement between tiles
if isfield(params,'beadmodel')
    % beadmodel = params.beadmodel;
    %based on median values, replace this with beadmodel @@ TODO @@
    matchparams.init(1,:)=[733 1.0214e-05 863];
    matchparams.init(2,:)=[465 -1.4153e-05 1451];
else
    matchparams.init_array=[]; % creates a array initialization based on stage displacements
    pvals_12 = [[733 1.02141e-05];[465 -1.4153e-05]];
    for it = 1:neigs_row_count
        for iadj = 1:2
            matchparams.init_array(iadj,:,it)=[pvals_12(iadj,:) pixshiftpertile(it,iadj,iadj)];
        end
    end
end

%% INITIALIZATION
% checkthese = [1 4 5 7]; % 0 - right - bottom - below
% indicies are 1 based,e.g. x = 1:dims(1), not 0:dims(1)-1
% xyz_umperpix = zeros(tile_count,3);
curvemodel = nan(3,3,neigs_row_count);
%medianResidualperTile = zeros(3,3,neigs_row_count);
paireddescriptor = cell(neigs_row_count,1);
% initialize
for ix = 1:neigs_row_count ,
    paireddescriptor{ix}.onx.valid = 0;
    paireddescriptor{ix}.onx.X = [];
    paireddescriptor{ix}.onx.Y = [];
    paireddescriptor{ix}.ony.valid = 0;
    paireddescriptor{ix}.ony.X = [];
    paireddescriptor{ix}.ony.Y = [];
    paireddescriptor{ix}.neigs = neigs(ix,:);
    paireddescriptor{ix}.count = [0 0];
end

% Make a template we'll use to initialize paired_descriptor_for_this_tile
% for each tile
paired_descriptor_for_this_tile_template=[];
paired_descriptor_for_this_tile_template{1}.valid = 0;
paired_descriptor_for_this_tile_template{1}.X = [];
paired_descriptor_for_this_tile_template{1}.Y = [];
paired_descriptor_for_this_tile_template{2}.valid = 0;
paired_descriptor_for_this_tile_template{2}.X = [];
paired_descriptor_for_this_tile_template{2}.Y = [];

%interiorTile_list = util.interior_tiles(scopeloc,1);

%%
try parfor_progress(0);catch;end
parfor_progress(neigs_row_count) ;

parfor neigs_row_index = 1:neigs_row_count ,
%tile_index_of_interest = find(strcmp('/2020-12-01/01/01916', scopeloc.relativepaths))
%for tile_index = tile_index_of_interest ,
    %% load descriptor pairs X (center) - Y (adjacent tile)
    idxcent = neigs(neigs_row_index,1);
    descent = descriptors{idxcent};
    if isempty(descent);continue;end
    descent = double(descent(:,1:3));
    if size(descent,1)<3;continue;end
    
    descent = util.correctTiles(descent,dims); % flip dimensions
       % ALT: I think this is because the fiducials are computed on the raw tile
       % stacks, which have to be flipped in x and y to get them into the proper
       % orientation for the rendered stack
    mout = nan(3,3);
    paired_descriptor_for_this_tile = paired_descriptor_for_this_tile_template;
    R_ = zeros(3); % median residual
    
    %%
    for iadj = 1:2 , %1:x-overlap, 2:y-overlap, 3:z-overlap
        %%
        % idaj : 1=right(+x), 2=bottom(+y), 3=below(+z)
        idxadj =  neigs(neigs_row_index,iadj+1);
        
        if isnan(idxadj);continue;end
        descadj = descriptors{idxadj};
        if isempty(descadj);continue;end
        
        descadj = double(descadj(:,1:3)); % descadj has x-y-z-w1-w2 format
        if size(descadj,1)<3;continue;end
        
        descadj = util.correctTiles(descadj,dims); % flip dimensions
        stgshift = 1000*(scopeloc.loc(idxadj,:)-scopeloc.loc(idxcent,:));
        pixshift = round(stgshift.*(dims-1)./(imsize_um));
        descadj = descadj + ones(size(descadj,1),1)*pixshift; % shift with initial guess based on stage coordinate
        
        %%
        nbound = [0 0];
        nbound(1) = max(pixshift(iadj),min(descadj(:,iadj))) - 15;
        nbound(2) = min(dims(iadj),max(descent(:,iadj)))+0 + 15;
        X = descent(descent(:,iadj)>nbound(1)&descent(:,iadj)<nbound(2),:);
        Y = descadj(descadj(:,iadj)>nbound(1)&descadj(:,iadj)<nbound(2),:);
        %%
        if size(X,1)<3 || size(Y,1)<3;continue;end
        % get descpair
        [X_,Y_] = match.descriptorMatch4XY(X,Y,matchparams);
        if size(X_,1)<3 || size(Y_,1)<3;continue;end
        Y_(:,iadj) = Y_(:,iadj)- pixshift(iadj);% move it back to original location after CDP
        %[X_e,Y_e,out_e,valid_e] = match.fcestimate(X_,Y_,iadj,matchparams,pinit_model);
        % [X_e2,Y_e2,out_e2,valid_e2] = match.fcestimate(X_e,Y_e,iadj,matchparams);
        %%
        if viz ,
            util.debug.vizMatchALT(scopeloc,neigs,descriptors,neigs_row_index,pixshift,iadj,X,Y,X_,Y_);  %#ok<UNRCH>
        end
        %%
        % get field curvature model
        if isfield(matchparams,'init_array') % overwrites any initialization with per tile values
            pinit_model = matchparams.init_array(:,:,neigs_row_index);
        elseif isfield(matchparams,'init')
            pinit_model = matchparams.init;
        else
            error('Unable to initialize pinit_model') ;
        end
%         pinit_model
%         [X_e,Y_e,out_e,valid_e] = match.fcestimate(X_,Y_,iadj,matchparams,pinit_model);
        
        [X_, Y_, model_for_this_pair, valid] = match.fcestimate(X_, Y_, iadj, matchparams, pinit_model) ;
        
        %%
        % flip back dimensions
        X_ = util.correctTiles(X_,dims);
        Y_ = util.correctTiles(Y_,dims);
        
        % store pairs
        mout(iadj,:) = model_for_this_pair;
        paired_descriptor_for_this_tile{iadj}.valid = valid;
        paired_descriptor_for_this_tile{iadj}.X = X_;
        paired_descriptor_for_this_tile{iadj}.Y = Y_;
        %R(:,iadj,ineig) = round(median(X_-Y_));
        R_(:,iadj) = round(median(X_-Y_));
    end
    %medianResidualperTile(:,:,neigs_row_index) = R_;
    curvemodel(:,:,neigs_row_index) = mout;
    
    paireddescriptor{neigs_row_index}.onx.valid = paired_descriptor_for_this_tile{1}.valid;
    paireddescriptor{neigs_row_index}.onx.X = paired_descriptor_for_this_tile{1}.X;
    paireddescriptor{neigs_row_index}.onx.Y = paired_descriptor_for_this_tile{1}.Y;
    
    paireddescriptor{neigs_row_index}.ony.valid = paired_descriptor_for_this_tile{2}.valid;
    paireddescriptor{neigs_row_index}.ony.X = paired_descriptor_for_this_tile{2}.X;
    paireddescriptor{neigs_row_index}.ony.Y = paired_descriptor_for_this_tile{2}.Y;
    paireddescriptor{neigs_row_index}.count = [size(paired_descriptor_for_this_tile{1}.X,1) size(paired_descriptor_for_this_tile{2}.X,1)];
    parfor_progress;
end
parfor_progress(0);
