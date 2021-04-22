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
do_show_visualizations = false ;
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
matchparams.viz = do_show_visualizations;
matchparams.fignum = fignum;
matchparams.opt.beta=2;
% matchparams.opt.method = 'nonrigid';

%%
pixshiftpertile = nan(neigs_row_count,3,2);
for ii=1:neigs_row_count
    for iadj = 1:2
        if isfinite(neigs(ii,iadj+1))
            um_shift = 1000*(scopeloc.loc(neigs(ii,iadj+1),:)-scopeloc.loc(neigs(ii,1),:));
            pixel_shift = round(um_shift.*(dims-1)./(imsize_um));
            pixshiftpertile(ii,:,iadj) = pixel_shift;
        end
    end
end
% if do_show_visualizations ,
%     these = isfinite(pixshiftpertile(:,1,1));onx = pixshiftpertile(these,1,1);medonx = median(onx);  %#ok<UNRCH>
%     these = isfinite(pixshiftpertile(:,2,2));ony = pixshiftpertile(these,2,2);medony = median(ony);
%     pixshiftpertile(isnan(pixshiftpertile(:,1,1)),:,1) = ones(sum(isnan(pixshiftpertile(:,1,1))),1)*[medonx,0,0];
%     pixshiftpertile(isnan(pixshiftpertile(:,2,2)),:,2) = ones(sum(isnan(pixshiftpertile(:,2,2))),1)*[0,medony,0];
%     figure(123), cla,subplot(121),hist(onx,100),hold on,subplot(122),hist(ony,100) %#ok<HIST>
% end
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
%for neigs_row_index = 1:neigs_row_count ,
%for neigs_row_index = round(neigs_row_count/2):neigs_row_count ,
%tile_index_of_interest = find(strcmp('/2020-12-01/01/01916', scopeloc.relativepaths))
%for tile_index = tile_index_of_interest ,
    %% load descriptor pairs X (center) - Y (adjacent tile)
    central_tile_index = neigs(neigs_row_index,1);
    central_tile_flipped_fiducials_and_descriptors = descriptors{central_tile_index};  %#ok<PFBNS>
    if isempty(central_tile_flipped_fiducials_and_descriptors);continue;end
    central_tile_flipped_fiducials = double(central_tile_flipped_fiducials_and_descriptors(:,1:3));
    if size(central_tile_flipped_fiducials,1)<3;continue;end
    
    central_tile_fiducials = util.correctTiles(central_tile_flipped_fiducials, dims) ;  % flip dimensions
       % ALT: I think this is because the fiducials are computed on the raw tile
       % stacks, which have to be flipped in x and y to get them into the proper
       % orientation for the rendered stack
    mout = nan(3,3) ;
    paired_descriptor_for_this_tile = paired_descriptor_for_this_tile_template;
    %R_ = zeros(3); % median residual
    
    %%
    for iadj = 1:2 , %1:x-overlap, 2:y-overlap, 3:z-overlap
        %%
        % idaj : 1=right(+x), 2=bottom(+y), 3=below(+z)
        other_tile_index = neigs(neigs_row_index,iadj+1);  %#ok<PFBNS>
        
        if isnan(other_tile_index);continue;end
        other_tile_raw_fiducials_and_descriptors = descriptors{other_tile_index};
        if isempty(other_tile_raw_fiducials_and_descriptors);continue;end
        
        other_tile_raw_fiducials = double(other_tile_raw_fiducials_and_descriptors(:,1:3)); % descadj has x-y-z-w1-w2 format
        if size(other_tile_raw_fiducials,1)<3;continue;end
        
        other_tile_fiducials = util.correctTiles(other_tile_raw_fiducials,dims);  % flip dimensions
        um_shift = 1000*(scopeloc.loc(other_tile_index,:)-scopeloc.loc(central_tile_index,:)) ;  %#ok<PFBNS>  % um
        pixel_shift = round(um_shift.*(dims-1)./(imsize_um));
        %other_tile_fiducials = other_tile_fiducials + ones(size(other_tile_fiducials,1),1)*pixel_shift ;  % shift with initial guess based on stage coordinate
        other_tile_shifted_fiducials = other_tile_fiducials + pixel_shift ;  % shift with initial guess based on stage coordinate
        
        %%
        nbound = [0 0];
        nbound(1) = max(pixel_shift(iadj),min(other_tile_shifted_fiducials(:,iadj))) - 15;
        nbound(2) = min(dims(iadj),max(central_tile_fiducials(:,iadj)))+0 + 15;
        central_tile_fiducials_near_overlap = central_tile_fiducials(central_tile_fiducials(:,iadj)>nbound(1)&central_tile_fiducials(:,iadj)<nbound(2),:);
        other_tile_shifted_fiducials_near_overlap = other_tile_shifted_fiducials(other_tile_shifted_fiducials(:,iadj)>nbound(1)&other_tile_shifted_fiducials(:,iadj)<nbound(2),:);
        %%
        if size(central_tile_fiducials_near_overlap,1)<3 || size(other_tile_shifted_fiducials_near_overlap,1)<3 ,
            continue
        end
        % get descpair
        [matched_central_tile_fiducials, matched_other_tile_shifted_fiducials] = ...
            match.descriptorMatch4XY(central_tile_fiducials_near_overlap, other_tile_shifted_fiducials_near_overlap, matchparams) ;
        if size(matched_central_tile_fiducials,1)<3 || size(matched_other_tile_shifted_fiducials,1)<3 ,
            continue
        end
        matched_other_tile_fiducials = matched_other_tile_shifted_fiducials ;
        matched_other_tile_fiducials(:,iadj) = matched_other_tile_shifted_fiducials(:,iadj) - pixel_shift(iadj);  % move it back to original location after CDP
        %[X_e,Y_e,out_e,valid_e] = match.fcestimate(matched_central_tile_fiducials,matched_other_tile_fiducials,iadj,matchparams,pinit_model);
        % [X_e2,Y_e2,out_e2,valid_e2] = match.fcestimate(X_e,Y_e,iadj,matchparams);
        %%
        if do_show_visualizations ,
            util.debug.vizMatchALT(scopeloc, ...
                                   neigs, ...
                                   descriptors, ...
                                   neigs_row_index, ...
                                   pixel_shift, ...
                                   iadj, ...
                                   central_tile_fiducials_near_overlap, ...
                                   other_tile_shifted_fiducials_near_overlap, ...
                                   matched_central_tile_fiducials, ...
                                   matched_other_tile_fiducials);  %#ok<UNRCH>
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
%         [X_e,Y_e,out_e,valid_e] = match.fcestimate(matched_central_tile_fiducials,matched_other_tile_fiducials,iadj,matchparams,pinit_model);
        
        [fced_matched_central_tile_fiducials, fced_matched_other_tile_fiducials, model_for_this_pair, is_valid] = ...
            match.fcestimate(matched_central_tile_fiducials, matched_other_tile_fiducials, iadj, matchparams, pinit_model) ;
        
        %%
        % flip back dimensions
        flipped_fced_matched_central_tile_fiducials = util.correctTiles(fced_matched_central_tile_fiducials, dims) ;
        flipped_fced_matched_other_tile_fiducials   = util.correctTiles(fced_matched_other_tile_fiducials  , dims) ;
        
        % store pairs
        mout(iadj,:) = model_for_this_pair;
        paired_descriptor_for_this_tile{iadj}.valid = is_valid;
        paired_descriptor_for_this_tile{iadj}.X = flipped_fced_matched_central_tile_fiducials;
        paired_descriptor_for_this_tile{iadj}.Y = flipped_fced_matched_other_tile_fiducials;
        %R(:,iadj,ineig) = round(median(X_-Y_));
        %R_(:,iadj) = round(median(flipped_fced_matched_central_tile_fiducials-flipped_fced_matched_other_tile_fiducials));
    end
    %medianResidualperTile(:,:,neigs_row_index) = R_;
    curvemodel(:,:,neigs_row_index) = mout ;
    
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
