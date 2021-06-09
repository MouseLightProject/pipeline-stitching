function [paired_descriptor_from_tile_index, curve_model_from_axis_index_from_tile_index] = ...
        xymatch(descriptors, tile_index_from_axis_index_from_tile_index, scopeloc, params, model)
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
    tile_count = size(tile_index_from_axis_index_from_tile_index,1) ;
    if ~isequal(tile_index_from_axis_index_from_tile_index(:,1), (1:tile_count)') ,
        error('Something has gone terribly wrong, error code 1138') ;
    end
    if nargin<5
        model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
    end
    
    debug = 0;
    %res = 0;
    do_run_in_debug_mode = params.do_run_in_debug_mode ;
    do_show_visualizations = logical(params.viz) ;
    fignum = 101;
    
    dims = params.imagesize;
    imsize_um = params.imsize_um;
    
    % slid = [[75 960];[0 dims(2)];[0 dims(3)]];
    % expensionshift = [0 0 20]; % HEURISTICS:: tissue expends, so overlap is bigger between tiles
    
    %%
    optimopts = statset('nlinfit');
    optimopts.RobustWgtFun = 'bisquare';
    
    cpd_options = struct() ;
    cpd_options.method='nonrigid';
    % opt.method='nonrigid_lowrank';
    cpd_options.beta=2;            % the width of Gaussian kernel (smoothness), higher numbers make transformation more stiff
    cpd_options.lambda=16;          % regularization weight
    cpd_options.viz=0;              % show every iteration
    cpd_options.outliers=0.9;       % use 0.7 noise weight
    cpd_options.fgt=0;              % do not use FGT (default)
    cpd_options.normalize=1;        % normalize to unit variance and zero mean before registering (default)
    cpd_options.corresp=1;          % compute correspondence vector at the end of registration (not being estimated by default)
    %     opt.max_it=100;         % max number of iterations
    %     opt.tol=1e-10;          % tolerance
    
    matchparams = struct() ;
    matchparams.model = model;
    matchparams.optimopts = optimopts;
    matchparams.opt = cpd_options ;
    matchparams.initial_distance_threshold = 40 ;  % distance between target and projected point has to be less than this number
    matchparams.registered_distance_threshold = 20 ;  % distance between target and projected point has to be less than this number
    matchparams.debug = debug;
    matchparams.viz = do_show_visualizations;
    matchparams.fignum = fignum;
    % matchparams.opt.method = 'nonrigid';
    
    %%
    pixshiftpertile = nan(tile_count,3,2);
    for ii=1:tile_count
        for axis_index = 1:2
            if isfinite(tile_index_from_axis_index_from_tile_index(ii,axis_index+1))
                um_shift_xyz = ...
                    1000*(scopeloc.loc(tile_index_from_axis_index_from_tile_index(ii,axis_index+1),:) - ...
                    scopeloc.loc(tile_index_from_axis_index_from_tile_index(ii,1),:)) ;
                pixel_shift_ijk = round(um_shift_xyz.*(dims-1)./(imsize_um)) ;
                pixshiftpertile(ii,:,axis_index) = pixel_shift_ijk ;
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
    for axis_index = 1:2
        % replace any nan rows with median
        these = isnan(pixshiftpertile(:,1,axis_index));
        pixshiftpertile(:,:,axis_index) = util.rowreplace(pixshiftpertile(:,:,axis_index),these,meds(axis_index,:));
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
        for it = 1:tile_count
            for axis_index = 1:2
                matchparams.init_array(axis_index,:,it)=[pvals_12(axis_index,:) pixshiftpertile(it,axis_index,axis_index)];
            end
        end
    end
    
    %% INITIALIZATION
    % checkthese = [1 4 5 7]; % 0 - right - bottom - below
    % indicies are 1 based,e.g. x = 1:dims(1), not 0:dims(1)-1
    % xyz_umperpix = zeros(tile_count,3);
    curve_model_from_axis_index_from_tile_index = nan(3,3,tile_count);
    %medianResidualperTile = zeros(3,3,neigs_row_count);
    paired_descriptor_from_tile_index = cell(tile_count,1);
    % initialize
    for ix = 1:tile_count ,
        paired_descriptor_from_tile_index{ix}.onx.valid = 0;
        paired_descriptor_from_tile_index{ix}.onx.X = [];
        paired_descriptor_from_tile_index{ix}.onx.Y = [];
        paired_descriptor_from_tile_index{ix}.ony.valid = 0;
        paired_descriptor_from_tile_index{ix}.ony.X = [];
        paired_descriptor_from_tile_index{ix}.ony.Y = [];
        paired_descriptor_from_tile_index{ix}.neigs = tile_index_from_axis_index_from_tile_index(ix,:);
        paired_descriptor_from_tile_index{ix}.count = [0 0];
    end
   
    
    %interiorTile_list = util.interior_tiles(scopeloc,1);
    
    %%
    try parfor_progress(0);catch;end
    parfor_progress(tile_count) ;

    if do_run_in_debug_mode ,
        for tile_index = 1:tile_count ,
            [did_succeed, curve_model_from_axis_index, paired_descriptor_from_axis_index] = ...
                xy_match_for_single_tile(tile_index, ...
                                         tile_index_from_axis_index_from_tile_index, ...
                                         descriptors, ...
                                         dims, ...
                                         scopeloc, ...
                                         imsize_um, ...
                                         matchparams, ...
                                         do_show_visualizations) ;
            if did_succeed ,
                curve_model_from_axis_index_from_tile_index(:,:,tile_index) = curve_model_from_axis_index ;

                paired_descriptor_from_tile_index{tile_index}.onx.valid = paired_descriptor_from_axis_index{1}.valid;
                paired_descriptor_from_tile_index{tile_index}.onx.X = paired_descriptor_from_axis_index{1}.X;
                paired_descriptor_from_tile_index{tile_index}.onx.Y = paired_descriptor_from_axis_index{1}.Y;

                paired_descriptor_from_tile_index{tile_index}.ony.valid = paired_descriptor_from_axis_index{2}.valid;
                paired_descriptor_from_tile_index{tile_index}.ony.X = paired_descriptor_from_axis_index{2}.X;
                paired_descriptor_from_tile_index{tile_index}.ony.Y = paired_descriptor_from_axis_index{2}.Y;

                paired_descriptor_from_tile_index{tile_index}.count = ...
                    [size(paired_descriptor_from_axis_index{1}.X,1) size(paired_descriptor_from_axis_index{2}.X,1)];
            end
            
            parfor_progress;
        end
    else
        parfor tile_index = 1:tile_count ,
            [did_succeed, curve_model_from_axis_index, paired_descriptor_from_axis_index] = ...
                xy_match_for_single_tile(tile_index, ...
                                         tile_index_from_axis_index_from_tile_index, ...
                                         descriptors, ...
                                         dims, ...
                                         scopeloc, ...
                                         imsize_um, ...
                                         matchparams, ...
                                         do_show_visualizations) ;
            if did_succeed ,
                curve_model_from_axis_index_from_tile_index(:,:,tile_index) = curve_model_from_axis_index ;

                paired_descriptor_from_tile_index{tile_index}.onx.valid = paired_descriptor_from_axis_index{1}.valid;
                paired_descriptor_from_tile_index{tile_index}.onx.X = paired_descriptor_from_axis_index{1}.X;
                paired_descriptor_from_tile_index{tile_index}.onx.Y = paired_descriptor_from_axis_index{1}.Y;

                paired_descriptor_from_tile_index{tile_index}.ony.valid = paired_descriptor_from_axis_index{2}.valid;
                paired_descriptor_from_tile_index{tile_index}.ony.X = paired_descriptor_from_axis_index{2}.X;
                paired_descriptor_from_tile_index{tile_index}.ony.Y = paired_descriptor_from_axis_index{2}.Y;

                paired_descriptor_from_tile_index{tile_index}.count = ...
                    [size(paired_descriptor_from_axis_index{1}.X,1) size(paired_descriptor_from_axis_index{2}.X,1)];
            end
            
            parfor_progress;
        end
    end
    parfor_progress(0);
end



function [did_succeed, curve_model_from_axis_index, paired_descriptor_from_axis_index] = ...
        xy_match_for_single_tile(tile_index, ...
                                 tile_index_from_axis_index_from_tile_index, ...
                                 descriptors, ...
                                 dims, ...
                                 scopeloc, ...
                                 imsize_um, ...
                                 matchparams, ...
                                 do_show_visualizations)
                             
    % Make a template we'll use to initialize paired_descriptor_for_this_tile
    % for each tile
    paired_descriptor_for_this_tile_template=[];
    paired_descriptor_for_this_tile_template{1}.valid = 0;
    paired_descriptor_for_this_tile_template{1}.X = [];
    paired_descriptor_for_this_tile_template{1}.Y = [];
    paired_descriptor_for_this_tile_template{2}.valid = 0;
    paired_descriptor_for_this_tile_template{2}.X = [];
    paired_descriptor_for_this_tile_template{2}.Y = [];
    
    %% load descriptor pairs X (center) - Y (adjacent tile)
    central_tile_index = tile_index_from_axis_index_from_tile_index(tile_index,1);
    central_tile_flipped_fiducials_and_descriptors = descriptors{central_tile_index}; 
    if isempty(central_tile_flipped_fiducials_and_descriptors) ,
        did_succeed = false ;
        curve_model_from_axis_index = [] ;
        paired_descriptor_from_axis_index = [] ;
        return
    end
    central_tile_flipped_fiducials = double(central_tile_flipped_fiducials_and_descriptors(:,1:3));
    if size(central_tile_flipped_fiducials,1) < 3 ,
        did_succeed = false ;
        curve_model_from_axis_index = [] ;
        paired_descriptor_from_axis_index = [] ;
        return
    end
    
    central_tile_fiducials = util.correctTiles(central_tile_flipped_fiducials, dims) ;  % flip dimensions
       % ALT: I think this is because the fiducials are computed on the raw tile
       % stacks, which have to be flipped in x and y to get them into the proper
       % orientation for the rendered stack
    curve_model_from_axis_index = nan(3,3) ;
    paired_descriptor_from_axis_index = paired_descriptor_for_this_tile_template;
    %R_ = zeros(3); % median residual
    
    %%
    for axis_index = 1:2 , %1:x-overlap, 2:y-overlap, 3:z-overlap
        %axis_index
        %%
        % idaj : 1=right(+x), 2=bottom(+y), 3=below(+z)
        other_tile_index = tile_index_from_axis_index_from_tile_index(tile_index,axis_index+1); 
        
        if isnan(other_tile_index);continue;end
        other_tile_raw_fiducials_and_descriptors = descriptors{other_tile_index};
        if isempty(other_tile_raw_fiducials_and_descriptors);continue;end
        
        other_tile_raw_fiducials = double(other_tile_raw_fiducials_and_descriptors(:,1:3)); % descadj has x-y-z-w1-w2 format
        if size(other_tile_raw_fiducials,1)<3;continue;end
        
        other_tile_fiducials = util.correctTiles(other_tile_raw_fiducials,dims);  % flip dimensions
        um_shift_xyz = 1000*(scopeloc.loc(other_tile_index,:)-scopeloc.loc(central_tile_index,:)) ;    % um
        pixel_shift_ijk = round(um_shift_xyz.*(dims-1)./(imsize_um)) ;
        %other_tile_fiducials = other_tile_fiducials + ones(size(other_tile_fiducials,1),1)*pixel_shift ;  % shift with initial guess based on stage coordinate
        other_tile_shifted_fiducials = other_tile_fiducials + pixel_shift_ijk ;  % shift with initial guess based on stage coordinate
        
        %%
        % We use "i" below to denote i, j, or k, depending on axis_index
        i_lower_bound = max(pixel_shift_ijk(axis_index),min(other_tile_shifted_fiducials(:,axis_index))) - 15;
        i_upper_bound = min(dims(axis_index),max(central_tile_fiducials(:,axis_index))) + 15 ;
        is_central_tile_fiducial_near_overlap = ...
            i_lower_bound<central_tile_fiducials(:,axis_index) & central_tile_fiducials(:,axis_index)<i_upper_bound ;
        central_tile_fiducials_near_overlap = ...
            central_tile_fiducials(is_central_tile_fiducial_near_overlap,:) ;
        is_other_tile_shifted_fiducial_near_overlap = ...
            i_lower_bound<other_tile_shifted_fiducials(:,axis_index) & other_tile_shifted_fiducials(:,axis_index)<i_upper_bound ;
        other_tile_shifted_fiducials_near_overlap = ...
            other_tile_shifted_fiducials(is_other_tile_shifted_fiducial_near_overlap,:) ;
        
        %%
        if size(central_tile_fiducials_near_overlap,1)<3 || size(other_tile_shifted_fiducials_near_overlap,1)<3 ,
            continue
        end
        
        % Run the matching algorithm
        [matched_central_tile_fiducials, matched_other_tile_shifted_fiducials] = ...
            match.descriptorMatch4XY(central_tile_fiducials_near_overlap, other_tile_shifted_fiducials_near_overlap, matchparams) ;
        match_count = size(matched_central_tile_fiducials,1) ; %#ok<NASGU>
        if size(matched_central_tile_fiducials,1)<3 || size(matched_other_tile_shifted_fiducials,1)<3 ,
            continue
        end
        matched_other_tile_fiducials = matched_other_tile_shifted_fiducials ;
        matched_other_tile_fiducials(:,axis_index) = ...
            matched_other_tile_shifted_fiducials(:,axis_index) - pixel_shift_ijk(axis_index);  % move it back to original location after CDP
        
        %%
        if do_show_visualizations ,
            util.debug.vizMatchALT(scopeloc, ...
                                   tile_index_from_axis_index_from_tile_index, ...
                                   descriptors, ...
                                   tile_index, ...
                                   pixel_shift_ijk, ...
                                   axis_index, ...
                                   central_tile_fiducials_near_overlap, ...
                                   other_tile_shifted_fiducials_near_overlap, ...
                                   matched_central_tile_fiducials, ...
                                   matched_other_tile_fiducials);  %#ok<UNRCH>
        end
        
        %%
        % get field curvature model
        if isfield(matchparams,'init_array') % overwrites any initialization with per tile values
            pinit_model = matchparams.init_array(:,:,tile_index);
        elseif isfield(matchparams,'init')
            pinit_model = matchparams.init;
        else
            error('Unable to initialize pinit_model') ;
        end        
        [fced_matched_central_tile_fiducials, fced_matched_other_tile_fiducials, curve_model, is_valid] = ...
            match.fcestimate(matched_central_tile_fiducials, matched_other_tile_fiducials, axis_index, matchparams, pinit_model, dims) ;
        
        %%
        % flip back dimensions
        flipped_fced_matched_central_tile_fiducials = util.correctTiles(fced_matched_central_tile_fiducials, dims) ;
        flipped_fced_matched_other_tile_fiducials   = util.correctTiles(fced_matched_other_tile_fiducials  , dims) ;
        
        % store pairs
        curve_model_from_axis_index(axis_index,:) = curve_model ;
        paired_descriptor_from_axis_index{axis_index}.valid = is_valid ;
        paired_descriptor_from_axis_index{axis_index}.X = flipped_fced_matched_central_tile_fiducials ;
        paired_descriptor_from_axis_index{axis_index}.Y = flipped_fced_matched_other_tile_fiducials ;
        %R(:,iadj,ineig) = round(median(X_-Y_));
        %R_(:,iadj) = round(median(flipped_fced_matched_central_tile_fiducials-flipped_fced_matched_other_tile_fiducials));
    end    
    did_succeed = true ;
end

