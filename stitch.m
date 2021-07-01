function stitch(tile_folder_path, ...
                pipeline_output_folder_path, ...
                stitching_output_folder_path, ...
                stitch_options) 
    %STICHING pipeline. Reads scope generated json file and returns a yml
    %configuration file that goes into renderer. Requires calling cluster jobs
    %to create subresults, i.e. descriptors. These functions can also run in
    %local machines with proper settings.
    %
    % Note that pipeline_output_folder_path points to files that are used as
    % *input* to this function.
    %
    % All outputs are stored in the folder stitching_output_folder_path.

    % NOTES:
    % directionMap = containers.Map({'-X','-Y','X','Y','-Z','Z'},[ 2, 3, 4, 5, 6, 7]);
    % directions = 'Z';

    % $Author: base $	$Date: 2016/09/21 11:52:40 $
    % Copyright: HHMI 2016

    % The recommended way to call this, in most cases, is with three arguments.
    % E.g.
    %     tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/2019-04-17/Tiling' ;
    %     pipeline_output_folder = '/nrs/mouselight/pipeline_output/2019-04-17' ;
    %     stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2019-04-17/stitching-output' ;
    %     stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;

    % Deal with optional arguments
    if ~exist('stitch_options', 'var') || isempty(stitch_options) ,
        stitch_options = struct() ;
    end    
    if ~isfield(stitch_options, 'do_force_computations') || isempty(stitch_options.do_force_computations) ,
        do_force_computations = false ;
    else
        do_force_computations = stitch_options.do_force_computations ;
    end
    if ~isfield(stitch_options, 'do_perform_field_correction') || isempty(stitch_options.do_perform_field_correction) ,
        do_perform_field_correction = true ;
    else
        do_perform_field_correction = stitch_options.do_perform_field_correction ;        
    end
    if ~isfield(stitch_options, 'do_run_in_debug_mode') || isempty(stitch_options.do_run_in_debug_mode) ,
        do_run_in_debug_mode = false ;
    else
        do_run_in_debug_mode = stitch_options.do_run_in_debug_mode ;        
    end
    if ~isfield(stitch_options, 'do_show_visualizations') || isempty(stitch_options.do_show_visualizations) ,
        do_show_visualizations = false ;
    else
        do_show_visualizations = stitch_options.do_show_visualizations ;        
    end
    
    % Define the paths to the different input files
    landmark_folder_path = fullfile(pipeline_output_folder_path,'stage_3_descriptor_output') ;
    match_folder_path = fullfile(pipeline_output_folder_path,'stage_4_point_match_output') ;

    % Use all the cores
    use_this_fraction_of_cores(1) ;

    % Want group-writable outputs
    system('umask u=rwx,g=rwx,o=rx') ;

    % Create the output folder if it doesn't exist
    if ~exist(stitching_output_folder_path, 'file') ,
        mkdir(stitching_output_folder_path) ;
    end

    % Get the path to the file storing the location of each tile according to the
    % scope hardware.
    scopeloc_file_path = fullfile(stitching_output_folder_path, 'scopeloc.mat') ;

    % We use the background channel for finding matches, because the lipofuscin
    % granules show up better on that channel.
    [neuron_channel_index, background_channel_index] = read_channel_semantics_file(tile_folder_path) ;  %#ok<ASGLU>
    desc_ch = { sprintf('%d', background_channel_index) } ;

    % Read .acquisition file for each tile, and populate scopeloc with the
    % stage positions of each tile.
    if exist(scopeloc_file_path, 'file') && ~do_force_computations ,  
        %load(scopefile, 'scopeloc', 'neighbors', 'experimentfolder') ;
        load(scopeloc_file_path, 'scopeloc') ;
    else
        is_sample_post_2016_04_04 = true ;  % set this to 1 for datasets acquired after 160404
        scopeloc = getScopeCoordinates(stitching_output_folder_path, tile_folder_path, is_sample_post_2016_04_04) ;  % parse from acqusition files
        neighbors = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
        experimentfolder = stitching_output_folder_path ;  % for backwards compatibility
        save(scopeloc_file_path,'scopeloc','neighbors','experimentfolder') ;
    end

    % Load z-plane matched landmarks
    regpts_file_path = fullfile(stitching_output_folder_path,'regpts.mat') ;
    if exist(regpts_file_path, 'file')  && ~do_force_computations ,
        load(regpts_file_path, 'regpts') ;
    else
        directions = 'Z';
        % load finished tile matches. find badly matched or missing tile pairs
        [regpts, featmap] = loadMatchedFeaturesNew(scopeloc, match_folder_path, background_channel_index, directions) ;    
        save(regpts_file_path, '-v7.3', 'regpts', 'featmap')
    end

    % Scope field curvature estimation
    % i) Finds matches on x&y
    % ii) Finds field curvature based on matched points
    % iii) Creates a 3D linear transform by jointly solving a linear system of
    %      equations
    scope_params_per_tile_file_path = fullfile(stitching_output_folder_path,'scopeparams_pertile.mat') ;
    if exist(scope_params_per_tile_file_path, 'file') && ~do_force_computations ,   
        load(scope_params_per_tile_file_path, 'scopeparams', 'curvemodel', 'params') ;
    else 
        fprintf('Running field curvature stage...\n') ;
        scopeacqparams = readScopeFile(fileparts(scopeloc.filepath{1}));        
        params = struct() ;
        params.scopeacqparams = scopeacqparams;
        params.imsize_um = [scopeacqparams.x_size_um scopeacqparams.y_size_um scopeacqparams.z_size_um];
        params.overlap_um = [scopeacqparams.x_overlap_um scopeacqparams.y_overlap_um scopeacqparams.z_overlap_um];
        params.imagesize = [1024 1536 251];
        params.do_run_in_debug_mode = do_run_in_debug_mode ;
        params.viz = do_show_visualizations ;
        params.debug = 0;
        params.Ndivs = 4;
        params.Nlayer = 4;
        params.htop = 5;
        params.expensionratio = 1;
        params.order = 1;
        params.applyFC = do_perform_field_correction ;
        params.singleTile = 1;
        [curvemodel, scopeparams] = ...
            estimate_field_curvature_for_all_tiles(scopeloc, landmark_folder_path, desc_ch, params) ;
        fprintf('Done running field curvature stage.\n') ;
        % scopeparams contains the per-tile *linear* transforms
        save(scope_params_per_tile_file_path, 'scopeparams', 'curvemodel', 'params', '-v7.3') ;
    end

    % Make a video
    descriptor_match_quality_video_file_path = fullfile(stitching_output_folder_path, 'descriptor-match-quality.avi') ;
    if ~exist(descriptor_match_quality_video_file_path, 'file') && ~do_force_computations ,
        fprintf('Making descriptorMatchQuality video...\n') ;
        %load(scope_params_per_tile_file_path,'scopeparams')
        descriptorMatchQuality(regpts,scopeparams,scopeloc,descriptor_match_quality_video_file_path)
        % createThumb(regpts,scopeparams,scopeloc,video_file_path)
        % descriptorMatchQualityHeatMap(regpts,scopeparams{end},scopeloc,video_file_path)
        % descriptorMatchQualityHeatMap_forPaper(regpts,scopeparams{end},scopeloc,video_file_path)
        fprintf('Done making descriptorMatchQuality video.\n') ;
    end

%     % Make thumbnails showing tile overlap, posibly ofther info
%     create_thumb_video_file_path = fullfile(stitching_output_folder_path, 'thumb.avi') ;
%     if ~exist(create_thumb_video_file_path, 'file') ,
%         fprintf('Making createThumb video...\n') ;
%         createThumb(regpts,scopeparams,scopeloc,create_thumb_video_file_path)
%         fprintf('Done making createThumb video.\n') ;
%     end

    % Make a heatmap of some sort
    descriptor_match_quality_heatmap_video_file_path = fullfile(stitching_output_folder_path, 'descriptor-match-quality-heatmap.avi') ;
    if ~exist(descriptor_match_quality_heatmap_video_file_path, 'file') && ~do_force_computations ,
        fprintf('Making descriptorMatchQualityHeatMap video...\n') ;
        descriptorMatchQualityHeatMap(regpts,scopeparams,scopeloc,descriptor_match_quality_heatmap_video_file_path) ;
        fprintf('Done making descriptorMatchQualityHeatMap video.\n') ;
    end

    % Compute the 3D vector field
    vecfield3D_file_path = fullfile(stitching_output_folder_path,'vecfield3D.mat') ;
    if exist(vecfield3D_file_path, 'file') && ~do_force_computations ,
        load(vecfield3D_file_path, 'vecfield3D', 'params') ;
    else
        fprintf('Running vectorField3D stage...\n') ;
        vecfield3D = vectorField3D(params, scopeloc, regpts, scopeparams, curvemodel, []) ;
        save(vecfield3D_file_path, 'vecfield3D', 'params') ;
    end

    % Finally, output the yaml file(s)
    big = 1;
    ymldims = [params.imagesize 2];  % [1024 1536 251 2]
    root = vecfield3D.root;
    targetidx = 1:size(scopeloc.gridix,1) ;
    outfile = fullfile(stitching_output_folder_path,sprintf('%s.control.yml',date));
    writeYML(outfile, targetidx, vecfield3D, big, ymldims, root) ;
    system(sprintf('cp %s %s',outfile,fullfile(stitching_output_folder_path,'tilebase.cache.yml'))) ;

    % Declare victory
    fprintf('Done with stitching.\n') ;
end
