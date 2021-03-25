function stitch(tile_folder_path, pipeline_output_folder_path, stitching_output_folder_path)
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

    % Define the paths to the different input files
    descoutput = fullfile(pipeline_output_folder_path,'stage_3_descriptor_output') ;
    matchoutput = fullfile(pipeline_output_folder_path,'stage_4_point_match_output') ;

    % Use all the cores
    use_this_fraction_of_cores(1) ;

    % Want group-writable outputs
    system('umask u=rwx,g=rwx,o=rx') ;

    % Create the output folder if it doesn't exist
    if ~exist(stitching_output_folder_path, 'file') ,
        mkdir(stitching_output_folder_path) ;
    end

    scopefile = fullfile(stitching_output_folder_path,'scopeloc.mat');
    descriptorfolder = descoutput;
    matchfolder = matchoutput;

    desc_ch = {'1'};
    descriptorfile = fullfile(stitching_output_folder_path,sprintf('descriptors_ch%s.mat',desc_ch{:})); % accumulated descriptor file


    %% 0: INTIALIZE
    % Read .acquisition file for each tile, and populate scopeloc with the
    % stage positions of each tile.
    if exist(scopefile, 'file') ,  
        load(scopefile, 'scopeloc', 'neighbors', 'experimentfolder') ;
    else
        is_sample_post_2016_04_04 = true ;  % set this to 1 for datasets acquired after 160404
        scopeloc = getScopeCoordinates(stitching_output_folder_path, tile_folder_path, is_sample_post_2016_04_04) ;  % parse from acqusition files
        neighbors = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
        experimentfolder = stitching_output_folder_path ;  % for backwards compatibility
        save(scopefile,'scopeloc','neighbors','experimentfolder')
    end

    %%
    % 1: LOAD MATCHED FEATS
    regpts_file_path = fullfile(stitching_output_folder_path,'regpts.mat') ;
    if exist(regpts_file_path, 'file') ,
        load(regpts_file_path, 'regpts', 'featmap') ;
    else
        %fprintf('Loading matched features stage...\n') ;
        %load(scopefile,'scopeloc');
        directions = 'Z';
        checkversion = 1; % 1: loads the version with "checkversion" extension and overwrites existing match if there are more matched points
        % load finished tile matches. find badly matched or missing tile pairs
        [regpts,featmap] = loadMatchedFeatures(scopeloc,matchfolder,directions,checkversion);    
        save(regpts_file_path, '-v7.3', 'regpts', 'featmap')
    end
    % if ~exist(fullfile(stitching_output_folder_path,'regpts_1stiter.mat'),'file') % faster to make a copy
    %     unix(sprintf('cp %s %s',regpts_file_path,fullfile(stitching_output_folder_path,'regpts_1stiter.mat')))
    % end

    %%
    % 2 scope params estimation.
    % i) finds matches on x&y
    % ii) finds field curvature based on matched points
    % iii) creates a 3D affine model by jointly solving a linear system of
    % equations

    scope_params_per_tile_file_path = fullfile(stitching_output_folder_path,'scopeparams_pertile.mat') ;
    if ~exist(descriptorfile, 'file') || ~exist(scope_params_per_tile_file_path, 'file') ,   
        %%
        fprintf('Running tileProcessor stage...\n') ;
        % paramater setting for descrtiptor match
        scopeacqparams = readScopeFile(fileparts(scopeloc.filepath{1}));
        params.scopeacqparams = scopeacqparams;
        params.imsize_um = [scopeacqparams.x_size_um scopeacqparams.y_size_um scopeacqparams.z_size_um];
        params.overlap_um = [scopeacqparams.x_overlap_um scopeacqparams.y_overlap_um scopeacqparams.z_overlap_um];
        params.imagesize = [1024 1536 251];

        params.viz = 0;
        params.debug = 0;
        params.Ndivs = 4;
        params.Nlayer = 4;
        params.htop = 5;
        params.expensionratio = 1;
        params.order = 1;
        params.applyFC = 1;
        params.beadparams = [];%PLACEHOLDER FOR BEADS, very unlikely to have it...
        params.singleTile = 1;

        [descriptors,paireddescriptor,curvemodel,scopeparams] = ...
            tileProcessor(scopeloc,descriptorfolder,desc_ch,params);
        save(descriptorfile,'descriptors','-v7.3')
        save(scope_params_per_tile_file_path,'paireddescriptor', ...
            'scopeparams', 'curvemodel','params','-v7.3')
    end

    %%
    video_file_path = fullfile(stitching_output_folder_path, sprintf('%s-1stiter-ch1-%s.avi','video',date())) ;
    if ~exist(video_file_path, 'file') ,
        %%
        fprintf('Running descriptorMatchQuality stage...\n') ;
        load(scope_params_per_tile_file_path,'scopeparams')
        descriptorMatchQuality(regpts,scopeparams{end},scopeloc,video_file_path)
        % createThumb(regpts,scopeparams,scopeloc,video_file_path)
        % descriptorMatchQualityHeatMap(regpts,scopeparams{end},scopeloc,video_file_path)
        % descriptorMatchQualityHeatMap_forPaper(regpts,scopeparams{end},scopeloc,video_file_path)
    end

    %%
    vecfield3D_file_path = fullfile(stitching_output_folder_path,'vecfield3D.mat') ;
    if ~exist(vecfield3D_file_path, 'file') ,
        fprintf('Running vectorField3D stage...\n') ;
        %load(scopefile,'scopeloc')
        %load(regpts_file_path,'regpts')
        load(scope_params_per_tile_file_path, 'scopeparams', 'curvemodel', 'params')

        vecfield3D = vectorField3D(params,scopeloc,regpts,scopeparams{end},curvemodel{end},[]);
        save(vecfield3D_file_path,'vecfield3D','params')
        save(fullfile(stitching_output_folder_path,sprintf('%s_%s',datestr(now,'mmddyyHHMMSS'),'vecfield3D')),'vecfield3D','params')
    end

    % Finally, output the yaml file(s)
    %%
    %load(scopefile,'scopeloc') ;
    load(vecfield3D_file_path,'vecfield3D','params') ;
    vecfield = vecfield3D ;

    %%
    % checkthese = [1 4 5 7]; % 0 - right - bottom - below
    % [neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
    params.big = 1;
    params.ymldims = [params.imagesize 2];%[1024 1536 251 2]
    params.root = vecfield.root;

    targetidx = 1:size(scopeloc.gridix,1);
    params.outfile = fullfile(stitching_output_folder_path,sprintf('%s.control.yml',date));
    writeYML(params, targetidx(:)', vecfield);
    unix(sprintf('cp %s %s',params.outfile,fullfile(stitching_output_folder_path,'tilebase.cache.yml'))) ;
    
    params.big=0 ;
    params.outfile = sprintf('%s/%s.old.control.yml',stitching_output_folder_path,date);
    writeYML(params, targetidx(:)', vecfield)
    unix(sprintf('cp %s %s',params.outfile,fullfile(stitching_output_folder_path,'tilebase.cache_old.yml'))) ;

    fprintf('Done with stitching.\n') ;
end
