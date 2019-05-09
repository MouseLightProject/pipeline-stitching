function stitch(sample_date_or_input_folder_name,pipelineoutputfolder,experimentfolder)
%STICHING pipeline. Reads scope generated json file and returns a yml
%configuration file that goes into renderer. Requires calling cluster jobs
%to create subresults, i.e. descriptors. These functions can also run in
%local machines with proper settings.
%
% [OUTPUTARGS] = STICHING(jsonfile)
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

% NOTES:
% directionMap = containers.Map({'-X','-Y','X','Y','-Z','Z'},[ 2, 3, 4, 5, 6, 7]);
% directions = 'Z';

% $Author: base $	$Date: 2016/09/21 11:52:40 $
% Copyright: HHMI 2016

% I've only ever called this like e.g. "stitch('2019-03-26')"  
%     --ALT, 2019-05-09

%% MAKE SURE PATHS etc are correct
%runfull = true;
if nargin==1
    %brain = '2018-08-01';
    sample_date = sample_date_or_input_folder_name;
    sample_date_or_input_folder_name = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/acquisition/%s',sample_date);
    pipelineoutputfolder = sprintf('/nrs/mouselight/pipeline_output/%s',sample_date);
    %arch = lower(computer('arch'));
    if ispc() ,
        error('windows machine, set the input using input arguments')
    else
        experimentfolder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s-%s', sample_date, getenv('USER'));
    end    
elseif nargin<1
    error('At least pass brain id')
end

addpath(genpath('./common'))
addpath(genpath('./functions'))
% classifierinput = inputfolder;
% raw input to descriptor generotion

piperun = true ;

if piperun
    if isequal(sample_date, '2017-09-25')
        classifierinput = sample_date_or_input_folder_name;
        descoutput ='/nrs/mouselight/cluster/classifierOutputs/2017-09-25/classifier_output'
        matchoutput = descoutput;
    elseif isequal(sample_date, '2018-08-15')
        descoutput = fullfile(pipelineoutputfolder,'stage_2_descriptor_output')
        matchinput = descoutput;
        matchoutput = fullfile(pipelineoutputfolder,'stage_3_point_match_output')
    else
        classifieroutput = fullfile(pipelineoutputfolder,'stage_2_classifier_output')
        descinput = classifieroutput;
        descoutput = fullfile(pipelineoutputfolder,'stage_3_descriptor_output')
        matchinput = descoutput;
        matchoutput = fullfile(pipelineoutputfolder,'stage_4_point_match_output')
    end
end


matfolder = fullfile(experimentfolder,'matfiles/');

if ~exist(experimentfolder, 'file') ,
    mkdir(experimentfolder) ;
end
unix(sprintf('chmod g+rxw %s',experimentfolder));
unix(sprintf('umask g+rxw %s',experimentfolder));
if ~exist(matfolder, 'file') ,
    mkdir(matfolder) ;
end


scopefile = fullfile(matfolder,'scopeloc.mat');
if piperun
    descriptorfolder = descoutput;
    matchfolder = matchoutput;
else
    descriptorfolder = fullfile(experimentfolder,'classifier_output');
    matchfolder = descriptorfolder;
end

desc_ch = {'0'};
descriptorfile = fullfile(matfolder,sprintf('descriptors_ch%s.mat',desc_ch{:})); % accumulated descriptor file
matchedfeatfile = fullfile(matfolder,sprintf('feats_ch%s.mat',desc_ch{:})); % accumulated descriptor file


%% 0: INTIALIZE
% read scope files and populate stage coordinates
if ~exist(scopefile, 'file') ,  
    newdash = 1; % set this to 1 for datasets acquired after 160404
    [scopeloc] = getScopeCoordinates(sample_date_or_input_folder_name,newdash);% parse from acqusition files
    [neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
    save(scopefile,'scopeloc','neighbors','experimentfolder')
end

% %% BULLSHIT CURATION STUFF
% % obsolute after pipeline, TODO: fix missing condition for tile runs
% % rather then channel logs
% if false ,  % this part doesn't seem to work if you give the function a single argument...
%     curationH5(classifierinput,classifieroutput)
%     % checkmissingProb(classifierinput,classifieroutput)
%     checkmissingDesc(descinput,descoutput)
%     checkmissingMatch(matchinput,matchoutput)
% end

%%
% 1: LOAD MATCHED FEATS
regpts_file_path = fullfile(matfolder,'regpts.mat') ;
if ~exist(regpts_file_path, 'file') ,
    fprintf('Loading matched features stage...\n') ;
    load(scopefile,'scopeloc');
    directions = 'Z';
    checkversion = 1; % 1: loads the version with "checkversion" extension and overwrites existing match if there are more matched points
    % load finished tile matches. find badly matched or missing tile pairs
    [regpts,featmap] = loadMatchedFeatures(scopeloc,matchfolder,directions,checkversion);    
    save(regpts_file_path, '-v7.3', 'regpts', 'featmap')
end
if ~exist(fullfile(matfolder,'regpts_1stiter.mat'),'file') % faster to make a copy
    unix(sprintf('cp %s %s',regpts_file_path,fullfile(matfolder,'regpts_1stiter.mat')))
end

% if false ,  % iterate on missing tiles (ANOTHER BULLSHIT)    
%     addpath(genpath('/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe'),'-end')
%     %pointmatch_task(brain,runlocal)
%     directions = 'Z';
%     ch=desc_ch{1};
%     [~,sample] = fileparts(experimentfolder);
%     runlocal=1;
%     pointmatch_task_local(sample,inputfolder,descriptorfolder,matchfolder,matfolder,directions,ch,runlocal)
%     rmpath(genpath('/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe'))
% end

%%
% 2 scope params estimation.
% i) finds matches on x&y
% ii) finds field curvature based on matched points
% iii) creates a 3D affine model by jointly solving a linear system of
% equations

scope_params_per_tile_file_path = fullfile(matfolder,'scopeparams_pertile.mat') ;
if ~exist(descriptorfile, 'file') || ~exist(scope_params_per_tile_file_path, 'file') ,   
    %%
    fprintf('Running tileProcessor stage...\n') ;
    load(scopefile,'scopeloc')
    % paramater setting for descrtiptor match
    scopeacqparams = readScopeFile(fileparts(scopeloc.filepath{1}));
    % params.sample = brain;
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

%     if 0
%         [descriptors,paireddescriptor,curvemodel,scopeparams] = ...
%             tileProcessor_debug(scopeloc,descriptorfolder,desc_ch,params);
%     else
    [descriptors,paireddescriptor,curvemodel,scopeparams] = ...
        tileProcessor(scopeloc,descriptorfolder,desc_ch,params);
    save(descriptorfile,'descriptors','-v7.3')
    save(scope_params_per_tile_file_path,'paireddescriptor', ...
        'scopeparams', 'curvemodel','params','-v7.3')
%     end
end

%%
video_file_path = fullfile(pwd(), 'videos', sprintf('%s-1stiter-ch1-%s.avi',sample_date,date())) ;
if ~exist(video_file_path, 'file') ,
    %%
    fprintf('Running descriptorMatchQuality stage...\n') ;
    load(scopefile,'scopeloc')
    load(scope_params_per_tile_file_path,'scopeparams')
    load(regpts_file_path,'regpts')
    if ~exist('./videos', 'file') ,
        mkdir('./videos')
    end
    descriptorMatchQuality(regpts,scopeparams{end},scopeloc,video_file_path)
    %     createThumb(regpts,scopeparams,scopeloc,video_file_path)
    % descriptorMatchQualityHeatMap(regpts,scopeparams{end},scopeloc,video_file_path)
%     descriptorMatchQualityHeatMap_forPaper(regpts,scopeparams{end},scopeloc,video_file_path)
end

%%
vecfield3D_file_path = fullfile(matfolder,'vecfield3D.mat') ;
if ~exist(vecfield3D_file_path, 'file') ,
    fprintf('Running vectorField3D stage...\n') ;
    load(scopefile,'scopeloc')
    load(regpts_file_path,'regpts')
    load(scope_params_per_tile_file_path, 'scopeparams', 'curvemodel', 'params')
    
    vecfield3D = vectorField3D(params,scopeloc,regpts,scopeparams{end},curvemodel{end},[]);
    save(vecfield3D_file_path,'vecfield3D','params')
    save(fullfile(matfolder,sprintf('%s_%s',datestr(now,'mmddyyHHMMSS'),'vecfield3D')),'vecfield3D','params')
end

% Finally, output the yaml file(s)
%%
% 4
load(scopefile,'scopeloc','experimentfolder')
load(vecfield3D_file_path,'vecfield3D','params')
vecfield = vecfield3D;

%%
% checkthese = [1 4 5 7]; % 0 - right - bottom - below
% [neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
params.big = 1;
params.ymldims = [params.imagesize 2];%[1024 1536 251 2]
sub = false ;
params.root = vecfield.root;

if sub
    targetidx = getTargetIDx(scopeloc,neighbors);
    copytiles2target('./test_copt',scopeloc,targetidx(1))
    params.outfile = fullfile(experimentfolder,sprintf('%s_sub.control.yml',date));
else
    targetidx = 1:size(scopeloc.gridix,1);
    params.outfile = fullfile(experimentfolder,sprintf('%s.control.yml',date));
end
writeYML(params, targetidx(:)', vecfield);
unix(sprintf('cp %s %s',params.outfile,fullfile(experimentfolder,'tilebase.cache.yml')));
%
if ~sub
    params.big=0 ;
    params.outfile = sprintf('%s/%s.old.control.yml',experimentfolder,date);
    writeYML(params, targetidx(:)', vecfield)
    unix(sprintf('cp %s %s',params.outfile,fullfile(experimentfolder,'tilebase.cache_old.yml')))
end
return


% %% 0.1: FLAT RUN
% % generate yml for without any optimization
% if false
%     %%
%     load(scopefile,'scopeloc','neighbors','imsize_um','experimentfolder','inputfolder')
%     if scope==1
%         scope1_beadparams = load('./beadparams/scope1_beadparams');
%         scopeparams = scope1_beadparams.scope1_beadparams;
%     else
%         scope2_beadparams = load('./beadparams/scope2_beadparams');
%         scopeparams = scope2_beadparams.scope2_beadparams;
%     end
%     vecfield = vectorField_flatrun(params,scopeloc,scopeparams,2);
%     
%     load ./matfiles/xypaireddescriptor paireddescriptor R curvemodel
%     [scopeparams,scopeparams_,paireddescriptor_,curvemodel_] = homographyPerTile6Neighbor(...
%         beadparams,neighbors,scopeloc,paireddescriptor,R,curvemodel,imsize_um);
%     vecfield3D_flat_4neig = vectorField_flatrun_pertile(params,scopeloc,scopeparams_,curvemodel_,[]);
%     save pertile_4neig scopeparams scopeparams_ paireddescriptor_ curvemodel_ vecfield3D_flat_4neig
% end
% %% stitching quality test
% if true
%     fprintf('Running stitching quality test...\n') ;
%     load(fullfile(matfolder,'scopeloc'),'scopeloc','imsize_um','experimentfolder','inputfolder')
%     load(fullfile(matfolder,'vecfield'),'vecfield','params')
%     %%
%     clc
%     params.big = 1;
%     params.dims = [params.imagesize 2]%[1024 1536 251 2]
%     sub = 0;
%     inds_ = inds(1)';
%     neigs = neighbors(inds_,checkthese);
%     targetidx = neigs([1 3])
%     params.root = vecfield.root;
%     %%
%     %
%     params.outfile = sprintf('%s%s-%d_%d_sub_1.tilebase.cache.yml',experimentfolder,date,targetidx);
%     params.outfile
%     vecfield_ = vecfield;
%     vecfield_.path{targetidx(1)} = '/00000';
%     writeYML(params, targetidx(:)', vecfield_)
%     paramoutput = '/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/2019-03-26/set_parameters_sub'
%     converter(params.outfile,paramoutput,'y1-ccx2-fixxed')
%     %%
%     
%     params.outfile = sprintf('%s%s-%d_%d_sub_2.tilebase.cache.yml',experimentfolder,date,targetidx);
%     params.outfile
%     vecfield_ = vecfield;
%     vecfield_.path{targetidx(2)} = '/00000';
%     writeYML(params, targetidx(:)', vecfield_)
%     paramoutput = '/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/2019-03-26/set_parameters_sub'
%     converter(params.outfile,paramoutput,'y2-ccx2-fixxed')
% end
end
