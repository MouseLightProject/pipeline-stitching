function [outputArgs] = main(inputfolder,pipelineoutputfolder,experimentfolder)
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

%% MAKE SURE PATHS etc are correct
runfull = false;
if nargin==1
    %brain = '2018-08-01';
    brain = inputfolder;
    inputfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/acquisition/%s',brain);
    pipelineoutputfolder = sprintf('/nrs/mouselight/pipeline_output/%s',brain);
    arch = lower(computer('arch'));
    if arch(1:2) == 'pc'
        error('windows machine, set the input using input arguments')
    else
        experimentfolder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s-%s',brain,getenv('USER'));
    end
    
elseif nargin<1
    error('At least pass brain id')
end

addpath(genpath('./common'))
addpath(genpath('./functions'))
% classifierinput = inputfolder;
% raw input to descriptor generotion

piperun = 1;

if piperun
    if brain=='2017-09-25'
        classifierinput = inputfolder;
        descoutput ='/nrs/mouselight/cluster/classifierOutputs/2017-09-25/classifier_output'
        matchoutput = descoutput;
    elseif brain=='2018-08-15'
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

mkdir(experimentfolder)
unix(sprintf('chmod g+rxw %s',experimentfolder));
unix(sprintf('umask g+rxw %s',experimentfolder));
mkdir(matfolder)


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
if runfull & 0
    newdash = 1; % set this to 1 for datasets acquired after 160404
    [scopeloc] = getScopeCoordinates(inputfolder,newdash);% parse from acqusition files
    [neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
    save(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
end

%% BULLSHIT CURATION STUFF
% obsolute after pipeline, TODO: fix missing condition for tile runs
% rather then channel logs
if 0
    curationH5(classifierinput,classifieroutput)
    % checkmissingProb(classifierinput,classifieroutput)
    checkmissingDesc(descinput,descoutput)
    checkmissingMatch(matchinput,matchoutput)
end

%%
% 1: LOAD MATCHED FEATS
if runfull
    load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder');
    directions = 'Z';
    checkversion = 1; % 1: loads the version with "checkversion" extension and overwrites existing match if there are more matched points
    % load finished tile matches. find badly matched or missing tile pairs
    [regpts,featmap] = loadMatchedFeatures(scopeloc,matchfolder,directions,checkversion);
    
    save(fullfile(matfolder,'regpts.mat'),'regpts','featmap')
    if ~exist(fullfile(matfolder,'regpts_1stiter.mat'),'file') % faster to make a copy
        unix(sprintf('cp %s %s',fullfile(matfolder,'regpts.mat'),fullfile(matfolder,'regpts_1stiter.mat')))
    end
end

if 0 % iterate on missing tiles (ANOTHER BULLSHIT)
    
    addpath(genpath('/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe'),'-end')
    %pointmatch_task(brain,runlocal)
    directions = 'Z';
    ch=desc_ch{1};
    [~,sample] = fileparts(experimentfolder);
    runlocal=1;
    pointmatch_task_local(sample,inputfolder,descriptorfolder,matchfolder,matfolder,directions,ch,runlocal)
    rmpath(genpath('/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe'))
end

%%
% 2 scope params estimation.
% i) finds matches on x&y
% ii) finds field curvature based on matched points
% iii) creates a 3D affine model by jointly solving a linear system of
% equations

if runfull
    
    %%
    load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
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

    if 0
        [descriptors,paireddescriptor,curvemodel,scopeparams] = ...
            tileProcessor_debug(scopeloc,descriptorfolder,desc_ch,params);
    else
        [descriptors,paireddescriptor,curvemodel,scopeparams] = ...
            tileProcessor(scopeloc,descriptorfolder,desc_ch,params);
        save(descriptorfile,'descriptors','-v7.3')
        save(fullfile(matfolder,'scopeparams_pertile'),'paireddescriptor', ...
            'scopeparams', 'curvemodel','params','-v7.3')
    end
end

%%
if runfull
    %%
    load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
    load(fullfile(matfolder,'scopeparams_pertile'),'scopeparams')
    load(fullfile(matfolder,'regpts'),'regpts')
    mkdir('./videos')
    videofile = sprintf('./videos/%s-1stiter-ch1-%s',brain,date)
    descriptorMatchQuality(regpts,scopeparams{end},scopeloc,videofile)
    %     createThumb(regpts,scopeparams,scopeloc,videofile)
    % descriptorMatchQualityHeatMap(regpts,scopeparams{end},scopeloc,videofile)
%     descriptorMatchQualityHeatMap_forPaper(regpts,scopeparams{end},scopeloc,videofile)
end

%%
if runfull
    load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
    load(fullfile(matfolder,'regpts'),'regpts')
    load(fullfile(matfolder,'scopeparams_pertile'),'paireddescriptor', ...
        'scopeparams', 'curvemodel','params')
    
    vecfield3D = vectorField3D(params,scopeloc,regpts,scopeparams{end},curvemodel{end},[]);
    if 1
        save(fullfile(matfolder,sprintf('%s_%s',datestr(now,'mmddyyHHMMSS'),'vecfield3D')),'vecfield3D','params')
        save(fullfile(matfolder,'vecfield3D'),'vecfield3D','params')
    end
end
%%
% 4
load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
load(fullfile(matfolder,'vecfield3D'),'vecfield3D','params')
vecfield = vecfield3D;

%%
% checkthese = [1 4 5 7]; % 0 - right - bottom - below
% [neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
params.big = 1;
params.ymldims = [params.imagesize 2];%[1024 1536 251 2]
sub = 0;
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
    params.big=0
    params.outfile = sprintf('%s/%s.old.control.yml',experimentfolder,date);
    writeYML(params, targetidx(:)', vecfield)
    unix(sprintf('cp %s %s',params.outfile,fullfile(experimentfolder,'tilebase.cache_old.yml')))
end
return
%% 0.1: FLAT RUN
% generate yml for without any optimization
if 0
    %%
    load(scopefile,'scopeloc','neighbors','imsize_um','experimentfolder','inputfolder')
    if scope==1
        scope1_beadparams = load('./beadparams/scope1_beadparams');
        scopeparams = scope1_beadparams.scope1_beadparams;
    else
        scope2_beadparams = load('./beadparams/scope2_beadparams');
        scopeparams = scope2_beadparams.scope2_beadparams;
    end
    vecfield = vectorField_flatrun(params,scopeloc,scopeparams,2);
    
    load ./matfiles/xypaireddescriptor paireddescriptor R curvemodel
    [scopeparams,scopeparams_,paireddescriptor_,curvemodel_] = homographyPerTile6Neighbor(...
        beadparams,neighbors,scopeloc,paireddescriptor,R,curvemodel,imsize_um);
    vecfield3D_flat_4neig = vectorField_flatrun_pertile(params,scopeloc,scopeparams_,curvemodel_,[]);
    save pertile_4neig scopeparams scopeparams_ paireddescriptor_ curvemodel_ vecfield3D_flat_4neig
end
%% stitching quality test
if 0
    load(fullfile(matfolder,'scopeloc'),'scopeloc','imsize_um','experimentfolder','inputfolder')
    load(fullfile(matfolder,'vecfield'),'vecfield','params')
    %%
    clc
    params.big = 1;
    params.dims = [params.imagesize 2]%[1024 1536 251 2]
    sub = 0;
    inds_ = inds(1)';
    neigs = neighbors(inds_,checkthese);
    targetidx = neigs([1 3])
    params.root = vecfield.root;
    %%
    %
    params.outfile = sprintf('%s%s-%d_%d_sub_1.tilebase.cache.yml',experimentfolder,date,targetidx);
    params.outfile
    vecfield_ = vecfield;
    vecfield_.path{targetidx(1)} = '/00000';
    writeYML(params, targetidx(:)', vecfield_)
    paramoutput = '/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/2017-02-13/set_parameters_sub'
    converter(params.outfile,paramoutput,'y1-ccx2-fixxed')
    %%
    
    params.outfile = sprintf('%s%s-%d_%d_sub_2.tilebase.cache.yml',experimentfolder,date,targetidx);
    params.outfile
    vecfield_ = vecfield;
    vecfield_.path{targetidx(2)} = '/00000';
    writeYML(params, targetidx(:)', vecfield_)
    paramoutput = '/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/2017-02-13/set_parameters_sub'
    converter(params.outfile,paramoutput,'y2-ccx2-fixxed')
end
end
