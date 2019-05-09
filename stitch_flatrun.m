function main_flatrun()
addpath(genpath('./common'))
addpath(genpath('./functions'))
brain = '2017-12-19';
tag='';

%% %%%%%%%%%%%%
directionMap = containers.Map({'-X','-Y','X','Y','-Z','Z'},[ 2, 3, 4, 5, 6, 7]);
directions = 'Z';

%% PATCH, fix this after dm11/tier2 merge
% raw input to descriptor generotion
inputfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling',brain);
experimentfolder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s',brain);
% experimentfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s%s/',brain,tag)

descriptorfolder = fullfile(experimentfolder,'classifier_output');
matfolder = fullfile(experimentfolder,'matfiles/');
scopefile = fullfile(matfolder,'scopeloc.mat');
desc_ch = {'1'};
descriptorfile = fullfile(matfolder,sprintf('descriptors_ch%s.mat',desc_ch{:})); % accumulated descriptor file
matchedfeatfile = fullfile(matfolder,sprintf('feats_ch%s.mat',desc_ch{:})); % accumulated descriptor file

mkdir(matfolder)
unix(sprintf('umask g+rxw %s',matfolder));
unix(sprintf('chmod g+rxw %s',matfolder));
%% 0: INTIALIZE
% read scope files and populate stage coordinates
if 1
    newdash = 1; % set this to 1 for datasets acquired after 160404
    [scopeloc] = getScopeCoordinates(inputfolder,newdash);% parse from acqusition files
    [neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
    save(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
end
%%
if 1
    %%
    load(scopefile,'scopeloc','neighbors','imsize_um','experimentfolder','inputfolder')
    scope = 0
    if scope==1
        scope1_beadparams = load('./beadparams/scope1_beadparams');
        scopeparams = scope1_beadparams.scope1_beadparams;
    elseif scope==2
        scope2_beadparams = load('./beadparams/scope2_beadparams');
        scopeparams = scope2_beadparams.scope2_beadparams;
    else
        scopeparams = [];
    end
    %%
    params.viz = 0;
    params.debug = 0;
    params.Ndivs = 4;
    params.Nlayer = 4;
    params.htop = 5;
    params.imagesize = [1024 1536 251];
    params.expensionratio = 1;
    
    vecfield = vectorField_flatrun(params,scopeloc,scopeparams,0);
end

if 1
    % checkthese = [1 4 5 7]; % 0 - right - bottom - below
    params.big = 1;
    params.ymldims = [params.imagesize 2];%[1024 1536 251 2]
    params.root = vecfield.root;
    
    targetidx = 1:size(scopeloc.gridix,1);
    params.outfile = fullfile(experimentfolder,sprintf('%s-flatrun.control.yml',date));
    
    writeYML(params, targetidx(:)', vecfield)
%     unix(sprintf('cp %s %s',params.outfile,fullfile(experimentfolder,'tilebase.cache.yml')))
end
