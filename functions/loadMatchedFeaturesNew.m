function [regpts,featmap] = loadMatchedFeaturesNew(scopeloc, descriptorfolder, background_channel_index, directions, featmap)
% get list file
numTiles = size(scopeloc.loc,1);

if ~exist('directions', 'var') || isempty(directions) ,
    directions='Z'; % only check z
end
if ~exist('featmap', 'var') || isempty(featmap) ,
    for idx=1:numTiles
        [featmap(idx).X,featmap(idx).Y,featmap(idx).Z] = deal([]);
    end
end

%%
tmp=cell(1,numTiles);
parfor idx = 1:numTiles
    match_file_name = sprintf('channel-%d-match-%s.mat', background_channel_index, directions) ;
    match_file_path = fullfile(descriptorfolder,scopeloc.relativepaths{idx},match_file_name);
    if exist(match_file_path,'file')
        % load descriptors
        tmp{idx} = load(match_file_path);
    else
        fprintf('Warning: Missing match file %s\n', match_file_path) ;
    end
%     if checkversion
%         matchfile2 = fullfile(descriptorfolder,scopeloc.relativepaths{idx},sprintf('match-%s-%d.mat',directions,checkversion));
%         if exist(matchfile2,'file')
%             % load descriptors
%             tmp2 = load(matchfile2);
%             if size(tmp{idx}.paireddescriptor.X,1)<size(tmp2.paireddescriptor.X,1)
%                 % overwrite
%                 tmp{idx} = tmp2;
%             end
%         end
%     end
end
%%
for idx = 1:numTiles
    featmap(idx).(genvarname(directions)) = tmp{idx};
end

% legacy variable
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
regpts = cell(1,length(featmap));for ii=1:length(regpts);if ~isfield(regpts{ii},'X');regpts{ii}.X=[];regpts{ii}.Y=[];regpts{ii}.matchrate=0;end;end
for ii=1:length(regpts)
    if isempty(featmap(ii).Z);continue;end
    regpts{ii} = featmap(ii).Z.paireddescriptor;
    regpts{ii}.neigs = [ii neighbors(ii,[4 5 7])];
end

