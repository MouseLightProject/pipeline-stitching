function [descriptors,paireddescriptor,curvemodel,scopeparams] = tileProcessor(scopeloc,descriptorfolder,desc_ch,params)

checkthese = [1 4 5 7]; % 0 - right - bottom - below
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
% accumulator for steps
[paireddescriptor,medianResidualperTile,curvemodel,scopeparams]=deal([]);
model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model

%% get tile descriptors
descriptors = getDescriptorsPerFolder(descriptorfolder,scopeloc,desc_ch);
sprintf('Loaded descriptors')
%% curvature estimation using desctiptor match
[paireddescriptor{end+1},medianResidualperTile{end+1},curvemodel{end+1}] = match.xymatch(descriptors,neighbors(:,checkthese),scopeloc,params,model);
sprintf('X&Y descriptor match done')

%% interpolate tiles with missing parameters from adjecent tiles
[paireddescriptor{end+1},curvemodel{end+1},unreliable,neigbors_used] = match.curvatureOutlierElimination(paireddescriptor{end},curvemodel{end},scopeloc,params,model);
sprintf('outlier elimination done')

%%
% tile base affine
if params.singleTile
    [scopeparams{1},scopeparams{2},paireddescriptor{end+1},curvemodel{end+1}] = homographyPerTile6Neighbor(...
        params,neighbors,scopeloc,paireddescriptor{end},curvemodel{end},unreliable,neigbors_used);
    sprintf('per-tile affine estimation')
else
    % joint affine estimation
    [scopeparams{end+1}] = match.estimatejointaffine(paireddescriptor{end},neighbors,scopeloc,params,curvemodel{end},0);
    [scopeparams{end+1}, paireddescriptor{end+1}, curvemodel{end+1}] = match.affineOutlierElimination( scopeloc, scopeparams{end}, paireddescriptor{end},curvemodel{end},unreliable );
    sprintf('joint affine estimation')
end


