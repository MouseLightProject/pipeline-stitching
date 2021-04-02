function [curvemodel, scopeparams] = tileProcessor(scopeloc,descriptorfolder,desc_ch,params)

    checkthese = [1 4 5 7]; % 0 - right - bottom - below
    neighbors = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format    
    % accumulator for steps
    model = @(p,y) p(3) - p(2).*((y-p(1)).^2) ;  % FC model

    % get tile descriptors
    descriptors = getDescriptorsPerFolder(descriptorfolder, scopeloc, desc_ch) ;
    sprintf('Loaded descriptors')
    
    % curvature estimation using descriptor match
    [paireddescriptor_from_xy_match, curvemodel_from_xy_match] = match.xymatch(descriptors, neighbors(:,checkthese), scopeloc, params, model) ;
    sprintf('X&Y descriptor match done')

    % interpolate tiles with missing parameters from adjecent tiles
    [paireddescriptor_from_COE, curvemodel_from_COE, unreliable, neigbors_used] = ...
        match.curvatureOutlierElimination(paireddescriptor_from_xy_match, curvemodel_from_xy_match, scopeloc, params, model) ;
    sprintf('outlier elimination done')

    % tile base affine
    if params.singleTile ,
        [scopeparams, curvemodel] = ...
            homographyPerTile6Neighbor(params, neighbors, scopeloc, paireddescriptor_from_COE, curvemodel_from_COE, unreliable, neigbors_used) ;
        sprintf('per-tile affine estimation')
    else
        % joint affine estimation
        scopeparams_from_EJA = ...
            match.estimatejointaffine(paireddescriptor_from_COE, neighbors, scopeloc, params, curvemodel_from_COE, 0) ;
        [scopeparams, curvemodel] = ...
            match.affineOutlierElimination( scopeloc, scopeparams_from_EJA, paireddescriptor_from_COE, curvemodel_from_COE, unreliable) ;
        sprintf('joint affine estimation')
    end

%     % Collect results
%     curvemodel = {curvemodel1, curvemodel2, curvemodel3} ;
%     scopeparams = {scopeparams1, scopeparams2} ;
end
