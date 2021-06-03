function [scopeparamsUpdated, curvemodel] = ...
    homographyPerTile6Neighbor(params,neighbors,scopeloc,paireddescriptor,curvemodel)
    %HOMOGRAPHYPERTILE Summary of this function goes here
    % 
    % [OUTPUTARGS] = HOMOGRAPHYPERTILE(INPUTARGS) Explain usage here
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

    % $Author: base $	$Date: 2017/04/26 10:53:18 $	$Revision: 0.1 $
    % Copyright: HHMI 2017
    checkthese = [1 4 5 7]; % 0 - right - bottom - below
    imsize_um = params.imsize_um;
    neigs = neighbors(:,checkthese);%[id -x -y +x +y -z +z] format
    Ntiles = size(neigs,1);
    stgdisp = NaN(Ntiles,3);
    for itile = 1:Ntiles
        idxcent = neigs(itile,1);
        for iadj = 1:3
            idxadj = neigs(itile,1+iadj);
            if isnan(idxadj)
            else
                stgdisp(itile,iadj) = 1000*(scopeloc.loc(idxadj,iadj)-scopeloc.loc(idxcent,iadj));
            end
        end
    end
    %%
    % pix resolution based on curve model
    if 0 % based on overall displacement
        xyz_umperpix_model = [stagedisplacement(:)*ones(1,size(curvemodel,3))]./squeeze(curvemodel(:,3,:));
    else %based on per tile
        xyz_umperpix_model = abs(stgdisp'./squeeze(curvemodel(:,3,:)));
    end

    %%
    params.order = 1;
    [scopeparams,validthis] = util.estimateaffine(paireddescriptor,neighbors,scopeloc,params,curvemodel,0);

    %% affine outlier: TODO, reject based on tile corner pixel shift magnitudes!!
    %mean transformation
    aff = mean(reshape([scopeparams(:).affineglFC],3,3,[]),3);

    reliable = true(Ntiles,1);
    for itile = 1:Ntiles
        if isempty(scopeparams(itile).affineglFC) 
            reliable(itile) = false ;
            disp(sprintf('WARNING: Empty affine estimate @tile: %d [%4.2f %4.2f %4.2f]',itile,(1000*scopeloc.loc(itile,:)))) %#ok<*DSPS>
        elseif norm(scopeparams(itile).affineglFC-aff)/norm(aff)*100>1
            disp(sprintf('WARNING: Bad affine estimate @tile: %d [%4.2f %4.2f %4.2f]',itile,(1000*scopeloc.loc(itile,:))))
            reliable(itile) = false ;
        end
    end
    
    scopeparamsUpdated = scopeparams ;
    if any(reliable) ,
        inliers = find(reliable) ;
        % for every tiles estimate an affine
        anchors = scopeloc.gridix(inliers,1:3);
        queries = scopeloc.gridix(:,1:3);
        IDX = knnsearch(anchors,queries,'K',1,'distance',@distfun);%W=[1 1 100000]
        % fill missing 
        for itile = 1:Ntiles
            ianch = inliers(IDX(itile));
            if itile == ianch; continue;end % skip if ancher is tile itself
            paireddescriptor{itile}.onx = paireddescriptor{ianch}.onx;
            paireddescriptor{itile}.ony = paireddescriptor{ianch}.ony;
            curvemodel(:,:,itile) = curvemodel(:,:,ianch);
            scopeparamsUpdated(itile) = scopeparams(ianch);
        end
    end
end
