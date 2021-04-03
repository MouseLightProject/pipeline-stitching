function [scopeparams,validthis] = estimateaffine(paireddescriptor,neighbors,scopeloc,params,curvemodel,old)
%ESTIMATEAFFINE Summary of this function goes here
%
% [OUTPUTARGS] = ESTIMATEAFFINE(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2017/04/26 14:29:14 $	$Revision: 0.1 $
% Copyright: HHMI 2017
%%
if nargin<6
    old = 1;
end
%beadreg = 0;
checkthese = [1 4 5 7]; % 0 - right - bottom - below
neigs4 = neighbors(:,[1 2 3 4 5]);% left - above - right - bottom
neigs = neighbors(:,checkthese);%[id -x -y +x +y -z +z] format
% beadparams_=beadparams;
dims = params.imagesize;
order = params.order; % power to weight shift
imsize_um = params.imsize_um;

% initialize empty scope params
% beadparamsZmatch_X = [];
% beadparamsZmatch_Y = [];
% beadparamsZmatchdispstage=[];
% TODO: fill below once we get bead images for each scope/objective
% scope = params.scope;
% beadparamsZmatch_X = beadparams.allX{3};
% beadparamsZmatch_Y = beadparams.allY{3};
% if scope==1 % patch
%     beadparams.dispstage{3}=beadparams_.dispstage{3}/150*144.5;
% end
% beadparamsZmatchdispstage = beadparams.dispstage{3};

xlocs = 1:dims(1);
ylocs = 1:dims(2);
[xy2,xy1] = ndgrid(ylocs(:),xlocs(:));
xy = [xy1(:),xy2(:)];

%%
Ntiles = size(neigs,1);
edges = cell(1,Ntiles);
for itile = 1:Ntiles
    edge = paireddescriptor{itile}.neigs(2:3);
    counts = paireddescriptor{itile}.count;
    edges{itile} = [[itile;itile], (edge(:)), counts(:)];
end
edges = cat(1,edges{:});
edges(any(isnan(edges),2),:)=[];
G = sparse(edges(:,1),edges(:,2),edges(:,3),Ntiles,Ntiles);
G = max(G,G');

%%
%parfor_progress(Ntiles)
if old
    skipinds = any(isnan(neigs4(:,[4 5])),2);
else
    % keeps tile if [(-x|+x) & (-y|+y)]: if there exists pairs on x&y
    skipinds = ~(any(isfinite(neigs4(:,[2 4])),2) & any(isfinite(neigs4(:,[3 5])),2)) ;
    % keeps tile only if there exists a consecutice pairs on x&y
    % skipinds = any(isnan(neigs4(:,2:3)),2)&any(isnan(neigs4(:,4:5)),2);
end
validthis = zeros(1,Ntiles);
scopeparams = struct_with_shape_and_fields([1 Ntiles], {'imsize_um', 'dims', 'affinegl', 'affineglFC'}) ;

%%
%try; parfor_progress(0);catch;end
parfor_progress(Ntiles) ;
parfor itile = 1:Ntiles ,
%for itile = Ntiles-2000:Ntiles ,
%for itile = 16618:16619
    %%
    %itile
    neigs_this_tile = neigs(itile,:) ;
    neigs4_this_tile = neigs4(itile,:) ;
    scopeparams(itile).imsize_um = imsize_um;
    scopeparams(itile).dims = dims;
    neiginds = find(G(itile,:));
    theseinds = setdiff(neiginds,paireddescriptor{itile}.neigs(2:3));

    if skipinds(itile)

    else
        allX = [] ;
        allY = [] ;
        sdisp = [] ;
        siz = zeros(1); % right/below + left/above
        stgdisp = zeros(3,1); % right/below + left/above
        % right adjacency
        if ~isnan(neigs_this_tile(2))
            siz(1) = size(paireddescriptor{itile}.onx.X,1);
            allX = [allX;paireddescriptor{itile}.onx.X];
            allY = [allY;paireddescriptor{itile}.onx.Y];
            stgdisp(:,1) = 1000*(scopeloc.loc(neigs_this_tile(2),:)-scopeloc.loc(neigs_this_tile(1),:)); %#ok<PFBNS>
            sdisp = [sdisp,1000*stgdisp(:,1)*ones(1,siz(1))];
            validthis(itile) = 1;
        end
        if ~isnan(neigs_this_tile(3))
            siz(end+1) = size(paireddescriptor{itile}.ony.X,1);  %#ok<AGROW>
            allX = [allX;paireddescriptor{itile}.ony.X];  
            allY = [allY;paireddescriptor{itile}.ony.Y];  
            stgdisp(:,end+1) = 1000*(scopeloc.loc(neigs_this_tile(3),:)-scopeloc.loc(neigs_this_tile(1),:));  %#ok<AGROW>
            sdisp = [sdisp,1000*stgdisp(:,end)*ones(1,siz(end))];  
            validthis(itile) = 1;
        end

        if ~isempty(theseinds) && ~old
            % check left/above
            for ii=1:length(theseinds)
                if find(paireddescriptor{theseinds(ii)}.neigs==itile)==2 ,  %#ok<PFBNS> % left
                    % left adjacency
                    ileft = neigs4_this_tile(2);
                    if ~isnan(ileft)
                        siz(end+1) = size(paireddescriptor{ileft}.onx.X,1);  %#ok<AGROW>
                        allX = [allX;paireddescriptor{ileft}.onx.Y];  
                          % flip X<->Y  
                        allY = [allY;paireddescriptor{ileft}.onx.X];  
                        stgdisp(:,end+1) = 1000*(scopeloc.loc(ileft,:)-scopeloc.loc(neigs_this_tile(1),:));  %#ok<AGROW>
                        sdisp = [sdisp,1000*stgdisp(:,end)*ones(1,siz(end))];  
                        validthis(itile) = 1;
                    end
                elseif find(paireddescriptor{theseinds(ii)}.neigs==itile)==3 % above
                    iabove = neigs4_this_tile(3);
                    siz(end+1) = size(paireddescriptor{iabove}.ony.X,1);  %#ok<AGROW>
                    allX = [allX;paireddescriptor{iabove}.ony.Y];
                      % flip X<->Y
                    allY = [allY;paireddescriptor{iabove}.ony.X];
                    stgdisp(:,end+1) = 1000*(scopeloc.loc(iabove,:)-scopeloc.loc(neigs_this_tile(1),:));  %#ok<AGROW>
                    sdisp = [sdisp,1000*stgdisp(:,end)*ones(1,siz(end))];
                    validthis(itile) = 1;
                else
                end
            end
        end

%         if isfield(params,'beadparams') && ~isempty(params.beadparams) ,
%             % append bead params if exists
%             beadparamsZmatch_X = params.beadparams.allX{3};
%             beadparamsZmatch_Y = params.beadparams.allY{3};
%             beadparamsZmatchdispstage = params.beadparams.dispstage{3};
%         else
        num = max(1,round(size(allX,1)/2));
        beadparamsZmatch_X=ones(num,1)*[0 0 250];
        beadparamsZmatch_Y=ones(num,1)*[0 0 0];
        beadparamsZmatchdispstage = (ones(num,1)*[0 0 250]*1e3)';
%         end
        %siz(end+1) = size(beadparamsZmatch_X,1);
        allX = [allX;beadparamsZmatch_X];
        allY = [allY;beadparamsZmatch_Y];
        sdisp = [sdisp,beadparamsZmatchdispstage];
        %allXFC = allX;
        %allYFC = allY;

        if validthis(itile) ,
            suballX = allX+1;
            suballY = allY+1;

            locs = util.fcshift(curvemodel(:,:,itile),order,xy,dims,suballX) ;
            allXFC = locs-1;
            locs = util.fcshift(curvemodel(:,:,itile),order,xy,dims,suballY) ;
            allYFC = locs-1;

            Dall = (allY-allX)';
            DallFC = (allYFC-allXFC)';
    %         inds=zeros(1,sum(siz(:)));
    %         siz12= [sum(siz(:,1:2),2) siz(:,3)]';
    %         idxsub = [0;cumsum(siz12(:))]';
    %         for ii=1:2:length(idxsub)-1
    %             inds([idxsub(ii)+1:idxsub(ii+1)]) = 1;
    %         end

%             if itile == 15843 || itile == 15844 ,
%                 keyboard
%             end

            glS_=sdisp/Dall;
            glSFC_=sdisp/DallFC;
            scopeparams(itile).affinegl = glS_;
            scopeparams(itile).affineglFC = glSFC_;
        end
    end
    parfor_progress() ;
end
parfor_progress(0) ;
    
end
