%%
if 0
    numtiles = size(scopeloc.loc,1);
    psh = nan(numtiles,3,2);
    for ii=1:numtiles
        for iadj = 1:2
            if isfinite(neigs(ii,iadj+1))
                stgshift = 1000*(scopeloc.loc(neigs(ii,iadj+1),:)-scopeloc.loc(neigs(ii,1),:));
                pixshift = round(stgshift.*(dims-1)./(imsize_um));
                psh(ii,:,iadj) = pixshift;
            end
        end
    end
    
    these = isfinite(psh(:,1,1));
    onx = psh(these,1,1);
    these = isfinite(psh(:,2,2));
    ony = psh(these,2,2);
    
    figure(123), cla
    subplot(121)
    hist(onx,100)
    hold on
    subplot(122)
    hist(ony,100)
end
ineig = 11719
% function iter_on_ineig(ineig)
idxcent = neigs(ineig,1);
descent = descriptors{idxcent};
descent = double(descent(:,1:3));
descent = util.correctTiles(descent,dims); % flip dimensions
mout = zeros(3,3);
paireddescriptor_ = paireddesctemp;
R_ = zeros(3); % median residual
%
iadj = 1%1:size(neigs,2)-2 %1:x-overlap, 2:y-overlap, 3:z-overlap
% idaj : 1=right(+x), 2=bottom(+y), 3=below(+z)
idxadj =  neigs(ineig,iadj+1);
descadj = descriptors{idxadj};
descadj = double(descadj(:,1:3)); % descadj has x-y-z-w1-w2 format
descadj = util.correctTiles(descadj,dims); % flip dimensions
stgshift = 1000*(scopeloc.loc(idxadj,:)-scopeloc.loc(idxcent,:));
pixshift = round(stgshift.*(dims-1)./(imsize_um));
descadj = descadj + ones(size(descadj,1),1)*pixshift; % shift with initial guess based on stage coordinate
nbound = [0 0];
nbound(1) = max(pixshift(iadj),min(descadj(:,iadj))) - 15;
nbound(2) = min(dims(iadj),max(descent(:,iadj)))+0 + 15;
X = descent(descent(:,iadj)>nbound(1)&descent(:,iadj)<nbound(2),:);
Y = descadj(descadj(:,iadj)>nbound(1)&descadj(:,iadj)<nbound(2),:);

% get descpair
[X_,Y_] = match.descriptorMatch4XY(X,Y,matchparams);
Y_(:,iadj) = Y_(:,iadj)- pixshift(iadj);% move it back to original location after CDP

% get field curvature model
if isfield(matchparams,'init_array') % overwrites any initialization with per tile values
    pinit_model = matchparams.init_array(:,:,ineig);
elseif isfield(matchparams,'init')
    pinit_model = matchparams.init;
end

% pinit_model
% [X_e,Y_e,out_e,valid_e] = match.fcestimate(X_,Y_,iadj,matchparams,pinit_model);
%%
% [X_e,Y_e,out_e,valid_e] = match.fcestimate(X_,Y_,iadj,matchparams);
matchparams.viz = 1
[X_e2,Y_e2,out_e2,valid_e2] = match.fcestimate(X_,Y_,iadj,matchparams,pinit_model)
matchparams.viz = 0
% [X_e2,Y_e2,out_e2,valid_e2] = match.fcestimate(X_e2,Y_e2,iadj,matchparams,pinit_model);

%%
if viz & 0;
    util.debug.vizMatch(scopeloc,neigs,descriptors,ineig,imsize_um,iadj);
    %%
    myplot3(X-1,{'bo','MarkerSize',12,'LineWidth',1})
    myplot3(Y-1,{'yo','MarkerSize',12,'LineWidth',1})
    % delete(findobj('Color','r'))
    hold on
    XX = [X_(:,1),Y_(:,1),nan*X_(:,1)]'-1;
    YY = [X_(:,2),Y_(:,2),nan*X_(:,2)]'-1;
    plot(XX,YY,'r')
end
if size(X_,1)<3 | size(Y_,1)<3;continue;end
%%
% get field curvature model
[X_,Y_,out,valid] = match.fcestimate(X_,Y_,iadj,matchparams);

%%
% flip back dimensions
X_ = util.correctTiles(X_,dims);
Y_ = util.correctTiles(Y_,dims);

% store pairs
mout(iadj,:) = out;
paireddescriptor_{iadj}.valid = valid;
paireddescriptor_{iadj}.X = X_;
paireddescriptor_{iadj}.Y = Y_;
%R(:,iadj,ineig) = round(median(X_-Y_));
R_(:,iadj) = round(median(X_-Y_));
