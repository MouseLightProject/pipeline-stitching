tottile = size(curvemodel{1},3);
rt1= zeros(tottile,3);
rt2= zeros(tottile,3);
siz1 = zeros(tottile,1);
siz2 = zeros(tottile,1);
for ii=1:tottile,rt1(ii,:)=curvemodel{1}(1,:,ii);siz1(ii)=size(paireddescriptor{1}{ii}.onx.X,1);end
for ii=1:tottile,rt2(ii,:)=curvemodel{1}(2,:,ii);siz2(ii)=size(paireddescriptor{1}{ii}.ony.X,1);end
%%
% close all
rt = rt2;
siz = siz2;
interior = util.interior_tiles(scopeloc,3);

% generate weighted histogram
xx = rt(:,3);
siz_uni = (siz>-1) .* (xx>0 & siz>0) .* (interior>0);
% replicate each entry siz many times
C = cell(1,tottile);
for ii=1:tottile
    C{ii} = repmat(xx(ii),siz_uni(ii),1);
end
C = cat(1,C{:});

figure, 
hist(C,100)
title('interior')
%%
set1 = abs(rt(:,3)-863)>3 & rt(:,1)>1;
set2 = abs(rt(:,3)-863)<3 & rt(:,1)>1;
rt(set1,:)
rt(set2,:)


%% %%%%%%%%%%%%%%%%%


% read tiles
iadj = 2; % 1 for x, 2 for y
idxcent = neigs(ineig,1);
idxadj = neigs(ineig,iadj+1);


[Ic,filec] = util.getTilefromId(scopeloc,idxcent);
[It, filet] = util.getTilefromId(scopeloc,idxadj);
Ic_m = rot90(max(Ic,[],3),2);
It_m = rot90(max(It,[],3),2);
dims = size(Ic);dims=dims([2 1 3]);
%%
idxcent = neigs(ineig,1);
descent = descriptors{idxcent};
descent = double(descent(:,1:3));
descent = util.correctTiles(descent,dims); % flip dimensions

descadj = descriptors{idxadj};
descadj = double(descadj(:,1:3)); % descadj has x-y-z-w1-w2 format
descadj = util.correctTiles(descadj,dims); % flip dimensions

stgshift = 1000*(scopeloc.loc(idxadj,:)-scopeloc.loc(idxcent,:));
pixshift = round(stgshift.*(dims-1)./(imsize_um));
descadj = descadj + ones(size(descadj,1),1)*pixshift; % shift with initial guess based on stage coordinate

%%
figure(43), cla
RA = imref2d(size(Ic_m),[1 dims(1)],[1 dims(2)]);
RB = imref2d(size(It_m),[1 dims(1)]+pixshift(1),[1 dims(2)]);
imshowpair(imadjust(Ic_m),RA,imadjust(It_m),RB,'falsecolor','Scaling','joint','ColorChannels','green-magenta')
hold on
myplot3(descent-1,{'bo','MarkerSize',6,'LineWidth',1})
myplot3(descadj-1,{'ro','MarkerSize',6,'LineWidth',1})

%%
figure(32)
showMatchedFeatures(Ic_m,It_m,X_(:,1:2),Y_(:,1:2))

%%
figure(13), clf, cla
imshow(Ic_m,[11 15]*1e3,'XData',[1 dims(1)],'YData',[1 dims(2)])
hold on
imshow(It_m,[11 15]*1e3,'XData',[1 dims(1)]+pixshift(1),'YData',[1 dims(2)])
xlim([1 dims(1)+pixshift(1)])
alpha(.5)
% 
% 
% hold on, myplot3(descent,'o')


% nbound = [0 0];
% nbound(1) = max(pixshift(iadj),min(descadj(:,iadj)));
% nbound(2) = min(dims(iadj),max(descent(:,iadj)))+0;
% X = descent(descent(:,iadj)>nbound(1)&descent(:,iadj)<nbound(2),:);
% Y = descadj(descadj(:,iadj)>nbound(1)&descadj(:,iadj)<nbound(2),:);

%%



load(fullfile(matfolder,'scopeparams_pertile_singletile'),'paireddescriptor', ...
    'scopeparams', 'R', 'curvemodel','scopeparams_', 'paireddescriptor_', ...
    'curvemodel_','params')
%%
[X_,Y_] = match.descriptorMatch(X,Y,matchparams);
Y_(:,iadj) = Y_(:,iadj)- pixshift(iadj);% move it back to original location after CDP
% get field curvature model
[X_,Y_,out] = match.fcestimate(X_,Y_,iadj,matchparams);



%%
figure
myplot3(X_,'.')
hold on
myplot3(Y_,'.')

%%
xlocs = 1:dims(1);
ylocs = 1:dims(2);
[xy2,xy1] = ndgrid(ylocs(:),xlocs(:));
xy = [xy1(:),xy2(:)];
order = params.order;
dims = params.imagesize;
%%
icent = 5136
ineig = 5163

% ineig=5136
kk=3

% icent = 5461
% ineig = 5488
% Acent = -scopeparams{2}(icent).affineglFC/1e3;
Acent = -glSFC_2/1e3
Aneig = Acent;%-scopeparams{2}(ineig).affineglFC/1e3;

Xcent = paireddescriptor{3}{icent}.ony.X;
Ycent = paireddescriptor{3}{icent}.ony.Y;
cur_cent=curvemodel{3}(:,:,icent);
cur_neig=curvemodel{3}(:,:,ineig);

% Xcent = paireddescriptor{icent}.ony.X;
% Ycent = paireddescriptor{icent}.ony.Y;
% cur_cent=curvemodel(:,:,icent);
% cur_neig=curvemodel(:,:,ineig);

% cur(2,2)= 0;
[Xcent_fc,xshift2D,yshift2D] = util.fcshift(cur_cent,order,xy,dims,Xcent+1);
Xcent_fc = Xcent_fc-1;
[Ycent_fc] = util.fcshift(cur_neig,order,xy,dims,Ycent+1);
Ycent_fc = Ycent_fc-1;


% sdisp = ones(size(Xcent,1),1)*1000*(scopeloc.loc(icent,:)-scopeloc.loc(ineig,:));
% Acent_ = sdisp'/[Ycent-Xcent]';
% res=Acent_*Ycent' -[ Acent_*Xcent'+(scopeloc.loc(icent,:)-scopeloc.loc(ineig,:))'*1000];
% norm(res(2,:))

% Acent_ = sdisp'/[Ycent_fc-Xcent_fc]';
% res=Acent_*(Ycent_fc - Xcent_fc)'-(scopeloc.loc(icent,:)-scopeloc.loc(ineig,:))'*1000;
% (abs(res(2,:)))>.35 % 1 pix
%
% figure
% myplot3(Xcent,'b.')
% %%
% glSFC_ = sdisp/DallFC
% res = [glSFC_*(Dall)-sdisp]/1e3;
% res = sqrt(sum(res.^2,1));
%
% glSFC_2 = sdisp(:,res<.5)/DallFC(:,res<.5)
%
% max(sqrt(sum(res.^2,1)))

%
shift = median(Xcent-Ycent);
shift_fc = median(Xcent_fc-Ycent_fc);
locX = [[Acent 1e3*scopeloc.loc(icent,:)']*[Xcent ones(size(Xcent,1),1)]']';
locY = [[Aneig 1e3*scopeloc.loc(ineig,:)']*[Ycent ones(size(Xcent,1),1)]']';
locX_fc = [[Acent 1e3*scopeloc.loc(icent,:)']*[Xcent_fc ones(size(Xcent,1),1)]']';
locY_fc = [[Aneig 1e3*scopeloc.loc(ineig,:)']*[Ycent_fc ones(size(Xcent,1),1)]']';

figure(icent+kk)
dcm = datacursormode( gcf );
dcm.UpdateFcn = @NewCallback;
clf
subplot(311)
hold on
myplot3(locX,'b.')
myplot3(locX(25,:),'r*')
myplot3(locY,'ro')
axis equal
set(gca,'YDir','reverse')
title(sprintf('%d/%d/%d/-',shift))
subplot(312)
hold on
myplot3(locX_fc,'b.')
myplot3(locY_fc,'ro')
axis equal
set(gca,'YDir','reverse')
title(sprintf('%d/%d/%d/fc',round(shift_fc)))

if 1
    figure(icent*10+kk)
    dcm = datacursormode( gcf );
    dcm.UpdateFcn = @NewCallback;
    clf
    subplot(311)
    hold on
    myplot3(Xcent,'b.')
    myplot3(Xcent(25,:),'b*')
    myplot3(Ycent+shift,'ro')
    set(gca,'YDir','reverse')
    axis equal
    title(sprintf('%d/%d/%d',shift))
    subplot(312)
    hold on
    myplot3(Xcent_fc,'b.')
    myplot3(Ycent_fc+shift_fc,'ro')
    set(gca,'YDir','reverse')
    axis equal
    title(sprintf('%d/%d/%d',round(shift_fc)))
end
%%
figure(icent)
dcm = datacursormode( gcf );
dcm.UpdateFcn = @NewCallback;
clf
subplot(311)
hold on
myplot3(Xcent,'b.')
myplot3(Ycent+shift,'ro')
axis equal
title(sprintf('%d/%d/%d',shift))
subplot(312)
hold on
myplot3(Xcent_fc,'b.')
myplot3(Ycent_fc+shift_fc,'ro')
axis equal
title(sprintf('%d/%d/%d',round(shift_fc)))
datacursormode(gcf)

%%



%%
% for debugging scripts
ineig = 5136
tifname = scopeloc.filepath{ineig}
filefold = fileparts(tifname);
[~,tifname] = fileparts(filefold);
tifpath = fullfile(filefold,sprintf('%s-ngc.0.tif',tifname));
Ic = mytiffread(tifpath);
%
try
    ineigx = neighbors(ineig,4)
    ineigy = neighbors(ineig,5)
catch
    ineigx = neigs(ineig,2)
    ineigy = neigs(ineig,3)
end

tifname = scopeloc.filepath{ineigy}
filefold = fileparts(tifname);
[~,tifname] = fileparts(filefold);
tifpath = fullfile(filefold,sprintf('%s-ngc.0.tif',tifname));
In = mytiffread(tifpath);

%%
ImC = max(Ic,[],3);
ImN = max(In,[],3);
ImX = ImC;
ImY = ImN;

ImX = rot90(ImC,2);
ImY = rot90(ImN,2);
%%
curvemodel(:,3,[5162 5163 5164])
%%
figure(7)
cla
imagesc(density,'Xdata',[gridx],'Ydata',[gridy])
axis tight
hold on,
plot(x,y,'m.')

%%
figure(11), cla
imshow(ImX,[11 15]*1e3)
hold on
myplot3(Xcent,'o')
%%
myplot3(descent,'d')
myplot3(X,'o')
myplot3(X_,'go')
% myplot3(X_(prob_inliers,:),'go')

%%
figure(21), cla
imshow(ImY,[11 15]*1e3)
hold on
myplot3(Y`cent,'o')
%%
myplot3(descadj-pixshift,'d')
myplot3(Y-pixshift,'o')
myplot3(Y_,'+')

%%
dispvec = X_-Y_;
figure,
plot(X_(:,1),dispvec(:,2),'.')
% xlim([min(dispvec(:,1)) max(dispvec(:,1))])
axis tight

%%
figure(1),
cla
imshow(ImC,[11 15]*1e3)
hold on
myplot3(paireddescriptor{ineig}.onx.X,'o')
% myplot3(descriptors{ineig}(:,1:3),'d')
% camroll(90)


%%
figure(2),
imshow(ImN,[11 15]*1e3)
hold on
myplot3(paireddescriptor{ineig}.ony.Y,'o')

%%
figure(2),
imshow(ImXp1,[])
hold on
myplot3(Y_,'o')
% myplot3(descent,'o')
%%
old=[-376.7212665974943206, -1.4150539397464144, 44.3289718325040667, 0.0, 71108753.59580696,
    2.3374766378038374, -297.0775556009028833, -13.8336760551395486, 0.0, 14667134.00371697,
    0.2542090768871866, 0.7845732492211573, 888.2112152288162861, 0.0, 35412268.12609629,
    0.0, 0.0, 0.0, 1.0, 0.0,
    -0.0, 0.0, 0.0, 0.0, 1.0]
new=[-376.2621493364757725, -1.4195567316836457, 44.7492455332472261, 0.0, 71108071.94444494,
    2.2485697966564264, -296.7169950905342830, -14.2602659659113709, 0.0, 14667817.50149242,
    0.0439235442969926, 0.5983385231528975, 887.8041013352451500, 0.0, 35413245.98378776,
    0.0, 0.0, 0.0, 1.0, 0.0,
    -0.0, 0.0, 0.0, 0.0, 1.0];

round(old*[1;1;1;0;1]/1000)'
round(new*[1;1;1;0;1]/1000)'
round(old*[1024;1536;1;0;1]/1000)'
round(new*[1024;1536;1;0;1]/1000)'
%%
err=[]
for itt = 33:35
    er1=scopeloc.loc(find(scopeloc.gridix(:,1)==220&scopeloc.gridix(:,2)==itt&scopeloc.gridix(:,3)==779),:)
    err(end+1,:) = er1(1,:);
end
%%
err3=[]
for ix = 218:220
    for iy = 33:35
        er1=scopeloc.loc(find(scopeloc.gridix(:,1)==ix&scopeloc.gridix(:,2)==iy&scopeloc.gridix(:,3)==781),:)
        err3(end+1,:) = er1(1,:);
    end
end
diff(err3)
%%
er2=scopeloc.loc(find(scopeloc.gridix(:,2)==33&scopeloc.gridix(:,3)==780),:)
er3=scopeloc.loc(find(scopeloc.gridix(:,2)==34&scopeloc.gridix(:,3)==780),:)

