function [Ic_m,It_m] = vizMatch(scopeloc,neigs,descriptors,ineig,imsize_um,iadj,X,Y,X_,Y_)
% read tiles
% iadj = 2; % 1 for x, 2 for y
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
if iadj==1
    RB = imref2d(size(It_m),[1 dims(1)]+pixshift(1),[1 dims(2)]);
else
    RB = imref2d(size(It_m),[1 dims(1)],[1 dims(2)]+pixshift(2));
end
imshowpair(imadjust(Ic_m),RA,imadjust(It_m),RB,'falsecolor','Scaling','joint','ColorChannels','green-magenta')
hold on
myplot3(descent-1,{'bo','MarkerSize',6,'LineWidth',1})
myplot3(descadj-1,{'yo','MarkerSize',6,'LineWidth',1})
%%
myplot3(X-1,{'bo','MarkerSize',12,'LineWidth',1})
myplot3(Y-1,{'yo','MarkerSize',12,'LineWidth',1})
% delete(findobj('Color','r'))
hold on
Y_2 = Y_;
Y_2(:,iadj) = Y_2(:,iadj) + pixshift(iadj);
XX = [X_(:,1),Y_2(:,1),nan*X_(:,1)]'-1;
YY = [X_(:,2),Y_2(:,2),nan*X_(:,2)]'-1;
plot(XX,YY,'r')
