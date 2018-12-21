% curvature statistics

clear tr fd
for ii=1:Nneig
    fd(ii) = size(paireddescriptor{ii}.onx.X,1);
    if fd(ii)
        tr(ii,:) = median(paireddescriptor{ii}.onx.Y-paireddescriptor{ii}.onx.X);
    end
end
disponx=tr(all(tr,2),:);

clear tr fd
for ii=1:Nneig
    fd(ii) = size(paireddescriptor{ii}.ony.X,1);
    if fd(ii)
        tr(ii,:) = median(paireddescriptor{ii}.ony.Y-paireddescriptor{ii}.ony.X);
    end
end
dispony=tr(all(tr,2),:);
figure(100), 
subplot(211)
histogram(disponx(:,1),'BinWidth',1)
title('median displacement on y based on descriptor match')
subplot(212)
histogram(dispony(:,2),'BinWidth',1)
title('median displacement on x based on descriptor match')
export_fig(fullfile('./visualization_figures/curvatureest/','median_displacement.png'),'-transparent')
export_fig(fullfile('./visualization_figures/curvatureest/','stage_displacement.png'),'-transparent')

%%
idir = 1;
str = {'onx','ony'};

model_on_idir = squeeze(curvemodel(idir,:,:))';
model_on_idir = model_on_idir(all(model_on_idir,2),:);
hf = figure(101);clf
subplot(311);
h_p1 = histogram(model_on_idir(:,1));
% xlim([-3 +3]+getmod(model_on_idir(:,1)))
title('p1')

subplot(312)
h_p2 = histogram(model_on_idir(:,2));
% xlim(sort([.999 1.001]*getmod(model_on_idir(:,2))))
title('p2')

subplot(313)
h_p3 = histogram(model_on_idir(:,3));
%     xlim([-3 +3]+getmod(model_on_idir(:,3)))
title('p3')

% leg.FontSize = 32;
h_p1.BinWidth = .1;
%     h_p2.BinWidth = .2;
h_p3.BinWidth = .2;
% export_fig(fullfile('./visualization_figures/curvatureest/',sprintf('pest-%s.png',str{idir})),'-transparent')



