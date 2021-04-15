function descriptorMatchQualityHeatMap(regpts, scopeparams, scopeloc, video_file_name)
    %DESCRIPTORMATCHQUALITY Summary of this function goes here
    %
    % [OUTPUTARGS] = DESCRIPTORMATCHQUALITY(INPUTARGS) Explain usage here
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

    % $Author: base $	$Date: 2016/11/15 11:13:34 $	$Revision: 0.1 $
    % Copyright: HHMI 2016

    % populate afftile
    tile_count = length(scopeparams) ;
    linear_transform_from_tile_index = zeros(3, 3, tile_count) ;
    for tile_index = 1 : tile_count ,
        if ~isempty(scopeparams(tile_index).affineglFC) ,
            linear_transform_from_tile_index(:,:,tile_index) = scopeparams(tile_index).affineglFC ;  % nm
        end
    end
    mflip = -eye(3) ; %mflip(3,3)=1;
    linear_transform = median(linear_transform_from_tile_index,3) * mflip ;
    affine_transform_from_tile_index = zeros(3, 4, tile_count) ;
    for tile_index = 1 : tile_count ,
        % form an affine matrix
        affine_transform = [linear_transform scopeloc.loc(tile_index,:)'*1e6] ;  % nm
        affine_transform_from_tile_index(:,:,tile_index) = affine_transform ;  % nm
    end

    %lattice_xyz_from_tile_index = scopeloc.gridix(:,1:3) ;
    %neighbors = buildNeighbor(lattice_xyz_from_tile_index) ;  %[id -x -y +x +y -z +z] format
    %checkthese = [1 4 5 7]; % 0 - below

    % % For each tile, compute the offset to the z+1 tile
    % dz_to_plus_z_tile = zeros(1, tile_count) ;  % mm
    % for tile_index = 1:tile_count ,
    %     if any(isnan(neighbors(tile_index,[1 7]))) ,
    %         dz_to_plus_z_tile(tile_index) = nan ;
    %     else
    %         dz_to_plus_z_tile(tile_index) = diff(scopeloc.loc(neighbors(tile_index,[1 7]),3));  % mm
    %     end
    % end

    lattice_z_from_tile_index = scopeloc.gridix(:,3) ;  % column
    lattice_z_from_z_layer_index = unique(lattice_z_from_tile_index) ;
    maximum_lattice_z = lattice_z_from_z_layer_index(end) ;
    z_layer_count = length(lattice_z_from_z_layer_index) ;

    max_xy_discrepancy_from_z_layer_index = zeros(1, z_layer_count) ;
    max_xyz_discrepancy_from_z_layer_index = zeros(1, z_layer_count) ;
    %xyz_t = cell(1, z_layer_count) ;
    %xyz_tp1 = cell(1, z_layer_count) ;
    self_fiducial_xyz_from_pair_index_from_z_layer_index = cell(1, z_layer_count) ;  % this is pair index within layer
    other_fiducial_xyz_from_pair_index_from_z_layer_index = cell(1, z_layer_count) ;  % this is pair index within layer
    for z_layer_index = 1 : z_layer_count ,
        lattice_z = lattice_z_from_z_layer_index(z_layer_index) ;

        fprintf('    Layer %d of %d\n', lattice_z, maximum_lattice_z) ;
        is_in_this_layer_from_tile_index = (lattice_z_from_tile_index'==lattice_z);

        tile_index_from_foo_index = find(is_in_this_layer_from_tile_index) ;  % a "foo" is a tile in this layer
        foo_count = length(tile_index_from_foo_index) ;
        %match_count_from_foo_index = zeros(1, foo_count) ;
        is_bar_from_foo_index = false(1, foo_count) ;
        for foo_index = 1 : foo_count ,
            tile_index = tile_index_from_foo_index(foo_index) ;
            regpts_for_this_tile = regpts{tile_index} ;
            is_bar_from_foo_index(foo_index) = (size(regpts_for_this_tile.X,1)>0) && ~isnan(regpts_for_this_tile.neigs(4)) ;
        end

        foo_index_from_bar_index = find(is_bar_from_foo_index) ;  % a "bar" is a tile in this layer with a non-zero match count & non-nan .neigs(4)
        bar_count = length(foo_index_from_bar_index) ;
        tile_index_from_bar_index = tile_index_from_foo_index(foo_index_from_bar_index) ;

        self_fiducial_xyz_from_pair_index_from_bar_index = cell(1, bar_count) ;  % this is pair index within tile
        other_fiducial_xyz_from_pair_index_from_bar_index = cell(1, bar_count) ;
        max_xy_discrepancy_from_bar_index = zeros(1, bar_count) ;
        max_xyz_discrepancy_from_bar_index = zeros(1, bar_count) ;
        for bar_index = 1 : bar_count ,
            tile_index = tile_index_from_bar_index(bar_index) ;
            regpts_for_this_tile = regpts{tile_index} ;
            neigs = regpts_for_this_tile.neigs ;
            % apply transforms to 
            affine_transform_for_this_tile = affine_transform_from_tile_index(:,:,neigs(1)) ;
            fiducial_xyz_for_this_tile_from_pair_index = add_ones_column(regpts_for_this_tile.X) * affine_transform_for_this_tile' ;  % nm
            affine_transform_for_other_tile = affine_transform_from_tile_index(:,:,neigs(4)) ;
            fiducial_xyz_for_other_tile_from_pair_index = add_ones_column(regpts_for_this_tile.Y) * affine_transform_for_other_tile' ;  % nm
            self_fiducial_xyz_from_pair_index_from_bar_index{bar_index} = fiducial_xyz_for_this_tile_from_pair_index ;  % nm
            other_fiducial_xyz_from_pair_index_from_bar_index{bar_index} = fiducial_xyz_for_other_tile_from_pair_index ;  % nm
            xyz_discrepancy_vector_for_this_tile_from_pair_index = ...
                fiducial_xyz_for_this_tile_from_pair_index - fiducial_xyz_for_other_tile_from_pair_index ;  % nm
            max_xy_discrepancy_from_bar_index(bar_index) = max(sqrt(sum(xyz_discrepancy_vector_for_this_tile_from_pair_index(:,1:2).^2,2)))/1e3 ;  % um
            max_xyz_discrepancy_from_bar_index(bar_index) = max(sqrt(sum(xyz_discrepancy_vector_for_this_tile_from_pair_index(:,1:3).^2,2)))/1e3 ;  % um
        end
        self_fiducial_xyz_from_pair_index_from_z_layer_index{z_layer_index} = cat(1,self_fiducial_xyz_from_pair_index_from_bar_index{:});  % nm
        other_fiducial_xyz_from_pair_index_from_z_layer_index{z_layer_index} = cat(1,other_fiducial_xyz_from_pair_index_from_bar_index{:});  % nm
        if bar_count > 0 ,
            max_xy_discrepancy_from_z_layer_index(z_layer_index) = max(max_xy_discrepancy_from_bar_index) ;  % um
            max_xyz_discrepancy_from_z_layer_index(z_layer_index) = max(max_xyz_discrepancy_from_bar_index) ;  % um
        end
    end

    centerPoint = 25;
    scalingIntensity = 4;
    max_xy_discrepancy = max(max_xy_discrepancy_from_z_layer_index) ;
    min_xy_discrepancy = min(max_xy_discrepancy_from_z_layer_index) ;
    newMap = newColMap(centerPoint, scalingIntensity, min_xy_discrepancy, max_xy_discrepancy) ;

    trmax = cellfun(@max,self_fiducial_xyz_from_pair_index_from_z_layer_index,'UniformOutput',false) ;
    trmin = cellfun(@min,self_fiducial_xyz_from_pair_index_from_z_layer_index,'UniformOutput',false) ;
    Rmax = max(cat(1,trmax{:}))+scopeparams(1).imsize_um*1e3.*[1 1 0];
    Rmin = min(cat(1,trmin{:}))-scopeparams(1).imsize_um*1e3.*[1 1 0];

    % plot descriptors
    myfig = 100;
    figure(myfig) ;
    cla() ;
    clf() ;
    hold on
    F(z_layer_count) = struct('cdata',[],'colormap',[]);

    for z_layer_index = 1 : z_layer_count ,
        lattice_z = lattice_z_from_z_layer_index(z_layer_index) ;

        is_in_this_layer_from_tile_index = (lattice_z_from_tile_index'==lattice_z) ;
        tile_count_in_this_layer = sum(is_in_this_layer_from_tile_index) ;
        x = scopeloc.loc(is_in_this_layer_from_tile_index,1)*1e6-scopeparams(1).imsize_um(1)*1e3;  % nm
        y = scopeloc.loc(is_in_this_layer_from_tile_index,2)*1e6-scopeparams(1).imsize_um(1)*1e3;  % nm
        w = scopeparams(1).imsize_um(1)*ones(tile_count_in_this_layer,1)*1e3;  % nm
        h = scopeparams(1).imsize_um(2)*ones(tile_count_in_this_layer,1)*1e3;  % nm

        self_fiducial_xyz_from_pair_index = self_fiducial_xyz_from_pair_index_from_z_layer_index{z_layer_index} ;  % pair index within layer
        other_fiducial_xyz_from_pair_index = other_fiducial_xyz_from_pair_index_from_z_layer_index{z_layer_index} ;  % pair index within layer        

        is_layer_empty = isempty(self_fiducial_xyz_from_pair_index) ;
        if is_layer_empty ,
            interpolated_xy_discrepancy_from_index_of_tile_within_layer = [] ;
            interpolated_xyz_discrepancy_from_index_of_tile_within_layer = [] ;
        else
            X = self_fiducial_xyz_from_pair_index(:,1) ;  % nm
            Y = self_fiducial_xyz_from_pair_index(:,2) ;  % nm
            Z = self_fiducial_xyz_from_pair_index(:,3) ;  % nm
            U = other_fiducial_xyz_from_pair_index(:,1) - X ;  % nm
            V = other_fiducial_xyz_from_pair_index(:,2) - Y ;  % nm
            W = other_fiducial_xyz_from_pair_index(:,3) - Z ;  % nm

            xyz_discrepancy_from_pair_index_within_layer = sqrt(U.^2+V.^2+W.^2) / 1e3 ;  % um
            xy_discrepancy_from_pair_index_within_layer = sqrt(U.^2+V.^2) / 1e3 ;  % um
            xy_discrepancy_from_xyz = ...
                scatteredInterpolant(other_fiducial_xyz_from_pair_index, xy_discrepancy_from_pair_index_within_layer, 'linear', 'nearest') ;
            xyz_discrepancy_from_xyz = ...
                scatteredInterpolant(other_fiducial_xyz_from_pair_index, xyz_discrepancy_from_pair_index_within_layer, 'linear', 'nearest') ;
            tile_xyz_from_index_of_tile_within_layer = scopeloc.loc(is_in_this_layer_from_tile_index,:)*1e6 ;  % nm
            interpolated_xy_discrepancy_from_index_of_tile_within_layer = xy_discrepancy_from_xyz(tile_xyz_from_index_of_tile_within_layer) ;
            interpolated_xyz_discrepancy_from_index_of_tile_within_layer = xyz_discrepancy_from_xyz(tile_xyz_from_index_of_tile_within_layer) ;
        end
        % initalize canvas
        figure(myfig), cla, clf
        h1=subaxis(2, 1, 1, 'sh', 0.03, 'sv', -0.03, 'padding', .04, 'margin', 0);
        cla, hold on
        fprintf('    Layer %d of %d\n', lattice_z, maximum_lattice_z) ;    
        %disp(['    Layer ' num2str(lattice_z) ' of ' num2str(max(lattice_z_from_tile_index))]);
        caxis([min_xy_discrepancy 0.9*max_xy_discrepancy])
        set(h1, 'XTick', []);
        set(h1, 'YTick', []);
        set(h1,'Color',[1 1 1]*.5)
        set(h1,'Ydir','reverse')
        set(h1,'Box','on')
        set(h1, ...
            'XColor'      , [.3 .3 .3], ...
            'YColor'      , [.3 .3 .3], ...
            'LineWidth'   , 10        );
    %     text(Rmin(1)+5e5,Rmin(2)+7e5,'Z',...
    %         'FontSize',40,'Color','k','HorizontalAlignment','left')
        view(0,90)
        colorbar
        colormap(newMap)

        for index_of_tile_within_layer = 1 : tile_count_in_this_layer ,
            rectangle('Position', [x(index_of_tile_within_layer) ...
                                   y(index_of_tile_within_layer) ...
                                   w(index_of_tile_within_layer) ...
                                   h(index_of_tile_within_layer)] , ...
                      'EdgeColor', 'r') ;
            if ~is_layer_empty ,
                xp = [x(index_of_tile_within_layer) ...
                      x(index_of_tile_within_layer) ...
                      x(index_of_tile_within_layer)+w(index_of_tile_within_layer) ...
                      x(index_of_tile_within_layer)+w(index_of_tile_within_layer)] ;
                yp = [y(index_of_tile_within_layer) ...
                      y(index_of_tile_within_layer)+h(index_of_tile_within_layer) ...
                      y(index_of_tile_within_layer)+h(index_of_tile_within_layer) ...
                      y(index_of_tile_within_layer)] ;
                interpolated_xyz_discrepancy = ...
                    interpolated_xyz_discrepancy_from_index_of_tile_within_layer(index_of_tile_within_layer) ;
                interpolated_xy_discrepancy = ...
                    interpolated_xy_discrepancy_from_index_of_tile_within_layer(index_of_tile_within_layer) ;
                interpolated_z_discrepancy = sqrt(interpolated_xyz_discrepancy^2 - interpolated_xy_discrepancy^2) ;
                patch(xp, yp, interpolated_z_discrepancy) ;
            end
        end
        if ~is_layer_empty ,
            scatter(self_fiducial_xyz_from_pair_index(:,1),self_fiducial_xyz_from_pair_index(:,2),6,'filled', ...
                    'MarkerFaceAlpha',.2,'MarkerFaceColor',[1 1 1]*.5)
        end
        axis equal
        xlim([Rmin(1) Rmax(1)])
        ylim([Rmin(2) Rmax(2)])
        %%
        h2=subaxis(2, 1, 2, 'sh', 0.03, 'sv', -0.03, 'padding', 0.04, 'margin', 0);
        cla, hold on
        set(h2, 'XTick', []);
        set(h2, 'YTick', []);
        set(h2,'Color',[1 1 1]*.5)
        set(h2,'Ydir','reverse')
        set(h2,'Box','on')
        set(h2, ...
            'XColor'      , [.3 .3 .3], ...
            'YColor'      , [.3 .3 .3], ...
            'LineWidth'   , 10        );
        view(0,90)
        axis equal
        colorbar
        caxis([min(max_xy_discrepancy_from_z_layer_index) max(max_xy_discrepancy_from_z_layer_index)*.9])
        for index_of_tile_within_layer = 1 : tile_count_in_this_layer ,
            rectangle('Position', [x(index_of_tile_within_layer) ...
                                   y(index_of_tile_within_layer) ...
                                   w(index_of_tile_within_layer) ...
                                   h(index_of_tile_within_layer)]) ;
            if is_layer_empty
                continue
            end

            xp = [x(index_of_tile_within_layer) ...
                  x(index_of_tile_within_layer) ...
                  x(index_of_tile_within_layer)+w(index_of_tile_within_layer) ...
                  x(index_of_tile_within_layer)+w(index_of_tile_within_layer)] ;
            yp = [y(index_of_tile_within_layer) ...
                  y(index_of_tile_within_layer)+h(index_of_tile_within_layer) ...
                  y(index_of_tile_within_layer)+h(index_of_tile_within_layer) ...
                  y(index_of_tile_within_layer)] ;
            interpolated_xy_discrepancy = interpolated_xy_discrepancy_from_index_of_tile_within_layer(index_of_tile_within_layer) ;
            patch(xp, yp, interpolated_xy_discrepancy) ;
            %         if length(vec{t})>=fix(ii)
            %             vv=vec{t}{fix(ii)};
            %             if ~isempty(vv)
            %                 vv=normr(vv(1:2));
            %                 %                 quiver(x(ii)+w(ii)/2,y(ii)+h(ii)/2,vv(1),vv(2),w(ii)/4,...
            %                 %                     'LineWidth',2,'Color','m','MaxHeadSize',4)
            %
            %                 headWidth = 6;
            %                 headLength = 3;
            %                 ah = annotation('arrow',...
            %                     'headStyle','plain',...
            %                     'HeadLength',headLength,'HeadWidth',headWidth,...
            %                     'Color','r');
            %                 set(ah,'parent',gca);
            %
            %                 x0 = x(ii)+w(ii)/2;
            %                 y0 = y(ii)+h(ii)/2;
            %                 vx0 = vv(1)*w(ii)/3;
            %                 vy0 = vv(2)*h(ii)/3;
            %
            %                 xn0 = (x0-Rmin(1))/(Rmax(1)-Rmin(1));
            %                 yn0 = (y0-Rmin(2))/(Rmax(2)-Rmin(2));
            %                 vnx0 = vx0/(Rmax(1)-Rmin(1));
            %                 vny0 = vy0/(Rmax(2)-Rmin(2));
            %                 set(ah,'position',[x(ii)+w(ii)/2,y(ii)+h(ii)/2,vv(1)*w(ii)/3,vv(2)*h(ii)/3]);
            %                 set(ah,'position',[xn0,yn0,vnx0,vny0]);
            %             end
            %         end
        end
        if ~is_layer_empty ,
            scatter(self_fiducial_xyz_from_pair_index(:,1), self_fiducial_xyz_from_pair_index(:,2), 6, 'filled', ...
                'MarkerFaceAlpha',.2,'MarkerFaceColor',[1 1 1]*.5)
        end

    %     text(Rmin(1)+5e5,Rmin(2)+7e5,'Lateral',...
    %         'FontSize',40,'Color','k','HorizontalAlignment','left')
        text(Rmin(1)+1e5, Rmax(2)-7e5, num2str(lattice_z), ...
            'FontSize',40, 'Color','k', 'HorizontalAlignment','left') ;

        xlim([Rmin(1) Rmax(1)]) ;
        ylim([Rmin(2) Rmax(2)]) ;
        drawnow() ;
        F(z_layer_index) = getframe(gcf) ;
    end
    
    % Make a video
    v = VideoWriter(video_file_name,'Motion JPEG AVI');
    % v.CompressionRatio = 3;
    v.Quality = 100;
    v.FrameRate=2;
    open(v)
    writeVideo(v,F)
    close(v)

    % Make a slide show (a folder of the frames, each a file)
    [video_parent_path, video_base_name] = fileparts(video_file_name) ;
    slide_show_folder_path = fullfile(video_parent_path, video_base_name) ;
    write_slide_show(slide_show_folder_path, F) ;
end

