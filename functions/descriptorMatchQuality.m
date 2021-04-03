function descriptorMatchQuality(regpts,scopeparams,scopeloc,videofile)
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
tile_count = size(regpts,2) ;
afftile = zeros(3, 4, tile_count);  % will hold an affine transform for each tile, constructed from the global um-from-voxels transform (3x3) and the translation from scopeloc (the nominal tile offset)
afftransform = scopeparams.affineglFC;  % 3x3, nm
mflip = eye(3);%mflip(3,3)=1;
afftransform = afftransform*mflip;
for tile_index = 1:tile_count ,
    %%
    % form an affine matrix
    tform = [afftransform scopeloc.loc(tile_index,:)'*1e6];  % 1e6 converts mm to nm
    afftile(:,:,tile_index) = tform;
end
%%
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3));  % tile_count x 7, gives tile index of neighbor in each direction, [id -x -y +x +y -z +z] order 
%checkthese = [1 4 5 7]; % 0 - below

% For each tile, compute the offset to the z+1 tile
dz_to_plus_z_tile = zeros(1, tile_count) ;  % mm
for tile_index = 1 : tile_count ,
    these_neighbors = neighbors(tile_index, :) ;
    if any(isnan(these_neighbors([1 7])))
        dz_to_plus_z_tile(tile_index) = nan;
    else
        dz_to_plus_z_tile(tile_index) = diff(scopeloc.loc(these_neighbors([1 7]),3)) ;  % mm
    end
end

%%
lattice_z_from_tile_index = scopeloc.gridix(:,3) ;  % column
lattice_z_from_z_layer_index = unique(lattice_z_from_tile_index) ;
%minimum_lattice_z = lattice_z_from_z_layer_index(1) ;
maximum_lattice_z = lattice_z_from_z_layer_index(end) ;
z_layer_count = length(lattice_z_from_z_layer_index) ;
XYZ_t = cell(1, z_layer_count) ;
for z_layer_index = 1 : z_layer_count ,
    %lattice_z = lattice_z_values(1:end-1)' ,
    lattice_z = lattice_z_from_z_layer_index(z_layer_index) ;
    fprintf('    Layer %d of %d\n', lattice_z, maximum_lattice_z) ;

    is_in_this_layer_from_tile_index = (lattice_z_from_tile_index'==lattice_z);
    %cnt = 0;
    %Npts = 0;
    %clear layer layerp1 poslayer_t poslayer_tp1
    %layer = [] ;
    tile_index_from_foo_index = find(is_in_this_layer_from_tile_index) ;  % a "foo" is a tile in this layer
    foo_count = length(tile_index_from_foo_index) ;
    match_count_from_foo_index = zeros(1, foo_count) ;
    for foo_index = 1 : foo_count ,
        tile_index = tile_index_from_foo_index(foo_index) ;
        regpts_for_this_tile = regpts{tile_index} ;
        if ~isfield(regpts_for_this_tile,'X') ,
            match_count = 0 ;
        else
            match_count = size(regpts_for_this_tile.X,1) ;
        end
        match_count_from_foo_index(foo_index) = match_count ;
    end    
    is_match_count_nonzero_from_foo_index = (match_count_from_foo_index>0) ;
    foo_index_from_bar_index = find(is_match_count_nonzero_from_foo_index) ;  % a "bar" is a tile in this layer with a non-zero match count
    bar_count = length(foo_index_from_bar_index) ;
    tile_index_from_bar_index = tile_index_from_foo_index(foo_index_from_bar_index) ;
    layer = cell(1, bar_count) ;
    poslayer_t = cell(1, bar_count) ;
    for bar_index = 1 : bar_count ,
        tile_index = tile_index_from_bar_index(bar_index) ;
        regpts_for_this_tile = regpts{tile_index} ;
        neigs = regpts_for_this_tile.neigs ;
        layer{bar_index} = regpts_for_this_tile.X ;
        %layerp1{bar_index} = regpts{id_ix}.Y;  %#ok<AGROW>
        % apply tforms
        poslayer_t_for_this_tile = ...
            [regpts_for_this_tile.X ones(size(regpts_for_this_tile.X,1),1)]*afftile(:,:,neigs(1))';  % map the z-paired fiducials from this tile to render coords
        %Layer_tp1 = [regpts_for_this_tile.Y ones(size(regpts_for_this_tile.Y,1),1)]*afftile(:,:,neigs(4))';
        %Npts = Npts+size(Layer_t,1);
        poslayer_t{bar_index} = poslayer_t_for_this_tile;
        %poslayer_tp1{bar_index} = Layer_tp1;  %#ok<AGROW>
    end
    if bar_count == 0 ,
        XYZ_t_for_this_layer = zeros(0, 3) ;
    else
        XYZ_t_for_this_layer = cat(1,poslayer_t{:}) ;  % collect the matched fiducials for all tiles in this layer
    end
    XYZ_t{z_layer_index} = XYZ_t_for_this_layer ;
end
trmax=cellfun(@max,XYZ_t,'UniformOutput',false);
trmin=cellfun(@min,XYZ_t,'UniformOutput',false);
if length(scopeparams) < 2 ,
    imsize_um = scopeparams.imsize_um;
else
    imsize_um = scopeparams(1).imsize_um;
end

Rmax = max(cat(1,trmax{:}))+imsize_um*1e3.*[1 1 0];
Rmin = min(cat(1,trmin{:}))-imsize_um*1e3.*[1 1 0];
% plot desctiptors
myfig = figure() ;
%
%figure(myfig), cla, clf
%hold(myaxes,'on') ;
%loops = latticeZRange(end)-latticeZRange(1);
%F(loops) = struct('cdata',[],'colormap',[]);
%frame_count = z_layer_count ;  % no frame for the highest-z layer
F(z_layer_count) = struct('cdata',[], 'colormap',[]) ;
for z_layer_index = 1 : z_layer_count  ,
    lattice_z = lattice_z_from_z_layer_index(z_layer_index) ;
    %%
%     idxtest = sliceinds(375)
    is_in_this_layer_from_tile_index = (lattice_z_from_tile_index'==lattice_z);
    sliceinds = find(is_in_this_layer_from_tile_index);
    if any(sliceinds)
    end
    myaxes = axes('Parent', myfig, ...
                  'Ydir', 'reverse', ...
                  'XLim', [Rmin(1) Rmax(1)], ...
                  'YLim', [Rmin(2) Rmax(2)], ...
                  'View', [0 90], ...
                  'Color',[1 1 1]*.8) ;
    %hold(myaxes, 'on') ;
    XYZ_t_this = XYZ_t{z_layer_index} ;
    %X = XYZ_t_this(:,1);
    %Y = XYZ_t_this(:,2);
    %Z = XYZ_t_this(:,3);
    %XYZ_tp1_this = XYZ_tp1{t} ;    
    %U = XYZ_tp1_this(:,1)-XYZ_t_this(:,1);
    %V = XYZ_tp1_this(:,1)-XYZ_t_this(:,2);
    %W = XYZ_tp1_this(:,1)-XYZ_t_this(:,3);
    % plot boxes
    %disp(['    Layer ' num2str(t) ' of ' num2str(max(lattice_z_from_tile_index))]);
    fprintf('    Layer %d of %d\n', lattice_z, maximum_lattice_z) ;
    %plot3(myaxes, XYZ_t_this(:,1),XYZ_t_this(:,2),XYZ_t_this(:,3),'.') ;
    line('Parent', myaxes, ...
         'XData', XYZ_t_this(:,1), ...
         'YData', XYZ_t_this(:,2), ...
         'LineStyle', 'none', ...
         'Color', 'b', ...
         'Marker', '.') ;
%             'ZData', XYZ_t_this(:,3), ...

    x = scopeloc.loc(is_in_this_layer_from_tile_index,1)*1e6;
    y = scopeloc.loc(is_in_this_layer_from_tile_index,2)*1e6;
    w = imsize_um(1)*ones(sum(is_in_this_layer_from_tile_index),1)*1e3;
    h = imsize_um(2)*ones(sum(is_in_this_layer_from_tile_index),1)*1e3;
    for tile_index=1:sum(is_in_this_layer_from_tile_index)
        rectangle(myaxes, 'Position', [x(tile_index) y(tile_index) w(tile_index) h(tile_index)]) ;
        mystr = sprintf('%05d\n(%d)\n%d:%d',sliceinds(tile_index),tile_index,scopeloc.gridix(sliceinds(tile_index),1),scopeloc.gridix(sliceinds(tile_index),2)) ;
        %text(myaxes, x(ii)+w(ii)/2,y(ii)+h(ii)/2,mystr,'Color','r','HorizontalAlignment','center') ;
    end
    title(myaxes, [num2str(lattice_z),' - ', num2str(max(dz_to_plus_z_tile(is_in_this_layer_from_tile_index)))]) ;
    %set(myaxes,'Ydir','reverse')
    %xlim(myaxes, [Rmin(1) Rmax(1)])
    %ylim(myaxes, [Rmin(2) Rmax(2)])
    text(myaxes, ...
         Rmin(1)+7e5, Rmax(2)-5e5, ...
         num2str(lattice_z), ...
        'FontSize',60, ...
        'Color','k', ...
        'HorizontalAlignment','center') ; 
    %set(myaxes, 'View', [0 90]) ;
    %set(myaxes, 'Color',[1 1 1]*.8) ;
    drawnow() ;
    F(z_layer_index) = getframe(myfig) ;
    delete(myaxes) ;
end
delete(myfig) ;
%%
v = VideoWriter(videofile,'Motion JPEG AVI');
% v.CompressionRatio = 3;
v.FrameRate=5;
v.open() ;
v.writeVideo(F) ;
v.close() ;
%%

% Write out a slide show version
[video_parent_path, video_base_name] = fileparts(videofile) ;
slide_show_folder_path = fullfile(video_parent_path, video_base_name) ;
write_slide_show(slide_show_folder_path, F) ;

end









