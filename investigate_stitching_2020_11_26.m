sample_date = '2020-11-26' ;
do_force_computation = false ;

script_folder_path = fileparts(mfilename('fullpath')) ;
memo_folder_path = fullfile(script_folder_path, 'memos', sample_date) ;

% Build an index of the paths to raw tiles
%raw_tile_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling', sample_date) ;
raw_tile_path = sprintf('/nearline/mouselight/data/RAW_archive/%s/Tiling', sample_date) ;
raw_tile_index = compute_or_read_from_memo(memo_folder_path, ...
                                           'raw-tile-index', ...
                                           @()(build_raw_tile_index(raw_tile_path)), ...
                                           do_force_computation) ;
ijk1_from_tile_index = raw_tile_index.ijk1_from_tile_index ;
xyz_from_tile_index = raw_tile_index.xyz_from_tile_index ;  % um
relative_path_from_tile_index = raw_tile_index.relative_path_from_tile_index ;
tile_index_from_tile_ijk1 = raw_tile_index.tile_index_from_tile_ijk1 ;

% Display some features of the raw tile index
tile_lattice_shape = size(tile_index_from_tile_ijk1)
tile_count = length(relative_path_from_tile_index) 

% There's a region in the render that looks like many tiles are "doubled".  What
% is going on here?  One is near [74835.9, 18946.8, 31796.0] in the render.
% There is a doubled process that is very clear.  What raw tile is this from? JW
% says 2020-12-01/01/01916/01916-ngc.0.tif.  What does that look like?

this_tile_relative_path = '2020-12-01/01/01916'
imagery_file_relative_path = imagery_file_relative_path_from_relative_path(this_tile_relative_path, 0, '.mj2') ;  % 0 is channel index
imagery_file_path = fullfile(raw_tile_path, imagery_file_relative_path)
raw_tile_stack_yxz_flipped = read_16bit_grayscale_mj2(imagery_file_path) ;
raw_tile_stack_yxz = flip(flip(raw_tile_stack_yxz_flipped, 1), 2) ;
raw_tile_stack_yxz_mip = max(raw_tile_stack_yxz, [], 3) ;
f = figure() ;
a = axes(f) ;
imshow(raw_tile_stack_yxz_mip, 'Parent', a) ;    
dijk = [ 0 0 0 ] ;
title_string = sprintf('[ %s ], %s', strtrim(sprintf('%g ', dijk)), this_tile_relative_path) ;
title(a, title_string) ;

% Find this tile in the index, look at the six tiles around it.  Do any of them
% have that big horizontal neurite?

%%
% Find this tile in the index
this_tile_index = find(strcmp(this_tile_relative_path, relative_path_from_tile_index))
this_tile_ijk1 = ijk1_from_tile_index(this_tile_index, :)
%dijk_from_neighbor_index = [-1 0 0 ; +1 0 0 ; 0 -1 0 ; 0 +1 0 ; 0 0 -1 ; 0 0 +1]
dijk_from_neighbor_index = [0 0 +1]
neighbor_count = size(dijk_from_neighbor_index, 1) ;
ijk1_from_neighbor_index = this_tile_ijk1 + dijk_from_neighbor_index

tile_index_from_neighbor_index = zeros(neighbor_count, 1) ;
for neighbor_index = 1 : neighbor_count ,
    ijk1 = ijk1_from_neighbor_index(neighbor_index, :) ;
    tile_index = find(all(ijk1==ijk1_from_tile_index, 2)) ;
    tile_index_from_neighbor_index(neighbor_index) = tile_index ;
end

relative_path_from_neighbor_index = relative_path_from_tile_index(tile_index_from_neighbor_index) ;

for neighbor_index = 1 : neighbor_count ,
    neighbor_relative_path = relative_path_from_neighbor_index{neighbor_index}
    neighbor_tile_index = tile_index_from_neighbor_index(neighbor_index)
    neighbor_tile_ijk1 = ijk1_from_neighbor_index(neighbor_index, :)
    imagery_file_relative_path = imagery_file_relative_path_from_relative_path(neighbor_relative_path, 0, '.mj2') ;  % 0 is channel index
    imagery_file_path = fullfile(raw_tile_path, imagery_file_relative_path) 
    raw_tile_stack_yxz_flipped = read_16bit_grayscale_mj2(imagery_file_path) ;
    raw_tile_stack_yxz = flip(flip(raw_tile_stack_yxz_flipped, 1), 2) ;
    raw_tile_stack_yxz_mip = max(raw_tile_stack_yxz, [], 3) ;
    f = figure() ;
    a = axes(f) ;  %#ok<LAXES>
    imshow(raw_tile_stack_yxz_mip, 'Parent', a) ;    
    dijk = dijk_from_neighbor_index(neighbor_index, :) ;
    title_string = sprintf('[ %s ] %s', strtrim(sprintf('%g ', dijk)), neighbor_relative_path) ;
    title(a, title_string) ;
end

% The z+1 stack MIP looks funny.  Has that big bright side-to-side process from
% the central tile (which is odd b/c that process was pretty much in the middle
% of that stack in z), but that process is shift down in y (higher y).  Maybe a
% half-cut was taken in this plane, and that's an issue?

% Where is that process in z inthe z+1 stack?
neighbor_index = 1 ;  % z+1
relative_path = relative_path_from_neighbor_index{neighbor_index} ;
imagery_file_relative_path = imagery_file_relative_path_from_relative_path(relative_path, 0, '.mj2') ;  % 0 is channel index
imagery_file_path = fullfile(raw_tile_path, imagery_file_relative_path) 
raw_tile_stack_yxz_flipped = read_16bit_grayscale_mj2(imagery_file_path) ;
raw_tile_stack_yxz = flip(flip(raw_tile_stack_yxz_flipped, 1), 2) ;
lapwing(raw_tile_stack_yxz)

% Looks like it's at very low-index elements of the stack.  (Those are indeed low
% z values in the rendered stack.)  

% Based on lining those two up, there about 154 z voxel offset between them.
% (Thus 154 um.)  Which seems high. 

% According to the .acquisition files for these two tiles, they're offset in z
% by ~159 um.  So that's consistent.  Also according to the .acquisition files,
% they're ideally supposed to have 75 um of overlap in z.  So the overlap here
% is about 251-154==97 um.  That's more than there should be, but not a ton
% more.

%%
% Let's look at the distibution of xyz offsets between tiles and their +1
% neighbors
dxyz_from_tile_index_from_axis_index = nan(3, tile_count, 3) ;  % last index is the offset axis dimension.  NB: We're putting the xyzs in the cols here
is_valid_from_from_tile_index_from_axis_index = false(tile_count, 3) ;
for axis_index = 1 : 3 ,
    dijk = ((1:3)==axis_index) ;
    for central_tile_index = 1 : tile_count ,
        central_tile_xyz = xyz_from_tile_index(central_tile_index, :) ;
        central_tile_ijk1 = ijk1_from_tile_index(central_tile_index,:) ;
        other_tile_ijk1 = central_tile_ijk1 + dijk ;
        if all(other_tile_ijk1 <= tile_lattice_shape) , 
            other_tile_index = index_using_rows(tile_index_from_tile_ijk1, other_tile_ijk1) ;
            if isfinite(other_tile_index) ,
                other_tile_xyz = xyz_from_tile_index(other_tile_index, :) ;
                dxyz = (other_tile_xyz - central_tile_xyz)' ;
                dxyz_from_tile_index_from_axis_index(:,central_tile_index,axis_index) = dxyz ;
                is_valid_from_from_tile_index_from_axis_index(central_tile_index,axis_index) = true ;
            end
        end        
    end
end

%%
% plot those
f = figure('color', 'w') ;
a = axes(f) ;
view(a, -160, 20) ;
axis(a, 'vis3d') ;
grid(a, 'on') ;
for axis_index = 1 : 3 ,
    dxyz_from_tile_index = dxyz_from_tile_index_from_axis_index(:,:,axis_index) ;
    is_valid_from_tile_index = is_valid_from_from_tile_index_from_axis_index(:,axis_index) ;
    dxyz_from_valid_tile_index = dxyz_from_tile_index(:,is_valid_from_tile_index) ;
    h = line('Parent', a, ...
             'XData', dxyz_from_valid_tile_index(1,:), ...
             'YData', dxyz_from_valid_tile_index(2,:), ...
             'ZData', dxyz_from_valid_tile_index(3,:)) ;
    h.LineStyle = 'none' ;
    color_tag = 'rgb' ;
    h.Color = color_tag(axis_index) ;
    h.Marker = '.' ;
end
%tag = 'xyz' ;
%xlim(a, [-1 400]) ;
%ylim(a, [-1 600]) ;
%zlim(a, [-1 400]) ;
xlabel(a, 'x (um)') ;
ylabel(a, 'y (um)') ;
zlabel(a, 'z (um)') ;
%a.DataAspectRatio = [1 1 1] ;
%title(a, tag(axis_index)) ;

% Looks like the z offset can vary a lot, the x and y offsets are only a few um

%%
axis_index = 3 ;
dxyz_from_tile_index = dxyz_from_tile_index_from_axis_index(:,:,axis_index) ;
is_valid_from_tile_index = is_valid_from_from_tile_index_from_axis_index(:,axis_index) ;
dxyz_from_valid_tile_index = dxyz_from_tile_index(:,is_valid_from_tile_index) ;
xyz_from_valid_tile_index = (xyz_from_tile_index(is_valid_from_tile_index,:))' ;  % moved xyzs into cols
z_from_valid_tile_index = xyz_from_valid_tile_index(axis_index,:) ;
dz_from_valid_tile_index = dxyz_from_valid_tile_index(axis_index,:) ;
f = figure('color', 'w') ;
a = axes(f) ;
plot(z_from_valid_tile_index-min(z_from_valid_tile_index), dz_from_valid_tile_index, '.k') ;
xlabel('z (relative to min, um)') ;
ylabel('dz (um)') ;
xlim([0 11000]) ;
ylim([0 220]) ;


%%
ijk1_from_valid_tile_index = (ijk1_from_tile_index(is_valid_from_tile_index, :))' ;  % moved ijk1s into cols
k1_from_valid_tile_index = ijk1_from_valid_tile_index(axis_index, :) ;
f = figure('color', 'w') ;
a = axes(f) ;
plot(k1_from_valid_tile_index, dz_from_valid_tile_index, '.k') ;
xlabel('lattice z') ;
ylabel('dz (um)') ;
ylim([0 220]) ;



