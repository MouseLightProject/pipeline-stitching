sample_date = '2020-11-26' ;
do_force_computation = false ;

script_folder_path = fileparts(mfilename('fullpath')) ;
memo_folder_path = fullfile(script_folder_path, 'memos') ;

% Build an index of the paths to raw tiles
raw_tile_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling', sample_date) ;
raw_tile_index = compute_or_read_from_memo(memo_folder_path, ...
                                           sprintf('raw-tile-index-%s', sample_date), ...
                                           @()(build_raw_tile_index(raw_tile_path)), ...
                                           do_force_computation) ;
ijk1_from_tile_index = raw_tile_index.ijk1_from_tile_index ;
xyz_from_tile_index = raw_tile_index.xyz_from_tile_index ;  % um
relative_path_from_tile_index = raw_tile_index.relative_path_from_tile_index ;
tile_index_from_ijk1 = raw_tile_index.tile_index_from_ijk1 ;

% Display some features of the raw tile index
tile_lattice_shape = size(tile_index_from_ijk1)
tile_count = length(relative_path_from_tile_index) 

% There's a region in the render that looks like many tiles are "doubled".  What
% is going on here?  One is near [74835.9, 18946.8, 31796.0] in the render.
% There is a doubled process that is very clear.  What raw tile is this from? JW
% says 2020-12-01/01/01916/01916-ngc.0.tif.  What does that look like?

% Find this tile in the index
tile_relative_path ='2020-12-01/01/01916'
this_tile_index = find(strcmp(tile_relative_path, relative_path_from_tile_index))
this_tile_ijk1 = ijk1_from_tile_index(this_tile_index, :)

raw_tile_stack_yxz_flipped = read_16bit_grayscale_tif(fullfile(raw_tile_path, '2020-12-01/01/01916/01916-ngc.0.tif')) ;
raw_tile_stack_yxz = flip(flip(raw_tile_stack_yxz_flipped, 1), 2) ;
raw_tile_stack_yxz_mip = max(raw_tile_stack_yxz, [], 3) ;
f = figure() ;
a = axes(f) ;
imshow(raw_tile_stack_yxz_mip, 'Parent', a) ;    
title_string = sprintf('[ %s ]', strtrim(sprintf('%g ', this_tile_ijk1))) ;
title(a, title_string) ;

% Find this tile in the index, look at the six tiles around it.  Do any of them
% have that big horizontal neurite?

%%

% The z+1 stack MIP looks funny.  Has that big bright side-to-side process from
% the central tile (which is odd b/c that process was pretty much in the middle
% of that stack in z), but that process is shift down in y (higher y).  Maybe a
% half-cut was taken in this plane, and that's an issue?

% Where is that process in z inthe z+1 stack?
other_tile_ijk1 = this_tile_ijk1 + [ 0 0 1 ] ;
tile_index = find(all(other_tile_ijk1==ijk1_from_tile_index, 2)) 
tile_relative_path = relative_path_from_tile_index{tile_index}  %#ok<FNDSB>
imagery_file_relative_path = imagery_file_relative_path_from_relative_path(tile_relative_path, 0) ;  % 0 is channel index
imagery_file_path = fullfile(raw_tile_path, imagery_file_relative_path) 
raw_tile_stack_yxz_flipped = read_16bit_grayscale_tif(imagery_file_path) ;
raw_tile_stack_yxz = flip(flip(raw_tile_stack_yxz_flipped, 1), 2) ;
raw_tile_stack_yxz_mip = max(raw_tile_stack_yxz, [], 3) ;
f = figure() ;
a = axes(f) ;
imshow(raw_tile_stack_yxz_mip, 'Parent', a) ;    
title_string = sprintf('[ %s ]', strtrim(sprintf('%g ', other_tile_ijk1))) ;
title(a, title_string) ;

% Theres a clear +y shift in the z+1 tile of about 200 voxels.  This may be the
% source of the troubles in this sample.

% So that central tile is in the z==22 layer, and the z+1 is in the z==23 layer,
% of course.

% Is 22-->23 always where the problem lies?  That would be nice...

% This next tile is near the transition in the dorsal part of the brain, near the
% midline.
tile_relative_path ='2020-12-01/01/01252'
this_tile_index = find(strcmp(tile_relative_path, relative_path_from_tile_index))
this_tile_ijk1 = ijk1_from_tile_index(this_tile_index, :)

% Show a MIP of the tile
[~,raw_tile_file_name_stem] = fileparts2(tile_relative_path) ;
raw_tile_file_name = [raw_tile_file_name_stem '-ngc.0.tif'] ;
raw_tile_stack_yxz_flipped = read_16bit_grayscale_tif(fullfile(raw_tile_path, tile_relative_path, raw_tile_file_name)) ;
raw_tile_stack_yxz = flip(flip(raw_tile_stack_yxz_flipped, 1), 2) ;
raw_tile_stack_yxz_mip = max(raw_tile_stack_yxz, [], 3) ;
f = figure() ;
a = axes(f) ;
imshow(raw_tile_stack_yxz_mip, 'Parent', a) ;    
title_string = sprintf('[ %s ]', strtrim(sprintf('%g ', this_tile_ijk1))) ;
title(a, title_string) ;
lapwing(raw_tile_stack_yxz) ;   % There's a little "hole" in this stack at slice z=166

% Get the z+1 tile
other_tile_ijk1 = this_tile_ijk1 + [0 0 1] 
tile_index = find(all(other_tile_ijk1==ijk1_from_tile_index, 2)) 
tile_relative_path = relative_path_from_tile_index{tile_index}  %#ok<FNDSB>

% Show a MIP of the tile
[~,raw_tile_file_name_stem] = fileparts2(tile_relative_path) ;
raw_tile_file_name = [raw_tile_file_name_stem '-ngc.0.tif'] ;
raw_tile_stack_yxz_flipped = read_16bit_grayscale_tif(fullfile(raw_tile_path, tile_relative_path, raw_tile_file_name)) ;
raw_tile_stack_yxz = flip(flip(raw_tile_stack_yxz_flipped, 1), 2) ;
raw_tile_stack_yxz_mip = max(raw_tile_stack_yxz, [], 3) ;
f = figure() ;
a = axes(f) ;
imshow(raw_tile_stack_yxz_mip, 'Parent', a) ;    
title_string = sprintf('[ %s ]', strtrim(sprintf('%g ', other_tile_ijk1))) ;
title(a, title_string) ;

% Get the z+2 tile
other_tile_ijk1 = this_tile_ijk1 + [0 0 2] 
tile_index = find(all(other_tile_ijk1==ijk1_from_tile_index, 2)) 
tile_relative_path = relative_path_from_tile_index{tile_index} 

% Show a MIP of the tile
[~,raw_tile_file_name_stem] = fileparts2(tile_relative_path) ;
raw_tile_file_name = [raw_tile_file_name_stem '-ngc.0.tif'] ;
raw_tile_stack_yxz_flipped = read_16bit_grayscale_tif(fullfile(raw_tile_path, tile_relative_path, raw_tile_file_name)) ;
raw_tile_stack_yxz = flip(flip(raw_tile_stack_yxz_flipped, 1), 2) ;
raw_tile_stack_yxz_mip = max(raw_tile_stack_yxz, [], 3) ;
f = figure() ;
a = axes(f) ;
imshow(raw_tile_stack_yxz_mip, 'Parent', a) ;    
title_string = sprintf('[ %s ]', strtrim(sprintf('%g ', other_tile_ijk1))) ;
title(a, title_string) ;




