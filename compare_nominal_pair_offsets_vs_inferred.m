raw_tile_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling', sample_date) ;
pipeline_output_folder_path = sprintf('/nrs/mouselight/pipeline_output/%s', sample_date)  %#ok<NOPTS> 
line_fixed_tile_path = fullfile(pipeline_output_folder_path, 'stage_1_line_fix_output')
landmark_folder_path = fullfile(pipeline_output_folder_path, 'stage_3_descriptor_output')
match_folder_path = fullfile(pipeline_output_folder_path, 'stage_4_point_match_output')
this_folder_path = fileparts(mfilename('fullpath')) ;
memo_folder_path = fullfile(this_folder_path, sprintf('memos-%s', sample_date)) ;


% Build an index of the paths to raw tiles
raw_tile_index = compute_or_read_from_memo(memo_folder_path, ...
                                           'raw-tile-index', ...
                                           @()(build_raw_tile_index(raw_tile_path)), ...
                                           do_force_computation) ;
ijk1_from_tile_index = raw_tile_index.ijk1_from_tile_index ;
nominal_xyz_from_tile_index = raw_tile_index.xyz_from_tile_index ;  % um
relative_path_from_tile_index = raw_tile_index.relative_path_from_tile_index ;
tile_index_from_ijk1 = raw_tile_index.tile_index_from_ijk1 ;

% Display some features of the raw tile index
tile_lattice_shape = size(tile_index_from_ijk1)
tile_count = length(relative_path_from_tile_index) 

% Read channel semantics
[neuron_channel_index, background_channel_index] = read_channel_semantics_file(raw_tile_path) ;
working_channel_index = background_channel_index ;

% Read in a single tile to get the tile shape
middle_tile_index = round(tile_count/2) ;
relative_path = relative_path_from_tile_index{middle_tile_index} ;
imagery_file_relative_path = imagery_file_relative_path_from_relative_path(relative_path, working_channel_index) ;
imagery_file_path = fullfile(line_fixed_tile_path, imagery_file_relative_path) 
raw_tile_stack_yxz_flipped = read_16bit_grayscale_tif(imagery_file_path) ;
tile_shape_jik = size(raw_tile_stack_yxz_flipped) ;
tile_shape_ijk = tile_shape_jik([2 1 3]) ;

% Collect the landmarks for the background channel
ijk0_from_landmark_index_from_tile_index = ...
    compute_or_read_from_memo(memo_folder_path, ...
                              sprintf('landmarks-channel-%d', working_channel_index), ...
                              @()(collect_landmarks(landmark_folder_path, relative_path_from_tile_index, working_channel_index, tile_shape_ijk)), ...
                              do_force_computation) ;

% Count the landmarks in each tile
landmark_count_from_tile_index = cellfun(@(a)(size(a,1)), ijk0_from_landmark_index_from_tile_index) ;
median_landmark_count = median(landmark_count_from_tile_index)

% % Find a tile with a landmark count near the median
% [~, tile_index] = min(abs(landmark_count_from_tile_index - median_landmark_count))

% Look at the tile
tile_index = 12000 ;
relative_path = relative_path_from_tile_index{tile_index} ;
imagery_file_relative_path = imagery_file_relative_path_from_relative_path(relative_path, working_channel_index) ;
imagery_file_path = fullfile(line_fixed_tile_path, imagery_file_relative_path) 
raw_tile_stack_yxz_flipped = read_16bit_grayscale_tif(imagery_file_path) ;
raw_tile_stack_yxz = flip(flip(raw_tile_stack_yxz_flipped, 1), 2) ;
raw_tile_stack_yxz_mip = max(raw_tile_stack_yxz, [], 3) ;
f = figure() ;
a = axes(f) ;
imshow(raw_tile_stack_yxz_mip, 'Parent', a) ;    
title_string = sprintf('Tile %d (%s)', tile_index, relative_path) ;
title(a, title_string) ;

% Add the landmark points to the MIP image
ijk0_from_landmark_index = ijk0_from_landmark_index_from_tile_index{tile_index} ;  % xyz order
ij0_from_landmark_index = ijk0_from_landmark_index(:,1:2) ;  % xy order
ij1_from_landmark_index = ij0_from_landmark_index + 1 ;
hold on ;
plot(ij1_from_landmark_index(:,1), ij1_from_landmark_index(:,2), 'r.') ;
hold off ;

% % Just look at a single plane
% k1 = 133 ;
% k0 = k1 - 1 ;
% raw_tile_plane_yx = raw_tile_stack_yxz(:,:,k1) ;
% ij0_from_plane_landmark_index = ij0_from_landmark_index(ijk0_from_landmark_index(:,3)==k0,:) ;
% f = figure() ;
% a = axes(f) ;
% imshow(raw_tile_plane_yx, 'Parent', a) ;    
% title_string = sprintf('Tile %d (%s), k1=%d', tile_index, relative_path, k1) ;
% title(a, title_string) ;
% hold on ;
% plot(ij0_from_plane_landmark_index(:,1)+1, ij0_from_plane_landmark_index(:,2)+1, 'r.') ;
% hold off ;

