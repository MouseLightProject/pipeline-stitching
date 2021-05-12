sample_date = '2020-11-26'  %#ok<NOPTS>
do_force_computation = false ;
path_kind_to_use_for_imagery = 'raw' ;  % this will end up using the line-fixed imagery, since the raw has been replaced with line-fixed

raw_tile_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling', sample_date) ;
pipeline_output_folder_path = sprintf('/nrs/mouselight/pipeline_output/%s', sample_date)  %#ok<NOPTS> 
line_fixed_tile_path = fullfile(pipeline_output_folder_path, 'stage_1_line_fix_output')
landmark_folder_path = fullfile(pipeline_output_folder_path, 'stage_3_descriptor_output')
match_folder_path = fullfile(pipeline_output_folder_path, 'stage_4_point_match_output')
this_folder_path = fileparts(mfilename('fullpath')) ;
memo_folder_path = fullfile(this_folder_path, sprintf('memos-%s', sample_date)) ;

if ~exist('path_kind_to_use_for_imagery', 'var') || isempty(path_kind_to_use_for_imagery) || strcmp(path_kind_to_use_for_imagery, 'line-fixed') ,
    imagery_tile_path = line_fixed_tile_path ;
elseif strcmp(path_kind_to_use_for_imagery, 'raw') ,
    imagery_tile_path = raw_tile_path ;
else
    error('Variable path_to_use_for_imagery is ''%s''.  Should be either ''line-fixed'' or ''raw''.', path_kind_to_use_for_imagery) ;
end

% Build an index of the paths to raw tiles
raw_tile_index = compute_or_read_from_memo(memo_folder_path, ...
                                           'raw-tile-index', ...
                                           @()(build_raw_tile_index(raw_tile_path)), ...
                                           do_force_computation) ;
tile_shape_ijk = raw_tile_index.tile_shape_ijk
spacing_um_xyz = raw_tile_index.spacing_um_xyz   % um
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
interesting_tile_index = 12000 ;  % seems to have matches not at landmarks
interesting_tile_index = 11000 ;
interesting_tile_index = find(strcmp('2020-12-01/01/01916', relative_path_from_tile_index))  % one we investigated earlier

relative_path = relative_path_from_tile_index{interesting_tile_index} ;
imagery_file_relative_path = imagery_file_relative_path_from_relative_path(relative_path, working_channel_index) ;
imagery_file_path = fullfile(imagery_tile_path, imagery_file_relative_path) 
raw_tile_stack_yxz_flipped = read_16bit_grayscale_tif(imagery_file_path) ;
raw_tile_stack_yxz = flip(flip(raw_tile_stack_yxz_flipped, 1), 2) ;
raw_tile_stack_yxz_mip = max(raw_tile_stack_yxz, [], 3) ;
f = figure('color', 'w') ;
a = axes(f) ;
imshow(raw_tile_stack_yxz_mip, ...
       'Parent', a, ...
       'DisplayRange', [min(min(raw_tile_stack_yxz_mip)) max(max(raw_tile_stack_yxz_mip))], ...
       'Colormap', parula(256)) ;
title_string = sprintf('Tile %d (%s)', interesting_tile_index, relative_path) ;
title(a, title_string) ;
colorbar(a) ;

% Add the landmark points to the MIP image
self_ijk0_from_match_index = ijk0_from_landmark_index_from_tile_index{interesting_tile_index} ;  % xyz order
self_ij0_from_match_index = self_ijk0_from_match_index(:,1:2) ;  % xy order
self_ij1_from_match_index = self_ij0_from_match_index + 1 ;
hold(a,'on') ;
plot(a, self_ij1_from_match_index(:,1), self_ij1_from_match_index(:,2), 'r.') ;
hold(a, 'off') ;

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

%%
%
% Matches
%

% Count the number of z-face pairs 
has_z_plus_1_tile_from_ijk1 = false(tile_lattice_shape) ;
has_z_plus_1_tile_from_ijk1(:,:,1:end-1) = isfinite(tile_index_from_ijk1(:,:,1:end-1)) & isfinite(tile_index_from_ijk1(:,:,2:end)) ;
pair_count = sum(sum(sum(has_z_plus_1_tile_from_ijk1))) 

% Want to know if has z+1 tile from tile_index
has_z_plus_1_tile_from_tile_index = false(tile_count,1) ;
for tile_index = 1 : tile_count ,
    ijk1 = ijk1_from_tile_index(tile_index,:) ;
    has_z_plus_1_tile = has_z_plus_1_tile_from_ijk1(ijk1(1), ijk1(2), ijk1(3)) ;
    has_z_plus_1_tile_from_tile_index(tile_index) = has_z_plus_1_tile ;
end
assert(sum(has_z_plus_1_tile_from_tile_index) == pair_count) ;
self_tile_index_from_pair_index = find(has_z_plus_1_tile_from_tile_index) ;

% Collect the z-face matches from disk
match_info = ...
    compute_or_read_from_memo(...
        memo_folder_path, ...
        sprintf('z-face-matches-%s-channel-%d', sample_date, background_channel_index), ...
        @()(collect_z_face_matches(match_folder_path, ...
                                   relative_path_from_tile_index, ...
                                   background_channel_index, ...
                                   has_z_plus_1_tile_from_tile_index, ...
                                   tile_shape_ijk)), ...
        do_force_computation) ;
self_ijk0_from_match_index_from_tile_index = match_info.self_ijk0_from_match_index_from_tile_index ;
neighbor_ijk0_from_match_index_from_tile_index = match_info.neighbor_ijk0_from_match_index_from_tile_index ;

% Count the matches per-tile
z_match_count_from_tile_index = cellfun(@(a)(size(a,1)), self_ijk0_from_match_index_from_tile_index) ;
check_z_match_count_from_tile_index = cellfun(@(a)(size(a,1)), neighbor_ijk0_from_match_index_from_tile_index) ;
assert(isequal(z_match_count_from_tile_index, check_z_match_count_from_tile_index)) ;
z_match_count_from_pair_index = z_match_count_from_tile_index(has_z_plus_1_tile_from_tile_index) ;

max_z_match_count = max(z_match_count_from_pair_index)
median_z_match_count = median(z_match_count_from_pair_index)

% % Find a tile with a match count near the median
% [~, pair_index] = min(abs(z_match_count_from_pair_index - median_z_match_count))
% tile_index = tile_index_from_pair_index(pair_index)

% % Look at the tile
% %tile_index = 12000 ;
% relative_path = relative_path_from_tile_index{tile_index} ;
% imagery_file_relative_path = imagery_file_relative_path_from_relative_path(relative_path, working_channel_index) ;
% imagery_file_path = fullfile(imagery_tile_path, imagery_file_relative_path) 
% raw_tile_stack_yxz_flipped = read_16bit_grayscale_tif(imagery_file_path) ;
% raw_tile_stack_yxz = flip(flip(raw_tile_stack_yxz_flipped, 1), 2) ;
% raw_tile_stack_yxz_mip = max(raw_tile_stack_yxz, [], 3) ;
% f = figure() ;
% a = axes(f) ;
% imshow(raw_tile_stack_yxz_mip, ...
%        'Parent', a, ...
%        'DisplayRange', [min(min(raw_tile_stack_yxz_mip)) max(max(raw_tile_stack_yxz_mip))], ...
%        'Colormap', parula(256)) ;
% title_string = sprintf('Tile %d (%s)', tile_index, relative_path) ;
% title(a, title_string) ;
% colorbar(a) ;

% Add the matched points to the MIP image
self_ijk0_from_match_index = self_ijk0_from_match_index_from_tile_index{interesting_tile_index} ;  % xyz order
self_ij0_from_match_index = self_ijk0_from_match_index(:,1:2) ;  % xy order
self_ij1_from_match_index = self_ij0_from_match_index + 1 ;
hold(a,'on') ;
plot(a, self_ij1_from_match_index(:,1), self_ij1_from_match_index(:,2), 'ro') ;
hold(a, 'off') ;

%
% Show the matched points in the neighbor tile
%
interesting_tile_ijk1 = ijk1_from_tile_index(interesting_tile_index, :) ;
neighbor_tile_ijk1 = interesting_tile_ijk1 + [0 0 1] ;
neighbor_tile_index = index_using_rows(tile_index_from_ijk1, neighbor_tile_ijk1)

% Look at the neighbor tile
relative_path = relative_path_from_tile_index{neighbor_tile_index} ;
imagery_file_relative_path = imagery_file_relative_path_from_relative_path(relative_path, working_channel_index) ;
imagery_file_path = fullfile(imagery_tile_path, imagery_file_relative_path) 
raw_tile_stack_yxz_flipped = read_16bit_grayscale_tif(imagery_file_path) ;
raw_tile_stack_yxz = flip(flip(raw_tile_stack_yxz_flipped, 1), 2) ;
raw_tile_stack_yxz_mip = max(raw_tile_stack_yxz, [], 3) ;
f = figure('color', 'w') ;
a = axes(f) ;
imshow(raw_tile_stack_yxz_mip, ...
       'Parent', a, ...
       'DisplayRange', [min(min(raw_tile_stack_yxz_mip)) max(max(raw_tile_stack_yxz_mip))], ...
       'Colormap', parula(256)) ;
title_string = sprintf('Tile %d (%s)', neighbor_tile_index, relative_path) ;
title(a, title_string) ;
colorbar(a) ;

% Add the landmark points to the MIP image
self_ijk0_from_landmark_index = ijk0_from_landmark_index_from_tile_index{neighbor_tile_index} ;  % xyz order
self_ij0_from_landmark_index = self_ijk0_from_landmark_index(:,1:2) ;  % xy order
self_ij1_from_landmark_index = self_ij0_from_landmark_index + 1 ;
hold(a, 'on') ;
plot(a, self_ij1_from_landmark_index(:,1), self_ij1_from_landmark_index(:,2), 'r.') ;
hold(a, 'off') ;

% Add the matched points *from this tile* to the image
neighbor_ijk0_from_match_index = neighbor_ijk0_from_match_index_from_tile_index{interesting_tile_index} ;  % xyz order
neighbor_ij0_from_match_index = neighbor_ijk0_from_match_index(:,1:2) ;  % xy order
neighbor_ij1_from_match_index = neighbor_ij0_from_match_index + 1 ;
hold(a, 'on') ;
plot(a, neighbor_ij1_from_match_index(:,1), neighbor_ij1_from_match_index(:,2), 'ro') ;
hold(a, 'off') ;

% OK, those look good

% For each tile pair, compute the tile index of the z+1 tile
% Already have the self-tile index for each pair (above)
self_tile_ijk1_from_pair_index = ijk1_from_tile_index(self_tile_index_from_pair_index,:) ;
neighbor_tile_ijk1_from_pair_index = self_tile_ijk1_from_pair_index + [0 0 1] ;
neighbor_tile_index_from_pair_index = index_using_rows(tile_index_from_ijk1, neighbor_tile_ijk1_from_pair_index) ;

% For each tile pair, compute the nominal offset between them.
self_tile_nominal_xyz_from_pair_index = nominal_xyz_from_tile_index(self_tile_index_from_pair_index, :) ;
neighbor_tile_nominal_xyz_from_pair_index = nominal_xyz_from_tile_index(neighbor_tile_index_from_pair_index, :) ;
nominal_dxyz_from_pair_index = neighbor_tile_nominal_xyz_from_pair_index - self_tile_nominal_xyz_from_pair_index

% For each tile pair, compute an offset based on the matches
empirical_dxyz_from_pair_index = zeros(pair_count, 3) ;
for pair_index = 1 : pair_count ,
    self_tile_index = self_tile_index_from_pair_index(pair_index, :) ;
    
    self_ijk0_from_match_index = self_ijk0_from_match_index_from_tile_index{self_tile_index} ;  % xyz order
    neighbor_ijk0_from_match_index = neighbor_ijk0_from_match_index_from_tile_index{self_tile_index} ;  % xyz order
    dijk_from_match_index = self_ijk0_from_match_index - neighbor_ijk0_from_match_index ;
    consensus_dijk = mean(dijk_from_match_index, 1) ;
    empirical_dxyz = spacing_um_xyz .* consensus_dijk ;
    empirical_dxyz_from_pair_index(pair_index,:) = empirical_dxyz ;    
end

empirical_from_nominal_dxyz_from_pair_index = empirical_dxyz_from_pair_index - nominal_dxyz_from_pair_index

f = figure('color', 'w') ;
a = axes(f) ;
plot3(a, ...
      nominal_dxyz_from_pair_index(:,1), ...
      nominal_dxyz_from_pair_index(:,2), ...
      nominal_dxyz_from_pair_index(:,3), ...
      'b.') ;
grid(a, 'on') ;
axis(a, 'vis3d') ;
xlabel('x (um)') ;
ylabel('y (um)') ;
zlabel('z (um)') ;

f = figure('color', 'w') ;
a = axes(f) ;
plot3(a, ...
      empirical_dxyz_from_pair_index(:,1), ...
      empirical_dxyz_from_pair_index(:,2), ...
      empirical_dxyz_from_pair_index(:,3), ...
      'r.') ;
grid(a, 'on') ;
axis(a, 'vis3d') ;
xlabel('x (um)') ;
ylabel('y (um)') ;
zlabel('z (um)') ;

f = figure('color', 'w') ;
a = axes(f) ;
plot3(a, ...
      empirical_from_nominal_dxyz_from_pair_index(:,1), ...
      empirical_from_nominal_dxyz_from_pair_index(:,2), ...
      empirical_from_nominal_dxyz_from_pair_index(:,3), ...
      '.') ;
grid(a, 'on') ;
axis(a, 'vis3d') ;
xlabel('x (um)') ;
ylabel('y (um)') ;
zlabel('z (um)') ;




% How big is the empirical shift in x-y plane?
empirical_dxy_from_pair_index = empirical_dxyz_from_pair_index(:,1:2) ;
empirical_dxy_norm_from_pair_index = vecnorm(empirical_dxy_from_pair_index, 2,  2) ;
median_empirical_dxy_norm = median(empirical_dxy_norm_from_pair_index, 1, 'omitnan')


% Make a histogram
max_empirical_dxy_norm = max(empirical_dxy_norm_from_pair_index, [], 'omitnan')

% Compute a histogram of those
bin_edges = 0 : 5 : 400 ;
bin_centers = ( bin_edges(1:end-1)+bin_edges(2:end) ) / 2 ;
counts_from_bin_index = histcounts(empirical_dxy_norm_from_pair_index, bin_edges) ;
%assert(sum(counts_from_bin_index) == tile_count) ;
fraction_from_bin_index = counts_from_bin_index / pair_count ;

% Plot the histogram
f = figure('color', 'w') ;
a = axes(f) ;
b = bar(a, bin_centers, 100*fraction_from_bin_index) ;
b.EdgeColor = 'none' ;
xlabel('Norm of empirical shift in xy plane') ;
ylabel('Fraction of pairs (%)') ;
%ylim([0 100]) ;
title(sample_date) ;

% What about our tile pair?
interesting_pair_index = find(self_tile_index_from_pair_index==interesting_tile_index)
empirical_dxy_norm_for_interesting_pair = empirical_dxy_norm_from_pair_index(interesting_pair_index)
fraction_above_50 = mean(empirical_dxy_norm_from_pair_index>50, 1)

% What's the max shift?
[max_dxy_norm, max_dxy_norm_pair_index] = max(empirical_dxy_norm_from_pair_index)

% plot the dxy norm vs the match count
f = figure('color', 'w') ;
a = axes(f) ;
plot(empirical_dxy_norm_from_pair_index, z_match_count_from_pair_index, '.') ;
xlabel('Norm of empirical shift in xy plane (um)') ;
ylabel('Matches per pair')

% Let's look at the ones that are high in both (which we don't see in 9-26) see
% if there's a spatial pattern.

is_problem_from_pair_index = (empirical_dxy_norm_from_pair_index > 50) & (z_match_count_from_pair_index > 20)
problem_pair_count = sum(is_problem_from_pair_index)
pair_index_from_problem_pair_index = find(is_problem_from_pair_index) ;
self_tile_ijk1_from_problem_pair_index = self_tile_ijk1_from_pair_index(pair_index_from_problem_pair_index,:) ;

f = figure('color', 'w') ;
a = axes(f) ;
%plot3(self_tile_ijk1_from_pair_index(:,1), self_tile_ijk1_from_pair_index(:,2), self_tile_ijk1_from_pair_index(:,3), 'b.') ;
%hold on ;
plot3(self_tile_ijk1_from_problem_pair_index(:,1), self_tile_ijk1_from_problem_pair_index(:,2), self_tile_ijk1_from_problem_pair_index(:,3), 'ro') ;
%hold off ;
xlabel('x') ;
ylabel('y') ;
zlabel('z') ;
grid(a, 'on')
axis(a, 'vis3d') ;

