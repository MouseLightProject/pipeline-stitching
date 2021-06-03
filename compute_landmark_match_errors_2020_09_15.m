sample_date = '2020-09-15'  %#ok<NOPTS>
path_kind_to_use_for_imagery = 'line-fixed' 
do_force_computation = false ;

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
raw_tile_ijk1_from_tile_index = raw_tile_index.raw_ijk1_from_tile_index ;
tile_ijk1_from_tile_index = raw_tile_index.ijk1_from_tile_index ;
nominal_xyz_from_tile_index = raw_tile_index.xyz_from_tile_index ;  % um
relative_path_from_tile_index = raw_tile_index.relative_path_from_tile_index ;
tile_index_from_tile_ijk1 = raw_tile_index.tile_index_from_tile_ijk1 ;
does_exist_from_tile_ijk1 = raw_tile_index.does_exist_from_tile_ijk1;

% Display some features of the raw tile index
tile_lattice_shape = size(tile_index_from_tile_ijk1)
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



%
% Matches
%

% Count the number of z-face pairs 
has_z_plus_1_tile_from_ijk1 = false(tile_lattice_shape) ;
has_z_plus_1_tile_from_ijk1(:,:,1:end-1) = isfinite(tile_index_from_tile_ijk1(:,:,1:end-1)) & isfinite(tile_index_from_tile_ijk1(:,:,2:end)) ;
pair_count = sum(sum(sum(has_z_plus_1_tile_from_ijk1))) 

% Want to know if has z+1 tile from tile_index
has_z_plus_1_tile_from_tile_index = false(tile_count,1) ;
for tile_index = 1 : tile_count ,
    ijk1 = tile_ijk1_from_tile_index(tile_index,:) ;
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

% Make a stack reflecting the z-matches for each tile
tile_ijk1_from_pair_index = tile_ijk1_from_tile_index(self_tile_index_from_pair_index, :) ;
z_match_count_from_tile_ijk1 = nan(size(tile_index_from_tile_ijk1)) ;
for pair_index = 1 : pair_count ,
    tile_ijk1 = tile_ijk1_from_pair_index(pair_index,:) ;
    z_match_count_from_tile_ijk1(tile_ijk1(1), tile_ijk1(2), tile_ijk1(3)) = z_match_count_from_pair_index(pair_index) ;    
end

% Make and show a montage of that
z_match_count_from_tile_ijk1_montage = montage_from_stack_ijk(z_match_count_from_tile_ijk1) ;
f = figure('color', 'w') ;
a = axes(f) ;
imagesc(z_match_count_from_tile_ijk1_montage, [0 600]) ;
colorbar(a) ;
title('Z-Match Count') 
drawnow



%
% Compute errors
%

% Collect information about neighbors of each tile
neighbor_count = 3 ;  % x+1, y+1, z+1
has_neighbor_from_neighbor_index_from_tile_index = false(neighbor_count, tile_count) ;
neighbor_tile_index_from_neighbor_index_from_tile_index = nan(neighbor_count, tile_count) ;
for tile_index = 1 : tile_count ,
    tile_ijk1 = tile_ijk1_from_tile_index(tile_index,:) ;
    for neighbor_index = 1 : neighbor_count ,
        neighbor_dijk1 = ((1:3)==neighbor_index) ;
        neighbor_ijk1 = tile_ijk1 + neighbor_dijk1 ;
        if all(neighbor_ijk1 <= tile_lattice_shape) ,
            does_neighbor_exist = index_using_rows(does_exist_from_tile_ijk1, neighbor_ijk1) ;
            has_neighbor_from_neighbor_index_from_tile_index(neighbor_index, tile_index) = does_neighbor_exist ;
            if does_neighbor_exist ,
                neighbor_tile_index = index_using_rows(tile_index_from_tile_ijk1, neighbor_ijk1) ;
                neighbor_tile_index_from_neighbor_index_from_tile_index(neighbor_index, tile_index) = neighbor_tile_index ;
            end
        end
    end
end
   
% Collect the matches
self_ijk0_from_match_index_from_neighbor_index_from_tile_index = cell(neighbor_count, tile_count) ;
self_ijk0_from_match_index_from_neighbor_index_from_tile_index(1:2,:) = {zeros(0,3)} ;
self_ijk0_from_match_index_from_neighbor_index_from_tile_index(3,:) = self_ijk0_from_match_index_from_tile_index ;
neighbor_ijk0_from_match_idx_from_neighbor_idx_from_tile_idx = cell(neighbor_count, tile_count) ;
neighbor_ijk0_from_match_idx_from_neighbor_idx_from_tile_idx(1:2,:) = {zeros(0,3)} ;
neighbor_ijk0_from_match_idx_from_neighbor_idx_from_tile_idx(3,:) = neighbor_ijk0_from_match_index_from_tile_index ;
match_count_from_neighbor_index_from_tile_index = zeros(neighbor_count, tile_count) ;
match_count_from_neighbor_index_from_tile_index(3,:) = z_match_count_from_tile_index ;

% Assemble the stage-based affine transforms
stage_affine_transform_from_tile_index = zeros(3, 4, tile_count) ;
stage_affine_transform_from_tile_index(1,1,:) = spacing_um_xyz(1) ;
stage_affine_transform_from_tile_index(2,2,:) = spacing_um_xyz(2) ;
stage_affine_transform_from_tile_index(3,3,:) = spacing_um_xyz(3) ;
stage_affine_transform_from_tile_index(:,4,:) = nominal_xyz_from_tile_index' ;

% Compute the match errors for each tile, neighbor
stage_match_sse_from_neighbor_index_from_tile_index = ...
        compute_affine_landmark_match_error(self_ijk0_from_match_index_from_neighbor_index_from_tile_index, ...
                                            has_neighbor_from_neighbor_index_from_tile_index, ...
                                            neighbor_tile_index_from_neighbor_index_from_tile_index, ...
                                            neighbor_ijk0_from_match_idx_from_neighbor_idx_from_tile_idx, ...
                                            stage_affine_transform_from_tile_index) ;

stage_z_match_sse_from_tile_index = reshape(stage_match_sse_from_neighbor_index_from_tile_index(3,:), [tile_count 1]) ;
stage_z_match_sse_from_pair_index = stage_z_match_sse_from_tile_index(has_z_plus_1_tile_from_tile_index) ;
stage_z_match_mse_from_pair_index = stage_z_match_sse_from_pair_index ./ z_match_count_from_pair_index ;

total_stage_z_match_sse = sum(stage_z_match_sse_from_pair_index)
stage_z_match_rmse_from_pair_index = sqrt(stage_z_match_mse_from_pair_index) ;
max_stage_z_match_rmse = max(stage_z_match_rmse_from_pair_index)
median_stage_z_match_rmse = median(stage_z_match_rmse_from_pair_index, 'omitnan')

% make a z-match RMSE stack
stage_z_match_rmse_from_tile_ijk1 = nan(size(tile_index_from_tile_ijk1)) ;
for pair_index = 1 : pair_count ,
    tile_ijk1 = tile_ijk1_from_pair_index(pair_index,:) ;
    stage_z_match_rmse_from_tile_ijk1(tile_ijk1(1), tile_ijk1(2), tile_ijk1(3)) = stage_z_match_rmse_from_pair_index(pair_index) ;    
end

stage_z_match_rmse_from_tile_ijk1_montage = montage_from_stack_ijk(stage_z_match_rmse_from_tile_ijk1) ;
f = figure('color', 'w') ;
a = axes(f) ;
imagesc(stage_z_match_rmse_from_tile_ijk1_montage, [0 100]) ;
colorbar(a) ;
title('Stage Z-Match RMSE (um)') 
drawnow

% scatter plot of RMSE and matches per pair
f = figure('color', 'w') ;
a = axes(f) ;
plot(a, z_match_count_from_pair_index, stage_z_match_rmse_from_pair_index, '.') ;
xlabel(a, 'Stage Z-match count per tile pair') ;
ylabel(a, 'Stage Z-match RMSE per tile pair') ;

% Load the per-tile affine transforms the stitcher outputs
stitching_output_folder_path = fullfile(memo_folder_path, 'stitching-output') ;
vecfield_file_path = fullfile(stitching_output_folder_path, 'vecfield3D.mat') ;
mat = load(vecfield_file_path, 'vecfield3D') ;
vecfield = mat.vecfield3D ;
raw_final_affine_transform_from_tile_index = vecfield.tform ;  % 5 x 5 x tile_count
raw_baseline_affine_transform_from_tile_index = vecfield.afftile ;  % 3 x 4 x tile_count
targets_from_tile_index = vecfield.control ;
cpg_i0_values = vecfield.xlim_cntrl ;
cpg_j0_values = vecfield.ylim_cntrl ;
cpg_k0_values_from_tile_index = vecfield.zlim_cntrl ;

% All these transforms assume the zero-based indexing, but they also annoyingly
% assume that the the tiles are flipped in x and y relative to how you'd think
% they would be.  (B/c that's how mouselight tiles are acquired.)
% So the control point grid (CPG) i0 and j0 values are in this 'flipped'
% coordinate system.  I don't really want deal with these flipped coordinates,
% so we need to unflip them.
%
% Also the final transform is a 5x5 matrix for each tile, so want to fix that.

rare_final_affine_transform_from_tile_index = raw_final_affine_transform_from_tile_index([1 2 3 5],[1 2 3 5],:) ;
medium_final_affine_transform_from_tile_index = rare_final_affine_transform_from_tile_index(1:3,:,:) ;
A_flipped_per_tile = medium_final_affine_transform_from_tile_index(:,1:3,:) ;
b_flipped_per_tile = medium_final_affine_transform_from_tile_index(:,4,:) ;
S = diag([-1 -1 +1]) ;
n_vector = [tile_shape_ijk(1)-1 tile_shape_ijk(2)-1 0]' ;
A_per_tile = pagemtimes(A_flipped_per_tile, S) ;
b_per_tile = pagemtimes(A_flipped_per_tile, n_vector) + b_flipped_per_tile ;
final_affine_transform_from_tile_index = 1e-3 * horzcat(A_per_tile, b_per_tile) ;  % nm->um

% Compute the match errors for each tile, neighbor
final_affine_match_sse_from_neighbor_index_from_tile_index = ...
        compute_affine_landmark_match_error(self_ijk0_from_match_index_from_neighbor_index_from_tile_index, ...
                                            has_neighbor_from_neighbor_index_from_tile_index, ...
                                            neighbor_tile_index_from_neighbor_index_from_tile_index, ...
                                            neighbor_ijk0_from_match_idx_from_neighbor_idx_from_tile_idx, ...
                                            final_affine_transform_from_tile_index) ;

final_affine_z_match_sse_from_tile_index = reshape(final_affine_match_sse_from_neighbor_index_from_tile_index(3,:), [tile_count 1]) ;
final_affine_z_match_sse_from_pair_index = final_affine_z_match_sse_from_tile_index(has_z_plus_1_tile_from_tile_index) ;
final_affine_z_match_mse_from_pair_index = final_affine_z_match_sse_from_pair_index ./ z_match_count_from_pair_index ;
final_affine_z_match_rmse_from_pair_index = sqrt(final_affine_z_match_mse_from_pair_index) ;

total_final_affine_z_match_sse = sum(final_affine_z_match_sse_from_pair_index)
max_final_affine_z_match_rmse = max(final_affine_z_match_rmse_from_pair_index)
median_final_affine_z_match_rmse = median(final_affine_z_match_rmse_from_pair_index, 'omitnan')

tile_ijk1_from_pair_index = tile_ijk1_from_tile_index(self_tile_index_from_pair_index, :) ;

% make a z-match RMSE stack
final_affine_z_match_rmse_from_tile_ijk1 = nan(size(tile_index_from_tile_ijk1)) ;
for pair_index = 1 : pair_count ,
    tile_ijk1 = tile_ijk1_from_pair_index(pair_index,:) ;
    final_affine_z_match_rmse_from_tile_ijk1(tile_ijk1(1), tile_ijk1(2), tile_ijk1(3)) = final_affine_z_match_rmse_from_pair_index(pair_index) ;    
end

final_affine_z_match_rmse_from_tile_ijk1_montage = montage_from_stack_ijk(final_affine_z_match_rmse_from_tile_ijk1) ;
f = figure('color', 'w') ;
a = axes(f) ;
imagesc(final_affine_z_match_rmse_from_tile_ijk1_montage, [0 100]) ;
colorbar(a) ;
title('Final affine Z-Match RMSE (um)') 
drawnow

% scatter plot of RMSE and matches per pair
f = figure('color', 'w') ;
a = axes(f) ;
plot(a, z_match_count_from_pair_index, final_affine_z_match_rmse_from_pair_index, '.') ;
xlabel(a, 'Final affine Z-match count per tile pair') ;
ylabel(a, 'Final affine Z-match RMSE per tile pair') ;

final_affine_sse_ratio = total_final_affine_z_match_sse/total_stage_z_match_sse


