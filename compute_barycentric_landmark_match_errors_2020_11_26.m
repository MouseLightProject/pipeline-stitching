sample_tag = '2020-11-26'  %#ok<NOPTS>
analysis_tag = 'production-classifier-z-match-count-threshold-50' 
path_kind_to_use_for_imagery = 'line-fixed' 
do_force_computation = true ;

raw_tile_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling', sample_tag) ;
pipeline_output_folder_path = sprintf('/nrs/mouselight/pipeline_output/%s', sample_tag)  %#ok<NOPTS> 
line_fixed_tile_path = fullfile(pipeline_output_folder_path, 'stage_1_line_fix_output')
landmark_folder_path = fullfile(pipeline_output_folder_path, 'stage_3_descriptor_output')
match_folder_path = fullfile(pipeline_output_folder_path, 'stage_4_point_match_output')
this_folder_path = fileparts(mfilename('fullpath')) ;
sample_memo_folder_path = fullfile(this_folder_path, 'memos', sample_tag) ;
analysis_memo_folder_path = fullfile(sample_memo_folder_path, analysis_tag) ;

if ~exist('path_kind_to_use_for_imagery', 'var') || isempty(path_kind_to_use_for_imagery) || strcmp(path_kind_to_use_for_imagery, 'line-fixed') ,
    imagery_tile_path = line_fixed_tile_path ;
elseif strcmp(path_kind_to_use_for_imagery, 'raw') ,
    imagery_tile_path = raw_tile_path ;
else
    error('Variable path_to_use_for_imagery is ''%s''.  Should be either ''line-fixed'' or ''raw''.', path_kind_to_use_for_imagery) ;
end

% Build an index of the paths to raw tiles
raw_tile_index = compute_or_read_from_memo(sample_memo_folder_path, ...
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
    compute_or_read_from_memo(sample_memo_folder_path, ...
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
        sample_memo_folder_path, ...
        sprintf('z-face-matches-channel-%d', background_channel_index), ...
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



%
% Compute things we'll need for all the error computations
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



%%
%
% Do the nominal afffine
%

% Assemble the stage-based affine transforms
stage_affine_transform_from_tile_index = zeros(3, 4, tile_count) ;
stage_affine_transform_from_tile_index(1,1,:) = spacing_um_xyz(1) ;
stage_affine_transform_from_tile_index(2,2,:) = spacing_um_xyz(2) ;
stage_affine_transform_from_tile_index(3,3,:) = spacing_um_xyz(3) ;
stage_affine_transform_from_tile_index(:,4,:) = nominal_xyz_from_tile_index' ;

% Compute the match errors for each tile, neighbor
stage_dist_from_match_idx_from_neighbor_idx_from_tile_idx = ...
        compute_affine_landmark_match_distance(self_ijk0_from_match_index_from_neighbor_index_from_tile_index, ...
                                               has_neighbor_from_neighbor_index_from_tile_index, ...
                                               neighbor_tile_index_from_neighbor_index_from_tile_index, ...
                                               neighbor_ijk0_from_match_idx_from_neighbor_idx_from_tile_idx, ...
                                               stage_affine_transform_from_tile_index) ;
stage_z_match_sse_from_neighbor_index_from_tile_index = cellfun(@(v)(sum(v.^2)), stage_dist_from_match_idx_from_neighbor_idx_from_tile_idx) ;
stage_z_match_sse_from_tile_index = reshape(stage_z_match_sse_from_neighbor_index_from_tile_index(3,:), [tile_count 1]) ;
stage_z_match_sse_from_pair_index = stage_z_match_sse_from_tile_index(has_z_plus_1_tile_from_tile_index) ;
stage_z_match_mse_from_pair_index = stage_z_match_sse_from_pair_index ./ z_match_count_from_pair_index ;

total_stage_z_match_sse = sum(stage_z_match_sse_from_pair_index)
stage_z_match_rmse_from_pair_index = sqrt(stage_z_match_mse_from_pair_index) 
max_stage_z_match_rmse = max(stage_z_match_rmse_from_pair_index)
median_stage_z_match_rmse = median(stage_z_match_rmse_from_pair_index, 'omitnan')

tile_ijk1_from_pair_index = tile_ijk1_from_tile_index(self_tile_index_from_pair_index, :) ;

% make a z-match RMSE stack
stage_z_match_count_from_tile_ijk1 = nan(size(tile_index_from_tile_ijk1)) ;
stage_z_match_rmse_from_tile_ijk1 = nan(size(tile_index_from_tile_ijk1)) ;
for pair_index = 1 : pair_count ,
    tile_ijk1 = tile_ijk1_from_pair_index(pair_index,:) ;
    stage_z_match_rmse_from_tile_ijk1(tile_ijk1(1), tile_ijk1(2), tile_ijk1(3)) = stage_z_match_rmse_from_pair_index(pair_index) ;    
    stage_z_match_count_from_tile_ijk1(tile_ijk1(1), tile_ijk1(2), tile_ijk1(3)) = z_match_count_from_pair_index(pair_index) ;    
end

stage_z_match_rmse_from_tile_ijk1_montage = montage_from_stack_ijk(stage_z_match_rmse_from_tile_ijk1) ;
f = figure('color', 'w') ;
a = axes(f) ;
imagesc(stage_z_match_rmse_from_tile_ijk1_montage, [0 100]) ;
colorbar(a) ;
title('Stage Z-Match RMSE (um)') 
drawnow

stage_z_match_count_from_tile_ijk1_montage = montage_from_stack_ijk(stage_z_match_count_from_tile_ijk1) ;
f = figure('color', 'w') ;
a = axes(f) ;
imagesc(stage_z_match_count_from_tile_ijk1_montage, [0 600]) ;
colorbar(a) ;
title('Stage Z-Match Count') 
drawnow

% scatter plot of RMSE and matches per pair
f = figure('color', 'w') ;
a = axes(f) ;
plot(a, z_match_count_from_pair_index, stage_z_match_rmse_from_pair_index, '.') ;
xlabel(a, 'Stage Z-match count per tile pair') ;
ylabel(a, 'Stage Z-match RMSE per tile pair') ;



%%
%
% Do the final affine
%

% Load the per-tile affine transforms the stitcher outputs
stitching_output_folder_path = fullfile(analysis_memo_folder_path, 'stitching-output') ;
vecfield_file_path = fullfile(stitching_output_folder_path, 'vecfield3D.mat') ;
mat = load(vecfield_file_path, 'vecfield3D') ;
vecfield = mat.vecfield3D ;
raw_final_affine_transform_from_tile_index = vecfield.tform ;  % 5 x 5 x tile_count
raw_baseline_affine_transform_from_tile_index = vecfield.afftile ;  % 3 x 4 x tile_count
flipped_cpg_i0_values = vecfield.xlim_cntrl ;
cpg_i0_values = tile_shape_ijk(1) - 1 - flipped_cpg_i0_values ;
flipped_cpg_j0_values = vecfield.ylim_cntrl ;
cpg_j0_values = tile_shape_ijk(2) - 1 - flipped_cpg_j0_values ;
cpg_k0_values_from_tile_index = vecfield.zlim_cntrl ;
cpg_shape_ijk = [ length(cpg_i0_values) length(cpg_j0_values) size(cpg_k0_values_from_tile_index,1) ] ;
targets_from_tile_index = 1e-3 * vecfield.control ;  % nm -> um
  % 100 x 3 x tile_count, 100 == 5*5*4, the number of grid points in x, y, z
  % These are in "z-major" order: x changes most quickly, y next most quickly, and
  % z changes most slowly.  It would probably make more sense to have this be 
  % 3 x (x grid count) x (y grid count) x (z grid count) x tile_count
  % Also note that they're in order of the "flipped" i0 and j0 coords.
  
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
final_affin_dist_from_match_idx_from_neighbor_idx_from_tile_idx = ...
        compute_affine_landmark_match_distance(self_ijk0_from_match_index_from_neighbor_index_from_tile_index, ...
                                               has_neighbor_from_neighbor_index_from_tile_index, ...
                                               neighbor_tile_index_from_neighbor_index_from_tile_index, ...
                                               neighbor_ijk0_from_match_idx_from_neighbor_idx_from_tile_idx, ...
                                               final_affine_transform_from_tile_index) ;

final_affine_match_sse_from_neighbor_index_from_tile_index = cellfun(@(v)(sum(v.^2)), final_affin_dist_from_match_idx_from_neighbor_idx_from_tile_idx) ;                                          
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



%%
%
% Use full barycentric transform
%

% Compute the match errors for each tile, neighbor
tic_id = tic() ;
[bary_distance_from_match_idx_from_neighbor_idx_from_tile_idx, is_within_cpg_from_match_idx_from_neighbor_idx_from_tile_idx] = ...
        compute_barycentric_landmark_match_distance(self_ijk0_from_match_index_from_neighbor_index_from_tile_index, ...
                                                    has_neighbor_from_neighbor_index_from_tile_index, ...
                                                    neighbor_tile_index_from_neighbor_index_from_tile_index, ...
                                                    neighbor_ijk0_from_match_idx_from_neighbor_idx_from_tile_idx, ...
                                                    cpg_i0_values , ...
                                                    cpg_j0_values , ...
                                                    cpg_k0_values_from_tile_index , ...
                                                    targets_from_tile_index) ;
elapsed_time = toc(tic_id) ;
fprintf('Time to compute barycentric landmark match errors: %g sec\n', elapsed_time) ;

% barycentric_z_match_sse_from_neighbor_index_from_tile_index = cellfun(@(v)(sum(v.^2)), bary_distance_from_match_idx_from_neighbor_idx_from_tile_idx) ;
% barycentric_match_count_from_neighbor_idx_from_tile_idx = cellfun(@(v)(sum(v)), is_within_cpg_from_match_idx_from_neighbor_idx_from_tile_idx) ;
% 
% barycentric_z_match_sse_from_tile_index = reshape(sum(barycentric_z_match_sse_from_neighbor_index_from_tile_index), [tile_count 1]) ;
% barycentric_match_count_from_tile_index = reshape(sum(barycentric_match_count_from_neighbor_idx_from_tile_idx), [tile_count 1]) ;
% barycentric_z_match_sse_from_pair_index = barycentric_z_match_sse_from_tile_index(has_z_plus_1_tile_from_tile_index) ;
% barycentric_z_match_count_from_pair_index = barycentric_match_count_from_tile_index(has_z_plus_1_tile_from_tile_index) ;
% barycentric_z_match_mse_from_pair_index = barycentric_z_match_sse_from_pair_index ./ barycentric_z_match_count_from_pair_index ;
% barycentric_z_match_rmse_from_pair_index = sqrt(barycentric_z_match_mse_from_pair_index) ;
% 
% total_barycentric_z_match_sse = sum(barycentric_z_match_sse_from_pair_index)
% max_barycentric_z_match_rmse = max(barycentric_z_match_rmse_from_pair_index)
% median_barycentric_z_match_rmse = median(barycentric_z_match_rmse_from_pair_index, 'omitnan')
% 
% tile_ijk1_from_pair_index = tile_ijk1_from_tile_index(self_tile_index_from_pair_index, :) ;
% 
% % make a z-match RMSE stack
% barycentric_z_match_rmse_from_tile_ijk1 = nan(size(tile_index_from_tile_ijk1)) ;
% for pair_index = 1 : pair_count ,
%     tile_ijk1 = tile_ijk1_from_pair_index(pair_index,:) ;
%     barycentric_z_match_rmse_from_tile_ijk1(tile_ijk1(1), tile_ijk1(2), tile_ijk1(3)) = barycentric_z_match_rmse_from_pair_index(pair_index) ;    
% end

[barycentric_in_cpg_total_ssd, barycentric_in_cpg_ssd_from_tile_index, barycentric_in_cpg_match_count_from_tile_index, ...
 barycentric_in_cpg_rmse_from_tile_index, barycentric_in_cpg_rmse_from_tile_ijk1] = ...
    compute_match_statistics(bary_distance_from_match_idx_from_neighbor_idx_from_tile_idx, ...
                             is_within_cpg_from_match_idx_from_neighbor_idx_from_tile_idx, ...
                             tile_ijk1_from_tile_index, ...
                             tile_index_from_tile_ijk1) ;


% barycentric_z_match_rmse_from_tile_ijk1_montage = montage_from_stack_ijk(barycentric_z_match_rmse_from_tile_ijk1) ;
% f = figure('color', 'w') ;
% a = axes(f) ;
% imagesc(barycentric_z_match_rmse_from_tile_ijk1_montage, [0 100]) ;
% colorbar(a) ;
% title('Final barycentric Z-Match RMSE (um)') 
% drawnow
% 
% % scatter plot of RMSE and matches per pair
% f = figure('color', 'w') ;
% a = axes(f) ;
% plot(a, barycentric_z_match_count_from_pair_index, barycentric_z_match_rmse_from_pair_index, '.') ;
% xlabel(a, 'Final barycentric Z-match count per tile pair') ;
% ylabel(a, 'Final barycentric Z-match RMSE per tile pair') ;

make_rmse_plots(barycentric_in_cpg_rmse_from_tile_index, ...
                barycentric_in_cpg_match_count_from_tile_index, ...
                barycentric_in_cpg_rmse_from_tile_ijk1, ...
                'Barycentric') ;
                         


%
% Want to re-compute the total SSE for stage transform, final affine, on just
% the in-CPG points.  This way we can compute how much the total squared
% distance has come down just on the matches for which we have a barycentric
% distance.
%
[stage_in_cpg_total_ssd, stage_in_cpg_ssd_from_tile_index, stage_in_cpg_match_count_from_tile_index, ...
 stage_in_cpg_rmse_from_tile_index, stage_in_cpg_rmse_from_tile_ijk1] = ...
    compute_match_statistics(stage_dist_from_match_idx_from_neighbor_idx_from_tile_idx, ...
                             is_within_cpg_from_match_idx_from_neighbor_idx_from_tile_idx, ...
                             tile_ijk1_from_tile_index, ...
                             tile_index_from_tile_ijk1) ;
[final_affine_in_cpg_total_ssd, final_affine_in_cpg_ssd_from_tile_index, final_affine_in_cpg_match_count_from_tile_index, ...
 final_affine_in_cpg_rmse_from_tile_index, final_affine_in_cpg_rmse_from_tile_ijk1] = ...
    compute_match_statistics(final_affin_dist_from_match_idx_from_neighbor_idx_from_tile_idx, ...
                             is_within_cpg_from_match_idx_from_neighbor_idx_from_tile_idx, ...
                             tile_ijk1_from_tile_index, ...
                             tile_index_from_tile_ijk1) ;
                                      
final_affine_in_cpg_ssd_ratio = final_affine_in_cpg_total_ssd/stage_in_cpg_total_ssd                                            
barycentric_ssd_ratio = total_barycentric_z_match_sse/stage_in_cpg_total_ssd

make_rmse_plots(stage_in_cpg_rmse_from_tile_index, ...
                stage_in_cpg_match_count_from_tile_index, ...
                stage_in_cpg_rmse_from_tile_ijk1, ...
                'In-CPG Stage') ;
            
make_rmse_plots(final_affine_in_cpg_rmse_from_tile_index, ...
                final_affine_in_cpg_match_count_from_tile_index, ...
                final_affine_in_cpg_rmse_from_tile_ijk1, ...
                'In-CPG Final affine') ;

% Make a montage of the ratio of the barycentric error to the final afffine
final_affine_to_barycentric_mse_ratio_from_tile_ijk1 = (final_affine_in_cpg_rmse_from_tile_ijk1.^2) ./ (barycentric_in_cpg_rmse_from_tile_ijk1.^2) ;
final_affine_to_barycentric_mse_ratio_from_tile_ijk1_montage = montage_from_stack_ijk(final_affine_to_barycentric_mse_ratio_from_tile_ijk1) ;
f1 = figure('color', 'w') ;
a1 = axes(f1) ;
imagesc(final_affine_to_barycentric_mse_ratio_from_tile_ijk1_montage, [0 10]) ;
colorbar(a1) ;
title(sprintf('Final affine MSE over barycentric')) ; 
drawnow() ;

            


