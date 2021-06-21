sample_tag = '2021-03-17-test-unscaled-region-1'  %#ok<NOPTS>
analysis_tag = 'adam-classifier' 
do_force_computation = false ;

raw_tile_path = '/groups/mousebrainmicro/mousebrainmicro/data/test-data/test-unscaled/2021-03-17-region1'
pipeline_output_folder_path = '/nrs/mouselight/pipeline_output/test-unscaled/2021-03-17-region-1-with-adam-classifier'
%line_fixed_tile_path = fullfile(pipeline_output_folder_path, 'stage_1_line_fix_output')
landmark_folder_path = fullfile(pipeline_output_folder_path, 'stage_3_descriptor_output')
match_folder_path = fullfile(pipeline_output_folder_path, 'stage_4_point_match_output')
this_folder_path = fileparts(mfilename('fullpath')) ;
sample_memo_folder_path = fullfile(this_folder_path, 'memos', sample_tag) ;
analysis_memo_folder_path = fullfile(sample_memo_folder_path, analysis_tag) ;

% Build an index of the paths to raw tiles
raw_tile_index = compute_or_read_from_memo(sample_memo_folder_path, ...
                                           'raw-tile-index', ...
                                           @()(build_raw_tile_index(raw_tile_path)), ...
                                           do_force_computation) ;
ijk1_from_tile_index = raw_tile_index.ijk1_from_tile_index ;
nominal_xyz_from_tile_index = raw_tile_index.xyz_from_tile_index ;  % um
relative_path_from_tile_index = raw_tile_index.relative_path_from_tile_index ;
tile_index_from_til_ijk1 = raw_tile_index.tile_index_from_tile_ijk1 ;

% Display some features of the raw tile index
tile_lattice_shape = size(tile_index_from_til_ijk1)
tile_count = length(relative_path_from_tile_index) 

% Read channel semantics
[neuron_channel_index, background_channel_index] = read_channel_semantics_file(raw_tile_path) ;
working_channel_index = background_channel_index ;

% Read in a single tile to get the tile shape
middle_tile_index = round(tile_count/2) ;
relative_path = relative_path_from_tile_index{middle_tile_index} ;
imagery_file_relative_path = imagery_file_relative_path_from_relative_path(relative_path, working_channel_index) ;
imagery_file_path = fullfile(raw_tile_path, imagery_file_relative_path) 
raw_tile_stack_yxz_flipped = read_16bit_grayscale_tif(imagery_file_path) ;
tile_shape_jik = size(raw_tile_stack_yxz_flipped) ;
tile_shape_ijk = tile_shape_jik([2 1 3]) ;

% Collect the landmarks for the background channel
ijk_from_landmark_index_from_tile_index = ...
    compute_or_read_from_memo(analysis_memo_folder_path, ...
                              sprintf('landmarks-channel-%d', working_channel_index), ...
                              @()(collect_landmarks(landmark_folder_path, relative_path_from_tile_index, working_channel_index, tile_shape_ijk)), ...
                              do_force_computation) ;

% Count the landmarks in each tile
landmark_count_from_tile_index = cellfun(@(a)(size(a,1)), ijk_from_landmark_index_from_tile_index) ;
max_landmark_count_from_tile_index = max(landmark_count_from_tile_index)

% Compute a histogram of those
bin_edges = 0 : 100 : 10e3 ;
bin_centers = ( bin_edges(1:end-1)+bin_edges(2:end) ) / 2 ;
counts_from_bin_index = histcounts(landmark_count_from_tile_index, bin_edges) ;
%assert(sum(counts_from_bin_index) == tile_count) ;
fraction_from_bin_index = counts_from_bin_index / tile_count ;

% Plot the histogram
f = figure('color', 'w') ;
a = axes(f) ;
b = bar(a, bin_centers, 100*fraction_from_bin_index) ;
b.EdgeColor = 'none' ;
xlabel('Landmarks per tile') ;
ylabel('Fraction of tiles (%)') ;
ylim([0 100]) ;
title(sample_tag) ;
f.Name = sprintf('%s-%s-landmarks-histogram', sample_tag, analysis_tag) ;

% Let's see a detail near zero
bin_edges = -0.5 : 1 : 100.5 ;
bin_centers = ( bin_edges(1:end-1)+bin_edges(2:end) ) / 2 ;
counts_from_bin_index = histcounts(landmark_count_from_tile_index, bin_edges) ;
%assert(sum(counts_from_bin_index) == tile_count) ;
fraction_from_bin_index = counts_from_bin_index / tile_count ;

% Plot the histogram
f = figure('color', 'w') ;
a = axes(f) ;
b = bar(a, bin_centers, 100*fraction_from_bin_index) ;
b.EdgeColor = 'none' ;
xlabel('Landmarks per tile') ;
ylabel('Fraction of tiles (%)') ;
ylim([0 1.4]) ;
title(sample_tag) ;
f.Name = sprintf('%s-%s-landmarks-histogram-detail', sample_tag, analysis_tag) ;



%%
%
% Matches
%

% Count the number of z-face pairs 
has_z_plus_1_tile_from_ijk1 = false(tile_lattice_shape) ;
has_z_plus_1_tile_from_ijk1(:,:,1:end-1) = isfinite(tile_index_from_til_ijk1(:,:,1:end-1)) & isfinite(tile_index_from_til_ijk1(:,:,2:end)) ;
pair_count = sum(sum(sum(has_z_plus_1_tile_from_ijk1))) 

% Want to know if has z+1 tile from tile_index
has_z_plus_1_tile_from_tile_index = false(tile_count,1) ;
for tile_index = 1 : tile_count ,
    ijk1 = ijk1_from_tile_index(tile_index,:) ;
    has_z_plus_1_tile = has_z_plus_1_tile_from_ijk1(ijk1(1), ijk1(2), ijk1(3)) ;
    has_z_plus_1_tile_from_tile_index(tile_index) = has_z_plus_1_tile ;
end
assert(sum(has_z_plus_1_tile_from_tile_index) == pair_count) ;

% Collect the z-face matches from disk
match_info = ...
    compute_or_read_from_memo(analysis_memo_folder_path, ...
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

% Compute a histogram of those
bin_edges = 0 : 50 : 2100 ;
bin_centers = ( bin_edges(1:end-1)+bin_edges(2:end) ) / 2 ;
counts_from_bin_index = histcounts(z_match_count_from_pair_index, bin_edges) ;
fraction_from_bin_index = counts_from_bin_index / pair_count ;

% Plot the histogram
f = figure('color', 'w') ;
a = axes(f) ;
b = bar(a, bin_centers, 100*fraction_from_bin_index) ;
b.EdgeColor = 'none' ;
xlabel('Matches per z-face pair') ;
ylabel('Fraction of pairs (%)') ;
ylim([0 100]) ;
title(sample_tag) ;
f.Name = sprintf('%s-%s-matches-histogram', sample_tag, analysis_tag) ;

% Let's see a detail near zero
bin_edges = -0.5 : 1 : 50.5 ;
bin_centers = ( bin_edges(1:end-1)+bin_edges(2:end) ) / 2 ;
counts_from_bin_index = histcounts(z_match_count_from_pair_index, bin_edges) ;
%assert(sum(counts_from_bin_index) == tile_count) ;
fraction_from_bin_index = counts_from_bin_index / pair_count ;

% Plot the histogram
f = figure('color', 'w') ;
a = axes(f) ;
b = bar(a, bin_centers, 100*fraction_from_bin_index) ;
b.EdgeColor = 'none' ;
xlabel('Matches per z-face pair') ;
ylabel('Fraction of pairs (%)') ;
ylim([0 100]) ;
title(sample_tag) ;
f.Name = sprintf('%s-%s-matches-histogram-detail', sample_tag, analysis_tag) ;

% Plot the number of matches vs the number of landmarks
landmark_count_from_pair_index = landmark_count_from_tile_index(has_z_plus_1_tile_from_tile_index) ;
f = figure('color', 'w') ;
a = axes(f) ;
plot(landmark_count_from_pair_index, z_match_count_from_pair_index, 'b.') ;
hold on ;
plot([0 10e3], [0 10e3], 'r') ;
xlabel('Landmarks per tile') ;
ylabel('Matches per tile (with z+1 tile)') ;
ylim([0 2000]) ;
xlim([0 10000]) ;
title(sample_tag) ;
f.Name = sprintf('%s-%s-matches-vs-landmarks', sample_tag, analysis_tag) ;

% What fraction of the landmarks end up as matches?
z_match_fraction_from_pair_index = z_match_count_from_pair_index ./ landmark_count_from_pair_index ;

% Compute a histogram of those
bin_edges = 0 : 0.01 : 1 ;
bin_centers = ( bin_edges(1:end-1)+bin_edges(2:end) ) / 2 ;
counts_from_bin_index = histcounts(z_match_fraction_from_pair_index, bin_edges) ;
fraction_from_bin_index = counts_from_bin_index / pair_count ;

% Plot the histogram
f = figure('color', 'w') ;
a = axes(f) ;
b = bar(a, bin_centers, 100*fraction_from_bin_index) ;
b.EdgeColor = 'none' ;
xlabel('Fraction of landmarks that become z-face matches per tile') ;
ylabel('Fraction of pairs (%)') ;
ylim([0 100]) ;
title(sample_tag) ;
f.Name = sprintf('%s-%s-match-fraction-histogram', sample_tag, analysis_tag) ;





