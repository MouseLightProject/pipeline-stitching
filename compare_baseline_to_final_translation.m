this_folder_path = fileparts(mfilename('fullpath')) ;
vecfield_file_path = fullfile(this_folder_path, '2020-09-15-stitching-output-verify/vecfield3D.mat') ;
mat = load(vecfield_file_path, 'vecfield3D', 'params') ;
vecfield = mat.vecfield3D ;

baseline_affine_transform_from_tile_index = vecfield.afftile ;  % 3 x 4 x tile_count
raw_final_affine_transform_from_tile_index = vecfield.tform ;  % 5 x 5 x tile_count
final_affine_transform_from_tile_index = raw_final_affine_transform_from_tile_index([1 2 3],[1 2 3 5],:) ;  % 3 x 4 x tile_count
tile_count = size(baseline_affine_transform_from_tile_index,3) ;

raw_baseline_translation_from_tile_index = 1e-3 * squeeze(baseline_affine_transform_from_tile_index(:,4,:)) ;  % 3 x tile_count, um
min_baseline_translation = min(raw_baseline_translation_from_tile_index,[],2) ;
baseline_translation_from_tile_index = raw_baseline_translation_from_tile_index - min_baseline_translation ;
final_translation_from_tile_index = 1e-3 * squeeze(final_affine_transform_from_tile_index(:,4,:)) - min_baseline_translation ;  % 3 x tile_count, um

translation_change_from_tile_index = final_translation_from_tile_index - baseline_translation_from_tile_index ;  % um

mean_translation_change = mean(translation_change_from_tile_index,2)
median_translation_change = median(translation_change_from_tile_index,2)

translation_change_mangnitude_from_tile_index = vecnorm(translation_change_from_tile_index) ;

median_translation_change_magnitude = median(translation_change_mangnitude_from_tile_index,2)
max_translation_change_magnitude = max(translation_change_mangnitude_from_tile_index)

bin_edges = 0 : 1 : (ceil(max_translation_change_magnitude)+1) ;
bin_centers = ( bin_edges(1:end-1)+bin_edges(2:end) ) / 2 ;

counts_from_bin_index = histcounts(translation_change_mangnitude_from_tile_index, bin_edges) ;
assert(sum(counts_from_bin_index) == tile_count) ;
fraction_from_bin_index = counts_from_bin_index / tile_count ;
f = figure('color', 'w') ;
a = axes(f) ;
b = bar(a, bin_centers, 100*fraction_from_bin_index) ;
b.EdgeColor = 'none' ;
xlabel('Translation change length (um)') ;
ylabel('Fraction of tiles (%)') ;


