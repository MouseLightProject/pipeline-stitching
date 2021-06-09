% Specify the input/output folders
sample_tag = '2021-03-17-test-unscaled'  %#ok<NOPTS>
analysis_tag = 'adam-classifier' 
do_force_computations = false
do_perform_field_correction = true
do_run_in_debug_mode = false
do_show_visualizations = true
this_folder_path = fileparts(mfilename('fullpath')) ;
tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/test-data/test-unscaled/2021-03-17'
pipeline_output_folder = '/nrs/mouselight/pipeline_output/test-unscaled/2021-03-17-adam-classifier'
sample_memo_folder_path = fullfile(this_folder_path, sprintf('memos-%s', sample_tag)) ;
analysis_memo_folder_path = fullfile(sample_memo_folder_path, analysis_tag) ;
stitching_output_folder_path = fullfile(analysis_memo_folder_path, 'stitching-output') ;

% % Call the function that does the real work
% stitch_options = struct('do_force_computations', do_force_computations, ...
%                         'do_perform_field_correction', do_perform_field_correction, ...
%                         'do_run_in_debug_mode', do_run_in_debug_mode, ...
%                         'do_show_visualizations', do_show_visualizations) ;
% stitch(tile_folder_path, ...
%        pipeline_output_folder, ...
%        stitching_output_folder_path, ...
%        stitch_options) ;

% Look at the between-tile shifts, as derived from the FC fits
scope_params_per_tile_file_path = fullfile(stitching_output_folder_path,'scopeparams_pertile.mat') ;
load(scope_params_per_tile_file_path, 'scopeparams', 'curvemodel', 'params') ;
scopeloc_file_path = fullfile(stitching_output_folder_path, 'scopeloc.mat') ;
load(scopeloc_file_path, 'scopeloc') ;

% Get the neighbor matrix
raw_tile_ijk_from_tile_index = scopeloc.gridix(:,1:3)
tile_ijk_from_tile_index = raw_tile_ijk_from_tile_index - min(raw_tile_ijk_from_tile_index, [], 1) + 1 
neighbors = buildNeighbor(tile_ijk_from_tile_index) %[id -x -y +x +y -z +z] format    
neighbor_tile_index_from_axis_index_from_tile_index = neighbors(:,[4 5 7])  % just the +1 tiles

% Compute the nominal tile offsets
nominal_tile_xyz_from_tile_index = 1000 * scopeloc.loc  % mm->um
tile_count = size(nominal_tile_xyz_from_tile_index, 1)
axis_count = 3 ;
nominal_tile_dxyz_from_axis_index_from_tile_index = nan(axis_count, axis_count, tile_count) ;
for tile_index = 1 : tile_count ,
    tile_xyz = nominal_tile_xyz_from_tile_index(tile_index,:) ;
    neighbor_tile_index_from_axis_index = neighbor_tile_index_from_axis_index_from_tile_index(tile_index,:) ;
    for axis_index = 1 : axis_count ,
        neighbor_tile_index = neighbor_tile_index_from_axis_index(axis_index) ;
        if isfinite(neighbor_tile_index) ,
            neighbor_tile_xyz = nominal_tile_xyz_from_tile_index(neighbor_tile_index,:) ;
            dxyz = neighbor_tile_xyz - tile_xyz ;  % um
            nominal_tile_dxyz_from_axis_index_from_tile_index(:, axis_index, tile_index) = dxyz' ;
        end
    end
end

spacing = params.imsize_um ./ (params.imagesize-1)

nominal_tile_dijk_from_axis_index_from_tile_index = nominal_tile_dxyz_from_axis_index_from_tile_index ./ spacing' 
