this_folder_path = fileparts(mfilename('fullpath')) ;
data_file_path = fullfile(this_folder_path, 'vectorfield3D-test-inputs-from-2020-09-15.mat') ;
load(data_file_path, 'params', 'scopeloc', 'regpts', 'scopeparams', 'curvemodel', 'tile_k_from_run_layer_index') ;
vecfield = vectorField3D(params, scopeloc, regpts, scopeparams, curvemodel, tile_k_from_run_layer_index) ;
