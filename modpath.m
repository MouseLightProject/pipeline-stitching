function modpath()
    path_to_this_file = mfilename('fullpath') ;
    path_to_this_folder = fileparts(path_to_this_file) ;
    addpath(genpath(fullfile(path_to_this_folder, 'common'))) ;
    addpath(genpath(fullfile(path_to_this_folder, 'functions'))) ;
    %addpath(fullfile(path_to_this_folder, 'mouselight_toolbox')) ;
    
    toolbox_path = fullfile(path_to_this_folder, 'tmt') ;
    modpath_path = fullfile(toolbox_path, 'modpath.m') ;
    run(modpath_path) ;    
end
