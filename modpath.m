function modpath()
    path_to_this_file = mfilename('fullpath') ;
    path_to_this_folder = fileparts(path_to_this_file) ;
    addpath(genpath(fullfile(path_to_this_folder, 'common'))) ;
    addpath(genpath(fullfile(path_to_this_folder, 'functions'))) ;
end


