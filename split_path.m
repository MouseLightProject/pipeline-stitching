function [parent, file_name] = split_path(path) 
    [parent, base, ext] = fileparts(path) ;
    file_name = horzcat(base, ext) ;
end
