function result = file_name_from_path(path) 
    [~, name, ext] = fileparts(path) ;
    result = horzcat(name, ext) ;
end
