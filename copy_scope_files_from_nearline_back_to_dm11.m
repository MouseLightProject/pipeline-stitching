function copy_scope_files_from_nearline_back_to_dm11(dest_raw_tiles_path, source_raw_tiles_path)
    
    function updated_state = my_callback(source_root_folder_absolute_path, current_folder_relative_path, folder_names, file_names, state)  %#ok<INUSL>
        source_folder_absolute_path = fullfile(source_root_folder_absolute_path, current_folder_relative_path) ;
        dest_folder_absolute_path = fullfile(dest_raw_tiles_path, current_folder_relative_path) ;
        ensure_folder_exists(dest_folder_absolute_path) ;
        
        % Filter out .tif, .mj2 files from the file names
        is_a_stack_from_file_name_index = cellfun(@is_a_stack, file_names) ;
        do_copy_from_file_name_index = ~is_a_stack_from_file_name_index ;
        file_names_to_copy = file_names(do_copy_from_file_name_index) ;
        
        for i = 1 : length(file_names_to_copy) ,
            file_name = file_names_to_copy{i} ;
            source_file_absolute_path = fullfile(source_folder_absolute_path, file_name) ;
            dest_file_absolute_path = fullfile(dest_folder_absolute_path, file_name) ;
            if ~exist(dest_file_absolute_path, 'file') ,
                %dest_file_absolute_path
                system_from_list_with_error_handling({'cp', '-T', source_file_absolute_path, dest_file_absolute_path}) ;
            end
        end
        
        updated_state = state ;  % not used
    end
   
    dirwalk(source_raw_tiles_path, @my_callback, []) ;    
end



function result = is_a_stack(file_name)
    result = contains(file_name, '.mj2') || contains(file_name, '.tif') ;
end
