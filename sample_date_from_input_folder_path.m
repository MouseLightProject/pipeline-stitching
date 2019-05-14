function result = sample_date_from_input_folder_path(input_folder_path) 
    result_or_empty = sample_date_from_input_folder_path_helper(input_folder_path) ;
    if isempty(result_or_empty) ,
        error('Unable to determine sample date for path %s', input_folder_path) ;
    else
        result = result_or_empty ;
    end
end


function result = sample_date_from_input_folder_path_helper(input_folder_path) 
    if isempty(input_folder_path) || isequal(input_folder_path, '/') ,
        % Down to the root, so give up        
        result = '' ;
    else
        [parent_path, leaf_name] = split_path(input_folder_path) ;
        if is_date(leaf_name) ,
            result = leaf_name ;
        else
            % Recurse, hopefully will find a date farther up the path
            % hierarchy
            result = sample_date_from_input_folder_path_helper(parent_path) ;
        end
    end
end


function result = is_date(str)
    parts = strsplit(str, '-') ;
    if length(parts)==3 ,
        year_as_string = parts{1} ;
        month_as_string = parts{2} ;
        day_as_string = parts{3} ;
        result = is_year(year_as_string) && is_month(month_as_string) && is_day(day_as_string) ;
    else
        result = false ;
    end
end


function result = is_year(str)
    if length(str)==4 ,
        year = str2double(str) ;
        result = ~isnan(year) && round(year)==year ;
    else
        result = false ;
    end
end


function result = is_month(str)
    if length(str)==2 ,
        month = str2double(str) ;
        result = ~isnan(month) && round(month)==month && 1<=month && month<=12 ;
    else
        result = false ;
    end
end


function result = is_day(str)
    if length(str)==2 ,
        day = str2double(str) ;
        result = ~isnan(day) && round(day)==day && 1<=day && day<=31 ;
    else
        result = false ;
    end
end
