function write_slide_show(slide_show_folder_path, F)
    if ~exist(slide_show_folder_path) ,
        mkdir(slide_show_folder_path) ;
    end
    if isstruct(F) ,
        frame_count = length(F) ;
        digit_count = ceil(log10(frame_count)) ;
        file_name_template = sprintf('frame-%%0%dd.png', digit_count) ;
        for i = 1 : frame_count ,
            frame = F(i) ;
            frame_file_name = sprintf(file_name_template, i) ;
            frame_file_path = fullfile(slide_show_folder_path, frame_file_name) ;
            imwrite(frame.cdata, frame_file_path, 'png') ;
        end        
    else
        error('Don''t know how to make a slide show of that') ;
    end
end
