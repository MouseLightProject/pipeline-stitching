function stitch_2019_04_17_test()    
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/2019-04-17/Tiling' %#ok<NOPRT>
    pipeline_output_folder =  '/nrs/mouselight/pipeline_output/2019-04-17' %#ok<NOPRT>
    stitching_output_folder_path = '/data/taylora/2019-04-17/stitching-output'  %#ok<NOPRT>  % typically run on canopus
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
