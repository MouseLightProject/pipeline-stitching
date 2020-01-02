function stitch_2019_05_27_test()    
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/2019-05-27/Tiling' %#ok<NOPRT>
    pipeline_output_folder =  '/nrs/mouselight/pipeline_output/2019-05-27' %#ok<NOPRT>
    stitching_output_folder_path = '/groups/mousebrainmicro/mousebrainmicro/Software/pipeline-stitching-shared/2019-05-27-test-stitching-output'  %#ok<NOPRT>  % typically run on canopus
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
