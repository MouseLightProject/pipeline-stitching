function stitch_2019_05_27()    
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/2019-05-27/Tiling' %#ok<NOPRT>
    pipeline_output_folder =  '/nrs/mouselight/pipeline_output/2019-05-27' %#ok<NOPRT>
    stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2019-05-27/stitching-output_TF20190726'  %#ok<NOPRT>
    %%stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
