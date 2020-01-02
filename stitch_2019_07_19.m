function stitch_2019_08_08()  
    addpath(genpath('/groups/mousebrainmicro/mousebrainmicro/Software/pipeline-stitching-shared/'))
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/acquisition/2019-08-08' %#ok<NOPRT>
    pipeline_output_folder =  '/nrs/mouselight/pipeline_output/2019-07-19' %#ok<NOPRT>
    stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2019-07-19/stitching-output_TF'  %#ok<NOPRT>
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
