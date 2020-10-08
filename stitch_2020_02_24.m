function stitch_2020_02_24()  
    addpath(genpath('/groups/mousebrainmicro/mousebrainmicro/Software/pipeline-stitching-shared/'))
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/acquisition/2020-04-15' %#ok<NOPRT>
    pipeline_output_folder = '/nrs/mouselight/pipeline_output/2020-04-15' %#ok<NOPRT>
    stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2020-04-15/stitching-output_TF1'  %#ok<NOPRT>
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
