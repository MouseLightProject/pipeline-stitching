function stitch_2020_01_23()  
    addpath(genpath('/groups/mousebrainmicro/mousebrainmicro/Software/pipeline-stitching-shared/'))
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/2020-01-23/Tiling' %#ok<NOPRT>
    pipeline_output_folder = '/nrs/mouselight/pipeline_output/2020-01-23' %#ok<NOPRT>
    stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2020-01-23/stitching-output_TF1'  %#ok<NOPRT>
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
